"""Run the Replica audit and write every artifact used by the dashboard/API."""

from __future__ import annotations

import json
from datetime import datetime, timezone

import joblib
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from scipy.stats import kruskal

from replica.config import PATHS, load_config
from replica.data import build_cohorts, expression_with_clinical
from replica.models import fit_nmf, module_scores, top_module_genes
from replica.stats import bh, by, npc_test
from replica.survival import (
    bootstrap_cindex,
    bootstrap_delta_cindex,
    cross_validated_risk,
    design_matrix,
    fit_cox,
    joint_permutation_logrank,
    km_curve,
    likelihood_ratio,
    logrank_statistic,
    permutation_logrank,
    risk_score,
)

SEED = 20261123
CONFIG = load_config()
B = int(CONFIG["permutation"]["n_primary"])


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.ndarray,)):
        return value.tolist()
    if isinstance(value, (pd.Timestamp,)):
        return value.isoformat()
    return value


def endpoint_cols(cohort: str, endpoint: str) -> tuple[str, str, str]:
    if endpoint == "os":
        return "os_time", "os_event", "OS"
    if cohort == "tcga":
        return "dfs_time", "dfs_event", "DFS"
    return "dfs_time", "dfs_event", "RFS"


def cohort_subset(frame: pd.DataFrame, cohort: str, endpoint: str) -> pd.DataFrame:
    flag = "analysis_os" if endpoint == "os" else "analysis_secondary"
    return frame[frame[flag]].copy()


def subtype_counts(frame: pd.DataFrame, cohort: str) -> pd.DataFrame:
    out = (
        frame.groupby("subtype4", dropna=False)
        .agg(n=("patient_id", "size"), events_os=("os_event", "sum"))
        .reset_index()
    )
    out["cohort"] = cohort
    return out


def unadjusted_scorecard(frame: pd.DataFrame, cohort: str, endpoint: str) -> dict:
    time_col, event_col, label = endpoint_cols(cohort, endpoint)
    data = cohort_subset(frame, cohort, endpoint)
    data = data.dropna(subset=[time_col, event_col, "subtype4"])
    data = data[data[time_col] > 0]
    if len(data) < 30:
        return {"cohort": cohort, "endpoint": label, "n": len(data), "error": "too few samples"}
    lr = permutation_logrank(
        data[time_col].to_numpy(),
        data[event_col].to_numpy(),
        data["subtype4"].to_numpy(),
        n_perm=B,
        seed=SEED,
    )
    fit = fit_cox(data, time_col, event_col, [], subtype_col="subtype4")
    risk = risk_score(fit, data)
    ci = bootstrap_cindex(data[time_col].to_numpy(), data[event_col].to_numpy(), risk.to_numpy(), n_boot=1000, seed=SEED)
    summary = fit.summary()
    hrs = []
    for term in summary.index:
        hrs.append(
            {
                "term": term,
                "hr": float(summary.loc[term, "hr"]),
                "lower": float(summary.loc[term, "hr_lower"]),
                "upper": float(summary.loc[term, "hr_upper"]),
                "p": float(summary.loc[term, "p"]),
            }
        )
    return {
        "cohort": cohort,
        "endpoint": label,
        "n": int(len(data)),
        "events": int(data[event_col].sum()),
        "logrank_stat": float(lr["statistic"]),
        "df": int(lr["df"]),
        "p_value": float(lr["p_value"]),
        "cindex": ci["cindex"],
        "cindex_ci_lower": ci["ci_lower"],
        "cindex_ci_upper": ci["ci_upper"],
        "hrs": hrs,
        "null": lr["null"],
    }


def incremental_value(frame: pd.DataFrame, cohort: str, endpoint: str, base_covariates: list[str]) -> dict:
    time_col, event_col, label = endpoint_cols(cohort, endpoint)
    data = cohort_subset(frame, cohort, endpoint).dropna(
        subset=[time_col, event_col, *base_covariates, "subtype4"]
    )
    data = data[data[time_col] > 0]
    base_fit = fit_cox(data, time_col, event_col, base_covariates, subtype_col=None)
    full_fit = fit_cox(data, time_col, event_col, base_covariates, subtype_col="subtype4")
    lr = likelihood_ratio(full_fit, base_fit)

    risk_base = cross_validated_risk(data, time_col, event_col, base_covariates, subtype_col=None, seed=SEED)
    risk_full = cross_validated_risk(data, time_col, event_col, base_covariates, subtype_col="subtype4", seed=SEED)
    mask = risk_base.notna() & risk_full.notna()
    y_time = data.loc[mask, time_col].to_numpy()
    y_event = data.loc[mask, event_col].to_numpy()
    ci_base = bootstrap_cindex(y_time, y_event, risk_base[mask].to_numpy(), n_boot=1000, seed=SEED)
    ci_full = bootstrap_cindex(y_time, y_event, risk_full[mask].to_numpy(), n_boot=1000, seed=SEED)
    delta = bootstrap_delta_cindex(
        y_time, y_event, risk_base[mask].to_numpy(), risk_full[mask].to_numpy(), n_boot=1000, seed=SEED
    )
    return {
        "cohort": cohort,
        "endpoint": label,
        "n": int(mask.sum()),
        "events": int(y_event.sum()),
        "cindex_base": ci_base["cindex"],
        "cindex_base_lower": ci_base["ci_lower"],
        "cindex_base_upper": ci_base["ci_upper"],
        "cindex_full": ci_full["cindex"],
        "cindex_full_lower": ci_full["ci_lower"],
        "cindex_full_upper": ci_full["ci_upper"],
        "delta_cindex": delta["delta_cindex"],
        "delta_lower": delta["ci_lower"],
        "delta_upper": delta["ci_upper"],
        "lr_stat": lr["lr_stat"],
        "lr_df": lr["df"],
        "lr_p": lr["p_value"],
        "base_cox": base_fit,
        "full_cox": full_fit,
    }


def calibration_table(fit, frame: pd.DataFrame, time_col: str, event_col: str, horizon: float = 60.0) -> pd.DataFrame:
    design = design_matrix(frame, fit.base_covariates, subtype_col=fit.subtype_col, reference=fit.reference)
    design = design.reindex(columns=fit.covariates).astype(float).fillna(0.0)
    sf = fit.model.predict_survival_function(design, times=[horizon])
    predicted = sf.iloc[0].to_numpy()
    out = frame[[time_col, event_col]].copy()
    out["predicted_survival"] = predicted
    out["risk_quartile"] = pd.qcut(out["predicted_survival"], 4, labels=False, duplicates="drop") + 1
    rows = []
    for quartile, subset in out.groupby("risk_quartile"):
        km = KaplanMeierFitter().fit(subset[time_col], subset[event_col])
        observed = float(km.predict(horizon))
        rows.append(
            {
                "quartile": int(quartile),
                "n": int(len(subset)),
                "predicted_survival": float(subset["predicted_survival"].mean()),
                "observed_survival": observed,
                "horizon_months": horizon,
            }
        )
    return pd.DataFrame(rows)


def transport(discovery: pd.DataFrame, validation: pd.DataFrame, base_covariates: list[str]) -> dict:
    d = cohort_subset(discovery, "tcga", "os").dropna(subset=["os_time", "os_event", *base_covariates, "subtype4"])
    v = cohort_subset(validation, "metabric", "os").dropna(subset=["os_time", "os_event", *base_covariates, "subtype4"])
    base_fit = fit_cox(d, "os_time", "os_event", base_covariates, subtype_col=None)
    full_fit = fit_cox(d, "os_time", "os_event", base_covariates, subtype_col="subtype4")
    risk_base = risk_score(base_fit, v).to_numpy()
    risk_full = risk_score(full_fit, v).to_numpy()
    y_time = v["os_time"].to_numpy()
    y_event = v["os_event"].to_numpy()
    ci_base = bootstrap_cindex(y_time, y_event, risk_base, n_boot=1000, seed=SEED)
    ci_full = bootstrap_cindex(y_time, y_event, risk_full, n_boot=1000, seed=SEED)
    delta = bootstrap_delta_cindex(y_time, y_event, risk_base, risk_full, n_boot=1000, seed=SEED)
    cal = calibration_table(full_fit, v, "os_time", "os_event", horizon=60.0)
    # Calibration slope: refit a Cox model with the locked linear predictor.
    vv = v.copy()
    vv["locked_lp"] = np.log(np.clip(risk_full, 1e-9, None))
    slope_fit = CoxPHFitter().fit(vv[["os_time", "os_event", "locked_lp"]], duration_col="os_time", event_col="os_event")
    return {
        "n_discovery": int(len(d)),
        "n_validation": int(len(v)),
        "cindex_base": ci_base,
        "cindex_full": ci_full,
        "delta": delta,
        "calibration": cal,
        "calibration_slope": float(slope_fit.params_["locked_lp"]),
        "full_fit": full_fit,
        "base_fit": base_fit,
    }


def intclust_benchmark(frame: pd.DataFrame, base_covariates: list[str]) -> tuple[pd.DataFrame, dict]:
    data = frame[frame["analysis_os"] & frame["intclust"].notna()].copy()
    dummies = pd.get_dummies(data["intclust"], prefix="intclust", dtype=float)
    if "intclust_1" in dummies.columns:
        dummies = dummies.drop(columns=["intclust_1"])
    data = pd.concat([data, dummies], axis=1)
    ic_cols = dummies.columns.tolist()
    rows = []
    risks = {}
    for name, covars, subtype_col in [
        ("clinical", base_covariates, None),
        ("clinical + PAM50", base_covariates, "subtype4"),
        ("clinical + PAM50 + IntClust", base_covariates + ic_cols, "subtype4"),
    ]:
        risk = cross_validated_risk(data, "os_time", "os_event", covars, subtype_col=subtype_col, seed=SEED)
        mask = risk.notna()
        ci = bootstrap_cindex(
            data.loc[mask, "os_time"].to_numpy(),
            data.loc[mask, "os_event"].to_numpy(),
            risk[mask].to_numpy(),
            n_boot=1000,
            seed=SEED,
        )
        rows.append(
            {
                "model": name,
                "n": int(mask.sum()),
                "cindex": ci["cindex"],
                "ci_lower": ci["ci_lower"],
                "ci_upper": ci["ci_upper"],
            }
        )
        risks[name] = risk
    out = pd.DataFrame(rows)
    out["delta_vs_pam50"] = out["cindex"] - out.loc[out["model"] == "clinical + PAM50", "cindex"].iloc[0]
    return out, risks


def module_analysis(tcga: pd.DataFrame, metabric: pd.DataFrame, base_covariates: list[str], config: dict) -> tuple[pd.DataFrame, dict, pd.DataFrame, dict]:
    expr_tcga, _ = expression_with_clinical("tcga", tcga)
    expr_metabric, _ = expression_with_clinical("metabric", metabric)
    nmf_result, comparison = fit_nmf(
        expr_tcga,
        ranks=config["nmf"]["ranks"],
        max_iter=config["nmf"]["max_iter"],
        seed=SEED,
    )
    modules = top_module_genes(nmf_result, n=config["nmf"]["top_module_genes"])
    panel = pd.read_csv(PATHS.raw / "panel_genes.csv")
    entrez_to_symbol = dict(zip(panel["entrez"].astype(str), panel["hugo"]))
    module_symbols = {
        module: [entrez_to_symbol.get(g, g) for g in genes] for module, genes in modules.items()
    }
    scores_tcga = module_scores(expr_tcga, modules)
    scores_metabric = module_scores(expr_metabric, modules)

    tcga_m = tcga.set_index("sample_id").join(scores_tcga, how="inner")
    metabric_m = metabric.set_index("sample_id").join(scores_metabric, how="inner")
    rows = []
    for module in scores_tcga.columns:
        for cohort, data in [("tcga", tcga_m), ("metabric", metabric_m)]:
            data = data.dropna(subset=[module, "subtype4", "os_time", "os_event", *base_covariates])
            groups = [g[module].to_numpy() for _, g in data.groupby("subtype4") if len(g) >= 5]
            kw_p = float(kruskal(*groups).pvalue) if len(groups) > 1 else np.nan
            # Adjusted Cox per SD of module score.
            d = data.copy()
            d["module_sd"] = (d[module] - d[module].mean()) / d[module].std(ddof=0)
            base = fit_cox(d, "os_time", "os_event", base_covariates, subtype_col=None)
            full = fit_cox(d, "os_time", "os_event", base_covariates + ["module_sd"], subtype_col=None)
            lr = likelihood_ratio(full, base)
            hr = float(np.exp(full.model.params_["module_sd"]))
            rows.append(
                {
                    "module": module,
                    "cohort": cohort,
                    "n": int(len(d)),
                    "events": int(d["os_event"].sum()),
                    "subtype_kw_p": kw_p,
                    "cox_hr_per_sd": hr,
                    "cox_lr_stat": lr["lr_stat"],
                    "cox_p": lr["p_value"],
                }
            )
    out = pd.DataFrame(rows)
    out["cox_q"] = np.nan
    for (cohort,), group in out.groupby(["cohort"]):
        out.loc[group.index, "cox_q"] = bh(group["cox_p"].to_numpy(), q=config["multiplicity"]["module_q"])
    out["subtype_q"] = bh(out["subtype_kw_p"].fillna(1).to_numpy(), q=config["multiplicity"]["module_q"])

    # Module-level NPC across modules for OS in METABRIC, using a median split
    # so that the permutation null can be computed with the log-rank statistic.
    npc_payload = {}
    data = metabric_m.dropna(subset=["os_time", "os_event", "subtype4"])
    module_cols = [m for m in scores_metabric.columns if m in data.columns]
    observed = []
    high_groups = {}
    for module in module_cols:
        med = data[module].median()
        high = (data[module] >= med).astype(int).to_numpy()
        high_groups[module] = high
        observed.append(logrank_statistic(data["os_time"].to_numpy(), data["os_event"].to_numpy(), high))
    rng = np.random.default_rng(SEED)
    null = np.empty((500, len(module_cols)))
    for b in range(500):
        order = rng.permutation(len(data))
        t = data["os_time"].to_numpy()[order]
        e = data["os_event"].to_numpy()[order]
        for j, module in enumerate(module_cols):
            null[b, j] = logrank_statistic(t, e, high_groups[module])
    npc_payload = npc_test(np.asarray(observed), null, method="liptak")
    npc_payload["modules"] = module_cols
    npc_payload["null"] = null
    return (
        out,
        {
            "selected_rank": nmf_result.rank,
            "comparison": comparison,
            "modules": modules,
            "module_symbols": module_symbols,
        },
        scores_tcga,
        npc_payload,
    )


def main() -> None:
    PATHS.ensure()
    config = load_config()
    cohorts = build_cohorts(config)
    tcga, metabric = cohorts["tcga"], cohorts["metabric"]
    base_covariates = config["covariates"]["primary"]

    scorecard_rows = []
    nulls = {}
    for cohort, frame in [("tcga", tcga), ("metabric", metabric)]:
        for endpoint in ["os", "secondary"]:
            row = unadjusted_scorecard(frame, cohort, endpoint)
            nulls[(cohort, row["endpoint"])] = row.pop("null", None)
            scorecard_rows.append(row)
    scorecard = pd.DataFrame(scorecard_rows)

    incremental_rows = []
    fits = {}
    for cohort, frame in [("tcga", tcga), ("metabric", metabric)]:
        for endpoint in ["os", "secondary"]:
            row = incremental_value(frame, cohort, endpoint, base_covariates)
            fits[(cohort, row["endpoint"])] = (row.pop("base_cox"), row.pop("full_cox"))
            incremental_rows.append(row)
    incremental = pd.DataFrame(incremental_rows)

    # Primary multiplicity family: subtype LR and incremental LR in the two OS cohorts.
    primary = []
    for cohort in ["tcga", "metabric"]:
        s = scorecard[(scorecard.cohort == cohort) & (scorecard.endpoint == "OS")].iloc[0]
        i = incremental[(incremental.cohort == cohort) & (incremental.endpoint == "OS")].iloc[0]
        primary.append({"cohort": cohort, "test": "subtype log-rank", "p": s.p_value})
        primary.append({"cohort": cohort, "test": "incremental LR", "p": i.lr_p})
    primary = pd.DataFrame(primary)
    primary["q_bh"] = bh(primary["p"].to_numpy(), q=config["multiplicity"]["primary_q"])
    primary["q_by"] = by(primary["p"].to_numpy())
    scorecard["p_bh"] = np.nan
    incremental["p_bh"] = np.nan
    for _, row in primary.iterrows():
        if row["test"] == "subtype log-rank":
            idx = scorecard[(scorecard.cohort == row["cohort"]) & (scorecard.endpoint == "OS")].index
            scorecard.loc[idx, "p_bh"] = row["q_bh"]
        else:
            idx = incremental[(incremental.cohort == row["cohort"]) & (incremental.endpoint == "OS")].index
            incremental.loc[idx, "p_bh"] = row["q_bh"]

    transport_result = transport(tcga, metabric, base_covariates)
    transport_table = pd.DataFrame(
        [
            {"model": "clinical", **transport_result["cindex_base"]},
            {"model": "clinical + PAM50", **transport_result["cindex_full"]},
        ]
    )
    transport_table["delta_cindex"] = transport_result["delta"]["delta_cindex"]
    transport_table["delta_lower"] = transport_result["delta"]["ci_lower"]
    transport_table["delta_upper"] = transport_result["delta"]["ci_upper"]
    cal = transport_result["calibration"]

    intclust_table, _ = intclust_benchmark(metabric, base_covariates)
    modules, nmf_meta, scores_tcga, module_npc = module_analysis(tcga, metabric, base_covariates, config)

    # Global NPC: four components (OS/DFS in TCGA; OS/RFS in METABRIC), with
    # within-cohort joint outcome permutations preserving endpoint correlation.
    tcga_joint = joint_permutation_logrank(
        tcga, [("os_time", "os_event"), ("dfs_time", "dfs_event")], n_perm=B, seed=SEED
    )
    metabric_joint = joint_permutation_logrank(
        metabric, [("os_time", "os_event"), ("dfs_time", "dfs_event")], n_perm=B, seed=SEED
    )
    global_null = np.hstack([tcga_joint["null"], metabric_joint["null"]])
    global_observed = np.concatenate([tcga_joint["observed"], metabric_joint["observed"]])
    global_npc = npc_test(global_observed, global_null, method="liptak")
    global_npc["components"] = ["TCGA OS", "TCGA DFS", "METABRIC OS", "METABRIC RFS"]
    global_npc["null"] = global_null
    global_npc["observed"] = global_observed

    tcga_npc = npc_test(tcga_joint["observed"], tcga_joint["null"], method="liptak")
    metabric_npc = npc_test(metabric_joint["observed"], metabric_joint["null"], method="liptak")

    # KM curves for the four primary cohort/endpoint combinations.
    km_frames = []
    for cohort, frame in [("tcga", tcga), ("metabric", metabric)]:
        for endpoint in ["os", "secondary"]:
            time_col, event_col, label = endpoint_cols(cohort, endpoint)
            data = cohort_subset(frame, cohort, endpoint)
            curve = km_curve(data, time_col, event_col, "subtype4")
            if not curve.empty:
                curve["cohort"] = cohort
                curve["endpoint"] = label
                km_frames.append(curve)
    km = pd.concat(km_frames, ignore_index=True) if km_frames else pd.DataFrame()

    # Forest plot table: unadjusted HRs from the scorecard.
    forest_rows = []
    for _, row in scorecard.iterrows():
        for hr in row["hrs"]:
            forest_rows.append(
                {
                    "cohort": row["cohort"],
                    "endpoint": row["endpoint"],
                    "term": hr["term"].replace("subtype_", ""),
                    "hr": hr["hr"],
                    "lower": hr["lower"],
                    "upper": hr["upper"],
                    "p": hr["p"],
                }
            )
    forest = pd.DataFrame(forest_rows)

    # Decision logic.
    def get_score(cohort, endpoint="OS"):
        return scorecard[(scorecard.cohort == cohort) & (scorecard.endpoint == endpoint)].iloc[0]

    def get_inc(cohort, endpoint="OS"):
        return incremental[(incremental.cohort == cohort) & (incremental.endpoint == endpoint)].iloc[0]

    tcga_s, met_s = get_score("tcga"), get_score("metabric")
    tcga_i, met_i = get_inc("tcga"), get_inc("metabric")
    basal_tcga = forest[(forest.cohort == "tcga") & (forest.term == "Basal")]["hr"].iloc[0]
    basal_met = forest[(forest.cohort == "metabric") & (forest.term == "Basal")]["hr"].iloc[0]
    replicated = bool(tcga_s.p_bh < 0.05 and met_s.p_bh < 0.05 and basal_tcga > 1 and basal_met > 1)
    incremental_present = bool(
        tcga_i.delta_lower > 0 and met_i.delta_lower > 0 and tcga_i.p_bh < 0.05 and met_i.p_bh < 0.05
    )
    if replicated and incremental_present:
        decision = "Replicated and incremental"
    elif replicated:
        decision = "Replicated but not incremental"
    elif (tcga_s.p_bh < 0.05) != (met_s.p_bh < 0.05):
        decision = "Cohort-specific"
    else:
        decision = "Not replicated"

    results = {
        "meta": {
            "project": "Replica",
            "version": config["version"],
            "seed": SEED,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "n_permutations_primary": B,
            "n_permutations_npc": B,
            "n_permutations_module": 500,
        },
        "cohorts": {
            "tcga": {
                "label": "TCGA-BRCA PanCancer Atlas",
                "n_raw": int(len(tcga)),
                "n_analysis_os": int(tcga.analysis_os.sum()),
                "n_analysis_secondary": int(tcga.analysis_secondary.sum()),
                "subtypes": subtype_counts(tcga, "tcga").to_dict(orient="records"),
            },
            "metabric": {
                "label": "METABRIC",
                "n_raw": int(len(metabric)),
                "n_analysis_os": int(metabric.analysis_os.sum()),
                "n_analysis_secondary": int(metabric.analysis_secondary.sum()),
                "subtypes": subtype_counts(metabric, "metabric").to_dict(orient="records"),
            },
        },
        "scorecard": scorecard.to_dict(orient="records"),
        "primary_family": primary.to_dict(orient="records"),
        "incremental": incremental.to_dict(orient="records"),
        "transport": {
            "n_discovery": transport_result["n_discovery"],
            "n_validation": transport_result["n_validation"],
            "calibration_slope": transport_result["calibration_slope"],
            "table": transport_table.to_dict(orient="records"),
            "calibration": cal.to_dict(orient="records"),
        },
        "intclust": intclust_table.to_dict(orient="records"),
        "nmf": {
            "selected_rank": int(nmf_meta["selected_rank"]),
            "comparison": nmf_meta["comparison"].to_dict(orient="records"),
            "modules": nmf_meta["modules"],
            "module_symbols": nmf_meta["module_symbols"],
            "module_npc": {k: v for k, v in module_npc.items() if k not in {"null"}},
        },
        "modules": modules.to_dict(orient="records"),
        "npc": {
            "global": {k: v for k, v in global_npc.items() if k not in {"null"}},
            "tcga": {k: v for k, v in tcga_npc.items() if k not in {"null"}},
            "metabric": {k: v for k, v in metabric_npc.items() if k not in {"null"}},
        },
        "decision": {
            "label": decision,
            "replicated": replicated,
            "incremental": incremental_present,
            "basal_hr_tcga": float(basal_tcga),
            "basal_hr_metabric": float(basal_met),
        },
    }

    # Write artifacts.
    (PATHS.artifacts / "results.json").write_text(json.dumps(jsonable(results), indent=2))
    scorecard.drop(columns=[c for c in scorecard.columns if c == "hrs"]).to_csv(PATHS.tables / "scorecard.csv", index=False)
    incremental.to_csv(PATHS.tables / "incremental.csv", index=False)
    primary.to_csv(PATHS.tables / "primary_family.csv", index=False)
    transport_table.to_csv(PATHS.tables / "transport.csv", index=False)
    cal.to_csv(PATHS.tables / "calibration.csv", index=False)
    intclust_table.to_csv(PATHS.tables / "intclust.csv", index=False)
    modules.to_csv(PATHS.tables / "modules.csv", index=False)
    km.to_csv(PATHS.tables / "km_curves.csv", index=False)
    forest.to_csv(PATHS.tables / "forest.csv", index=False)
    pd.concat([subtype_counts(tcga, "tcga"), subtype_counts(metabric, "metabric")]).to_csv(
        PATHS.tables / "subtype_counts.csv", index=False
    )
    pd.DataFrame(
        [
            {"cohort": "tcga", "model": "clinical", **transport_result["cindex_base"]},
            {"cohort": "tcga", "model": "clinical+PAM50", **transport_result["cindex_full"]},
        ]
    ).to_csv(PATHS.tables / "transport_cindex.csv", index=False)
    pd.DataFrame({"component": global_npc["components"], "statistic": global_observed}).to_csv(
        PATHS.tables / "npc_components.csv", index=False
    )
    pd.DataFrame({"null_combined": global_npc["null_combined"]}).to_csv(
        PATHS.tables / "npc_null.csv", index=False
    )
    pd.DataFrame({"module": module_npc["modules"], "statistic": module_npc["observed_components"],
                  "p": module_npc["component_pvalues"]}).to_csv(PATHS.tables / "module_npc_components.csv", index=False)

    # Locked models and coefficients.
    joblib.dump(transport_result["full_fit"], PATHS.models / "discovery_full_cox.joblib")
    joblib.dump(transport_result["base_fit"], PATHS.models / "discovery_base_cox.joblib")
    joblib.dump(nmf_meta, PATHS.models / "nmf_meta.joblib")
    coef_rows = []
    for label, fit in [("discovery_clinical", transport_result["base_fit"]),
                       ("discovery_clinical_pam50", transport_result["full_fit"])]:
        for term, row in fit.summary().iterrows():
            coef_rows.append({"model": label, "term": term, "coef": row["coef"],
                              "hr": row["hr"], "lower": row["hr_lower"], "upper": row["hr_upper"],
                              "p": row["p"]})
    pd.DataFrame(coef_rows).to_csv(PATHS.tables / "model_coefficients.csv", index=False)

    # Attrition and data dictionary.
    pd.concat([pd.read_csv(PATHS.tables / f"attrition_{c}.csv").assign(cohort=c) for c in ["tcga", "metabric"]]).to_csv(
        PATHS.tables / "attrition.csv", index=False
    )
    print(json.dumps(jsonable(results["decision"]), indent=2))
    print(scorecard[["cohort", "endpoint", "n", "events", "p_value", "cindex", "p_bh"]].to_string(index=False))
    print(incremental[["cohort", "endpoint", "cindex_base", "cindex_full", "delta_cindex", "delta_lower", "delta_upper", "lr_p"]].to_string(index=False))
    print(transport_table.to_string(index=False))
    print(intclust_table.to_string(index=False))
    print(modules.to_string(index=False))
    print("global NPC p =", global_npc["combined_pvalue"], "components", global_npc["component_pvalues"])


if __name__ == "__main__":
    main()
