"""Censoring-aware survival utilities.

The permutation test permutes the *bivariate outcome* (time, event), which keeps
the censoring distribution intact.  For adjusted tests the permutation is
stratified by clinical covariates.  The log-rank statistic is implemented
directly so that thousands of permutations remain practical.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from sklearn.model_selection import StratifiedKFold
from scipy.stats import chi2



@dataclass
class CoxFit:
    model: CoxPHFitter
    design: pd.DataFrame
    data: pd.DataFrame
    duration_col: str
    event_col: str
    covariates: list[str]
    base_covariates: list[str]
    subtype_col: str | None
    reference: str

    @property
    def log_likelihood(self) -> float:
        return float(self.model.log_likelihood_)

    def summary(self) -> pd.DataFrame:
        table = self.model.summary.copy()
        table["hr"] = np.exp(table["coef"])
        table["hr_lower"] = np.exp(table["coef lower 95%"])
        table["hr_upper"] = np.exp(table["coef upper 95%"])
        return table


def design_matrix(
    frame: pd.DataFrame,
    covariates: Iterable[str],
    subtype_col: str | None = "subtype4",
    reference: str = "LumA",
) -> pd.DataFrame:
    cols = list(covariates)
    design = frame[cols].astype(float).copy()
    if subtype_col and subtype_col in frame.columns:
        dummies = pd.get_dummies(frame[subtype_col], prefix="subtype", dtype=float)
        keep = [c for c in dummies.columns if c != f"subtype_{reference}"]
        design = pd.concat([design, dummies[keep]], axis=1)
    return design


def fit_cox(
    frame: pd.DataFrame,
    duration_col: str,
    event_col: str,
    covariates: Iterable[str],
    subtype_col: str | None = "subtype4",
    reference: str = "LumA",
) -> CoxFit:
    design = design_matrix(frame, covariates, subtype_col=subtype_col, reference=reference)
    data = pd.concat(
        [frame[[duration_col, event_col]].astype(float), design], axis=1
    ).dropna()
    data = data[data[duration_col] > 0]
    if data[event_col].sum() < 5:
        raise ValueError("Too few events to fit a Cox model")
    model = CoxPHFitter()
    model.fit(data, duration_col=duration_col, event_col=event_col, show_progress=False)
    return CoxFit(
        model,
        design,
        data,
        duration_col,
        event_col,
        list(design.columns),
        list(covariates),
        subtype_col,
        reference,
    )


def likelihood_ratio(full: CoxFit, reduced: CoxFit) -> dict[str, float]:
    """LR test for the extra terms in ``full`` relative to ``reduced``."""
    # Refit both on the same rows to make the likelihoods comparable.
    common = full.data.index.intersection(reduced.data.index)
    if len(common) < 20:
        raise ValueError("Not enough common rows for a likelihood-ratio test")
    full_data = full.data.loc[common]
    reduced_data = reduced.data.loc[common]
    f = CoxPHFitter().fit(full_data, duration_col=full.duration_col, event_col=full.event_col, show_progress=False)
    r = CoxPHFitter().fit(reduced_data, duration_col=reduced.duration_col, event_col=reduced.event_col, show_progress=False)
    lr = 2.0 * (f.log_likelihood_ - r.log_likelihood_)
    df = max(1, full_data.shape[1] - reduced_data.shape[1])
    return {"lr_stat": float(lr), "df": int(df), "p_value": float(chi2.sf(lr, df))}


def risk_score(fit: CoxFit, frame: pd.DataFrame) -> pd.Series:
    design = design_matrix(
        frame, fit.base_covariates, subtype_col=fit.subtype_col, reference=fit.reference
    )
    design = design.reindex(columns=fit.covariates)
    design = design.astype(float).fillna(0.0)
    return fit.model.predict_partial_hazard(design)


def cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    """Harrell's C-index, vectorized over comparable later-time pairs.

    Pairs with equal observed times are not treated as comparable.  This is the
    conventional Harrell definition and matches lifelines to numerical
    precision on the validation tests.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    risk = np.asarray(risk, dtype=float)
    t = time[:, None]
    r = risk[:, None]
    e = event[:, None]
    comparable = (t < t.T) & (e == 1)
    concordant = float(np.sum(comparable & (r > r.T)))
    discordant = float(np.sum(comparable & (r < r.T)))
    tied = float(np.sum(comparable & (r == r.T)))
    denominator = concordant + discordant + tied
    if denominator == 0:
        return float("nan")
    return float((concordant + 0.5 * tied) / denominator)


def cross_validated_risk(
    frame: pd.DataFrame,
    duration_col: str,
    event_col: str,
    covariates: list[str],
    subtype_col: str | None = "subtype4",
    folds: int = 5,
    seed: int = 20261123,
) -> pd.Series:
    """Out-of-fold partial hazards, using the same folds for nested models."""
    data = frame.dropna(subset=[duration_col, event_col, *covariates, "subtype4"]).copy()
    data = data[data[duration_col] > 0]
    design = design_matrix(data, covariates, subtype_col=subtype_col)
    y = data[event_col].astype(int).to_numpy()
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    risks = pd.Series(np.nan, index=data.index, dtype=float)
    for train_idx, test_idx in splitter.split(data, y):
        train = data.iloc[train_idx]
        # Fit on the design columns; a fold can lack a rare subtype, so drop
        # all-zero columns and restore them as zero for prediction.
        train_design = design.iloc[train_idx].copy()
        keep = [c for c in train_design.columns if train_design[c].nunique() > 1]
        if not keep:
            continue
        fit_data = pd.concat([train[[duration_col, event_col]].astype(float), train_design[keep]], axis=1)
        model = CoxPHFitter()
        try:
            model.fit(fit_data, duration_col=duration_col, event_col=event_col, show_progress=False)
        except Exception:
            continue
        test_design = design.iloc[test_idx][keep].astype(float)
        risks.iloc[test_idx] = model.predict_partial_hazard(test_design).to_numpy()
    return risks


def bootstrap_cindex(
    time: np.ndarray,
    event: np.ndarray,
    risk: np.ndarray,
    n_boot: int = 1000,
    seed: int = 20261123,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    risk = np.asarray(risk, dtype=float)
    values = []
    n = len(time)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if event[idx].sum() < 5:
            continue
        values.append(cindex(time[idx], event[idx], risk[idx]))
    values = np.asarray(values)
    return {
        "cindex": cindex(time, event, risk),
        "ci_lower": float(np.quantile(values, 0.025)),
        "ci_upper": float(np.quantile(values, 0.975)),
        "n_boot": int(values.size),
    }


def bootstrap_delta_cindex(
    time: np.ndarray,
    event: np.ndarray,
    risk_base: np.ndarray,
    risk_full: np.ndarray,
    n_boot: int = 1000,
    seed: int = 20261123,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    risk_base = np.asarray(risk_base, dtype=float)
    risk_full = np.asarray(risk_full, dtype=float)
    n = len(time)
    deltas = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if event[idx].sum() < 5:
            continue
        deltas.append(
            cindex(time[idx], event[idx], risk_full[idx])
            - cindex(time[idx], event[idx], risk_base[idx])
        )
    deltas = np.asarray(deltas)
    observed = cindex(time, event, risk_full) - cindex(time, event, risk_base)
    return {
        "delta_cindex": float(observed),
        "ci_lower": float(np.quantile(deltas, 0.025)),
        "ci_upper": float(np.quantile(deltas, 0.975)),
        "n_boot": int(deltas.size),
    }


def km_curve(frame: pd.DataFrame, time_col: str, event_col: str, group_col: str = "subtype4") -> pd.DataFrame:
    rows = []
    for group, subset in frame.groupby(group_col):
        if len(subset) < 5:
            continue
        fitter = KaplanMeierFitter()
        fitter.fit(subset[time_col], subset[event_col], label=str(group))
        timeline = fitter.survival_function_.index.to_numpy()
        sf = fitter.survival_function_.iloc[:, 0].to_numpy()
        ci = fitter.confidence_interval_
        rows.append(
            pd.DataFrame(
                {
                    "group": str(group),
                    "time": timeline,
                    "survival": sf,
                    "lower": ci.iloc[:, 0].to_numpy(),
                    "upper": ci.iloc[:, 1].to_numpy(),
                }
            )
        )
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def logrank_statistic(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> float:
    """Multi-group log-rank chi-square statistic (G - 1 df)."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    group = np.asarray(group)
    _, group = np.unique(group, return_inverse=True)
    n_groups = int(group.max()) + 1
    event_times = np.unique(time[event == 1])
    observed = np.zeros(n_groups)
    expected = np.zeros(n_groups)
    covariance = np.zeros((n_groups, n_groups))
    for t in event_times:
        at_risk = time >= t
        n = int(at_risk.sum())
        if n == 0:
            continue
        events_at_t = (time == t) & (event == 1)
        d = int(events_at_t.sum())
        if d == 0:
            continue
        n_risk_g = np.bincount(group[at_risk], minlength=n_groups).astype(float)
        d_g = np.bincount(group[events_at_t], minlength=n_groups).astype(float)
        observed += d_g
        expected += d * n_risk_g / n
        if n > 1:
            p = n_risk_g / n
            covariance += (d * (n - d) / (n - 1)) * (np.diag(p) - np.outer(p, p))
    z = (observed - expected)[: n_groups - 1]
    v = covariance[: n_groups - 1, : n_groups - 1]
    if not np.isfinite(v).all() or np.linalg.matrix_rank(v) == 0:
        return 0.0
    return float(z @ np.linalg.pinv(v) @ z)


def permutation_logrank(
    time: np.ndarray,
    event: np.ndarray,
    group: np.ndarray,
    n_perm: int = 10000,
    seed: int = 20261123,
    strata: np.ndarray | None = None,
) -> dict[str, object]:
    """Permute the bivariate outcome, optionally within clinical strata.

    Permuting (time, event) pairs keeps the censoring distribution intact.  The
    test is conditional on the observed group vector.
    """
    rng = np.random.default_rng(seed)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    group = np.asarray(group)
    observed = logrank_statistic(time, event, group)
    null = np.empty(n_perm, dtype=float)
    if strata is None:
        for b in range(n_perm):
            order = rng.permutation(len(time))
            null[b] = logrank_statistic(time[order], event[order], group)
    else:
        strata = np.asarray(strata)
        for b in range(n_perm):
            perm_time = time.copy()
            perm_event = event.copy()
            for stratum in np.unique(strata):
                idx = np.flatnonzero(strata == stratum)
                if idx.size < 2:
                    continue
                order = rng.permutation(idx)
                perm_time[idx] = time[order]
                perm_event[idx] = event[order]
            null[b] = logrank_statistic(perm_time, perm_event, group)
    p_value = float((1.0 + np.sum(null >= observed)) / (n_perm + 1.0))
    return {
        "statistic": observed,
        "p_value": p_value,
        "null": null,
        "n_permutations": n_perm,
        "df": int(np.unique(group).size - 1),
    }


def joint_permutation_logrank(
    frame: pd.DataFrame,
    endpoints: list[tuple[str, str]],
    group_col: str = "subtype4",
    n_perm: int = 2000,
    seed: int = 20261123,
) -> dict[str, object]:
    """Permute complete outcome rows jointly across several endpoints.

    The same row permutation is applied to every endpoint, so the within-person
    correlation between OS and DFS/RFS is preserved in the null distribution.
    """
    cols = [c for endpoint in endpoints for c in endpoint]
    data = frame.dropna(subset=[*cols, group_col]).copy()
    for time_col, _ in endpoints:
        data = data[data[time_col] > 0]
    group = data[group_col].to_numpy()
    observed = np.array(
        [logrank_statistic(data[t].to_numpy(), data[e].to_numpy(), group) for t, e in endpoints]
    )
    rng = np.random.default_rng(seed)
    null = np.empty((n_perm, len(endpoints)), dtype=float)
    for b in range(n_perm):
        order = rng.permutation(len(data))
        for j, (t, e) in enumerate(endpoints):
            null[b, j] = logrank_statistic(
                data[t].to_numpy()[order], data[e].to_numpy()[order], group
            )
    return {
        "observed": observed,
        "null": null,
        "n": int(len(data)),
        "endpoints": endpoints,
        "n_permutations": n_perm,
    }


def stratified_permutation_lr(
    frame: pd.DataFrame,
    duration_col: str,
    event_col: str,
    base_covariates: list[str],
    n_perm: int = 500,
    seed: int = 20261123,
) -> dict[str, float]:
    """Adjusted test: permute outcomes within age x stage x ER strata and refit."""
    rng = np.random.default_rng(seed)
    data = frame.dropna(subset=[duration_col, event_col, *base_covariates, "subtype4"]).copy()
    data = data[data[duration_col] > 0]
    base = fit_cox(data, duration_col, event_col, base_covariates, subtype_col=None)
    full = fit_cox(data, duration_col, event_col, base_covariates, subtype_col="subtype4")
    observed = likelihood_ratio(full, base)["lr_stat"]
    strata = (
        pd.qcut(data["age10"], 3, labels=False, duplicates="drop").astype(str)
        + "|" + data["stage_high"].astype(int).astype(str)
        + "|" + data["er_pos"].astype(int).astype(str)
    ).to_numpy()
    null = []
    for _ in range(n_perm):
        perm = data.copy()
        for stratum in np.unique(strata):
            idx = np.flatnonzero(strata == stratum)
            if idx.size < 10:
                continue
            order = rng.permutation(idx)
            perm.iloc[idx, perm.columns.get_loc(duration_col)] = data.iloc[order][duration_col].to_numpy()
            perm.iloc[idx, perm.columns.get_loc(event_col)] = data.iloc[order][event_col].to_numpy()
        try:
            b = fit_cox(perm, duration_col, event_col, base_covariates, subtype_col=None)
            f = fit_cox(perm, duration_col, event_col, base_covariates, subtype_col="subtype4")
            null.append(likelihood_ratio(f, b)["lr_stat"])
        except Exception:
            continue
    null = np.asarray(null)
    return {
        "statistic": float(observed),
        "p_value": float((1.0 + np.sum(null >= observed)) / (null.size + 1.0)),
        "n_permutations": int(null.size),
    }
