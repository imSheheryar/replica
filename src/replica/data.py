"""Cohort construction and harmonisation.

All decisions here are intentionally explicit and counted.  The output is a
single analysis-ready table per cohort plus a machine-readable attrition table.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .config import PATHS, load_config


def _read(name: str) -> pd.DataFrame:
    path = PATHS.raw / name
    if not path.exists():
        raise FileNotFoundError(f"Missing raw file: {path}. Run scripts/fetch_clinical.py first.")
    return pd.read_csv(path, dtype=str)


def parse_event(series: pd.Series) -> pd.Series:
    """Parse cBioPortal event strings such as '1:DECEASED' and '0:LIVING'."""
    text = series.astype("string").str.strip().str.upper()
    out = pd.Series(np.nan, index=series.index, dtype="float")
    out[text.str.startswith("1:", na=False)] = 1.0
    out[text.str.startswith("0:", na=False)] = 0.0
    return out


def parse_time(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def parse_stage_tcga(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.upper().str.replace(r"\s+", " ", regex=True)
    text = text.str.replace("STAGE ", "", regex=False).str.strip()
    mapping = {
        "I": 1, "IA": 1, "IB": 1,
        "II": 2, "IIA": 2, "IIB": 2,
        "III": 3, "IIIA": 3, "IIIB": 3, "IIIC": 3,
        "IV": 4,
    }
    return text.map(mapping).astype("float")


def parse_stage_metabric(series: pd.Series) -> pd.Series:
    value = pd.to_numeric(series, errors="coerce")
    # Stage 0 is in-situ disease; this is an invasive cohort, so treat as unknown.
    return value.where(value.between(1, 4))


def parse_receptor(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip().str.upper()
    out = pd.Series(np.nan, index=series.index, dtype="float")
    out[text.isin(["POSITIVE", "POS", "YES"])] = 1.0
    out[text.isin(["NEGATIVE", "NEG", "NO"])] = 0.0
    return out


def parse_grade(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.extract(r"([123])", expand=False)
    return pd.to_numeric(text, errors="coerce")


def _first_sample_map(short: str) -> pd.DataFrame:
    frame = _read(f"{short}_sample_map.csv")
    frame = frame.dropna().drop_duplicates()
    frame = frame.sort_values(["patientId", "sampleId"]).groupby("patientId", as_index=False).first()
    return frame


def _tcga_clinical() -> pd.DataFrame:
    pan = _read("tcga_pan_clinical.csv")
    legacy = _read("tcga_legacy_clinical.csv")
    keep = ["patientId", "ER_STATUS_BY_IHC", "PR_STATUS_BY_IHC", "HER2_IHC_SCORE", "IHC_HER2"]
    legacy = legacy[[c for c in keep if c in legacy.columns]].drop_duplicates("patientId")
    clinical = pan.merge(legacy, on="patientId", how="left", validate="one_to_one")
    clinical = clinical.merge(_first_sample_map("tcga_pan"), on="patientId", how="left")
    return clinical


def _metabric_clinical() -> pd.DataFrame:
    patient = _read("metabric_clinical.csv")
    sample_long = _read("metabric_sample_clinical_long.csv")
    sample_cols = [
        "ER_STATUS", "PR_STATUS", "HER2_STATUS", "GRADE", "TUMOR_STAGE", "TUMOR_SIZE",
    ]
    sample_long = sample_long[sample_long["clinicalAttributeId"].isin(sample_cols)]
    sample = (
        sample_long.pivot_table(
            index=["patientId", "sampleId"],
            columns="clinicalAttributeId",
            values="value",
            aggfunc="first",
        )
        .reset_index()
        .rename(columns={"ER_STATUS": "ER_STATUS_SAMPLE", "PR_STATUS": "PR_STATUS_SAMPLE",
                         "HER2_STATUS": "HER2_STATUS_SAMPLE", "GRADE": "GRADE_SAMPLE",
                         "TUMOR_STAGE": "TUMOR_STAGE_SAMPLE", "TUMOR_SIZE": "TUMOR_SIZE_SAMPLE"})
    )
    clinical = patient.merge(sample, on="patientId", how="left", validate="one_to_one")
    return clinical


def harmonize_tcga() -> pd.DataFrame:
    clinical = _tcga_clinical()
    out = pd.DataFrame({
        "cohort": "tcga",
        "cohort_label": "TCGA-BRCA PanCancer Atlas",
        "patient_id": clinical["patientId"],
        "sample_id": clinical["sampleId"],
        "subtype_raw": clinical["SUBTYPE"],
        "age": pd.to_numeric(clinical["AGE"], errors="coerce"),
        "stage_raw": clinical["AJCC_PATHOLOGIC_TUMOR_STAGE"],
        "stage": parse_stage_tcga(clinical["AJCC_PATHOLOGIC_TUMOR_STAGE"]),
        "er": parse_receptor(clinical["ER_STATUS_BY_IHC"]),
        "pr": parse_receptor(clinical["PR_STATUS_BY_IHC"]),
        "grade": np.nan,
        "her2": np.nan,
        "os_time": parse_time(clinical["OS_MONTHS"]),
        "os_event": parse_event(clinical["OS_STATUS"]),
        "dfs_time": parse_time(clinical["DFS_MONTHS"]),
        "dfs_event": parse_event(clinical["DFS_STATUS"]),
        "secondary_endpoint": "DFS",
        "chemotherapy": np.nan,
        "hormone_therapy": np.nan,
        "radio_therapy": np.nan,
    })
    return out


def harmonize_metabric() -> pd.DataFrame:
    clinical = _metabric_clinical()
    out = pd.DataFrame({
        "cohort": "metabric",
        "cohort_label": "METABRIC",
        "patient_id": clinical["patientId"],
        "sample_id": clinical["sampleId"],
        "subtype_raw": clinical["CLAUDIN_SUBTYPE"],
        "intclust": clinical["INTCLUST"],
        "threegene": clinical["THREEGENE"],
        "age": pd.to_numeric(clinical["AGE_AT_DIAGNOSIS"], errors="coerce"),
        "stage_raw": clinical["TUMOR_STAGE_SAMPLE"],
        "stage": parse_stage_metabric(clinical["TUMOR_STAGE_SAMPLE"]),
        "er": parse_receptor(clinical["ER_STATUS_SAMPLE"]),
        "pr": parse_receptor(clinical["PR_STATUS_SAMPLE"]),
        "grade": parse_grade(clinical["GRADE_SAMPLE"]),
        "her2": parse_receptor(clinical["HER2_STATUS_SAMPLE"]),
        "os_time": parse_time(clinical["OS_MONTHS"]),
        "os_event": parse_event(clinical["OS_STATUS"]),
        "dfs_time": parse_time(clinical["RFS_MONTHS"]),
        "dfs_event": parse_event(clinical["RFS_STATUS"]),
        "secondary_endpoint": "RFS",
        "chemotherapy": clinical["CHEMOTHERAPY"],
        "hormone_therapy": clinical["HORMONE_THERAPY"],
        "radio_therapy": clinical["RADIO_THERAPY"],
    })
    return out


def apply_subtype_mapping(frame: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    mapping = config["subtypes"]["mapping"]
    frame = frame.copy()
    frame["subtype"] = frame["subtype_raw"].map(mapping)
    primary = set(config["subtypes"]["primary"])
    frame["subtype4"] = frame["subtype"].where(frame["subtype"].isin(primary))
    frame["age10"] = frame["age"] / 10.0
    frame["stage_high"] = frame["stage"].map({1: 0.0, 2: 0.0, 3: 1.0, 4: 1.0})
    frame["er_pos"] = frame["er"]
    frame["pr_pos"] = frame["pr"]
    frame["grade3"] = frame["grade"].map({1: 0.0, 2: 0.0, 3: 1.0})
    frame["her2_pos"] = frame["her2"]
    return frame


def attrition(frame: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    primary = frame["subtype4"].notna()
    os_ok = frame["os_time"].gt(0) & frame["os_event"].notna()
    covars = ["age10", "stage_high", "er_pos", "pr_pos"]
    complete = frame[covars].notna().all(axis=1)
    rows = [
        ("raw patients", len(frame)),
        ("primary four-class subtype", int(primary.sum())),
        ("+ valid OS", int((primary & os_ok).sum())),
        ("+ complete clinical covariates", int((primary & os_ok & complete).sum())),
        ("secondary endpoint available", int((primary & frame["dfs_time"].gt(0) & frame["dfs_event"].notna()).sum())),
    ]
    return pd.DataFrame(rows, columns=["step", "n"])


def build_cohorts(config: dict[str, Any] | None = None, save: bool = True) -> dict[str, pd.DataFrame]:
    config = config or load_config()
    PATHS.ensure()
    frames = {
        "tcga": apply_subtype_mapping(harmonize_tcga(), config),
        "metabric": apply_subtype_mapping(harmonize_metabric(), config),
    }
    for key, frame in frames.items():
        frame["analysis_os"] = (
            frame["subtype4"].notna()
            & frame["os_time"].gt(0)
            & frame["os_event"].notna()
            & frame[["age10", "stage_high", "er_pos", "pr_pos"]].notna().all(axis=1)
        )
        frame["analysis_secondary"] = (
            frame["subtype4"].notna()
            & frame["dfs_time"].gt(0)
            & frame["dfs_event"].notna()
            & frame[["age10", "stage_high", "er_pos", "pr_pos"]].notna().all(axis=1)
        )
        if save:
            frame.to_csv(PATHS.processed / f"{key}_cohort.csv", index=False)
            attrition(frame, config).to_csv(PATHS.tables / f"attrition_{key}.csv", index=False)
    if save:
        pd.concat([attrition(f, config).assign(cohort=k) for k, f in frames.items()]).to_csv(
            PATHS.tables / "attrition.csv", index=False
        )
    return frames


def load_expression(cohort: str, genes: list[str] | None = None) -> pd.DataFrame:
    filename = {"tcga": "tcga_pan_expression.csv", "metabric": "metabric_expression.csv"}.get(
        cohort, f"{cohort}_expression.csv"
    )
    path = PATHS.raw / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing expression file: {path}. Run scripts/fetch_expression.py first.")
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if genes is not None:
        keep = [g for g in genes if g in frame.columns]
        frame = frame[keep]
    # TCGA RSEM values are linear; METABRIC microarray values are already log2.
    if frame.to_numpy().max() > 50:
        frame = np.log2(frame.clip(lower=0) + 1)
    return frame


def zscore_genes(frame: pd.DataFrame) -> pd.DataFrame:
    mean = frame.mean(axis=0)
    sd = frame.std(axis=0, ddof=0).replace(0, np.nan)
    return (frame - mean) / sd


def pathway_scores(cohort: str) -> pd.DataFrame:
    panel = pd.read_csv(PATHS.raw / "panel_genes.csv")
    panel["entrez"] = panel["entrez"].astype(int).astype(str)
    expr = zscore_genes(load_expression(cohort, panel["entrez"].unique().tolist()))
    scores: dict[str, pd.Series] = {}
    for pathway, group in panel.groupby("pathway"):
        genes = [g for g in group["entrez"].astype(str) if g in expr.columns]
        if genes:
            scores[pathway.replace("HALLMARK_", "")] = expr[genes].mean(axis=1)
    return pd.DataFrame(scores).sort_index()


def expression_with_clinical(cohort: str, frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (gene matrix, pathway scores) aligned to the cohort's sample IDs."""
    panel = pd.read_csv(PATHS.raw / "panel_genes.csv")
    panel["entrez"] = panel["entrez"].astype(int).astype(str)
    expr = zscore_genes(load_expression(cohort, panel["entrez"].unique().tolist()))
    expr = expr.reindex(frame["sample_id"].astype(str)).dropna(how="all")
    scores = pathway_scores(cohort).reindex(expr.index)
    return expr, scores
