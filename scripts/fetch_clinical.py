"""Download and flatten cBioPortal clinical data for the Replica cohorts.

The cBioPortal clinical-data endpoint returns long-format records.  This script
keeps the raw JSON response and writes a wide CSV plus a sample map so that every
downstream artifact can be traced back to an immutable download.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
API = "https://www.cbioportal.org/api"
STUDIES = {
    "tcga_pan": "brca_tcga_pan_can_atlas_2018",
    "tcga_legacy": "brca_tcga",
    "metabric": "brca_metabric",
}


def fetch(study: str, kind: str) -> list[dict]:
    url = f"{API}/studies/{study}/clinical-data"
    params = {"clinicalDataType": kind, "projection": "SUMMARY"}
    response = requests.get(url, params=params, timeout=120)
    response.raise_for_status()
    payload = response.json()
    (RAW / f"{study}_{kind.lower()}_clinical.json").write_text(json.dumps(payload))
    return payload


def to_wide(records: list[dict]) -> pd.DataFrame:
    frame = pd.DataFrame(records)
    if frame.empty:
        return frame
    frame = frame[["patientId", "clinicalAttributeId", "value"]].copy()
    frame["value"] = frame["value"].astype(str)
    wide = frame.pivot_table(
        index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first"
    ).reset_index()
    wide.columns.name = None
    return wide


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []
    for short, study in STUDIES.items():
        for kind in ("PATIENT", "SAMPLE"):
            records = fetch(study, kind)
            manifest.append(
                {
                    "study": study,
                    "short": short,
                    "kind": kind,
                    "records": len(records),
                    "unique_patients": len({r["patientId"] for r in records}),
                    "unique_samples": len({r.get("sampleId") for r in records if r.get("sampleId")}),
                }
            )
            if kind == "PATIENT":
                wide = to_wide(records)
                wide.to_csv(RAW / f"{short}_clinical.csv", index=False)
            else:
                frame = pd.DataFrame(records)
                if not frame.empty:
                    frame.to_csv(RAW / f"{short}_sample_clinical_long.csv", index=False)
                    sample_map = (
                        frame[["patientId", "sampleId"]]
                        .dropna()
                        .drop_duplicates()
                        .sort_values(["patientId", "sampleId"])
                    )
                    sample_map.to_csv(RAW / f"{short}_sample_map.csv", index=False)
    pd.DataFrame(manifest).to_csv(RAW / "clinical_manifest.csv", index=False)
    print(pd.DataFrame(manifest).to_string(index=False))


if __name__ == "__main__":
    main()
