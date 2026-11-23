"""Download a compact, pre-specified Hallmark marker panel from cBioPortal.

The full Hallmark union is ~4,400 genes.  Downloading every gene for every sample
would be several gigabytes of JSON, so the panel is fixed *before* outcome
analysis: for each Hallmark pathway take the 20 genes that occur in the fewest
Hallmark sets (ties broken alphabetically).  This keeps pathway-specific signal
and yields a 998-gene panel that is practical to reproduce.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
API = "https://www.cbioportal.org/api"
COHORTS = {
    "tcga_pan": {
        "study": "brca_tcga_pan_can_atlas_2018",
        "profile": "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "samples": "tcga_pan_sample_map.csv",
    },
    "metabric": {
        "study": "brca_metabric",
        "profile": "brca_metabric_mrna",
        "samples": "metabric_sample_map.csv",
    },
}
BATCH = 50


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def post_with_retry(url: str, payload: dict, attempts: int = 4) -> list[dict]:
    for attempt in range(attempts):
        try:
            response = requests.post(url, json=payload, timeout=300)
            if response.status_code == 200:
                return response.json()
            print(f"  HTTP {response.status_code}; retry {attempt + 1}/{attempts}")
        except requests.RequestException as exc:
            print(f"  {type(exc).__name__}; retry {attempt + 1}/{attempts}")
        time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url}")


def main() -> None:
    panel = pd.read_csv(RAW / "panel_gene_selection.csv")
    genes = json.loads((RAW / "genes.json").read_text())
    sym2id = {
        g["hugoGeneSymbol"]: g["entrezGeneId"]
        for g in genes
        if g.get("type") == "protein-coding"
    }
    panel["entrez"] = panel["hugo"].map(sym2id)
    panel = panel.dropna(subset=["entrez"]).astype({"entrez": int})
    panel.to_csv(RAW / "panel_genes.csv", index=False)
    entrez = [int(x) for x in sorted(panel["entrez"].unique())]
    print(f"panel: {panel['hugo'].nunique()} genes, {len(entrez)} Entrez IDs")

    manifest = []
    for short, meta in COHORTS.items():
        sample_map = pd.read_csv(RAW / meta["samples"])
        sample_ids = sorted(sample_map["sampleId"].unique())
        url = f"{API}/molecular-profiles/{meta['profile']}/molecular-data/fetch?projection=SUMMARY"
        frames = []
        started = time.time()
        for start in range(0, len(entrez), BATCH):
            batch = entrez[start : start + BATCH]
            records = post_with_retry(url, {"sampleIds": sample_ids, "entrezGeneIds": batch})
            if records:
                frames.append(pd.DataFrame(records)[["sampleId", "entrezGeneId", "value"]])
            print(
                f"{short}: batch {start // BATCH + 1}/{(len(entrez) - 1) // BATCH + 1} "
                f"({len(records)} records, {time.time() - started:.0f}s)"
            )
        long = pd.concat(frames, ignore_index=True)
        wide = long.pivot_table(
            index="sampleId", columns="entrezGeneId", values="value", aggfunc="first"
        )
        wide.columns = [str(int(c)) for c in wide.columns]
        wide = wide.sort_index()
        out = RAW / f"{short}_expression.csv"
        wide.to_csv(out)
        manifest.append(
            {
                "cohort": short,
                "study": meta["study"],
                "profile": meta["profile"],
                "samples_requested": len(sample_ids),
                "samples_returned": len(wide),
                "genes_requested": len(entrez),
                "genes_returned": wide.shape[1],
                "records": len(long),
                "seconds": round(time.time() - started, 1),
                "sha256": sha256(out),
            }
        )
        print(f"{short}: {wide.shape} -> {out}")
    pd.DataFrame(manifest).to_csv(RAW / "expression_manifest.csv", index=False)


if __name__ == "__main__":
    main()
