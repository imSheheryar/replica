"""Release gates for the frozen Replica artifacts.

These are deliberately simple, hard failures.  A dashboard should never be
published from a half-written or statistically impossible artifact.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts"
TABLES = ARTIFACTS / "tables"

REQUIRED = {
    "meta",
    "cohorts",
    "scorecard",
    "primary_family",
    "incremental",
    "transport",
    "intclust",
    "nmf",
    "modules",
    "npc",
    "decision",
}


def fail(message: str) -> None:
    print(f"FAIL: {message}")
    sys.exit(1)


def main() -> None:
    path = ARTIFACTS / "results.json"
    if not path.exists():
        fail("artifacts/results.json is missing; run the pipeline first")
    results = json.loads(path.read_text())
    missing = REQUIRED - set(results)
    if missing:
        fail(f"results.json missing keys: {sorted(missing)}")

    scorecard = pd.DataFrame(results["scorecard"])
    if scorecard.empty:
        fail("scorecard is empty")
    for column in ["cohort", "endpoint", "n", "events", "p_value", "cindex"]:
        if column not in scorecard:
            fail(f"scorecard missing column {column}")
    if not scorecard["p_value"].between(0, 1).all():
        fail("scorecard contains a p-value outside [0, 1]")
    if not scorecard["cindex"].between(0.5, 1.0).all():
        fail("scorecard contains a C-index below 0.5 or above 1")

    primary = pd.DataFrame(results["primary_family"])
    if len(primary) != 4:
        fail(f"primary family should have 4 tests, got {len(primary)}")
    if not primary["q_bh"].between(0, 1).all():
        fail("primary family has a BH q outside [0, 1]")

    npc = results["npc"]["global"]
    if not 0 <= npc["combined_pvalue"] <= 1:
        fail("global NPC p-value outside [0, 1]")
    if npc["n_permutations"] < 100:
        fail("too few NPC permutations for a release artifact")

    for name in ["scorecard.csv", "incremental.csv", "forest.csv", "km_curves.csv", "modules.csv"]:
        if not (TABLES / name).exists():
            fail(f"missing table {name}")

    site = ROOT / "site" / "index.html"
    if not site.exists():
        fail("site/index.html is missing; run scripts/build_site.py")
    html = site.read_text()
    if "fig-forest" not in html or "fig-km" not in html:
        fail("dashboard is missing required figures")

    print("OK: artifacts pass schema and sanity gates")
    print(f"  scorecard rows: {len(scorecard)}")
    print(f"  primary family rows: {len(primary)}")
    print(f"  global NPC p: {npc['combined_pvalue']:.4f}")
    print(f"  decision: {results['decision']['label']}")


if __name__ == "__main__":
    main()
