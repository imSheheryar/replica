# MLOps and release gates

Replica is a batch statistical product, not a continuously retrained model. The MLOps goal is to make every number traceable and every release reproducible.

## Stages

| Stage | Command | Output | Gate |
|---|---|---|---|
| Ingest | `make clinical`, `make expression` | raw JSON/CSV + manifests + SHA-256 | source counts match expected study sizes |
| Harmonise | `make cohorts` | `data/processed/*_cohort.csv` | four-class subtype and valid OS present |
| Analyse | `make analysis` | `artifacts/results.json`, tables, locked models | schema, p-value range, C-index range, primary family size |
| Present | `make site` | `site/index.html` + assets | required figures and sections present |
| Validate | `make validate` | exit code | all gates pass |
| Serve | `make api` / Docker | FastAPI + static dashboard | `/health` returns 200 |
| Test | `make test` | pytest | unit tests pass |
| Release | Git tag + CI | GitHub release / Pages | CI green, artifact hash recorded |

## Data gates

The ingestion scripts write `clinical_manifest.csv` and `expression_manifest.csv`. Each row records the study, profile, requested/returned counts, records, elapsed time, and SHA-256. A release artifact is invalid if the raw files do not match the manifest.

The harmonisation step writes `attrition.csv` with the count at each exclusion. A cohort with fewer than 100 patients in the primary OS analysis fails the gate; this is a sanity check, not a scientific threshold.

## Model gates

The discovery model is fitted once and serialised to `artifacts/models/discovery_full_cox.joblib`. The model card records the feature set. No validation outcome is used for feature selection. The NMF gene lists are serialised in `results.json`; the METABRIC scores are a projection, not a re-fit.

## Statistical gates

`scripts/validate_artifacts.py` checks:

- `results.json` has every required top-level key;
- the scorecard is non-empty and every p-value is in [0, 1];
- every C-index is in [0.5, 1] (a C-index below 0.5 indicates a sign error or a pathological model);
- the primary family has exactly four tests and every BH q is in [0, 1];
- the global NPC p-value is in [0, 1] and uses at least 100 permutations;
- the dashboard contains the forest and KM figures.

These gates are not a substitute for scientific review. They catch broken releases, not bad science.

## Serving gates

The FastAPI service is read-only. It does not accept patient data, does not fit models, and does not expose row-level records. `/health` is used by Docker and the hosting platform. The static dashboard is served from `site/`; GitHub Pages publishes the same files.

## Reproducibility

Every run uses seed 20261123. The config version is written into `results.json`. The package versions are pinned in `pyproject.toml` lower bounds; the Docker image pins the Python version. For a stricter release, build the Docker image and archive its digest alongside the artifact hash.

## Rollback

If a release fails a gate after publication, revert the `site/` and `artifacts/` commit, redeploy Pages, and open an audit-log entry. Do not edit a published p-value in place. Re-run the pipeline from the frozen config and record the new artifact hash.
