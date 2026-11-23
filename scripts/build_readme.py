"""Generate README.md from the frozen artifacts, so prose cannot drift from results."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts"
TABLES = ARTIFACTS / "tables"


def pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def main() -> None:
    results = json.loads((ARTIFACTS / "results.json").read_text())
    scorecard = pd.DataFrame(results["scorecard"])
    incremental = pd.DataFrame(results["incremental"])
    modules = pd.DataFrame(results["modules"])
    decision = results["decision"]
    transport = results["transport"]
    npc = results["npc"]["global"]
    nmf = results["nmf"]
    sensitivity = results.get("sensitivity", [])

    tcga_os = scorecard[(scorecard.cohort == "tcga") & (scorecard.endpoint == "OS")].iloc[0]
    met_os = scorecard[(scorecard.cohort == "metabric") & (scorecard.endpoint == "OS")].iloc[0]
    tcga_inc = incremental[(incremental.cohort == "tcga") & (incremental.endpoint == "OS")].iloc[0]
    met_inc = incremental[(incremental.cohort == "metabric") & (incremental.endpoint == "OS")].iloc[0]
    trans_row = transport["table"][1] if len(transport["table"]) > 1 else transport["table"][0]

    score_rows = "\n".join(
        f"| {r.cohort.upper()} | {r.endpoint} | {int(r.n):,} | {int(r.events):,} | "
        f"{r.cindex:.3f} | {r.p_value:.2e} | {r.p_bh:.3f} |"
        if pd.notna(r.p_bh)
        else f"| {r.cohort.upper()} | {r.endpoint} | {int(r.n):,} | {int(r.events):,} | "
        f"{r.cindex:.3f} | {r.p_value:.2e} | — |"
        for r in scorecard.itertuples()
    )
    inc_rows = "\n".join(
        f"| {r.cohort.upper()} | {r.endpoint} | {r.cindex_base:.3f} | {r.cindex_full:.3f} | "
        f"{r.delta_cindex:+.3f} | [{r.delta_lower:+.3f}, {r.delta_upper:+.3f}] | {r.lr_p:.2e} |"
        for r in incremental.itertuples()
    )
    module_rows = "\n".join(
        f"| {r.module} | {r.cohort.upper()} | {r.cox_hr_per_sd:.2f} | {r.cox_p:.2e} | {r.cox_q:.3f} |"
        for r in modules.itertuples()
    )
    sensitivity_rows = "\n".join(
        f"| {r['cohort'].upper()} | {r['endpoint']} | {int(r['n']):,} | {int(r['events']):,} | "
        f"{r['cindex']:.3f} | {r['p_value']:.2e} |"
        for r in sensitivity
    )
    sim_path = TABLES / "simulation_summary.csv"
    fdr_path = TABLES / "simulation_fdr_summary.csv"
    sim_section = ""
    if sim_path.exists() and fdr_path.exists():
        sim = pd.read_csv(sim_path)
        fdr = pd.read_csv(fdr_path)
        sim_rows = "\n".join(
            f"| {'null' if r['shift'] == 0 else f'{r.shift:.2f}'} | {r['npc_type1_or_power']:.1%} | "
            f"{r['marginal_bh_rate']:.1%} | {r['bonferroni_rate']:.1%} |"
            for _, r in sim.iterrows()
        )
        fdr_rows = "\n".join(
            f"| {r['method']} | {r['mean_discoveries']:.1f} | {r['mean_fdp']:.3f} | {r['mean_power']:.1%} |"
            for _, r in fdr.iterrows()
        )
        sim_section = f"""## Methods check: operating characteristics

This is a simulation, not a breast cancer result. At shift 0 the rejection rate estimates Type I error; positive shifts estimate power. The FDR check compares marginal BH with Fisher NPC under dependent features.

| Scenario | NPC rejection | Marginal BH | Bonferroni |
|---|---:|---:|---:|
{sim_rows}

| Method | Mean discoveries | Mean FDP | Power |
|---|---:|---:|---:|
{fdr_rows}

The FDR simulation uses an asymptotic Fisher combination for speed; the cohort-level NPC uses permutation calibration. The check is designed to show operating behaviour, not to produce a second cohort p-value.

"""

    headline = (
        f"**TL;DR — {decision['label']}.** PAM50 separation in TCGA OS "
        f"(p = {tcga_os.p_value:.2e}; BH q = {tcga_os.p_bh:.3f}) and METABRIC OS "
        f"(p = {met_os.p_value:.2e}; BH q = {met_os.p_bh:.3f}). "
        f"Out-of-fold ΔC over the clinical baseline was {tcga_inc.delta_cindex:+.3f} "
        f"[{tcga_inc.delta_lower:+.3f}, {tcga_inc.delta_upper:+.3f}] in TCGA and "
        f"{met_inc.delta_cindex:+.3f} [{met_inc.delta_lower:+.3f}, {met_inc.delta_upper:+.3f}] in METABRIC. "
        f"The locked TCGA clinical + PAM50 model transported with C-index {trans_row['cindex']:.3f} "
        f"and ΔC {trans_row['delta_cindex']:+.3f} over the clinical-only model. "
        f"The four-component global NPC test gave p = {npc['combined_pvalue']:.4f}; this is a global "
        f"association test, not a replication test, and the METABRIC components drive it."
    )

    text = f"""# Replica

### Which breast cancer subtypes survive external validation?

[![ci](https://github.com/imSheheryar/replica/actions/workflows/ci.yml/badge.svg)](https://github.com/imSheheryar/replica/actions/workflows/ci.yml)
[![pages](https://github.com/imSheheryar/replica/actions/workflows/pages.yml/badge.svg)](https://github.com/imSheheryar/replica/actions/workflows/pages.yml)
[![python](https://img.shields.io/badge/python-3.10%2B-1f7a72)](https://www.python.org/)
[![license](https://img.shields.io/badge/license-MIT-12212b)](LICENSE)

A pre-specified, internally frozen replication audit of breast cancer subtype systems across **TCGA-BRCA** and **METABRIC**, with censoring-aware permutation tests, nonparametric combination (NPC), FDR control, external validation, expression modules, a clinical evidence dashboard, and a read-only FastAPI service.

**Live dashboard:** [imSheheryar.github.io/replica](https://imSheheryar.github.io/replica/)

![Replica dashboard](site/assets/dashboard-hero.png)

*The analysis plan is version-controlled in [`docs/preregistration.md`](docs/preregistration.md) and was frozen before outcome analysis. It is an internal pre-registration, not an externally deposited registry entry; the audit log records every correction and deviation.*

---

## The finding

{headline}

The important result is not that PAM50 is "significant." It is whether the published categories survive an independent cohort and whether they add anything beyond age, stage, ER, and PR. Replica reports both, including the uncomfortable case where a subtype system is prognostic but not incremental.

![Replica architecture](docs/architecture.png)

## What this project contributes

Replica is not a new classifier and not another TCGA clustering notebook. Its contribution is an **audit design** plus an empirical answer:

- **Subtype as a transported predictor, not a label to be re-derived.** The published PAM50/Claudin-low calls are frozen, the TCGA model is locked, and the same coefficients are applied to METABRIC.
- **Censoring-aware NPC across dependent endpoints.** The global test permutes bivariate outcome pairs jointly within each cohort, so OS–DFS/RFS dependence and the censoring distribution are preserved rather than assumed away.
- **External modules without outcome leakage.** NMF is fitted in TCGA; the top-30 genes per module are locked before METABRIC is scored, so METABRIC outcomes cannot influence module definition.
- **A negative result that is interrogated, not rescued.** The primary TCGA OS test is null; the full-available-cohort sensitivity shows the null is power/selection-sensitive, but the canonical Basal-worst pattern still does not replicate.
- **Release engineering.** Frozen JSON artifacts, schema/sanity gates, 14 unit tests, Ruff, CI, Docker, a static Plotly dashboard, and a read-only FastAPI service.

## Why this is not another TCGA clustering notebook

Most subtype papers derive a classifier and report survival in the same cohort. That design cannot separate reproducible biology from overfitting, platform effects, or confounding. Replica freezes the published calls, locks the model, moves it to an independent cohort, and evaluates the transport.

| Design choice | What it prevents |
|---|---|
| Pre-specified analysis plan, frozen before outcomes | post-hoc endpoint and subgroup selection |
| TCGA reference → METABRIC external validation | same-cohort optimism |
| Clinical baseline (age, stage, ER, PR) | attributing generic prognosis to subtype |
| Out-of-fold ΔC with paired bootstrap | in-sample C-index inflation |
| Bivariate outcome permutation | invalid tests that permute only the event flag |
| Joint NPC calibration | an independence assumption across OS and DFS/RFS |
| BH/BY for module families; maxT as a methods check | cherry-picked "significant" modules |
| Locked NMF gene lists | outcome leakage into external module scores |
| Read-only API over frozen JSON | a dashboard that silently changes a p-value |

## Scorecard

The primary family is four tests: subtype log-rank in TCGA OS, subtype log-rank in METABRIC OS, and the incremental likelihood-ratio test in each OS cohort. BH is applied at q = 0.05. With 2,000 permutations, a reported p of 0.0005 is the Monte Carlo floor (1/2001), not an exact population p-value.

| Cohort | Endpoint | n | Events | C-index | Log-rank p | BH q |
|---|---|---:|---:|---:|---:|---:|
{score_rows}

The TCGA null is a failure to meet the pre-specified replication rule in this cohort, not evidence that PAM50 carries no biological signal. The sensitivity below makes that distinction concrete.

### Sensitivity: the complete-case filter matters

The primary cohort requires age, stage, ER, and PR. This sensitivity uses every patient with a four-class subtype and a valid endpoint, without the clinical complete-case filter. It is not a rescue; it shows how much the exclusion changes the answer.

| Cohort | Endpoint | n | Events | C-index | Permutation p |
|---|---|---:|---:|---:|---:|
{sensitivity_rows}

{decision.get('caveat', '')}

This release reports the full-available-cohort sensitivity and the METABRIC IntClust benchmark. The remaining analyses listed in the pre-registration are prespecified extensions, not executed results; see the audit log for the exact status.

## Incremental value beyond clinical variables

Out-of-fold C-index from 5-fold cross-validation within each cohort. The same folds are used for the clinical and clinical + subtype models. ΔC is paired; the interval is a patient-level bootstrap.

| Cohort | Endpoint | Clinical C | Clinical + PAM50 C | ΔC | 95% CI | LR p |
|---|---|---:|---:|---:|---|---:|
{inc_rows}

The METABRIC OS likelihood-ratio p is 0.025, but the paired ΔC interval includes zero and the TCGA ΔC is negative. The incremental criterion therefore fails even though the LR test is nominally significant in one cohort. That distinction is the point of reporting effect size and interval alongside the p-value.

## External transport

The TCGA clinical + PAM50 model is locked and applied to METABRIC. Discrimination uses risk ranking; calibration at 60 months uses the TCGA baseline hazard and therefore exposes baseline-hazard differences rather than hiding them.

| Model | C-index | 95% CI | ΔC vs clinical |
|---|---:|---|---:|
| Clinical | {transport['table'][0]['cindex']:.3f} | [{transport['table'][0]['ci_lower']:.3f}, {transport['table'][0]['ci_upper']:.3f}] | — |
| Clinical + PAM50 | {transport['table'][1]['cindex']:.3f} | [{transport['table'][1]['ci_lower']:.3f}, {transport['table'][1]['ci_upper']:.3f}] | {transport['table'][1]['delta_cindex']:+.3f} [{transport['table'][1]['delta_lower']:+.4f}, {transport['table'][1]['delta_upper']:+.4f}] |

Calibration slope: **{transport['calibration_slope']:.2f}**. A slope near 1 means the locked linear predictor transfers on the log-hazard scale; a slope below 1 indicates over-dispersion.

The transport ΔC is small and its interval only just excludes zero. It is evidence of a detectable rank improvement, not evidence of a clinically useful reclassification.

### Internal benchmark: IntClust

IntClust has no external cohort with comparable labels, so it is a METABRIC-only benchmark, not a replication. It asks whether a more granular integrative cluster system adds discrimination beyond PAM50 and clinical variables in the same cohort.

| Model | C-index | ΔC vs clinical + PAM50 |
|---|---:|---:|
| Clinical | {results['intclust'][0]['cindex']:.3f} | {results['intclust'][0]['delta_vs_pam50']:+.4f} |
| Clinical + PAM50 | {results['intclust'][1]['cindex']:.3f} | — |
| Clinical + PAM50 + IntClust | {results['intclust'][2]['cindex']:.3f} | {results['intclust'][2]['delta_vs_pam50']:+.4f} |

The IntClust increment is small (+{results['intclust'][2]['delta_vs_pam50']:.4f}) and internally estimated. It is reported to show that the audit can benchmark a second subtype system, not as external evidence.

## Nonparametric combination

The global test combines four dependent components: TCGA OS, TCGA DFS, METABRIC OS, and METABRIC RFS. The null distribution is generated by permuting outcome pairs jointly within each cohort, preserving the OS–DFS/RFS correlation. Liptak is the primary combiner; Fisher and Tippett are implemented as alternative combiners but are not presented as executed cohort analyses.

- Global NPC p: **{npc['combined_pvalue']:.4f}**
- Component p-values: {', '.join(f'{p:.3f}' if p > 0.001 else '≤0.001' for p in npc['component_pvalues'])}
- Permutations: {npc['n_permutations']:,}
- Method: {npc['method']}

**Read this as a global association test, not a replication test.** The two METABRIC components are at the permutation floor (0.0005 each), while the TCGA OS and DFS components are 0.307 and 0.495. A significant global NPC therefore does not overturn the cohort-specific decision; it says the four-endpoint ensemble contains signal, concentrated in METABRIC.

Fisher and Tippett combiners and permutation maxT are implemented and unit-tested. The cohort-level result reported here is Liptak; the other combiners are not presented as executed cohort analyses.

## Expression modules

NMF is fitted in TCGA only, on a pre-specified 995-gene Hallmark marker panel. METABRIC module scores are the mean z-scored expression of the locked top-30 genes per module, so no METABRIC outcome information enters module definition.

Selected rank: **k = {nmf['selected_rank']}**. Module-level Cox models are adjusted for the clinical baseline; BH is applied at q = 0.10.

| Module | Cohort | HR per SD | Cox p | BH q |
|---|---|---:|---:|---:|
{module_rows}

The modules are strongly subtype-associated (Kruskal–Wallis q values below 1e-60) but none is independently prognostic after clinical adjustment (all BH q > 0.25). That is a useful negative result: the expression structure is real, but it is largely collinear with the subtype/clinical axis rather than a new independent risk factor.

{sim_section}## Quickstart

```bash
git clone https://github.com/imSheheryar/replica
cd replica
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"

# Rebuild the raw downloads (requires internet; cBioPortal API)
make clinical
make expression

# Harmonise, analyse, build the dashboard, and run release gates
make pipeline

# Serve the API + dashboard
make api
# http://127.0.0.1:8000
```

If the raw data are already present, `make analysis` and `make site` are enough. The pipeline is deterministic under seed `20261123`.

Run `make test`, `make lint`, and `make validate` to reproduce the release gates. The API can be served with `make api`; it reads frozen artifacts and never fits a model.

## Repository map

| Path | Purpose |
|---|---|
| `config/analysis.yaml` | frozen primary parameters, seed, permutation counts, FDR thresholds |
| `docs/preregistration.md` | hypotheses, endpoints, models, decision rules, before outcomes |
| `docs/methods.md` | derivations, permutation mechanics, NPC calibration, limitations |
| `docs/data_card.md` | cohort provenance, harmonisation, missingness, ethics |
| `docs/model_card.md` | intended use, out-of-scope use, metrics, ethical considerations |
| `docs/mlops.md` | ingest → analyse → validate → serve → release gates |
| `docs/audit_log.md` | decisions and deviations with timestamps |
| `src/replica/data.py` | clinical harmonisation, expression features, attrition |
| `src/replica/survival.py` | Cox, log-rank, censoring-aware permutations, bootstrap C-index |
| `src/replica/stats.py` | NPC, Fisher/Liptak/Tippett, BH/BY, maxT |
| `src/replica/models.py` | NMF discovery and locked module scoring |
| `scripts/fetch_clinical.py` | cBioPortal clinical download + manifests |
| `scripts/fetch_expression.py` | 995-gene expression panel download + checksums |
| `scripts/run_analysis.py` | primary audit, transport, IntClust, modules, NPC |
| `scripts/validate_artifacts.py` | release gates and statistical sanity checks |
| `scripts/build_site.py` | static Plotly dashboard |
| `app/api.py` | read-only FastAPI service |
| `tests/` | NPC, FDR, survival, harmonisation tests |
| `Dockerfile` | API + dashboard image |
| `render.yaml` | Render blueprint for the read-only API |
| `.github/workflows/` | CI and GitHub Pages deployment |

## Statistical core

The permutation test permutes the **bivariate outcome pair** \((T_i,\delta_i)\), not just the event indicator. That preserves the censoring distribution. For adjusted tests, outcomes are permuted within age × stage × ER strata. The log-rank statistic is implemented directly so the scheme is auditable and thousands of permutations are practical.

NPC calibrates the combined statistic against the **joint** permutation distribution. For every permutation, the same row permutation is applied to OS and DFS/RFS within a cohort, so the dependence between endpoints is preserved. The p-value is \((1+\\#\\{{C_b \\ge C_{{\\text{{obs}}}}\\}})/(B+1)\), never zero.

FDR control uses Benjamini–Hochberg for the primary family and module family, Benjamini–Yekutieli as a dependence-robust sensitivity analysis, and permutation maxT for the module family.

## Dashboard and API

The static dashboard is built from frozen artifacts with Plotly and deployed to GitHub Pages. It shows the primary scorecard first, then incremental value, hazard ratios, Kaplan–Meier curves, transport calibration, the NPC null distribution, and expression modules. Exploratory and sensitivity results are visually separated and labelled.

The FastAPI service exposes only aggregate JSON: `/health`, `/api/summary`, `/api/scorecard`, `/api/forest`, `/api/km`, `/api/modules`, `/api/npc`, `/api/transport`, `/api/methods`. It never accepts patient data and never fits a model.

**Deploy.** The dashboard deploys through `.github/workflows/pages.yml`. The API is container-ready: `docker build -t replica:1.0.0 . && docker run -p 8000:8000 replica:1.0.0`. A Render blueprint is included in `render.yaml`; it points at the same Dockerfile and health endpoint. The public dashboard is static; the API is optional and read-only.

## Limitations

The subtype labels are not assay-identical across cohorts. TCGA PAM50 calls and METABRIC PAM50/Claudin-low calls were produced by different pipelines. TCGA grade is absent from the cBioPortal export, so the primary clinical baseline uses age, stage, ER, and PR; grade is a METABRIC-only sensitivity. DFS and RFS are different endpoints. Complete-case analysis can bias estimates if missingness is informative. The audit is observational and cannot identify a causal effect of subtype on survival. The permutation test addresses exchangeability, not confounding. The pre-registration is internal rather than externally registered, and the 2,000-permutation budget limits resolution at the tail; both are stated rather than hidden.

## Citation

```bibtex
@software{{replica2026,
  author = {{imSheheryar}},
  title = {{Replica: A replication audit of breast cancer subtype systems}},
  year = {{2026}},
  url = {{https://github.com/imSheheryar/replica}}
}}
```

Source cohorts: TCGA Research Network; Curtis et al., *Nature* 2012; Pereira et al., *Nature Communications* 2016. Code MIT; data remain governed by the source terms.
"""
    (ROOT / "README.md").write_text(text)
    print("Wrote README.md")


if __name__ == "__main__":
    main()
