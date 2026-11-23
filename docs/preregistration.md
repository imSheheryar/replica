# Pre-registration — Replica

**Version:** 1.0, frozen 2026-11-23, before accessing any survival outcome by subtype in the analysis code.
**Status:** internal, version-controlled analysis plan. It is not deposited in an external registry (for example, OSF or ClinicalTrials.gov). The commit history and audit log provide the timestamp; readers should treat this as a pre-specified audit plan, not as an externally registered trial.
**Project:** Replica — a replication audit of breast cancer subtype systems.
**Primary question:** Which published breast cancer subtype systems replicate across independent cohorts, and do they add prognostic value beyond routine clinical variables?

This document is deliberately written before the primary analysis. It fixes the cohorts, endpoint definitions, model specifications, permutation schemes, multiplicity families, and decision rules. Exploratory analyses are labelled as such and cannot be promoted to primary evidence after seeing results.

## 1. Why this is a replication audit, not another clustering paper

Breast cancer intrinsic subtypes are widely used as if they were portable biological categories. Most publications derive a classifier in one cohort and report survival associations in that same cohort. That design cannot distinguish a reproducible biological axis from cohort-specific overfitting, platform effects, or confounding by age, stage, grade, and receptor status.

Replica asks three questions that a single-cohort analysis cannot answer:

1. **Transport:** does a subtype system fitted in a discovery cohort retain discrimination in an independent cohort?
2. **Increment:** does it add information beyond cheap, routinely available clinical variables?
3. **Robustness:** does the answer survive a censoring-aware permutation test, multiplicity control, and reasonable changes in endpoint and cohort construction?

The contribution is the audit framework and the empirical answer, not a new cluster label.

## 2. Cohorts and eligibility

| Role | Cohort | Source | Inclusion |
|---|---|---|---|
| Reference / discovery | TCGA-BRCA PanCancer Atlas | cBioPortal `brca_tcga_pan_can_atlas_2018` | Female or male patients with a primary breast tumour, non-missing PAM50 `SUBTYPE`, OS time > 0, and OS status |
| External validation | METABRIC | cBioPortal `brca_metabric` | Patients with non-missing `CLAUDIN_SUBTYPE`, OS time > 0, and OS status |
| Clinical supplement | TCGA-BRCA Firehose Legacy | cBioPortal `brca_tcga` | ER/PR by IHC joined to PanCancer patients by `patientId`; used only for covariates, never for outcomes |
| Internal benchmark | METABRIC | cBioPortal `brca_metabric` | Patients with non-missing `INTCLUST` and OS/RFS |

Only one primary tumour sample per patient is retained. If multiple samples exist, the sample with the lexicographically smallest sample ID is used, without reference to outcome. The TCGA PanCancer cohort is the reference/discovery cohort: it has published PAM50 calls and supplies the locked model and NMF modules. METABRIC is external because its PAM50/Claudin-low calls were generated independently.

The TCGA Firehose cohort overlaps the PanCancer cohort at the patient level (1,084 of 1,084 PanCancer patients are present in Firehose). ER/PR are therefore a within-patient clinical annotation, not an independent cohort.

## 3. Subtype harmonisation

The primary four-class comparison is **LumA, LumB, Basal, HER2**. TCGA labels `BRCA_LumA`, `BRCA_LumB`, `BRCA_Basal`, `BRCA_Her2` map to these classes. METABRIC `CLAUDIN_SUBTYPE` values `LumA`, `LumB`, `Basal`, `Her2` map directly.

The following labels are analysed separately and never silently merged into the four-class primary analysis:

- `Normal` (TCGA `BRCA_Normal`; METABRIC `Normal`): a low-tumour-content / normal-like class whose interpretation is contested.
- `claudin-low` (METABRIC only): a published extension of PAM50, not one of the original four intrinsic classes.
- `NC` (METABRIC only): not classified; excluded from subtype models.

Sensitivity analyses repeat the primary analysis with (a) Normal included as a fifth class and (b) claudin-low included as a fifth class in METABRIC. These are labelled sensitivity, not primary.

## 4. Endpoints

| Endpoint | TCGA definition | METABRIC definition | Role |
|---|---|---|---|
| Overall survival (OS) | `OS_MONTHS`, `OS_STATUS` (`1:DECEASED`) | `OS_MONTHS`, `OS_STATUS` (`1:DECEASED`) | Primary, both cohorts |
| Disease-free / relapse-free survival | TCGA `DFS_MONTHS`, `DFS_STATUS` (`1:Recurred/Progressed`) | METABRIC `RFS_MONTHS`, `RFS_STATUS` (`1:Recurred`) | Secondary, cohort-specific; combined only through NPC with a prespecified caveat |

OS is the primary endpoint because its definition is comparable. DFS and RFS are not the same endpoint: TCGA DFS includes recurrence, progression, and death; METABRIC RFS is relapse-free survival. They are therefore analysed separately and combined only as a prespecified sensitivity analysis. A difference between the OS and RFS/DFS conclusions is informative, not a failure to be hidden.

Times are converted to months. Zero or negative times are excluded. Patients with missing time or status are excluded from that endpoint, not imputed.

## 5. Covariates

The clinical baseline is deliberately restricted to variables that a clinic could obtain without a gene-expression assay:

- age at diagnosis (continuous, per 10 years);
- pathologic stage (I/II vs III/IV; missing as a separate category in sensitivity only);
- ER status (positive vs negative/indeterminate);
- PR status (positive vs negative/indeterminate).

Histologic grade is recorded in METABRIC but is absent from the TCGA PanCancer and TCGA Firehose clinical exports available through cBioPortal. Including it only in METABRIC would make the discovery and validation models non-comparable, so grade is excluded from the primary baseline and used only in a METABRIC-only sensitivity analysis. This decision was made before any survival model was fitted and is recorded in the audit log.

HER2 is not available in the TCGA PanCancer clinical file and is not available at comparable quality in TCGA Firehose. It is therefore **not** part of the primary clinical baseline. METABRIC HER2 is used in a prespecified sensitivity analysis, and the primary METABRIC model is also reported without HER2 so that the discovery and validation models use the same covariate set.

Missing covariates are handled by complete-case analysis in the primary model. A sensitivity analysis uses median/mode imputation plus missingness indicators. The number excluded at each step is reported.

## 6. Primary estimands and models

### 6.1 Replication of subtype prognostic separation

For each cohort and each endpoint, fit a Cox proportional hazards model with subtype as the only covariate. Report:

- likelihood-ratio test for the subtype term (3 df for four classes);
- hazard ratios for LumB, Basal, and HER2 relative to LumA, with 95% confidence intervals;
- Kaplan-Meier curves by subtype;
- concordance index (C-index) with 95% bootstrap CI.

A subtype system is considered to **replicate** if the subtype likelihood-ratio test is significant at the prespecified alpha of 0.05 in both cohorts **and** the direction of the Basal-vs-LumA hazard ratio is the same in both cohorts with a meta-analytic 95% CI excluding 1. The second condition prevents a statistically significant but biologically incoherent reversal from being called replication.

### 6.2 Incremental value beyond clinical variables

For each cohort and endpoint, compare:

- **Clinical model:** age, stage, grade, ER, PR.
- **Clinical + subtype model:** the same covariates plus subtype.

The primary incremental estimand is the change in C-index (ΔC) from the clinical model to the clinical + subtype model, with a paired bootstrap 95% CI. A likelihood-ratio test is reported as a secondary criterion. Increment is considered **present** if the lower bound of the 95% CI for ΔC is greater than 0 and the LR test p-value is < 0.05. Because ΔC can be small even when the LR test is significant, both are reported and neither is hidden.

### 6.3 External transport

The clinical + subtype Cox model is fitted on TCGA (discovery) and applied to METABRIC (validation) using the fitted coefficients. Because baseline hazards are cohort-specific, the primary transport metric is the C-index for risk ranking. Calibration at 60 months is reported as a secondary metric using the discovery baseline hazard; systematic miscalibration is expected and is interpreted as a cohort/baseline-hazard difference, not as evidence against discrimination.

A **cohort-recalibrated** transport analysis refits only the baseline hazard in METABRIC and reports calibration again. This separates discrimination from baseline transport.

### 6.4 IntClust benchmark

Within METABRIC, compare clinical, clinical + PAM50, and clinical + PAM50 + IntClust models using 5-fold cross-validated C-index. IntClust is not available in TCGA, so this is an internal benchmark, not an external replication. The question is whether a 10-class integrative cluster adds anything beyond PAM50 and clinical variables in the cohort where it was derived.

### 6.5 Expression modules

A 998-gene Hallmark marker panel is fixed before outcome analysis (20 most pathway-specific genes per Hallmark pathway; ties alphabetical). Within TCGA, genes are log2-transformed where needed, z-scored, and the top 500 variable genes are used for non-negative matrix factorisation (NMF). NMF ranks k = 3,...,8 are compared by silhouette on the sample loading matrix and reconstruction error; the selected k is the smallest rank within 0.02 of the best silhouette, favouring interpretability. Module scores in METABRIC are the mean z-scored expression of the top 30 genes by NMF loading, so no METABRIC outcome information enters module definition.

For each module, test:

- association with subtype (Kruskal-Wallis);
- association with OS (Cox per SD, likelihood-ratio test).

The module-level null is controlled with Benjamini-Hochberg FDR at q = 0.10 across modules and endpoints. A secondary NPC analysis combines the module-level tests, with the joint null distribution obtained by permutation.

## 7. Nonparametric combination (NPC)

### 7.1 Why NPC

The audit has multiple partially dependent endpoints (OS, RFS/DFS) and multiple subtype contrasts. A Bonferroni correction is conservative under positive dependence; a parametric multivariate model imposes distributional assumptions that are questionable for censored, skewed, high-dimensional data. NPC combines marginal test statistics without assuming a joint parametric distribution and calibrates the combination by permutation, so the dependence among components is preserved empirically.

### 7.2 Combining functions

For marginal p-values \(p_1,\ldots,p_K\), the implemented combining functions are:

- **Fisher:** \(C_F = -2\sum_{k=1}^{K}\log p_k\), null \(\chi^2_{2K}\) under independence (used only as a reference; the permutation calibration does not rely on it).
- **Liptak:** \(C_L = \sum_{k=1}^{K}\Phi^{-1}(1-p_k)\).
- **Tippett:** \(C_T = \max_{k}(1-p_k)\).
- **Mahalanobis-style directional NPC (sensitivity):** combine signed z-scores with a permutation-estimated covariance matrix.

The primary combination is Liptak on marginal permutation p-values, because it is sensitive to consistent moderate effects across components without being dominated by a single endpoint. Fisher and Tippett are implemented as alternative combiners; the released cohort-level result is Liptak (see execution status below).

### 7.3 Permutation scheme

The primary unadjusted test permutes the **bivariate outcome pair** \((T_i,\delta_i)\) across subjects. Permuting the pair, rather than the event indicator alone, preserves the observed censoring distribution and the dependence between time and event status. The test is conditional on the observed covariate/ subtype vector and is valid under the sharp null of no association between subtype and the survival outcome, with the usual assumption that censoring is independent of subtype. This is a censoring-aware randomisation test; it does not assume exponential or proportional hazards.

For adjusted tests, outcomes are permuted **within strata** of the clinical covariates (age tertile × stage I/II vs III/IV × ER positive vs negative). This preserves the clinical covariate-outcome relationship under the null and tests whether subtype adds information beyond those strata. Strata with fewer than 10 subjects are pooled into an adjacent stratum before permutation. The exact stratification is recorded in the artifact metadata.

The frozen execution budget is 2,000 permutations for primary and exploratory cohort tests. The simulation study uses 80 permutations per replicate because it averages over many replicates; its operating-characteristic estimates are not used for cohort inference. A fixed seed (20261123) is used. The permutation p-value is \((1 + \#\{T_b \ge T_{obs}\})/(B+1)\), never 0. With 2,000 permutations the Monte Carlo floor is \(1/2001 = 0.0005\); a reported p of 0.0005 therefore means "at or below the floor," not an exact population p-value. An earlier draft specified 10,000 permutations; the executed pipeline and the released configuration use 2,000. This is recorded as a documentation correction in the audit log and changes Monte Carlo resolution, not the decision rule.

For NPC, each component statistic is computed on every permutation. Marginal permutation p-values are obtained from the same permutation ensemble, and the combined statistic is calibrated against the permutation distribution of the combined statistic. This preserves the dependence among components and avoids an independence assumption.

### 7.4 What the permutation test does not solve

Permutation tests do not fix confounding, measurement error, selection bias, or non-random missingness. They test exchangeability under a specified null. The audit therefore reports covariate balance, missingness, and sensitivity analyses rather than treating a small permutation p-value as causal evidence.

## 8. Multiplicity control

The primary family is the four tests in the replication scorecard: the subtype log-rank test in TCGA-OS, the subtype log-rank test in METABRIC-OS, the incremental likelihood-ratio test in TCGA-OS, and the incremental likelihood-ratio test in METABRIC-OS. ΔC is the paired effect size for the incremental tests; its bootstrap interval is reported alongside the LR p-value, not as a separate p-value in the family. Benjamini-Hochberg FDR at q = 0.05 is applied within this family. The scorecard reports raw and BH-adjusted p-values.

Secondary families:

- endpoint-specific models (OS, DFS/RFS): BH within endpoint;
- module tests: BH at q = 0.10 across modules and endpoints;
- sensitivity analyses: reported raw, labelled sensitivity, and not used to rescue a failed primary result.

The Benjamini-Yekutieli procedure is reported as a sensitivity analysis because the tests are positively dependent. Permutation maxT adjusted p-values are also reported for the module family.

## 9. Decision rules

The audit has four possible conclusions, and all four are publishable:

1. **Replicated and incremental:** subtype log-rank significant in both cohorts, Basal HR direction concordant, the incremental LR q < 0.05 in both cohorts, and the paired ΔC CI excludes 0 in both cohorts.
2. **Replicated but not incremental:** subtype separation replicates, but ΔC CI includes 0 in at least one cohort. This is the most likely outcome and is the central empirical contribution: the subtype is prognostic, but not additive to routine clinical variables.
3. **Cohort-specific:** significant in one cohort but not the other, or a direction reversal. Report transport failure and investigate platform, case-mix, and endpoint differences.
4. **Not replicated:** no significant separation in the validation cohort. Report the null result and the power analysis; do not reframe it as a discovery.

A negative or mixed result is not a failure of the project. The project's value is the audit, the framework, and the honest evidence.

## 10. Sensitivity analyses (prespecified)

- Include Normal as a fifth class.
- Include claudin-low as a fifth class in METABRIC.
- Impute missing clinical covariates with median/mode plus missingness indicators.
- Use DFS in TCGA and RFS in METABRIC separately; combine with NPC only as a sensitivity analysis.
- Exclude patients with stage IV disease.
- Exclude patients with zero follow-up time.
- Use an alternative NMF rank (k ± 1).
- Use Fisher and Tippett NPC combinations.
- Recalibrate the baseline hazard in METABRIC.
- Include HER2 in the METABRIC clinical model only.

**Execution status in this release.** The full-available-cohort complete-case sensitivity and the METABRIC IntClust benchmark were executed and are reported in `artifacts/` and the dashboard. The module family is controlled with BH at q = 0.10; the permutation maxT implementation is unit-tested but is not reported as a cohort-level result in this release. Fisher and Tippett NPC combiners are implemented and unit-tested; the cohort-level NPC result is Liptak, with Fisher/Tippett deferred. The Normal/claudin-low expansions, imputation sensitivity, stage-IV exclusion, alternative NMF rank, METABRIC baseline recalibration, and HER2-adjusted model are documented as prespecified extensions; they are not reported as results in this release. The audit log records this status so that a planned analysis is not mistaken for a completed one.

## 11. Power and simulation

A simulation study evaluates the operating characteristics of the framework:

- correlated multivariate outcomes with non-Gaussian errors;
- exponential censoring with 20–50% censoring;
- null and alternative effect sizes;
- compare marginal BH, Bonferroni, and NPC Fisher/Liptak;
- report type I error, power, and FDR at nominal 0.05.

The simulation is a methods check, not a substitute for the cohort analysis. Its code and seed are released.

## 12. Reporting and reproducibility

Every artifact records the input checksums, code version, random seed, and timestamp. The repository contains a data card, model card, methods document, audit log, tests, CI, Dockerfile, and a one-command pipeline. The dashboard is read-only and displays the prespecified primary results first; exploratory and sensitivity results are visually separated and labelled.

## 13. Deviations

Any deviation from this pre-registration is recorded in `docs/audit_log.md` with the date, reason, and whether it was made before or after seeing the relevant result. The analysis code reads the frozen configuration in `config/analysis.yaml`; changing a primary parameter requires a new configuration version.
