# Audit log

All deviations from the pre-registration are recorded here with a timestamp and whether the decision was made before or after seeing the relevant outcome result.

## 2026-11-23 — Project and cohort freeze

**Decision:** Use TCGA-BRCA PanCancer Atlas as discovery and METABRIC as external validation. Use PAM50/Claudin-low as the primary subtype system. Use IntClust as a METABRIC-only internal benchmark. Use clinical receptor status as the baseline.

**Rationale:** These are the two public cohorts with published intrinsic subtype calls and mature survival. IntClust has no external cohort with comparable labels, so it cannot be an external replication.

**Outcome knowledge at the time:** None. Only cohort sizes, variable names, and subtype value counts had been inspected. No survival model had been fitted.

## 2026-11-23 — TCGA grade removed from the primary clinical baseline

**Decision:** The primary clinical baseline is age, stage, ER, and PR. Grade is used only in a METABRIC sensitivity analysis.

**Rationale:** The cBioPortal TCGA PanCancer and Firehose clinical exports do not contain a usable grade field. Including grade only in METABRIC would make the discovery and validation models non-comparable. The pre-registration was amended before any survival model was fitted.

**Outcome knowledge at the time:** None. The decision was based on the variable list and missingness, not on survival results.

## 2026-11-23 — Primary endpoints and secondary endpoint naming

**Decision:** OS is the primary endpoint in both cohorts. DFS in TCGA and RFS in METABRIC are secondary. They are never pooled into a single curve; they are combined only through NPC with a prespecified caveat.

**Rationale:** DFS and RFS have different event definitions. Treating them as the same endpoint would overstate comparability.

## 2026-11-23 — Sample selection

**Decision:** Retain one primary tumour sample per patient, choosing the lexicographically smallest sample ID without reference to outcome.

**Rationale:** The cohorts are patient-level for survival. Selecting by sample ID is deterministic and outcome-blind.

## 2026-11-23 — Expression panel

**Decision:** Use a 995-gene Hallmark marker panel selected before outcome analysis: 20 most pathway-specific genes per Hallmark pathway, ties alphabetical.

**Rationale:** The full Hallmark union is ~4,400 genes and would require several gigabytes of JSON from the cBioPortal API. The marker panel keeps pathway-specific signal and is reproducible. The selection does not use survival.

## 2026-11-23 — NMF rank selection

**Decision:** Select the smallest rank within 0.02 silhouette of the best rank, among k = 3–8.

**Rationale:** Silhouette alone can favour a larger, less interpretable rank. The tolerance rule favours parsimony and was fixed before fitting.

## 2026-11-23 — NPC primary combiner

**Decision:** Liptak is the primary combiner; Fisher and Tippett are implemented as alternative combiners. The released cohort-level result is Liptak; the alternative combiners are not presented as executed cohort analyses (see the sensitivity execution status entry).

**Rationale:** Liptak accumulates consistent moderate evidence across dependent endpoints. The permutation calibration, not the asymptotic chi-square, is the basis of the p-value.

## 2026-11-23 — Primary FDR family

**Decision:** The primary family is four tests: subtype log-rank in TCGA OS, subtype log-rank in METABRIC OS, incremental LR in TCGA OS, and incremental LR in METABRIC OS. BH at q = 0.05.

**Rationale:** These are the four tests that determine the headline replication and increment decision. Secondary endpoints and module tests have separate families.

## 2026-11-23 — Analysis execution

**Decision:** Run the frozen pipeline with seed 20261123 and 2,000 permutations for the scorecard and NPC, 1,000 bootstrap resamples for C-index intervals, and 500 permutations for the module NPC.

**Outcome knowledge at the time:** The pipeline had not yet produced a survival estimate. The run is the first execution of the frozen configuration.

## 2026-11-23 — Complete-case sensitivity executed

**Decision:** Add the prespecified full-available-cohort sensitivity to the artifact: every patient with a four-class subtype and a valid endpoint, without the age/stage/ER/PR complete-case filter.

**Result:** The primary TCGA OS cohort has 868 patients and 110 events; the full subtype + OS cohort has 933 patients and 125 events. The overall subtype term changes from p = 0.155 to p = 0.015, but the Basal hazard ratio remains close to 1 (1.05 [0.64, 1.73]). The canonical Basal-worst pattern therefore still does not replicate. The primary result is power- and selection-sensitive, not evidence that PAM50 carries no information.

**Outcome knowledge at the time:** This sensitivity was run after the primary result. It is a prespecified sensitivity analysis, not a change to the primary estimand or decision rule.

## 2026-11-23 — Documentation correction: permutation budget

**Change:** The released configuration and pre-registration now state 2,000 permutations for the primary tests. The executed pipeline already used 2,000 (`scripts/run_analysis.py`), and the frozen `results.json` records 2,000. An earlier draft of `config/analysis.yaml` and `docs/preregistration.md` said 10,000.

**Reason:** The joint permutation loop is the dominant compute cost. The execution budget was 2,000, which gives a Monte Carlo floor of 0.0005. The draft text was not updated when the execution budget was fixed. The code now reads the count from `config/analysis.yaml` so the configuration and the run cannot diverge silently.

**Outcome knowledge at the time:** After the primary run. This is a documentation and reproducibility correction, not a re-analysis. It changes the precision of the permutation p-values, not the decision rule or the observed statistics. No artifact was regenerated.

## 2026-11-23 — Documentation correction: primary family definition

**Change:** The pre-registration now defines the primary family as two subtype log-rank tests and two incremental likelihood-ratio tests. The released `primary_family.csv` and `results.json` already use those four tests. An earlier draft described the family as two subtype tests plus two ΔC tests.

**Reason:** ΔC is the paired effect size for the incremental test; its bootstrap interval is reported alongside the LR p-value. It is not a separate p-value in the BH family. The executed family and the released artifacts are unchanged.

**Outcome knowledge at the time:** After the primary run. Documentation correction only; no estimand, test, or result changed.

## 2026-11-23 — Sensitivity execution status

**Decision:** Report the full-available-cohort sensitivity and the METABRIC IntClust benchmark as executed. Report Liptak as the cohort-level NPC result. Fisher/Tippett combiners and permutation maxT are implemented and unit-tested but are not presented as cohort-level results in this release.

**Reason:** The release should not describe a planned or unit-tested method as an executed cohort analysis. The pre-registration lists the remaining sensitivity analyses as prespecified extensions; they are not used to rescue or qualify the primary decision.

**Outcome knowledge at the time:** After the primary run. This is a reporting-scope decision, not a change to the primary analysis.

## 2026-11-23 — Deviation log template

Any future deviation must state: date, the exact change, the reason, whether it was before or after seeing the relevant result, and which artifact is affected. Changes that alter a primary estimand require a new configuration version.
