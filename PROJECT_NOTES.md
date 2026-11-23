# Replica — project notes

## Scope
A replication audit of breast cancer subtype systems across independent cohorts, with a censoring-aware nonparametric combination (NPC) framework, FDR control, external validation and a clinical evidence dashboard.

**Primary question:** Which published breast cancer subtypes replicate across cohorts, and do they add prognostic value beyond routine clinical variables?

**Project name:** Replica
**Repo:** imSheheryar/replica (public release target)
**Tagline:** Which breast cancer subtypes survive external validation?

## Cohorts
- TCGA-BRCA PanCancer Atlas: `brca_tcga_pan_can_atlas_2018` (1,084 samples). Clinical `SUBTYPE` values: BRCA_LumA 499, BRCA_LumB 197, BRCA_Basal 171, BRCA_Her2 78, BRCA_Normal 36. OS 1,084; DFS 942; PFS 1,083.
- TCGA-BRCA Firehose Legacy: `brca_tcga` (1,108 samples). Has ER_STATUS_BY_IHC (1,048), PR_STATUS_BY_IHC (1,047), stage, OS, DFS. No HER2, no PAM50.
- METABRIC: `brca_metabric` (2,509 samples). Patient-level `CLAUDIN_SUBTYPE` (PAM50 + claudin-low): LumA 700, LumB 475, Her2 224, claudin-low 218, Basal 209, Normal 148, NC 6. `INTCLUST` 1,980; `THREEGENE` 1,764; OS 1,981; RFS 2,388. Sample-level ER/PR/HER2, stage, grade. Patient-level treatment: chemo 1,980, hormone 1,980, radio 1,980.

## cBioPortal API
Base: `https://www.cbioportal.org/api`
- Clinical attributes: `/studies/{studyId}/clinical-attributes`
- Clinical data: `/studies/{studyId}/clinical-data?clinicalDataType=PATIENT|SAMPLE&projection=SUMMARY`
- Molecular profiles: `/studies/{studyId}/molecular-profiles`
- Expression fetch: `POST /molecular-profiles/{profileId}/molecular-data/fetch?projection=SUMMARY` with body `{"sampleIds":[...],"entrezGeneIds":[...]}` works (tested METABRIC).
- TCGA expression profile: `brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna`
- METABRIC expression profile: `brca_metabric_mrna`

## Data decisions
- Primary subtype system: PAM50/Claudin-low harmonized to LumA, LumB, Basal, HER2 for cross-cohort comparison; Normal/Claudin-low analyzed separately as sensitivity.
- Secondary: IntClust within METABRIC (internal benchmark; no TCGA labels). THREEGENE within METABRIC as another benchmark.
- Clinical covariates common to both: age, stage, grade, ER, PR. HER2 available in METABRIC but not TCGA PanCancer; use in METABRIC sensitivity only.
- Endpoints: OS in both; RFS in METABRIC, DFS in TCGA. Combine OS + RFS/DFS via NPC with caveat; report separately.
- Expression panel: MSigDB Hallmark 2023.2.Hs symbols GMT downloaded to `data/raw/hallmark.gmt` (50 sets, 4,384 unique genes). Need map symbols to Entrez IDs for cBioPortal.
- cBioPortal gene mapping endpoint test failed: `POST /genes/fetch` with `{"hugoGeneSymbols":[...]}` returned 400. Try `GET /genes/{symbol}` or `GET /genes?hugoGeneSymbol=...` or POST a bare list.

## Datahub
- `https://cbioportal-datahub.s3.amazonaws.com/{study}.tar.gz` returns 403.
- UCSC Xena S3 downloads return 403.
- Use cBioPortal API for clinical + expression subsets.

## Planned analyses
1. Cohort harmonization and data cards.
2. PAM50 survival replication: KM, Cox with clinical adjustment, forest plots.
3. External transport: lock Cox coefficients on TCGA, apply to METABRIC; C-index, calibration, ΔC-index vs clinical-only; bootstrap CIs and permutation p-values.
4. IntClust internal incremental value in METABRIC beyond PAM50 + clinical.
5. NMF modules from Hallmark expression panel in TCGA; module scores in METABRIC; NPC/permutation association with subtype and survival; BH/BY FDR.
6. NPC framework: combine OS + RFS/DFS partial tests; Fisher/Liptak/Tippett; permutation engine; FDR.
7. Simulation: correlated non-Gaussian censored data; type I error, power, FDR vs parametric baselines.
8. Static dashboard + FastAPI + Docker + CI + model/data cards + pre-registration + audit log.

## Environment
- Project root: `/home/ubuntu/projects/replica`
- Python environment: `.venv` (Python 3.12 locally; CI uses 3.11).

## Release checklist
- Frozen artifacts and dashboard rebuilt from `config/analysis.yaml`.
- Unit tests, Ruff, and artifact validation pass.
- README generated from `artifacts/results.json`.
- Public repository and GitHub Pages deployment verified.
