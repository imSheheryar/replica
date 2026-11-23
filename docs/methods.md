# Methods — Replica

This document is the technical companion to the pre-registration. It gives the estimands, the exact permutation mechanics, the NPC calibration, and the failure modes.

## 1. Notation

For patient \(i=1,\ldots,n\), let \(T_i\) be the observed follow-up time, \(\delta_i\in\{0,1\}\) the event indicator, \(G_i\in\{\text{LumA},\text{LumB},\text{Basal},\text{Her2}\}\) the subtype, and \(Z_i\) the clinical covariate vector. The primary endpoint is overall survival (OS); the secondary endpoint is DFS in TCGA and RFS in METABRIC.

The Cox model is

\[
\lambda(t\mid G_i,Z_i)=\lambda_0(t)\exp\left(\beta^\top G_i+\gamma^\top Z_i\right),
\]

with LumA as the reference category. The subtype null is \(H_0:\beta=0\). The incremental null is that adding \(G\) does not improve the fit over \(Z\) alone.

## 2. Log-rank statistic

At each distinct event time \(t_j\), let \(n_j\) be the number at risk, \(d_j\) the number of events, \(n_{gj}\) the number at risk in group \(g\), and \(d_{gj}\) the events in group \(g\). The observed and expected event counts are

\[
O_g=\sum_j d_{gj},\qquad E_g=\sum_j d_j\frac{n_{gj}}{n_j}.
\]

The covariance contribution is

\[
V_j=\frac{d_j(n_j-d_j)}{n_j-1}\left(\operatorname{diag}(p_j)-p_jp_j^\top\right),
\quad p_{gj}=\frac{n_{gj}}{n_j}.
\]

The multi-group statistic is

\[
X^2=(O-E)_{1:G-1}^\top
\left(\sum_j V_j\right)_{1:G-1,1:G-1}^{-1}
(O-E)_{1:G-1},
\]

which is asymptotically \(\chi^2_{G-1}\) under the null and independent censoring. Replica uses the statistic as a permutation statistic, not only as an asymptotic test.

## 3. Censoring-aware permutation

Permuting the event indicator alone changes the censoring distribution and can invalidate the test. Replica permutes the **bivariate outcome pair** \((T_i,\delta_i)\). The pair carries both the event time and the censoring information, so the marginal censoring distribution is preserved.

For the unadjusted test:

1. Compute \(X^2_{\text{obs}}\) on the observed data.
2. For \(b=1,\ldots,B\), draw a permutation \(\pi_b\) of the rows, form \((T_{\pi_b(i)},\delta_{\pi_b(i)})\), and compute \(X^2_b\).
3. Report \(p=(1+\#\{b:X^2_b\ge X^2_{\text{obs}}\})/(B+1)\).

This is a conditional randomisation test. It is valid under the sharp null that subtype is independent of the survival outcome, with the usual independent-censoring assumption. It does not require proportional hazards, exponential tails, or a parametric copula.

For adjusted tests, the pair is permuted **within strata** of age tertile × stage (I/II vs III/IV) × ER status. This keeps the clinical covariate–outcome relationship fixed under the null and asks whether subtype adds information beyond those strata. Strata with fewer than 10 observations are not permuted independently; they are pooled with the nearest stratum before the test. The realised strata are recorded in the artifact metadata.

The permutation test does not address confounding, measurement error, selection bias, or non-random missingness. It tests exchangeability, not causation.

## 4. Nonparametric combination

Let \(T_1,\ldots,T_K\) be \(K\) test statistics for partially dependent hypotheses, all oriented so that larger values mean more evidence. Examples: OS and DFS/RFS log-rank statistics in two cohorts; module-level Cox statistics.

For each component \(k\), the marginal permutation p-value is

\[
p_k=\frac{1+\sum_{b=1}^{B}\mathbf{1}\{T_{k,b}\ge T_{k,\text{obs}}\}}{B+1}.
\]

Replica implements three combining functions:

\[
C_{\text{Fisher}}=-2\sum_{k=1}^{K}\log p_k,
\]

\[
C_{\text{Liptak}}=\sum_{k=1}^{K}\Phi^{-1}(1-p_k),
\]

\[
C_{\text{Tippett}}=\max_{k}(1-p_k).
\]

The naive null distribution of \(C\) assumes independence across components. That assumption is false here: OS and DFS/RFS are measured on the same patients, and the two cohorts share subtype biology. Replica therefore calibrates \(C\) against the **joint permutation distribution**. For every permutation \(b\), the same row permutation is applied to all endpoints within a cohort, and the component statistics are recomputed. The resulting \(B\times K\) matrix preserves the within-person and within-cohort dependence. Marginal p-values for each permutation are obtained by ranking each column, and the combined statistic is recomputed for every row. The NPC p-value is

\[
p_{\text{NPC}}=\frac{1+\#\{b:C_b\ge C_{\text{obs}}\}}{B+1}.
\]

The global test combines four components: TCGA OS, TCGA DFS, METABRIC OS, METABRIC RFS. Within-cohort permutations are applied jointly to OS and DFS/RFS, so the correlation between the two endpoints is preserved. Across cohorts, the rows are independent and are paired arbitrarily.

The primary combiner is Liptak because it accumulates consistent moderate evidence. Fisher and Tippett are implemented as alternative combiners and unit-tested; the released cohort-level NPC result is Liptak.

## 5. Multiplicity

The primary family has four tests: subtype log-rank in TCGA OS, subtype log-rank in METABRIC OS, incremental likelihood-ratio test in TCGA OS, and incremental likelihood-ratio test in METABRIC OS. Benjamini–Hochberg at \(q=0.05\) is applied to this family. Benjamini–Yekutieli is reported as a dependence-robust sensitivity analysis.

For the module family, BH is applied at \(q=0.10\) across modules and endpoints. A permutation maxT implementation is included and unit-tested, using the maximum statistic across modules in each permutation; the released module family is controlled with BH.

BH is

\[
q_{(i)}=\min_{j\ge i}\left(\frac{n}{j}p_{(j)}\right),
\]

with the usual monotone enforcement. BY multiplies by the harmonic number \(H_n=\sum_{j=1}^n 1/j\).

## 6. Incremental value

The clinical baseline is age per 10 years, stage III/IV vs I/II, ER positive vs negative, and PR positive vs negative. TCGA grade is not available in the cBioPortal export, so grade is not in the primary model; it is a METABRIC-only sensitivity analysis.

Incremental value is measured with **out-of-fold** predictions from 5-fold cross-validation, stratified by event status. The same folds are used for the clinical and clinical + subtype models. The estimand is

\[
\Delta C = C_{\text{full}}^{\text{OOF}} - C_{\text{clinical}}^{\text{OOF}}.
\]

The 95% interval is a paired bootstrap over patients, resampling the out-of-fold predictions. The likelihood-ratio test on the full data is reported as a secondary criterion. A positive ΔC with a CI excluding 0 and a BH-adjusted LR p < 0.05 is required for "incremental".

C-index is Harrell's concordance: the probability that, for a randomly chosen comparable pair, the patient with the shorter survival has the higher predicted risk. Ties in predictions count 0.5.

## 7. External transport

The clinical + subtype Cox model is fitted on TCGA and frozen. The fitted coefficients are applied to METABRIC, and the C-index is computed on the validation risk scores. The clinical-only model is transported in the same way, and the difference is bootstrapped.

Calibration at 60 months uses the TCGA baseline hazard:

\[
\hat S_i(60)=\hat S_0(60)^{\exp(\hat\eta_i)}.
\]

Patients are grouped into quartiles of predicted survival, and the observed KM survival at 60 months is compared with the mean prediction. A calibration slope is obtained by fitting a univariable Cox model with the locked linear predictor in METABRIC; a slope near 1 indicates that the log-hazard scale transports, while a slope below 1 indicates over-dispersion.

A cohort-specific baseline hazard is not re-estimated in the primary transport analysis, because that would change the model. The calibration plot is therefore expected to show baseline-hazard differences. That is a finding about transport, not a bug.

## 8. IntClust benchmark

IntClust is available only in METABRIC. Three models are compared with 5-fold cross-validated C-index:

1. clinical covariates;
2. clinical + PAM50;
3. clinical + PAM50 + IntClust dummies.

The question is whether a 10-class integrative cluster adds information beyond PAM50 and clinical variables in the cohort where it was derived. This is an internal benchmark, not external validation.

## 9. Expression modules

The expression panel is fixed before outcome analysis. For each of the 50 Hallmark pathways, Replica takes the 20 genes that occur in the fewest Hallmark sets, breaking ties alphabetically. This yields 995 mapped genes. The selection is pathway-specific and does not use survival.

Within TCGA, expression is log2-transformed where needed and genes are z-scored. NMF is fitted to the top 500 variable genes:

\[
X\approx WH,\quad W\ge0,\ H\ge0,
\]

with \(X\) the samples × genes matrix, \(W\) the sample loadings, and \(H\) the gene loadings. Ranks \(k=3,\ldots,8\) are compared by mean silhouette of the argmax module assignment on \(W\) and by reconstruction error. The selected rank is the smallest \(k\) within 0.02 silhouette of the best, favouring interpretability.

Module scores in METABRIC are the mean z-scored expression of the top 30 genes by \(H\) for each module. No METABRIC outcome information enters the gene list or the score. This is a locked projection, not a re-fitting.

Module association with subtype is tested by Kruskal–Wallis. Module association with OS is tested by a Cox model per standard deviation of the module score, adjusted for the clinical baseline, with a likelihood-ratio test. BH at \(q=0.10\) controls the module family.

## 10. Simulation study

The simulation evaluates operating characteristics under correlated, non-Gaussian, censored outcomes:

- \(n=160\), \(p=5\) endpoints, correlation \(\rho=0.45\);
- multivariate \(t_5\) errors;
- exponential censoring with 20–50% censoring;
- null and alternative shift sizes;
- 500 replicates per scenario;
- compare marginal BH, Bonferroni, and NPC Fisher/Liptak.

The simulation is a methods check, not evidence about breast cancer. It answers whether the NPC calibration controls type I error when endpoints are dependent and non-normal.

## 11. Software and reproducibility

Python 3.10+, pandas, NumPy, SciPy, scikit-learn, lifelines, statsmodels, Plotly, FastAPI. The pipeline is a sequence of scripts with a Makefile; every random draw uses seed 20261123. Artifacts include input checksums, package versions, and a timestamp. Tests cover the NPC calibration, BH/BY, harmonisation, Cox fitting, and permutation behaviour.

## 12. Known limitations

The subtype labels are not assay-identical across cohorts. TCGA labels come from a published PAM50 classifier; METABRIC labels come from a different pipeline. The clinical covariates are not measured identically. TCGA grade is missing from the public export. DFS and RFS are different endpoints. Complete-case analysis can bias estimates if missingness is informative. The audit is observational and cannot identify a causal effect of subtype on survival.

The complete-case filter matters empirically. The primary TCGA OS cohort has 868 patients and 110 events. Using every patient with a four-class subtype and valid OS raises this to 933 patients and 125 events; the overall subtype term becomes significant (permutation p = 0.015), but the Basal hazard ratio is 1.05 [0.64, 1.73]. The canonical Basal-worst pattern therefore still does not replicate. The primary negative result should be read as power- and selection-sensitive, not as proof that PAM50 carries no prognostic information.
