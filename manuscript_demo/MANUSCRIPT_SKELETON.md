# [Working title] Congenial multiple imputation of left- and interval-censored chemical exposures for single-pollutant and mixture health-effect models: a demonstration with prenatal PFAS and childhood outcomes in JECS

*Manuscript skeleton — a demonstration/methods paper for the generic MI + brms
pipeline with the `leftcens` censored-exposure block-FCS. Results tables/figures
are produced by the scripts in this folder (`analyze_imputation.R`,
`run_scenarios.R`, `analyze_bkmr.R`); the demonstration data are fully synthetic
and grounded in three JECS papers (see `PROVENANCE.md`).*

---

## Abstract (structured)

- **Background.** Environmental chemical exposures are routinely left-/interval-
  censored below reporting limits. The common fix — substituting a fraction of the
  limit (e.g. MDL/2) — is known to bias exposure–response estimates, and it is still
  widely used (incl. the reference JECS analyses). When the censored analyte is a
  focal *exposure*, imputing it *without the outcome* is uncongenial and biases the
  effect estimate.
- **Objective.** Demonstrate a reusable pipeline that imputes censored exposures
  **congenially** (outcome-aware, skew-robust, respecting the censoring interval)
  and feeds both single-pollutant and mixture (BKMR) models, and quantify when it
  matters.
- **Methods.** Synthetic cohort grounded in JECS (7 prenatal PFAS → Kawasaki
  disease [rare] and childhood wheeze/asthma [common]). Four methods compared
  against a complete-data oracle: half-limit substitution, and outcome-aware vs
  outcome-blind congenial imputation. Robustness across censoring rate and
  inter-analyte correlation.
- **Results.** [Fill from results.] Imputation ≫ substitution throughout; the
  outcome-inclusion (congeniality) benefit is *conditional* — greatest under low
  analyte correlation and heavy censoring.
- **Conclusion.** [Fill.] Practical guidance on when outcome-aware censored
  imputation is essential.

**Keywords:** left-censoring, interval-censoring, multiple imputation, congeniality,
BKMR, PFAS, exposure mixtures, JECS.

---

## 1. Introduction

- The censoring problem in exposure epidemiology (LOD/MDL/LCMRL; ND and DNQ tiers).
- Substitution methods (MDL/2, MDL/√2) and their documented bias (Lubin 2004;
  Helsel). Note the reference JECS papers used substitution (Iwata) / a Gibbs
  imputation not outcome-aware (Atagi).
- Multiple imputation for censored data; **congeniality** (Meng 1994): a focal
  censored *exposure* must be imputed using the outcome, or the exposure–response
  coefficient is biased (attenuated). Contrast with the outcome-vs-predictor
  asymmetry (von Hippel 2007).
- Mixtures: correlated exposures, WQS/BKMR (Bobb 2015); censoring feeds these too.
- **Contribution / gap.** A single reusable pipeline that (i) imputes censored
  exposures congenially and skew-robustly, (ii) scales to large cohorts, (iii) feeds
  single-pollutant *and* mixture models — and a quantification of *when* the
  outcome-aware step is essential. [Cite the pipeline + `leftcens`.]

## 2. Methods

### 2.1 Congenial censored-exposure imputation
- The block-FCS: `miceRanger` for MAR covariates + `leftcens::impute_censored_conditional()`
  for the censored exposures — a per-analyte, interval-censored, **sinh-arcsinh
  (skew-aware)** draw on the latent scale, conditioning on the outcome and covariates,
  proper MI. Multiple-imputation-then-deletion for missing outcomes.
- Reference the pipeline (`00_censored_exposure.R`) and the `leftcens` package.

### 2.2 Methods compared (vs a complete-data oracle)
1. **Oracle** — uncensored exposures (upper bound on performance).
2. **Substitution** — half-limit fill (the reference-paper method).
3. **Impute, outcome-blind** — congenial machinery but *without* the outcome.
4. **Impute, outcome-aware** — the pipeline.

### 2.3 Estimands and analysis models
- **Single-pollutant:** logistic odds ratio per log2 ng/mL, adjusted, one PFAS at a
  time; pooled by Rubin's rules.
- **Mixture:** BKMR overall mixture effect (q25→q75), `bkmr::kmbayes` per imputation,
  pooled; on an Iwata-style subsample (BKMR is MCMC-heavy — note the accelerated /
  A-BKMR option for large n). WQS optional.
- **Scoring:** bias, RMSE, and interval coverage **relative to the oracle analysis**.

### 2.4 Demonstration data (synthetic, JECS-grounded)
- 7 prenatal PFAS (medians, MDL/LCMRL, correlations from Lai 2025); three-tier
  interval censoring (ND / DNQ / quantified); Z→X determinants (parity, age, region).
- Two outcomes: Kawasaki disease (rare; Iwata 2024) and childhood wheeze/asthma
  (common; Atagi 2024). See `PROVENANCE.md` (Table S1) for the full parameter map.
- **Fully synthetic; not JECS participant data.** Known truth enables oracle scoring.

### 2.5 Robustness scenarios
- Censoring rate (real vs increased, via `limit_scale`) × inter-analyte correlation
  (low / real / high) × outcome (rare / common). Reps for Monte-Carlo error.

## 3. Results

### 3.1 Single-pollutant odds ratios (Table 2 / Fig 2)
- Imputation vs substitution vs oracle across PFAS. *[from `results/single_pollutant_*.csv`.]*
- **Y-inclusion contrast** (outcome-aware vs outcome-blind): conditional on
  correlation/censoring. *[from `results/scenario_grid.csv`.]*

### 3.2 Mixture / BKMR (Table 3 / Fig 3)
- Overall PFAS-mixture effect on KD under each method vs oracle. *[from
  `results/mixture_*.csv`.]* Single-exposure BKMR effects / PIPs; comparison with
  the substitution-fed BKMR (the §7.7 non-linear-congeniality question on real-grounded data).

### 3.3 Robustness surface (Fig 4)
- Method error vs oracle across the censoring × correlation grid; identify the regime
  where outcome-aware imputation is essential.

## 4. Discussion
- Headline: substitution is dominated by congenial imputation for both estimands.
- **When does outcome-inclusion matter?** Weak analyte correlation + heavy censoring;
  when analytes are strongly correlated they rescue each other (reconciles with the
  Gaussian-copula rationale, Lai 2025). Practical guidance.
- Rare vs common outcomes: information in the outcome bounds the achievable benefit.
- Mixtures: does congenial *linear* imputation feed BKMR adequately, or is a
  substantive-model-compatible draw needed? [§7.7.]
- Limitations: synthetic demonstration (though grounded); BKMR compute; single
  imputation model family; partial identification of the sub-limit exposure–response.

## 5. Conclusion
- [Fill.] A reusable, scalable route to congenial censored-exposure analysis for
  both single-pollutant and mixture models, with evidence on when it is essential.

## Reproducibility
- Data: `R/make_demo_data.R` (synthetic; seed-fixed). Analyses: `analyze_imputation.R`,
  `run_scenarios.R`, `analyze_bkmr.R`. Pipeline example: `examples/jecs_pfas_kd/`.
  Provenance of every DGP parameter: `PROVENANCE.md`. Software: the MI + brms
  pipeline; `leftcens` (>= 0.9.0); `bkmr`.

## References (to complete)
- Lai et al. 2025, *Ecotox. Environ. Saf.* 294:118107 (PFAS determinants).
- Iwata et al. 2024, *Environ. Int.* 183:108321 (PFAS & Kawasaki disease; BKMR/WQS).
- Atagi et al. 2024, *Environ. Res.* 240:117499 (PFAS & wheeze/asthma).
- Rubin 1987; Meng 1994 (congeniality); von Hippel 2007 (MID); Bartlett 2015 (SMC-FCS);
  Bobb et al. 2015 (BKMR); Lubin et al. 2004 (substitution bias); Helsel (NADA).
