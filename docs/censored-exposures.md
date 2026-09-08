# Censored (below-detection-limit) exposures

Imputing exposures reported only as "below the limit of detection" jointly with the outcome, using the `censored_exposure_block_fcs` strategy.

> **Read the scope limit first.** This is validated for additive / per-analyte exposure–response functions and **must not be used for mixture / BKMR analyses**.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

When a focal **exposure** is left- or interval-censored below an analytical
reporting limit (LOD / MDL / LCMRL), imputing it *without the outcome* biases the
exposure–response estimate (an uncongenial imputation). The pipeline provides an
opt-in strategy that imputes such exposures **congenially** — outcome-aware,
skew-robust, and respecting the censoring interval — via a two-engine block-FCS:
`miceRanger` for the MAR covariates and `leftcens::impute_censored_conditional()`
for the censored exposures. It handles multiple exposures and the three-tier
(non-detect / detected-not-quantified / quantified) interval structure.

- **Enable it** with `analysis_spec$imputation$strategy <- "censored_exposure_block_fcs"`
  plus an `analysis_spec$imputation$censored_exposure` block. See
  [`00_censored_exposure.R`](../00_censored_exposure.R) and
  [`examples/censored_exposure/`](../examples/censored_exposure/).
- **Input convention** (leftcens interval columns): each censored exposure `X`
  carries per-row bounds `X_lo` / `X_hi` — `X_lo == X_hi` where quantified;
  `X_lo <= 0`, `X_hi = LOD` for a left-censored non-detect; both finite for a
  detected-not-quantified interval. In the dictionary mark each exposure
  `role = exposure`, `impute_target = FALSE`, `use_in_model = TRUE`.
- **Requires** `leftcens` (>= 0.9.0), which exports `impute_censored_conditional()`.
  Not on CRAN: `remotes::install_github("fabregithub/leftcens@v0.9.0")`.
- **Config knobs** (`censored_exposure`): `exposure_vars`; a reduced `predictors`
  set (condition on the *other* analytes + outcome + key determinants, not the
  full covariate set — the X-block cost is ~`k^2`); `outer_sweeps`; `margin`
  (`"shash"` skew-aware / `"gaussian"` tobit); `log_scale`; `mid_delete_imputed_y`
  (multiple-imputation-then-deletion for missing outcomes); and `n_cores`
  (parallelises the `m` completed datasets).
- Produces the standard `imputed_###.rds` + manifest, so **Steps 4–12 are
  unchanged**. Method background, the validation study, and the design rationale are
  in [`validation/`](../validation/) — start at
  [`validation/README.md`](../validation/README.md), the status index for what has been
  validated, what is open, and what is next.

**Validation status.** On a Monte-Carlo study with a known exposure–response function
(300 replications per cell, `m` = 30, `n` = 800), this strategy recovers the focal
censored-exposure coefficient with **relative bias under 1% and 95% interval coverage
of 0.957** at both 20% and 40% non-detects — statistically indistinguishable from an
oracle fitted on the uncensored data. Imputing the same exposure *without* the outcome
(the conventional pre-step) attenuates that coefficient by **−5.9% at 20% ND and −14.1%
at 40% ND**, with coverage falling to 0.82. The strategy is also markedly more efficient
than complete-case analysis (RMSE 0.049 vs 0.072 at 40% ND). These figures are for the
shipped pipeline code, not a prototype of it — full tables and caveats in
[`validation/phase1/FINDINGS_v1.md`](../validation/phase1/FINDINGS_v1.md).

> ## ⚠️ Scope: additive exposure–response functions only. **Do not use this strategy for
> mixture / BKMR analyses.**
>
> The claim above holds for **additive / per-analyte** exposure–response functions.
>
> For **mixture / BKMR-style surfaces** with interactions and curvature, this strategy is
> **not merely inadequate — it is worse than doing nothing sophisticated.** Measured on the
> estimands a BKMR analysis actually reports, against analytic truth, at 200 replications
> per cell:
>
> | | pipeline | LOD/√2 substitution | complete case |
> |---|---|---|---|
> | mean abs. error vs the oracle, 21 cells | **14.6%** | **6.4%** | 10.5% |
> | pure curvature, 40% non-detects | **−57%** | −1.9% *(n.s.)* | −18% |
>
> It **destroys 57% of the curvature** at 40% non-detects, and on that estimand LOD/√2
> substitution is statistically indistinguishable from the oracle while this strategy is
> not. The reason is structural: the X block draws the censored exposure from a conditional
> that is **linear in the predictors**, so it imposes linearity on precisely the rows where
> curvature would appear. Imputing with the wrong functional form is worse than not
> imputing.
>
> Interaction estimands are more robust — they survive censoring of a single exposure
> (+2.6% to +5.7%) and break down (−12.5%) only when a second exposure is censored too.
>
> **Coverage will not warn you.** In these runs a −42% point bias sat behind 0.945–0.960
> coverage, because the pooled intervals ran 1.5–1.6× wider than the estimates' true
> sampling spread.
>
> Mixture analyses need **substantive-model-compatible imputation**, which this pipeline
> does not implement. Evidence:
> [`validation/phase1/FINDINGS_v3.md`](../validation/phase1/FINDINGS_v3.md); background in
> [`validation/PLAN_leftcensored_exposure_integration.md`](../validation/PLAN_leftcensored_exposure_integration.md) §7.7.

A follow-up robustness sweep (10 scenarios × 300 replications) confirmed the bias result
holds well outside the tested corner — with **every** exposure censored, at **n = 80 000**,
under right-skew, and at exposure correlation 0.8, relative bias never exceeded 2.7%.

**These results also hold under MAR, not only MCAR.** Every figure above was originally
measured with covariates missing *completely* at random. A 12-scenario, 1000-replication
sweep with covariate missingness driven by the outcome and an observed exposure — genuine
MAR, verified conditionally independent of the missing values themselves — found the engine
**slightly less** biased than under MCAR at matched missing fractions (−0.17% against
−1.67% at 40% missing), with coverage 0.948–0.967. That held across a doubling of the
mechanism's strength, with non-linear covariates, and with 20% missing outcomes on top.
The MCAR-based numbers above are therefore conservative rather than optimistic. See
[`validation/phase1/FINDINGS_v10.md`](../validation/phase1/FINDINGS_v10.md). Non-ignorable
(MNAR) missingness is outside what any imputation engine can address and is a matter of
study design, as noted above.

> **The figures above assume the outcome is linear in the covariates.** Every validation
> cell up to V10 generated `Y` as a linear function of `(log X, Z)`. When the outcome is
> instead curved in a covariate, **bias degrades by about 3 percentage points regardless of
> which imputer is used** — the shipped `bart` default moves from 0.53% to 3.10% mean
> absolute bias, and the other three imputers degrade by a statistically indistinguishable
> amount (spread 0.44 pp). This is a property of the pipeline rather than of the imputer
> choice, and it is not captured by the figures above.
>
> The mechanism is general, not specific to censored exposures: the Z block imputes
> covariates *conditioning on the outcome*, so when the outcome is a non-linear function of
> a covariate, that conditional is harder for every imputer. Any analysis with missing
> covariates is affected.
>
> **Coverage does not warn you** — it was unchanged at 0.953, because at `n` = 800 a 3%
> relative bias is only ≈0.28 standard errors. Point estimates move; intervals do not
> notice. This is the third consecutive validation track where calibration diagnostics
> missed a real bias finding.
>
> Measured with one curvature form and a **correctly specified** analysis model, so this is
> imputation-induced bias alone; omitting the non-linear term from your own model would add
> ordinary misspecification on top. Evidence:
> [`validation/phase1/FINDINGS_v11.md`](../validation/phase1/FINDINGS_v11.md).

> **Caveat on interval width — applies to the whole pipeline, not just this strategy.**
> The sweep found that 95% credible intervals become **anti-conservative when a
> substantial fraction of the data is imputed**: coverage was 0.95–0.97 with complete
> covariates, 0.93 at 40% MCAR covariates, and **0.89** with 20% MCAR covariates plus 20%
> missing outcomes. Point estimates stayed sound throughout (bias ≤2.7%) — the intervals
> are too narrow, not mis-centred.
>
> A follow-up study **attributed the cause**: the `miceRanger` imputation step is
> *improper* multiple imputation in Rubin's sense — it does not draw imputation-model
> parameters from a posterior — so the between-imputation variance is understated by
> **8–59%**, scaling with how much is imputed. Replacing that step with a proper Bayesian
> draw restores calibration. This is a property of the pipeline's **general** row-level
> imputation path, so it applies to any analysis with substantial missingness, not only to
> censored exposures. Raising `m` does not help, and the
> multiple-imputation-then-deletion step is not implicated (disabling it makes matters
> worse). Details and numbers:
> [`validation/phase1/FINDINGS_v4.md`](../validation/phase1/FINDINGS_v4.md).
>
> **Corrected in v1.5.0.** `analysis_spec$imputation$z_imputer` defaults to `"bart"`:
> the Z block draws imputations from a BART posterior, which is proper multiple imputation
> (parameters drawn from a posterior) *and* flexible enough for non-linear covariates.
> Across 7 scenarios × 300 replications this gave the best bias of any variant tested
> (≤1.3% everywhere) with **coverage 0.937–0.963**, and it costs ~4% of pipeline runtime.
>
> **⚠️ Added 2026-09-08 — a separate and larger limit.** `leftcens` draws each below-LOD
> value from a conditional **linear in its predictors** (the skew-aware part is the *margin*,
> not the mean). If a covariate has a genuinely **non-linear** relationship with a censored
> exposure, that cannot be represented: measured at **~20% bias in the exposure coefficient,
> flat from n = 800 to 12,800** — asymptotic, and not coverage-detectable. Check for curved or
> saturating covariate–exposure relationships before trusting this path, and pass linearising
> terms explicitly via `censored_exposure$predictors` if you find one. Detail and remedies:
> [Which covariates to adjust for](covariate-roles.md#a-non-linear-covariateexposure-relationship-breaks-the-censored-exposure-draw).
>
> **Scope, added 2026-09-07.** Every one of those 7 scenarios had a *causally inert*
> covariate — independent of the exposures, or generated from them. Where a covariate is a
> **confounder or a mediator** with substantial missingness, the same default carries **+5%
> to +10% bias** at n = 800 (falling to ~2–4% at n = 12,800) and **coverage does not reveal
> it**. The ≤1.3% figure does not transfer to that case. See
> [Which covariates to adjust for](covariate-roles.md).
>
> **BART is nonetheless the right default, and more clearly than before**: it is the only
> Z-block imputer measured whose bias *shrinks* with sample size. Requires the `dbarts`
> package; **if absent the pipeline warns and falls back to `forest_boot`**, which under a
> confounder carries **−23% bias that does not decay, with coverage reaching zero**. Install
> `dbarts` before a real analysis and check the log's `Z-block imputer:` line to confirm
> which imputer ran. Setting `z_imputer = "forest_boot"` or `"forest"` reproduces a v1.4.0 or
> pre-v1.4.0 analysis exactly — a **behaviour change**, so cite the pipeline version used,
> and do not expect the covariate path to be unbiased while you do.
>
> **What is and is not established.** Coverage is 0.937–0.963 across every tested
> condition. Interval *width* readings ran 0–8% above the ideal, but those deviations are
> within the diagnostic's own Monte-Carlo error at 300 replications — none is statistically
> significant. So the honest statement is that interval width is **consistent with
> calibration, bounded to roughly ±8%**, rather than proven exact. Full analysis, including
> two refuted alternative explanations, in
> [`validation/phase1/FINDINGS_v7.md`](../validation/phase1/FINDINGS_v7.md).
>
> **Recommended `m`:** 30, which is the shipped default. Interval width stabilises there
> and the fraction-of-missing-information rule agrees.

<details>
<summary>Earlier history of this caveat (v1.4.0)</summary>

> **Largely corrected in v1.4.0.** `analysis_spec$imputation$proper_draw` now defaults to
> `TRUE`, bootstrapping the imputation model's training data once per imputation so that
> the uncertainty reaches the intervals. This closes **66%** of the shortfall at 40%
> covariate missingness (coverage 0.933 → 0.947) and improves every case tested, without
> ever over-correcting. Set `proper_draw = FALSE` to reproduce a pre-v1.4.0 analysis
> exactly — note that this is a **behaviour change**, so cite the pipeline version used.
>
> **Residual limitation.** The correction is partial. Under *heavy* missingness in both
> covariates and the outcome, coverage reaches about **0.91** rather than 0.95, because a
> random forest is already a bagged estimator and an outer bootstrap injects less
> parameter uncertainty than a parametric posterior draw would. So in that regime: point
> estimates and their ordering are reliable, but treat interval widths as a mild lower
> bound and be cautious about borderline "significant" findings. Full analysis in
> [`validation/phase1/FINDINGS_v5.md`](../validation/phase1/FINDINGS_v5.md).

</details>

---

---

*[← Back to the main README](../README.md)*
