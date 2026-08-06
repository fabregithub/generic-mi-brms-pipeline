# Plan: congenial integration of left-censored exposure imputation into the MI + brms pipeline

Design/plan document for `generic-mi-brms-pipeline`. Captures the decision and the
rationale so the project can restart cold. Status: **agreed direction, not yet
implemented.**

---

## 1. Motivation & scope

The pipeline imputes MAR/MCAR data with **`miceRanger`**, fits **`brms`**, and pools
with Rubin's rules. It requires MNAR to be handled *before* entry. Left-censored
exposures (values below a lab reporting limit / LOD) are a special, benign MNAR:
the value exists, the mechanism (below LOD) is known, and the arbitrary lab cut
point is exactly what we want to correct for. `leftcens`/`gsimp` is our engine for
that.

**Problem discovered.** When the censored analyte is the **exposure** (a focal
*predictor* of a health outcome), imputing it as a **pre-step with `leftcens`
alone — i.e., without the outcome `Y`** — biases the exposure–response coefficient.
This must be fixed for the exposure–outcome estimand (the ERF). It does *not*
affect using `leftcens` for **descriptive** exposure statistics.

Two cases (only the first needs this rewrite):
- **Exposure as predictor** (e.g. chemical exposure → health outcome): needs
  congenial, joint handling. **This document.**
- **Censored analyte as outcome** (determinants of exposure): prefer modelling the
  censoring directly in `brms` (`resp | cens(lod) ~ ...`); no imputation of the
  outcome needed.

---

## 2. Theory — why this is necessary (not just efficiency)

1. **Congeniality.** Valid MI draws the missing predictor from
   `p(X | X_obs, Z, Y) ∝ p(Y | X, Z) · p(X | Z)`. `Y` appears *because* `X → Y`
   (Bayes inversion of the causal likelihood — using `Y` to impute `X` is **not**
   asserting `Y → X`). Dropping `Y` ⟺ assuming `X ⟂ Y | Z` ⟺ assuming `β = 0`,
   the very null under test. Result: `β̂` is **inconsistent** (bias direction is
   DGP-specific; for MCAR it is attenuation `(1−π)β`, for left-censoring it can go
   the other way — a theorem about inconsistency, not a sign).

2. **Marginal vs joint.** `leftcens` (no `Y`) reconstructs the **marginal**
   `p(X | analytes)` correctly — right for descriptives — but sets
   `Cov(X_imp, Y | Z) = 0` in the censored cells, destroying the **joint**
   `p(X, Y)`. The ERF (`β = Cov(X,Y|Z)/Var(X|Z)`) is a **joint** functional.
   Correct-for-the-marginal ≠ correct-for-the-ERF.

3. **The low-dose ERF.** The only data signal about the low-dose ERF is the
   censored subjects' **outcomes** `Y` (their `X` is unobserved). Imputing `X`
   without `Y` — or complete-case — discards that signal; the low-dose ERF then
   comes purely from the modelled tail. Using `Y` is what lets the low-dose
   outcomes speak.

4. **Outcome vs predictor asymmetry** (reconciles the pipeline's current defaults):
   - Missing **outcome** `Y` → do **not** impute it; estimate `β` from observed-`Y`
     rows, predict missing `Y` post-hoc (von Hippel 2007). The pipeline is right.
   - Missing **focal predictor** `X` → **must** use `Y` in its imputation. This is
     the one case the pipeline's "don't put `Y` in the MI" default does not cover.

5. **Irreducible limit.** Even done correctly, the ERF **shape** below the LOD is
   only *partially identified* (true low-dose `X` never observed); the modelled
   tail fills the gap. Report transparently. No method escapes this.

References: Rubin (1987); Meng (1994, congeniality); Moons (2006) & Sterne (2009,
include the outcome); von Hippel (2007, missing Ys; 2020, how many imputations);
Bartlett (2015, substantive-model-compatible FCS).

---

## 3. Architecture decision: two-engine **block-FCS**

`mice` (per-variable model refits) is too slow at ~80k × ~300; **`miceRanger`
(ranger backend) is required for speed** and is kept. But a single RF cannot
impute a censored `X` (it ignores the LOD bound), and `Y` must inform `X`. So run
a **block fully-conditional-specification** (a standard MICE extension), each outer
sweep:

- **Z block** (≈300 MAR covariates): `miceRanger`, conditioning on current `X`, `Y`.
  Fast RF, unchanged.
- **X block** (censored exposure): a **censored, skew-aware conditional draw**,
  conditioning on current `Z`, `Y`.
- Alternate for a few outer sweeps until stable. Both engines are fast, so 3–5
  sweeps is cheap even at scale.

`Y`: predictor (`use_in_model = TRUE`), not imputed. Handle its missingness by
**MID** (multiple-imputation-then-deletion): impute `Y` so it is complete for the
FCS predictor step, then delete the imputed-`Y` rows before the `brms` analysis
(they carry no information for `β`).

Why not the matrix `gsimp_impute()`: it conditions the censored analyte on the
**other analytes only** (continuous, positive, log scale). Our predictors are
**mixed-type** (`Y`, categorical/continuous `Z`) and cannot be bolted on as
"analyte columns" (the log-Gaussian copula assumes continuous positive margins).
We need the single-variable, general-predictor form below.

---

## 4. The one component to build: general-predictor interval-censored, skew-aware conditional draw

A compact function — the reusable core, usable both as a `mice.impute.leftcens`
method and standalone in the block loop. Contract like a `mice.impute.*` method.

**Input:** the censored variable `X` with per-cell interval bounds `(lo, hi)` (LOD
info), and a **design matrix** of predictors (`Y`, `Z`, other analytes — mixed
types handled by a formula, so binary `Y` is just a 0/1 covariate; no copula-column
problem).

**Draw:**
1. Fit `X`'s margin by **interval-censored MLE** (sinh-arcsinh, per `leftcens`'s
   `fit_shash_margin`) to capture skew; transform observed `X` and the bounds to a
   latent normal scale.
2. `survival::survreg(Surv(latent_lo, latent_hi, type="interval2") ~ predictors,
   dist="gaussian")` — an interval-censored regression on the general design matrix.
3. Draw the latent value **truncated below the latent LOD** (`rnorm_trunc`).
4. Back-transform through the margin to the working scale.

**Proper MI:** draw the regression coefficients (and margin parameters) from their
posterior/asymptotic distribution per imputation, so between-imputation variance is
honest (as in `leftcens`'s `margin_draw` hybrid). Fall back to Gaussian/tobit when
skew is negligible.

**Home:** best exported from `leftcens` as `mice.impute.leftcens()` (reusable), with
a thin standalone wrapper for the block loop. *(Package-side TODO in `leftcens`.)*

---

## 5. Orchestration (block-FCS skeleton)

```
initialise X (bounds midpoint or marginal censored draw), Z (simple fill)
for t in 1..T_outer:
    Z  <- miceRanger(Z | current X, Y, other Z)      # few inner iters
    X  <- leftcens_censored_draw(X | current Z, Y, analytes)   # component in §4
# repeat the whole thing for m completed datasets (different seeds)
# fit brms on each completed dataset; pool with Rubin's rules (MID: drop imputed-Y rows)
```

Method routing in the pipeline config: mark the exposure to use the censored
method (a per-variable **method** field), keep `Y` as `use_in_model = TRUE`
(predictor), `impute_target = FALSE`. No new *role* flag needed — the gap was
routing `X` through a `Y`-aware method, not a missing role.

---

## 6. Estimands & number of imputations (two-m recipe)

FMI is **estimand-specific** and ~**n-invariant** (empirically confirmed in the
`leftcens` E8 study); required `m` tracks FMI, not sample size.

- **Descriptive exposure statistics** (marginal means/quantiles; `Y` irrelevant):
  `leftcens::gsimp_mi(m ≈ 30)` — cheap, no `brms`. One `leftcens` pool; use a
  prefix for the pipeline if desired.
- **ERF / `brms` coefficients:** the block-FCS pipeline. `m` is set by the **brms
  estimand's FMI**. Because the focal censored exposure is now in the loop, do
  **not** assume the MAR-only `m ≈ 8`; **pilot once at `m ≈ 20`, read the FMI from
  Rubin's pooling, and use the von Hippel two-stage rule to certify/raise `m`** —
  one expensive run, not a sweep. If `Y` truly barely depends on the censored range
  the FMI (and `m`) stay small; if the exposure is the focal predictor, expect it
  to rise.

---

## 7. Validation — simulation study (quantify the bias)

A Monte Carlo study with a **known ERF** that measures the exposure-response bias
each imputation strategy incurs. Because everything is simulated, the truth is
known and bias is exact. Test **two ERF structures**: an **additive / per-analyte**
model and a **mixture-surface** model (BKMR-style) — the latter stresses the
low-dose region and interactions, which is exactly where censored exposures live.

### 7.1 Data-generating process (complete data, no missingness)
- **`ZZ`** (n × q covariates): mixed types (continuous + categorical/binary),
  possible correlation.
- **`XX`** (n × p exposures): correlated, **right-skewed** (log-normal mixture — the
  realistic case where the censored-imputation *shape* matters), optionally
  confounded by `ZZ`.
- **`Y`** from a fixed, known ERF, in two forms:
  - **(A) Additive:** `Y = α + Σ_j β_j g_j(X_j) + γᵀ ZZ + ε`. Estimand: the
    per-exposure effects `β_j` (or exposure-response curves `g_j`).
  - **(B) Mixture surface:** `Y = α + h(XX) + γᵀ ZZ + ε`, `h` non-linear with
    **interactions** (the BKMR target). Estimands: the overall mixture effect (all
    exposures q25→q75), single-exposure effects (others held fixed), pairwise
    interactions, and `h` at representative exposure profiles.
- Fix true `β`/`h`/`γ` so the ERF is known exactly.

### 7.2 Induce missingness
- **Left-censor each `X_j`** below a per-analyte LOD; vary the ND fraction
  (e.g. {20%, 40%}). MNAR-by-censoring.
- **MCAR (or MAR-on-`ZZ`) a fraction of `ZZ`**.

### 7.3 Procedures compared (fit the matched ERF model to each completed dataset)
1. **Oracle** — complete data, no missingness. Best-case reference.
2. **Complete-case** — drop rows with any censored `X` / missing `Z`. Unbiased for
   the additive slope (selection on a predictor) but **discards the low-exposure
   region** (loses the low-dose ERF); inefficient.
3. **Independent `leftcens` pre-step** — `leftcens` imputes `XX` *without* `Y`,
   `miceRanger` imputes `ZZ`, then fit the ERF. The current "MNAR-before-pipeline"
   architecture. **Hypothesised biased.**
4. **Pipeline as-is** — current practice (e.g. substitution / MICE without proper
   censoring handling). Status-quo reference.
5. **Combined block-FCS (proposed)** — censored, `Y`-aware conditional draw for
   `XX` inside the FCS; `miceRanger` for `ZZ`; alternate, condition on `Y`.
   **Hypothesised low bias.**
6. **Joint `brms` (gold standard, subset)** — `bf(Y ~ h(mi(XX)) + ZZ) +
   bf(XX | mi() + cens(lod) ~ ZZ)`; congenial by construction. Reference on a
   tractable subsample.

### 7.4 Metrics (vs the known ERF)
- **Bias, RMSE, and CI coverage** of the ERF estimand:
  - Additive: per-`β_j` (and curve bias).
  - Mixture: overall mixture effect, single-exposure effects, **interaction terms**,
    and `h` at representative exposure percentiles.
- Report per procedure × ERF form × ND fraction. Enough reps for small MC error.

### 7.5 Design factors
- ERF form (additive vs mixture surface); ND fraction of `XX`; MCAR fraction of
  `ZZ`; exposure skewness; exposure correlation. Sample size `n`: use moderate `n`
  for speed and confirm at one large `n` — the imputation FMI is ~n-invariant
  (`leftcens` E8), so `n` is a robustness axis, not a driver of required `m`.

### 7.6 Hypotheses the study settles
- **H1** — the independent `leftcens` pre-step (no `Y`) biases the ERF; the combined
  block-FCS removes most of it.
- **H2** — the bias is **larger for the mixture surface than the additive model**,
  because interactions and the low-dose surface depend most on the censored region.
- **H3** — complete-case is unbiased for additive slopes but loses the low-dose ERF;
  the combined method recovers it (using the censored subjects' `Y`).

### 7.7 Open challenge to flag (and quantify)
**Congeniality with a *non-linear* surface.** The block-FCS censored draw conditions
`X` on the predictors *linearly* (as does `leftcens`'s copula: `X` on other analytes
linearly). That may be **uncongenial with a non-linear `h(XX)`** (interactions,
curvature), leaving residual bias in the **mixture** case even for the "combined"
method. The principled fix is the joint model with the surface *in the imputation*
(substantive-model-compatible). A key output of this study is **whether the simpler
block-FCS is "good enough" for BKMR, or whether the non-linear-congeniality gap is
material** — a genuine research question, not a settled one.

### 7.8 Also verify
- **Convergence** of the outer block-FCS (trace exposure summaries across sweeps).
- Reproduce the pattern with the `leftcens` `validation/` machinery style so it is
  auditable.

---

## 8. Caveats / open issues

- Two-engine block-FCS + MID on `Y` is bespoke; §7 validation is mandatory before
  manuscript use.
- Low-dose ERF shape is partially identified (§2.5) — report the modelled tail
  assumption explicitly.
- Binary outcomes are fine **as predictors** in the censored regression (§4); the
  copula-column problem only arose in the (rejected) "add `Y` as an analyte" hack.
- MID bookkeeping: ensure imputed-`Y` rows are deleted before the `brms` fit.

## 9. Fallback

If the block-FCS proves fragile, the `brms` **joint model** (§7.3) is the
principled alternative for the exposure case — congenial by construction — at the
cost of scaling worse than `miceRanger` on 80k × 300. Keep it as the validation
reference regardless.

---

## 10. Concrete next steps (on restart)

1. **Build the §7 simulation harness first** — known-ERF generator (additive +
   mixture surface), censoring/MCAR injection, the six procedures, and the
   bias/coverage metrics. This quantifies the problem *and* becomes the acceptance
   test for the fix (and settles §7.7 for the BKMR case).
2. Implement the §4 component (ideally as `leftcens::mice.impute.leftcens()` +
   standalone wrapper).
3. Add per-variable method routing in the pipeline config; confirm
   `use_in_model = TRUE` places `Y` in the predictor set (and MID handling).
4. Build the §5 block-FCS orchestration around `miceRanger`.
5. Run the §7 study across procedures × ERF forms × censoring; tabulate bias.
6. Pilot `m` per §6; document the FMI and the chosen `m`.
