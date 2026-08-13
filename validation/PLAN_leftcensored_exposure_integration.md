# Plan: congenial integration of left-censored exposure imputation into the MI + brms pipeline

Design/plan document for `generic-mi-brms-pipeline`. Captures the decision and the
rationale so the project can restart cold. Status: **Phase 1 validation complete
(2026-08-13) — H1 confirmed, §7.7 gap shown material for mixtures; see §10 and
[`validation/phase1/FINDINGS.md`](phase1/FINDINGS.md). Phases 3–5 not yet
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

> **Emphasis note (see §9–§10).** The block-FCS below is the *production-scale*
> engine, but it is bespoke and — per §7.7 — is not guaranteed congenial with a
> non-linear mixture surface. The `brms` **joint model** (§9) *is* congenial by
> construction and needs no new package code. The recommended path therefore
> **validates the problem and the joint-model fix first (Phase 1, §10)** and only
> builds this block-FCS if the joint model cannot meet the scale actually required.
> Read §3–§5 as the design for that engine, not as an instruction to build it before
> the Phase-1 gate.

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

> **Shipped & validated (2026-08-13).** This exact recipe is now exported from
> `leftcens` (>= 0.9.0) as **`impute_censored_conditional()`** — a single-variable,
> general-predictor, interval-censored, skew-aware, proper-MI draw — plus the
> supporting margin API (`fit_shash_margin`, `draw_margin`, `x_to_z`/`z_to_x`). The
> Phase-1 harness's `cens_mi_y_shash` now *calls* that public function (no longer
> `:::` internals) and tracks the oracle on the additive DGP at both skew levels
> (bias ≤3%, coverage ~0.95), where the Gaussian version fails under skew. So §4's
> core is **done and de-risked**; the remaining package work is the thin
> `mice.impute.leftcens()` wrapper (deferred) and generalising to multiple censored
> analytes / the block loop (Phase 4). See [`validation/phase1/`](phase1/).

> **Engine clarification — `impute_censored_conditional()` vs the copula engine.**
> `leftcens`'s `gsimp_impute()` / `gsimp_mi()` / `preflight_*()` default to
> `imp_model = "copula"`, which conditions an analyte on the **other analytes only**
> — so it *cannot* include the outcome `Y`. The exposure case *requires* `Y`
> (congeniality, §2), so the pipeline and the reference use
> **`impute_censored_conditional()` with `margin = "shash"`** instead: it conditions
> on a **general** predictor set (`Y`, `Z`, analytes) while reusing the **same
> skew-robust sinh-arcsinh margin the copula engine is built on**. Net: *same skew
> handling as copula, different (`Y`-aware) conditioning; the copula engine itself is
> not on the production path.* `copula` (via `gsimp_mi`) appears in this project only
> in the **validation** harness's `leftcens_prestep` arm — the deliberately-biased
> "impute X without Y" comparison, not production.

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

> **Do NOT route 80k × 300 through the `mice` engine.** The slowness seen with
> plain `mice` is the **FCS engine refitting a model for all ~300 variables every
> iteration**, not the censored-imputation *method*. The block loop above sidesteps
> this: `miceRanger` (fast RF) handles the ~300-variable Z block, and the X block is
> a **single** call to the standalone `leftcens::impute_censored_conditional()` per
> sweep (one `survreg` fit, cheap even at n = 80k). The pipeline calls that function
> **directly** in this loop — it never enters the `mice` engine. `leftcens`'s
> `mice.impute.leftcens()` wrapper (same algorithm, mice `(y, ry, x, …)` interface,
> bounds via `blots`) is a **moderate-scale ecosystem convenience only**; using it
> *through* `mice()` on the full data reintroduces the 300-variable engine cost that
> `miceRanger` exists to avoid.
>
> Scale caveat for the X block: conditioning the censored `X` on all ~300 `Z` in one
> parametric Tobit is a single fit but its cost is ~`k²` in predictor width `k`.
> **Measured at `n` = 80k** (`run_scale_test.R`): ~1.4s/sweep at `k`=50 → **~2 min**
> for a 20×5-sweep X-block, vs 23s/sweep → **~39 min** at `k`=300. So condition on a
> **reduced predictor set** (confounders + `Y` + correlated analytes, or the
> ranger-imputed `Z`), not the full 300 — a Phase-4 design choice that keeps the
> X-block negligible next to `brms`. (`survreg` stayed numerically stable to `k`=300.)

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

> **Phase-1 result (2026-08-13) — the gap is MATERIAL for mixtures.** A congenial,
> `Y`-aware *linear* censored imputation (`cens_mi_y`, the §4-preview reference)
> tracks the oracle almost exactly on the **additive** DGP (bias ≤0.4%, coverage
> ~0.96 at 20% and 40% non-detects) — linear congenial imputation is sufficient
> there. On the **mixture surface** the *same* linear method is itself badly biased
> (**+7.7%** at 20% ND → **+19.6%** at 40% ND; coverage collapsing to **0.62**),
> because it is uncongenial with the non-linear analysis model (interactions +
> curvature). Conclusion: **the simple linear block-FCS draw is not good enough for
> BKMR-style mixtures** — the mixture case needs substantive-model-compatible
> imputation (the surface in the imputation model) or the joint-model-with-surface.
> Full table and reading: [`validation/phase1/FINDINGS.md`](phase1/FINDINGS.md).
> *Caveat:* the Phase-1 mixture estimand is a single local main-effect coefficient,
> a scaffold simplification; the definitive BKMR estimands (overall mixture effect,
> interactions, `h` at profiles) still need the estimand extension before this is
> manuscript-final — but the mechanism is demonstrated.

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

## 9. The joint model is the reference, not merely a fallback

The `brms` **joint model** (§7.3.6) — `bf(Y ~ h(mi(XX)) + ZZ) + bf(XX | mi() +
cens(lod) ~ ZZ)` — is **congenial by construction** and requires **no new package
code**. Two consequences the plan now leans on:

1. **It is the gold-standard reference for validation regardless of the production
   choice** — the yardstick every other procedure in §7 is measured against.
2. **It may be sufficient on its own.** The *only* reason to prefer the bespoke
   block-FCS (§3–§5) is scaling on 80k × 300. That justification is irrelevant at
   the moderate `n` used for validation, and §7.7 warns the block-FCS may stay
   *uncongenial* with a non-linear surface where the joint model would not. So the
   joint model is a genuine candidate for the **primary** exposure path, not a
   consolation prize.

**Decision gate (Phase 2, §10):** build the block-FCS only if Phase 1 shows the
joint model both (a) fixes the ERF bias and (b) cannot meet the scale actually
required in production. If it meets the scale, the block-FCS need not be built at
all; if it does not, Phase 1 has produced the validated reference the block-FCS
must match.

> **Phase-1 implementation note (2026-08-13).** The *naive* brms encoding
> `bf(X | mi() + cens(lod) ~ Z)` does **not** work as-is in the installed brms:
> when censored rows carry the LOD value (not `NA`), brms creates no latent
> per-observation parameters, so `mi(X)` in the outcome model silently uses the
> LOD-substituted values — it collapses to **LOD-substitution** and is biased
> (verified: it equals `lm` on LOD-filled `X`, ~+12% at 20% ND). The working
> congenial reference used in Phase 1 is instead a single-variable, `Y`-aware
> censored MI. **Two variants, and the skew sweep settled which is valid:**
> - `cens_mi_y` — a Gaussian (tobit) `survreg` of `X` on `(Y, Z, other X)`. Tracks
>   the oracle on Gaussian log-exposures but **breaks under right-skew** (bias to
>   −14%, coverage to **0.47** at 40% ND) — the tobit-under-skew failure that
>   motivated the `leftcens` copula default, now seen on the ERF coefficient.
> - **`cens_mi_y_shash`** (the default reference) — fits `X`'s margin by interval-
>   censored **sinh-arcsinh** MLE (`leftcens::fit_shash_margin`), maps `X` and the
>   LOD to a latent normal scale, runs the `Y`-aware censored regression *there*,
>   draws truncated below the latent LOD, back-transforms; proper-MI draws of both
>   margin and regression parameters. **Tracks the oracle at both skew levels**
>   (bias ≤3%, coverage ~0.95) and degrades gracefully to Gaussian. This *is* the
>   §4 component prototyped (single variable, general predictors, `Y`-aware,
>   skew-aware). A correct brms joint model remains possible but needs a careful
>   latent-variable encoding, not the naive form above.

---

## 10. Concrete next steps (on restart) — cheap-first, gated phasing

The ordering principle: **confirm the problem and the principled (joint-model) fix
before investing in the bespoke engine.** The most novel, highest-risk work — the §4
component and the §5 block-FCS — comes *after* a decision gate, not before. Each
phase is gated by the phase before it.

**Phase 1 — Confirm the problem and the congenial fix (no new package code). ✅ DONE
(2026-08-13).** Built the §7 harness in [`validation/phase1/`](phase1/) (known-ERF
generator: additive + mixture; censoring injection; bias/coverage metrics) with the
procedures that need no new code: (1) oracle, (2) complete-case, (3) independent
`leftcens` pre-step via `gsimp_mi`, (4) congenial `Y`-aware censored MI reference
(`cens_mi_y` — see §9 note; the naive brms joint encoding did not work and was
replaced). **Results (300 reps/cell):**
- **H1 confirmed (additive):** the no-`Y` pre-step attenuates the focal coefficient,
  −5.6% at 20% ND → **−14.4%** at 40% ND, coverage 0.92 → 0.82; `cens_mi_y` removes
  it (bias ≤0.4%, coverage ~0.96) and is more efficient than complete-case.
- **H3 confirmed:** complete-case is ~unbiased for the slope but inefficient (rmse up
  to ~4.5× oracle at 40% ND on the mixture surface).
- **§7.7 settled — gap material for mixtures:** the linear congenial method is itself
  biased on the mixture surface (+7.7% → **+19.6%**, coverage → 0.62).
- **§7.5 skew settled — the reference must be skew-aware:** a Gaussian/tobit `Y`-aware
  reference tracks the oracle on Gaussian log-exposures but **breaks under right-skew**
  (bias to −14%, coverage to **0.47**); a **sinh-arcsinh (shash) margin** version
  (`cens_mi_y_shash`, the §4 recipe) restores it (bias ≤3%, coverage ~0.95) at both skew
  levels. Independent confirmation of the `leftcens` copula switch on the ERF coefficient.
  The §4 component must therefore be **both** skew-aware and `Y`-aware — neither existing
  tool is both.
- Full table + reading: [`validation/phase1/FINDINGS.md`](phase1/FINDINGS.md).

**Phase 2 — Decision gate (read-out from Phase 1).**
- **Additive / per-analyte ERFs:** a `Y`-aware *linear* congenial imputation is
  sufficient (it tracks the oracle). Building the §4 linear component is justified and
  will work; no surface-in-imputation needed here.
- **Mixture / BKMR ERFs:** linear congenial imputation is **not** sufficient (§7.7 gap
  is material). The mixture path needs substantive-model-compatible imputation (the
  surface in the imputation model) or the joint-model-with-surface — the simple linear
  block-FCS draw should not be assumed adequate.
- **Scale question — RESOLVED (2026-08-13):** at `n` = 80k the X-block draw
  (`impute_censored_conditional`) is cheap and super-linear in predictor width `k`
  (~`k²`), not in `n`: ~0.9s/sweep at `k`=25, 1.4s at `k`=50, 3.2s at `k`=100,
  rising to 23s at `k`=300. Extrapolated to 20 datasets × 5 outer sweeps the whole
  X-block is **~1.6–2.4 min at `k` ≤ 50** (a rounding error next to `miceRanger` and
  `brms`), vs ~39 min if conditioning on all ~300 `Z`. `survreg` converged cleanly at
  every width. So a congenial approach **meets production scale** with a **reduced
  predictor set** — no speed reason to abandon it; the block-FCS architecture is
  viable. Timings: [`validation/phase1/run_scale_test.R`](phase1/run_scale_test.R).
- **Reference caveat:** the definitive mixture verdict needs the richer BKMR estimands
  (overall mixture effect, interactions) — a Phase-1 estimand extension — before it is
  manuscript-final; the mechanism, however, is demonstrated.

**Phase 2b — Estimand extension (recommended before the mixture verdict is final).**
Add the true BKMR estimands to the harness (overall mixture effect q25→q75, pairwise
interactions, `h` at profiles) and a mixture-appropriate reference, so the §7.7
conclusion rests on the estimands that actually matter, not the scaffold's single
local coefficient. Cheap relative to Phases 3–4.

**Phase 3 — Build the §4 component in `leftcens`** *(gate says build for the additive/
per-analyte case)*. **Core DONE (leftcens 0.9.0, 2026-08-13):** shipped the
single-variable, general-predictor, interval-censored, skew-aware, proper-MI draw as
**`leftcens::impute_censored_conditional()`**, plus the exported margin API
(`fit_shash_margin`, `draw_margin`, `x_to_z`/`z_to_x`); tests, docs, `R CMD check`
clean; the Phase-1 harness's `cens_mi_y_shash` now consumes it and still tracks the
oracle at both skew levels. **Remaining Phase-3 work:**
- The thin `leftcens::mice.impute.leftcens()` wrapper (moderate-scale ecosystem
  convenience; bounds via mice `blots`). *Deferred by scope choice.*
- Generalise the core to **multiple censored analytes** and the block-loop calling
  convention (for §5).
- **For mixtures, scope a substantive-model-compatible variant** (surface in the
  imputation model) — Phase 1 shows the plain linear draw is insufficient there.

> **Margin-API export task — DONE (leftcens 0.9.0).** The four pieces the component
> needs — `fit_shash_margin()` (interval-censored sinh-arcsinh MLE), `draw_margin()`
> (proper-MI parameter draw), and the `x_to_z()` / `z_to_x()` latent-scale transforms
> — were promoted from `:::` internals to **exported, documented, tested API** (plain
> functions, per the design decision). Downstream code (this harness included) no
> longer depends on unexported symbols.

**Phase 4 — Wire block-FCS into the pipeline.** **Core BUILT (2026-08-13):**
- New **`00_censored_exposure.R`** — `run_censored_exposure_block_fcs()`: the §5
  outer loop (miceRanger Z-block via the *reused* `run_row_level_imputation()`
  helper + `leftcens::impute_censored_conditional()` X-block, Y-aware), log-scale
  handling from the dictionary `scale`, and **MID** (Y imputed so it is complete for
  the X-block predictor step, then imputed-Y rows deleted). Returns the standard
  `imputed_list`, so `03_impute.R`'s existing file/manifest logic is reused untouched.
- **Method routing:** new opt-in `strategy = "censored_exposure_block_fcs"` +
  `imputation$censored_exposure` config block; input convention = `leftcens` interval
  columns `X_lo`/`X_hi`; exposures flagged `impute_target = FALSE, use_in_model = TRUE`.
- **Dispatch:** one `else if` branch in `03_impute.R` (main stream otherwise intact).
- **Validated (synthetic):** runs end-to-end, produces clean completed datasets,
  drops imputed-Y rows, and **recovers the exposure-response coefficient** (bias
  −1.7%, coverage 0.92 vs oracle +0.8%, true slope 0.80).

**End-to-end acceptance test — PASSED (2026-08-13).** Built a censored demo (n=500,
35% non-detects, MAR covariate, missing outcomes) with `expo`/`expo_lo`/`expo_hi`
columns and ran **Steps 1→4** in an isolated dir with `strategy =
"censored_exposure_block_fcs"`, `m=6`, `custom_formula = y ~ log(expo) + z1 + z2`:
- Step 1 validation accepted `log(expo)`; Steps 1–2 carried the `_lo`/`_hi` columns
  through (`prepare_raw_data` keeps all columns).
- Step 3 produced 6 clean completed datasets (470 rows after MID dropped the 30
  imputed-Y rows; `expo` complete/positive; `_lo`/`_hi` dropped).
- Step 4 (`brms`, cmdstanr) consumed them and fit all 6 with **zero errors**; the
  pooled `log(expo)` coefficient (0.66, 95% CrI [0.54, 0.77]) covers the
  dataset's oracle (0.71). Rigorous unbiasedness is from the 25-rep MC (−1.7%,
  coverage 0.92), not this single realization.
- Two bugs fixed while wiring: (a) the per-sweep NA-reset must **exclude the
  exposures** (else it blanks what the X-block imputes); (b) imputation log-scale
  is now a `censored_exposure$log_scale` config flag, decoupled from the dictionary
  `scale` to avoid double-transform with a `log()` in the model formula.

**Remaining Phase-4 work:**
- Add procedures 5–6 of §7 to the harness (the pipeline block-FCS itself) and re-run.
- Multiple censored exposures / interval (DNQ) cases beyond the single-exposure path.
- Assert the **MID invariant** in a test, not just in code.

**Phase 5 — Pilot `m`** per §6 on the chosen path; document the FMI and the final `m`.

*Note vs. earlier draft:* the previous §10 built the §4 component and block-FCS
(steps 2 & 4) before knowing whether the simpler joint model already suffices. The
phasing above turns that into an explicit gate so the bespoke, higher-risk engine is
built only when Phase 1 proves it is needed.
