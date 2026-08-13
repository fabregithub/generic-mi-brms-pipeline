# Phase 1 findings — left-censored exposure, focal-predictor bias

**Run:** 2026-08-13. Full grid, 300 reps/cell, `M`=30 imputations, `n`=800, parallel
(`NCORES`=14). Procedures: `oracle`, `complete_case`, `leftcens_prestep` (no-`Y`),
`cens_mi_y` (congenial, `Y`-aware). Estimand: `b_logX1`, the focal censored-exposure
coefficient (true value **0.40**). Reproduce with `./run_phase1.sh`; regenerate this
table from `results/latest.rds`.

## Results (all cells n_ok = 300)

| ERF | ND | procedure | est | rel. bias | rmse | coverage |
|---|---|---|---|---|---|---|
| additive | 0.2 | oracle           | 0.398 | −0.4% | 0.039 | 0.95 |
| additive | 0.2 | complete_case    | 0.398 | −0.5% | 0.060 | 0.94 |
| additive | 0.2 | **leftcens_prestep** | 0.377 | **−5.6%** | 0.046 | 0.92 |
| additive | 0.2 | **cens_mi_y**    | 0.399 | −0.3% | 0.041 | 0.96 |
| additive | 0.4 | oracle           | 0.400 | −0.0% | 0.038 | 0.95 |
| additive | 0.4 | complete_case    | 0.399 | −0.4% | 0.072 | 0.97 |
| additive | 0.4 | **leftcens_prestep** | 0.342 | **−14.4%** | 0.072 | 0.82 |
| additive | 0.4 | **cens_mi_y**    | 0.402 | +0.4% | 0.043 | 0.96 |
| mixture  | 0.2 | oracle           | 0.401 | +0.1% | 0.041 | 0.95 |
| mixture  | 0.2 | complete_case    | 0.401 | +0.3% | 0.083 | 0.97 |
| mixture  | 0.2 | leftcens_prestep | 0.415 | +3.8% | 0.044 | 0.95 |
| mixture  | 0.2 | **cens_mi_y**    | 0.431 | **+7.7%** | 0.053 | 0.90 |
| mixture  | 0.4 | oracle           | 0.400 | +0.0% | 0.040 | 0.95 |
| mixture  | 0.4 | complete_case    | 0.402 | +0.4% | **0.178** | 0.94 |
| mixture  | 0.4 | leftcens_prestep | 0.416 | +4.0% | 0.042 | 0.97 |
| mixture  | 0.4 | **cens_mi_y**    | 0.479 | **+19.6%** | 0.090 | **0.62** |

## Reading against the hypotheses

**H1 — the no-`Y` pre-step biases the ERF; the congenial method removes it. CONFIRMED
(additive).** `leftcens_prestep` attenuates the focal coefficient toward zero, worsening
with censoring (−5.6% → −14.4%) as coverage falls (0.92 → 0.82). `cens_mi_y` (congenial,
includes `Y`) sits on the oracle at both censoring levels (bias ≤0.4%, coverage ~0.96)
and is more efficient than complete-case (rmse 0.043 vs 0.072 at 40% ND). This is the
plan's central claim, demonstrated end-to-end.

**H3 — complete-case unbiased for the slope but inefficient. CONFIRMED.** Complete-case
bias is negligible everywhere, but its rmse balloons as censoring rises — up to 0.178 on
the mixture surface at 40% ND (~4.5× the oracle). The congenial method recovers that
efficiency.

**H2 / §7.7 — the non-linear-congeniality gap is MATERIAL for mixtures. NEW, decisive.**
On the mixture surface the *linear* congenial imputation (`cens_mi_y`) is itself badly
biased (+7.7% → +19.6%, coverage collapsing to 0.62), because a linear imputation model
is uncongenial with a non-linear analysis model (interactions + curvature). So §7.7's
open question — "is the simpler linear approach good enough for BKMR?" — resolves toward
**no** for mixtures.

## What it means for the plan (Phase-2 gate)

- **Additive / per-analyte ERFs:** a `Y`-aware *linear* congenial imputation is
  sufficient. The §4 component as specified will work; `cens_mi_y` is that component in
  miniature and is the reference implementation.
- **Mixture / BKMR ERFs:** linear congenial imputation leaves material bias → the mixture
  path needs substantive-model-compatible imputation (surface in the imputation model) or
  the joint-model-with-surface. Do not assume the simple linear block-FCS draw suffices.

## Skew robustness — the reference must be skew-aware (§7.5)

The main run above used `skew = 0` (Gaussian log-exposures), where a Gaussian/tobit
conditional is correctly specified. A follow-up sweep (`run_skew_sweep.R`, additive,
150 reps/cell) added right-skew and exposed that **the choice of margin is first-order**:

| skew | ND | `leftcens_prestep` (copula, no Y) | `cens_mi_y` (Gaussian+Y) | `cens_mi_y_shash` (shash+Y) |
|---|---|---|---|---|
| 0.00 | 0.2 | −4.3%, cov 0.92 | +1.1%, cov 0.95 | +1.0%, cov 0.97 |
| 0.00 | 0.4 | −15.7%, cov 0.81 | −0.6%, cov 0.95 | −1.1%, cov 0.99 |
| 0.75 | 0.2 | −1.2%, cov 0.93 | **−8.6%, cov 0.77** | **+1.2%, cov 0.93** |
| 0.75 | 0.4 | −5.7%, cov 0.97 | **−14.2%, cov 0.47** | **+2.7%, cov 0.95** |

*(oracle unbiased with ~0.95 coverage throughout; complete-case ~unbiased.)*

- **The Gaussian/tobit reference `cens_mi_y` breaks under right-skew** — from tracking the
  oracle at skew 0 to the worst arm at skew 0.75 (bias to −14%, coverage to **0.47**). Its
  Gaussian conditional draws non-detects from a too-heavy left tail — the exact tobit
  failure that motivated `leftcens`'s copula default. **A tobit reference is only a valid
  gold standard when the margin is Gaussian.**
- **The skew-aware `cens_mi_y_shash` holds at both skew levels** (bias ≤2.7%, coverage
  ~0.95) and degrades gracefully to Gaussian at skew 0. It fits X1's margin with an
  interval-censored **sinh-arcsinh** MLE (`leftcens::fit_shash_margin`), maps to a latent
  normal scale, runs the `Y`-aware censored regression there, and back-transforms — the
  §4 recipe in miniature. This is now the **default reference** in the harness.
- **The copula pre-step is margin-robust but still omits `Y`.** Independent confirmation of
  the `leftcens` copula switch on a new estimand (the ERF coefficient, not marginal
  recovery). Note its omitted-`Y` attenuation is smaller under this skewed DGP.

**Design lesson:** the §4 component must be **both** skew-aware **and** `Y`-aware. Neither
existing tool has both — `leftcens` copula is skew-aware but omits `Y`; the Gaussian
`cens_mi_y` includes `Y` but breaks under skew. `cens_mi_y_shash` combines them and is the
validated prototype.

## Method note — why not brms

The plan's named gold standard `bf(Y ~ mi(X)) + bf(X | mi() + cens(lod) ~ Z)` does not
work as-is in the installed brms: with censored rows carrying the LOD value (not `NA`),
brms creates no latent per-observation parameters, so `mi(X)` uses the LOD-substituted
values — it collapses to **LOD-substitution** and is biased (verified: brms estimate ==
`lm` on LOD-filled `X`, ~+12% at 20% ND). The working reference is `cens_mi_y`: a
left-censored Gaussian `survreg` of `X` on `(Y, Z, other X)`, proper-MI draws of the
coefficients + scale, imputation truncated below the LOD (`leftcens::rnorm_trunc`),
Rubin-pooled. `proc_brms_joint` is retained in the code but deprecated and out of the
default set.

## Caveats (scope of this run)

1. **Mixture estimand is a scaffold simplification** — a single local main-effect
   coefficient, not the true BKMR estimands (overall mixture effect, interactions, `h` at
   profiles). The mixture rows *demonstrate the §7.7 mechanism* but the definitive verdict
   needs the estimand extension (plan §10, Phase 2b).
2. **Only `X1` is censored** (`censor_which = 1L`); MCAR-on-`Z` is off. Skew has been swept
   (see above); the remaining robustness axes (censor all exposures, MCAR covariates,
   correlation, large `n`) are not yet swept.
3. **Scale untested** — moderate `n` only; the 80k × 300 production-scale question (the
   other half of the gate) is open.
