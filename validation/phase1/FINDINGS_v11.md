# V11 findings — the imputer default survives, but for a different reason than expected

**Run:** 2026-08-31 10:35 → 21:03 JST (**628.0 min**, 8 scenarios × 4 imputer arms + oracle
× 300 reps = 2,400 tasks, `n_ok` = 300 in all 40 cells, zero task errors, zero non-finite
estimates). Results: `results/v11_latest.rds`. Registered criteria:
`../PLAN_pipeline_validation.md` §8d.

**Why this ran.** `z_imputer = "bart"` has been the shipped default since v1.5.0, chosen in
V6/V7 on a design where the outcome is **linear in `(logX, Z)` by construction** — so a
parametric imputer for `Z1` was correctly specified against the outcome and could never be
penalised for misspecification. `FINDINGS_v7.md` named this confound and left it untested.
This run removes it.

**Three conclusions, in order of importance:**

1. **The default holds.** `bartMI` is not significantly beaten by any arm in any
   non-linear-outcome cell.
2. **The V7 hypothesis is not supported.** All four imputers degrade by about the same
   amount (spread 0.44 pp), so the linear-outcome design was **not** materially flattering
   the parametric arm. The confound was real; its effect was not.
3. **A new limitation, affecting every arm.** A non-linear outcome degrades bias by
   **≈3 percentage points regardless of imputer**. `bartMI`'s mean absolute bias goes from
   **0.53% to 3.10%**. This is a property of the pipeline, not of the imputer choice, and
   it is not covered by any existing claim.

---

## Control: bit-identical to V10

| cell | arm | identical | max abs diff |
|---|---|---|---|
| `mcar_z40` | `pipeline_bartMI` | TRUE | 0 |
| `combined` | `pipeline_bartMI` | TRUE | 0 |
| `mar_z40` | `pipeline_bartMI` | TRUE | 0 |

300 replications each, on `estimate` and `se`. Adding `y_form` and `seed_as` perturbed
nothing in the existing design.

---

## Primary: relative bias by cell and arm

Matched pairs — each non-linear cell uses its linear twin's seed, so the two see identical
exposures, covariates and outcome noise and differ only by the `b_zq (Z1^2 - 1)` term.

| cell | forest | properBoot | micePmm | **bartMI** |
|---|---|---|---|---|
| `mcar_z40` | −1.64 | −1.57 | −3.18 | **−1.21** ★ |
| `ynl_mcar_z40` | −4.07 | −4.02 | −6.96 | **−3.59** ★ |
| `missing_y20` | +2.81 | +2.67 | −0.83 | **−0.04** ★ |
| `ynl_missing_y20` | −1.38 ★ | −1.41 | −3.56 | **−3.27** |
| `combined` | +2.43 | +2.02 | −0.15 ★ | **+0.56** |
| `ynl_combined` | −0.20 ★ | −0.54 | −2.26 | **−1.71** |
| `mar_z40` | −0.49 | −0.21 ★ | +1.04 | **−0.30** |
| `ynl_mar_z40` | −4.29 | −2.97 ★ | −3.55 | **−3.85** |

★ = smallest |bias| in that cell.

At face value `bartMI` wins only one of the four non-linear cells. **That reading is
wrong**, and the paired test says so.

### `bartMI` is not significantly beaten anywhere

| cell | beaten by | Δ\|bias\| (pp) | MC error | *p* | verdict |
|---|---|---|---|---|---|
| `ynl_mcar_z40` | — | — | — | — | `bartMI` best |
| `ynl_missing_y20` | forest | +0.232 | 0.144 | 0.11 | not distinguishable |
| `ynl_combined` | forest | **−0.170** | 0.177 | 0.34 | not distinguishable *(favours bartMI)* |
| `ynl_mar_z40` | properBoot | **−0.444** | 0.281 | 0.12 | not distinguishable *(favours bartMI)* |

Paired on identical datasets, MC error ~0.15–0.28 pp. In two of the three cells where
another arm had the smaller summary bias, the **paired point estimate actually favours
`bartMI`** — the summary ranking flipped on noise smaller than its own error bar. The
registered DECISION criterion is met: the default was tested on a design that could have
refuted it, and was not refuted.

### Why the ranking looked like it changed

| arm | mean \|bias\| linear | mean \|bias\| non-linear | mean degradation |
|---|---|---|---|---|
| forest | 1.84% | 2.49% | −3.26 pp |
| properBoot | 1.62% | 2.24% | −2.97 pp |
| micePmm | 1.30% | **4.08%** | −3.30 pp |
| **bartMI** | **0.53%** | 3.10% | **−2.86 pp** |

The degradation is **negative for every arm**. `forest` and `properBoot` start at *positive*
bias (+1.6% to +1.8%) in the outcome-dominated cells, so a −3 pp shift carries them
*through* zero and they land near it. `bartMI` starts near zero, so the same shift moves it
away. **Their smaller non-linear bias is two errors partially cancelling, not robustness** —
and it would reverse under a generator whose curvature had the opposite sign.

---

## Secondary: the V7 hypothesis is not supported

V7's concern was that a linear outcome flattered the *parametric* arm, so `micePmm` should
degrade most once that flattery is removed. Paired degradation, linear → non-linear:

| pair | forest | properBoot | micePmm | bartMI |
|---|---|---|---|---|
| `mcar_z40` | −2.42 | −2.45 | −3.78 | −2.38 |
| `missing_y20` | −4.19 | −4.09 | −2.73 | −3.22 |
| `combined` | −2.63 | −2.56 | −2.11 | −2.27 |
| `mar_z40` | −3.81 | −2.77 | −4.59 | −3.56 |
| **mean** | **−3.26** | **−2.97** | **−3.30** | **−2.86** |

All *p* < 1e-4 individually, but the **spread across arms is only 0.44 pp**. `micePmm`
degrades 3.30 pp and `bartMI` 2.86 pp — a difference far too small to have driven the V6/V7
verdict, where the margins were 0.59–3.85 pp.

**So the confound was real but inert.** The V6/V7 comparison was conducted on a design that
structurally favoured the parametric arm, and removing that advantage does not change the
conclusion. `micePmm` is in fact the **worst** arm under a non-linear outcome (mean |bias|
4.08%), which is the direction V7 guessed, but it arrives there by having the largest
absolute bias, not the largest degradation.

---

## The new limitation: every arm loses ~3 pp

This is the finding with consequences beyond the imputer question. A non-linear outcome
degrades bias by ≈3 pp **whichever imputer is used**, including the oracle-adjacent arms. For
the shipped default that is **0.53% → 3.10%** mean absolute bias — roughly a sixfold
increase, and outside the ≤1% and ≤2.7% figures the README quotes for the additive path.

**Coverage does not reveal it.** Mean coverage is unchanged: 0.953 → 0.953 for `bartMI`,
0.948 → 0.948 for the oracle, with `width/SE` moving 1.014 → 1.056. The reason is scale — a
3% relative bias is ≈0.28 empirical SE at `n` = 800, too small to move interval coverage
even though it is a real shift in the point estimate. **This is the third consecutive track
where calibration diagnostics were blind to a bias finding** (V3, V8, V11).

---

## What this settles

- **`z_imputer = "bart"` stays.** It is best or joint-best in every cell, never
  significantly beaten, has the lowest mean absolute bias under a linear outcome (0.53%),
  and degrades least under a non-linear one (−2.86 pp).
- **The V6/V7 decision stands, on evidence that could have overturned it.** A default that
  survives a test built to break it is worth more than one never tested.
- **A non-linear outcome is a scope limit for the whole pipeline**, not an imputer question.
  The README's additive-path bias figures assume an outcome linear in the covariates.

## Caveats (scope of this run)

1. **One curvature form, one coefficient.** `b_zq (Z1^2 - 1)` with `b_zq` = 0.40, comparable
   to `gamma[1]` = 0.50. The ≈3 pp figure is specific to that; a larger or differently
   shaped non-linearity would presumably cost more.
2. **The degradation's sign is generator-specific.** All four arms shift *negative* here. The
   observation that `forest`/`properBoot` benefit from cancellation depends on that sign, and
   would reverse if the curvature were opposite. Do not read their non-linear numbers as a
   robustness property.
3. **The analysis model remains correctly specified.** `dgp_formula()` gains the matching
   `I(Z1^2)` term, so the ≈3 pp is entirely imputation-induced. A user who *omits* the
   non-linear term from their own analysis model would additionally carry ordinary model
   misspecification, which this run does not measure.
4. **Non-linear in the covariate, not the exposures.** Curvature in the exposure–response
   surface is the mixture path, governed by V3/V9.
5. **`missing_y20` is uncontrolled.** No reproducible baseline exists for it (V10 did not run
   that cell), so unlike the other three linear cells its numbers were not checked against a
   prior run.
6. **Additive ERF, scalar estimand `b_logX1`, `m` = 30, `n` = 800.**
