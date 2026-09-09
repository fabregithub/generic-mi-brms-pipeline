# V23 findings — the derived law mostly holds, its exponent test has no power, and item 07 is confirmed catastrophic

**Run:** 2026-09-09 07:56 → 13:00 JST (**303.7 min** against 4.8 h predicted; 2 `n` levels ×
7 cells × 4 arms + oracle × 500 reps = 7,000 tasks, `n_ok` = 500 in all 14 cells, zero task
errors, zero non-finite estimates). Results: `results/v23_n{800,3200}_latest.rds`. Verdict
arithmetic: [`analyze_v23.R`](analyze_v23.R). Registered in
[`run_v23_curvature.sh`](run_v23_curvature.sh)'s header, committed before the run.

**The first track designed from a derivation rather than a surprise.** V22 found ~20%
asymptotic bias in the censored-exposure draw by accident, in a secondary cell. Before
building anything this time the mechanism was derived — `leftcens` fits a conditional linear
in its predictors, so what governs the damage is `u`, the density-weighted residual of the
covariate arrow `g` after its best linear fit below the LOD — with a reason for `u` to enter
**squared** (it is an L²-projection residual, so the first-order term in the bias expansion
vanishes). Registered law: **bias% = 300·u²**, calibrated on one point.

**Three results.** The law's *level* predictions are good — four of six within 1.5 pp — but it
**over-predicts at large `u`** and misses the biggest cell by 8.8 pp. The **exponent test has
no power**: 1.64 with a CI containing both 1 and 2, so it neither confirms nor refutes the
orthogonality argument. And **item 07's grid exposure draw is confirmed catastrophic**, at up
to **+142%** with coverage 0.000.

---

## The law

`z21_exact` isolates the X block (exact covariate draw, shipped `leftcens` exposure draw).
`n` = 800. Reported in `bias/SE` per PLAN §11.

| cell | `u` | predicted | measured | miss | `bias/SE` | direction |
|---|---|---|---|---|---|---|
| `cv000` | 0.0000 | 0.00% | **−0.80%** | −0.80 | 0.07 | — *(the null)* |
| `cv040` | 0.0399 | 0.48% | **+0.89%** | +0.41 | 0.07 | away |
| `cv119` | 0.1190 | 4.25% | **+4.70%** | +0.46 | 0.34 | away |
| `cv153` | 0.1531 | 7.04% | **+3.37%** | −3.66 | 0.29 | away |
| `cv182` | 0.1818 | 9.91% | **+11.42%** | +1.50 | 0.79 | away |
| `cv262` | 0.2621 | 20.60% | +20.31% | −0.29 | 1.19 | away *(calibration)* |
| `cv364` | 0.3636 | 39.66% | **+30.82%** | **−8.84** | 1.48 | away |

**The null holds** — `cv000` sets `g ≡ 0`, and with nothing for a linear draw to miss the bias
is −0.80% (bias/SE 0.07). The oracle is unbiased in every cell (worst 0.63%).

**Four of six predictions land within 1.5 pp**, from a *single* calibration point. For a
one-parameter law derived rather than fitted, that is a strong showing.

**But the registered ±5 pp gate fails at `cv364`, by 8.8 pp, and the sign of the miss is
systematic.** The two largest-`u` cells are both over-predicted, and a freely fitted constant
with the exponent held at 2 gives **k = 249** rather than 300. **The law saturates at large
`u`** — a one-parameter `u²` form is not enough, exactly as this run was registered to
determine.

**The bias is asymptotic, confirmed cleanly.** Between `n` = 800 and 3,200 every non-null cell
shifts by **≤0.27 pp** (`cv262`: +20.31% → +20.37%). That is the derivation's central
structural claim and it is not in doubt.

**The direction is stable**: away from the null in every cell with `u` > 0.05. Anti-conservative
throughout, but at least predictable — unlike V15's curvature estimate, which inverts.

---

## The exponent test had no power, and the gate hid that

Registered: reject the quadratic form if the log–log slope's 95% CI excludes 2.

> slope **+1.638**, 95% CI **[+0.994, +2.282]**

The CI contains 2, so the gate passed. **It also contains 1, 1.5, and everything between** —
width 1.29 on six points. The gate could only ever have failed if the truth were far from
quadratic; it cannot distinguish quadratic from linear, which is the comparison that matters.

**So the orthogonality argument is not refuted and not confirmed.** The honest position is
that it survives a test too weak to have threatened it — and the point estimate of 1.64 leans
*below* 2, which is what the saturation at large `u` would produce.

**This is the second consecutive track where a registered gate resolved for the wrong
reason** — V22's ≥3% bar on BART missed while the conclusion it licensed was contradicted by
the neighbouring numbers. Both were single-threshold gates on a quantity whose informativeness
depends on precision nobody checked in advance. **The convention that follows: before
registering a CI-based gate, state what CI width would make it discriminating, and check the
design can deliver it.** Six points at 500 reps could not.

Distinguishing 1.64 from 2 needs more `u` levels at the low end, where the bias is small and MC
error dominates — a materially more expensive run than this one, and worth doing only if the
orthogonality argument is load-bearing for something else.

---

## The practical threshold, in the new currency

This is the first track graded in `bias/SE` rather than relative bias, and it makes the
guidance sharper than the percentages do:

| `u` | `bias/SE` | verdict against the 0.3 gate |
|---|---|---|
| ≤ 0.04 | 0.07 | comfortable |
| 0.119 | 0.34 | **at the gate** |
| 0.153 | 0.29 | at the gate |
| 0.182 | 0.79 | fails — half a detectable effect |
| 0.262 | 1.19 | fails |
| 0.364 | 1.48 | fails badly |

**The threshold sits at `u` ≈ 0.12.** Below it the defect is tolerable; above it, it consumes
a quarter of a detectable effect and rises fast. Three of seven cells fail the new default
gate — expected, since they were built to be severe, but it is the number a user needs.

Combined with the censoring-rate calculation (`u` grows ~51× between 5% and 70% non-detects
for a fixed arrow): **at 5–10% non-detects almost no realistic curvature reaches `u` = 0.12;
at 40%+ a moderately curved covariate relationship does.**

---

## Item 07's grid exposure draw: confirmed catastrophic, and it scales with `u`

The 22-rep probe suggested an inversion. At 500 reps it is not marginal:

| cell | `z21_exact` | `pipeline_bartMI` | `pipeline_micePmm` | **`smc_xgrid`** |
|---|---|---|---|---|
| `cv000` | −0.80% | −2.03% | −3.40% | **−0.92%** |
| `cv040` | +0.89% | +3.44% | +4.08% | **+16.80%** |
| `cv119` | +4.70% | +9.25% | +8.10% | **+43.83%** |
| `cv153` | +3.37% | +3.34% | −2.67% | **+8.57%** |
| `cv182` | +11.42% | +16.13% | +10.52% | **+52.58%** |
| `cv262` | +20.31% | +28.84% | +18.05% | **+94.06%** |
| `cv364` | +30.82% | +46.22% | +29.43% | **+142.06%** |

**`smc_xgrid` is fine when `g ≡ 0` (−0.92%) and degrades monotonically with `u`**, reaching
three times the shipped path's bias with coverage 0.000. The mechanism is the one predicted:
it proposes from the exposure prior and reweights by `p(Y | x, rest)`, so it never uses the
`Z1 = g(logX1)` information at all — discarding the linear-but-partly-right covariate
conditioning `leftcens` does use.

**Item 07 has been the only open build track for weeks. It must not ship without this cell
class in its acceptance set.** It remains the right fix for the mixture failure (V3/V9) and
the non-linear-outcome penalty (V13/V14), both of which are about the *outcome* likelihood —
but a covariate related non-linearly to the censored exposure is a case it makes several times
worse.

**`micePmm` tracks the exact-Z arm closely** (+29.43% against +30.82% at `cv364`), while
`bartMI` runs ~15 pp above both. Consistent with V21/V22: this is an X-block defect, so the
covariate imputer barely moves it, and BART adds its own separate covariate-draw bias on top.

---

## What this changes

| claim | status |
|---|---|
| The defect is derivable from what `leftcens` conditions on | **Confirmed** — the null cell is zero, the level predictions land, the bias is asymptotic |
| `u` is the governing quantity | **Confirmed** — four of six predictions within 1.5 pp from one calibration point |
| bias = 300·u², one parameter | **Insufficient** — over-predicts at large `u`; free constant 249; `cv364` misses by 8.8 pp |
| The orthogonality argument (exponent exactly 2) | **Untested** — CI [0.99, 2.28] contains 1 and 2 alike |
| The bias is asymptotic | **Confirmed** — ≤0.27 pp shift over a 4× range of `n` |
| Direction is stable and away from the null | **Confirmed** |
| Item 07's grid draw inverts under a non-linear covariate arrow | **Confirmed at 500 reps**, up to +142%, scaling with `u` |

---

## Caveats (scope of this run)

- **One arrow family** (`a·tanh(1.8x) + c·(x²−1)`), one estimand, one covariate. `u` is
  computed for a Gaussian exposure marginal; a skewed exposure would move the density weights
  and hence `u` itself.
- **Two `n` levels.** Enough for flatness, not for an `n`-exponent.
- **The exponent is unresolved and the design could not have resolved it.** See above.
- **`u` is derived under the assumption that the exposure conditional's *only* non-linearity
  is `g`.** With several curved covariates the residual would have to be taken jointly, which
  this does not test.
- **The `cv153` miss (−3.66 pp) is inside the gate but not inside MC error**, and it is the
  quadratic-only arrow — the one whose `g` is monotone over the censored region. That may be a
  second regime rather than noise; the sweep is too coarse to say.
