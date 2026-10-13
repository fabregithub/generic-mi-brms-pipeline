# V29 — `s²` is the driver, and the old refutation was a level artefact

**Question.** V28 left one live hypothesis: the chain's `derived/measured` ratio trends with the
arrow in the **fork** sweep (0.944 → 1.112) but is flat in the **pipe** sweep (~0.95). The arrow
moves the exposure conditional's `Z1` coefficient in both roles but moves
`Var(logX1) = sd_x² + delta_zx²` only in the fork — so the flat pipe acquits the coefficient and
implicates `Var(logX1)`. Does `s²` drive it, on its own?

**Outcome: P1 — yes, in BOTH roles, and a single curve in `s²` describes them.**

## 1. Design

`sd_x` moves `s²` with the arrow held at 0.60, and it does so in *both* roles — which the arrow
cannot do. Six cells, `n = 12800`, 200 reps, arms paired on byte-identical data
(`pipeline_properZ` vs `pipeline_properZ_ds`, differing only by the pre-v1.6.0 scale defect).
Registered before the run in `run_v29_sdx.sh` and `R/robustness.R`:

> **P1** ratio rises with `sd_x` in both roles → `s²` is the driver.
> **P2** fork only → something fork-specific survives.
> **P3** flat in both → `s²` is not the driver.

The pipe's `s²` range (0.49 → 1.69) is **wider** than the fork arrow sweep's (1.16 → 1.64), so a
flat pipe could not have been dismissed as too small a lever.

## 2. Result

| cell | role | `sd_x` | `s²` | measured | derived | ratio |
|---|---|---|---|---|---|---|
| `sxf07` | fork | 0.7 | 0.850 | +8.45 ± 0.09 | +7.81 | **0.924 ± 0.019** |
| `sxf10` | fork | 1.0 | 1.356 | +5.03 ± 0.04 | +5.02 | **0.998 ± 0.017** |
| `sxf13` | fork | 1.3 | 2.041 | +2.54 ± 0.03 | +2.74 | **1.077 ± 0.023** |
| `sxp07` | pipe | 0.7 | 0.489 | +3.54 ± 0.05 | +2.51 | **0.709 ± 0.053** |
| `sxp10` | pipe | 1.0 | 0.999 | +5.97 ± 0.05 | +5.58 | **0.935 ± 0.028** |
| `sxp13` | pipe | 1.3 | 1.688 | +6.60 ± 0.05 | +6.98 | **1.057 ± 0.024** |

Both roles rise: fork **+0.128 ± 0.010 per unit `s²`** (t = 13.0), pipe **+0.242 ± 0.074**
(t = 3.3). **P1.**

`sxf10` reproduces the canonical `zr_fork` geometry (`s²` = 1.356) as a within-run anchor and
reads 0.998 — consistent with V28's independent `fkd06` at 0.992.

## 3. The strong form: one curve, no role

Fitting `ratio ~ s²` on all six cells (inverse-variance weighted) and asking whether the role
adds anything:

| model | added term | p |
|---|---|---|
| `ratio ~ s²` | — | slope **+0.161 ± 0.033**, p = 0.008 |
| `+ role` | role offset | **0.79** |
| `+ s² × role` | role-specific slope | **0.21** |

**The role adds nothing.** Five of six cells sit within **0.023** of the pooled line; the
exception is `sxp07` (−0.142), the lowest-`s²` cell and the only one outside the fork's swept
range, so it cannot be told apart from curvature in `s²` at the low end. The line crosses 1 at
**`s²` ≈ 1.42**.

**This explains V28 completely.** The fork's arrow moves `s²` (`= sd_x² + delta_zx²`) and the
pipe's does not, so a sweep in the arrow shows a trend in the fork only — exactly what V28 saw.
There is no fork-specific mechanism. There never was a "fork residual"; there is an `s²`
dependence, seen through a parameterisation that confounded the two in one role.

## 4. Why the 2026-09-11 refutation read the way it did

`THEORY.md` recorded the `exp`-linearisation hypothesis as **refuted**: it predicts agreement
*degrading* as `s²` grows, and the two-block sweep over `s²` = 0.52 → 2.05 gave
**0.841 → 0.879 → 0.900 → 0.908** — read as moving steadily *toward* 1.

Those numbers rise with `s²`. **So does V29's.** The two sweeps agree on direction and disagree
only on **level** — and the level is what decides whether a monotone rise reads as "converging to
1" or as "crossing 1". The two-block instrument sits entirely below 1 because it misses the
pipeline's own level by 31%, a limitation `THEORY.md` states in the caution directly beneath the
refutation. Viewed from a baseline shifted down, crossing looks like converging.

**The refutation was an artefact of the instrument's level offset, not a finding.** The caution
attached to it was the right one, and it was load-bearing: *"can refute a hypothesis about the
closed form's shape but cannot adjudicate the pipeline's level."* The level was the whole
question.

## 5. What this does NOT establish

**`s²` governs the chain's error; the `exp` linearisation is not shown to be the mechanism.** A
linearisation error must vanish as the thing being linearised stops mattering — `ratio → 1` as
`s² → 0`. V29's ratio does the opposite: it falls to **0.709** at `s²` = 0.489 and crosses 1 at
`s²` ≈ 1.42. A pure linearisation error cannot produce a ratio that is *below* 1 at small `s²`.

So the honest statement is: **the chain's error is a monotone function of `s²` that passes
through zero near `s²` ≈ 1.4, and whatever produces it is not a single linearisation term.** At
least two `s²`-dependent effects of opposite sign are in play, or one link is mis-scaled in a way
that depends on `s²` — the truncation factor `Var(L|L≤c)/σ_x²` is the obvious suspect, since the
LOD is a quantile of `L` and therefore moves with `s²` by construction.

The practical reading is unchanged and now better supported: **the chain is accurate near
`s²` ≈ 1.4 and degrades in both directions away from it**, by about 0.16 in ratio per unit `s²`.

## 6. Reproduce

```bash
bash run_v29_sdx.sh        # 6 cells x 200 reps, n = 12800; 51 min measured
bash run_status.sh v29
Rscript analyze_v29.R
```

Cells `sxf07/10/13`, `sxp07/10/13` in `R/robustness.R` (axis `"sdx"`); raw results in
`results/v29_n12800_raw.csv`. The `sd_x` scenario field was added in V29 and defaults to 1.0, so
no pre-V29 cell moves.
