# V17 findings — the covariate's causal role decides everything, and the shipped default carries **10% bias** in a cell nobody had built

**Run:** 2026-09-03 12:28 → 13:51 JST (**82.6 min**, stages 2 and 3 of `run_v16_v17.sh`;
8 + 2 cells, 500 reps, `n_ok` = 500 in all 10 cells, zero task errors, zero non-finite
estimates). Results: `results/v17a_latest.rds`, `results/v17b_latest.rds`. Verdict
arithmetic: [`analyze_v16_v17.R`](analyze_v16_v17.R). Criteria:
`../PLAN_pipeline_validation.md` §8j.

**The question.** In every scenario this harness had ever run, `Z` was either independent of
the exposures or generated *from* them. **No track had a confounder, a mediator or a
collider.** So every covariate result on record — V4's variance attribution, V5–V7's imputer
comparison, V10's MAR cells, V12's draw shape, V13's block decomposition — was measured in
the one structure where the covariate is statistically real and causally inert.

**Three answers.** The structural prediction **held** (4 of 4 gates). The mechanism
prediction **failed** — MAR adds almost nothing, not the +4 to +6 pp derived. And the run
turned up something nobody registered a bar on: **the shipped default, `Y` in both blocks,
is biased +5% to +10.4% wherever the covariate is a confounder or a mediator**, against
≤2.3% in the structure every earlier track used.

---

## The headline: the structural 2×2, confirmed

Paired shift of `pipeline_noYz` against `pipeline_bartMI`, in percentage points of
`b₁ = 0.40`. MC error ≈ 0.22 pp.

| role | `Z` adjusted for? | `Z` on an X–Y path? | MCAR | MAR | predicted (MCAR) |
|---|---|---|---|---|---|
| precision *(from V16)* | yes | **no** | **+1.17** | +0.93 | ≈ 0 |
| **fork** (confounder) | yes | yes | **+25.47** | +26.18 | +26.7 |
| **pipe** (mediator) | yes | yes | **+28.17** | +28.90 | +30.4 |
| **collider** | **no** | yes | **−0.23** | +2.08 | −0.7 |
| **mixed** (fork + pipe) | yes | yes | **+22.58** | +22.01 | +19.5 |

**All four structural gates passed**: `|noYz|` > 10 pp in fork, pipe and mixed under both
mechanisms, and < 3 pp in precision under both.

**The Z block's `Y`-omission penalty is large exactly when the covariate is adjusted for
*and* lies on an open X–Y path** — confounding or mediating alike. It is ≈ zero when the
covariate is adjusted for but off any path, and ≈ zero when it is on a path the analysis
correctly omits.

**This is a factor of ~22 in the same arm.** `noYz` is +1.2 pp in the structure V12, V13 and
V16 all measured, and +25 to +28 pp one arrow away. Nothing about the imputation changed —
only what the covariate *is*.

**The pipe was the discriminating case and it is why the earlier qualifier was wrong.** An
earlier version of §8j said the condition binds a covariate's draw "through that covariate's
association with the exposure". A mediator has no confounding path, yet produces the largest
penalty of any role. **Association is not the criterion; adjusted-for-and-on-a-path is.**
The pipe's +28.2 pp also sits inside its pre-registered upper bound of +75 pp — the distance
from the direct effect (0.40) to the total (0.70) — which is the check that the mechanism is
the one claimed and not something else.

---

## The mechanism prediction failed: MAR adds almost nothing

Registered: MAR-on-`Y` should add **+4 to +6 pp** in every on-path cell, because omitting `Y`
there drops the variable the missingness mechanism depends on — a validity failure on top of
the congeniality one.

| role | MAR − MCAR, measured | predicted |
|---|---|---|
| fork | **+0.70 ± 0.29** | +5.7 |
| pipe | **+0.73 ± 0.31** | +5.7 |
| mixed | **−0.58 ± 0.31** | +4.4 |
| collider | **+2.31 ± 0.12** | +3.8 |

**Three of four are within ±1 pp of zero, and `mixed` has the wrong sign.** The registered
"positive in all four on-path cells" gate is missed, and `zrmar_fork` (+26.2 against +32.4)
and `zrmar_pipe` (+28.9 against +36.1) missed their magnitude bars as a direct consequence —
they were predicted as MCAR plus a mechanism term that is not there.

**The likely reason, offered as a conjecture and not a result.** The Z block conditions on
the exposures and on the other covariate, which are themselves associated with `Y`. A
flexible draw (BART) can recover most of the `Y`-relevant information from those proxies
without `Y` itself; the linear stand-in in `predict_roles.R` cannot, so it over-attributed
the loss to `Y`'s absence. That is testable — a `mice pmm` Z block should show a larger MAR
term than BART if this is right.

**What survives is the stronger version of the claim.** "The path structure dominates the
missingness mechanism" was registered qualitatively and is confirmed more sharply than
derived: the mechanism term is not merely smaller than the structural one, it is **within
noise of zero** except in the collider. **The DAG decides; MCAR-versus-MAR barely
registers.**

**The collider is the one cell where the mechanism matters** (+2.31 ± 0.12, the only clean
positive). `Z₁` there is largely a function of `Y`, so a `Y`-less draw is badly wrong, and
although the analysis omits `Z₁`, the X block still conditions on it. Its **leak is the
largest of any role** — −6.20 and −6.42 ± 0.12 — so the collider is also the strongest single
case for these blocks not being separable.

---

## Not predicted, and the most useful thing in the run

The gates above are all **paired** contrasts, in which the reference arm's own bias cancels.
That is what makes them clean, and it is also what hides this:

| cell | `oracle` bias | **`pipeline_bartMI` bias** | coverage | width/SE | FMI | bias / emp. SE |
|---|---|---|---|---|---|---|
| `zr_collider` | +0.7% | +0.2% | 0.954 | 0.98 | 0.33 | 0.02 |
| `zr_mixed` | +1.1% | +0.6% | 0.954 | 1.06 | 0.45 | 0.04 |
| **`zr_fork`** | −0.1% | **+4.96%** | 0.944 | 1.02 | 0.44 | **0.37** |
| **`zr_pipe`** | +0.2% | **+5.39%** | 0.942 | 1.03 | 0.41 | **0.35** |
| `zrmar_mixed` | +1.1% | +4.39% | 0.968 | 1.13 | 0.54 | 0.30 |
| **`zrmar_fork`** | −0.1% | **+8.52%** | 0.932 | 1.07 | 0.52 | **0.62** |
| **`zrmar_pipe`** | +0.2% | **+10.40%** | 0.938 | 1.12 | 0.51 | **0.67** |

`pipeline_bartMI` is the **shipped configuration** — `Y` in both blocks, BART Z block, the
recommended setup. The oracle is unbiased in every role, so the DGP and the analysis model
are correctly specified and this is the imputation's bias, not the design's.

**Against ≤2.3% in every precision cell (V16), a confounder or a mediator with missingness
takes the shipped pipeline to +5% and, under MAR, to +10.4%.** That is larger than anything
V4–V14 measured, and it is in the configuration the documentation recommends.

**Coverage does not detect it — the fourth independent instance** (V3, V8, V15, now V17).
Coverage sits at 0.93–0.97 and width/SE at 0.98–1.13, i.e. the intervals are honestly sized.
They cover because at `n` = 800 the bias is only **0.35–0.67 empirical SE**.

**And that is a sample-size accident.** Relative bias is roughly constant in `n` while the SE
falls as `1/√n`. Taking `zrmar_pipe` at face value — bias 0.0416, empirical SE 0.0624 —
`n` = 8,000 would give SE ≈ 0.0197 and bias ≈ 2.1 SE, at which point coverage collapses
toward 50%. **This is a projection from one cell, not a measurement**, and it is the obvious
next thing to check: the same cells at larger `n`. A defect that hides at `n` = 800 and
surfaces at `n` = 8,000 is the worst possible shape for a pipeline built for large-N data.

---

## V17b: `use_as_auxiliary` does not reach the X block — real, and small

`docs/variable-dictionary.md` documents `FALSE / FALSE / TRUE` as *"use as imputation
predictor only"*. The Z block honours it; `00_censored_exposure.R`'s `auto_preds` reads
`use_in_model` only, so an auxiliary is silently dropped from the **censored exposure draw**.
`auxZ_shipped` runs the real `auto_preds` path — **exercised here for the first time in any
track**, and it ran without error in 1,000 tasks. `auxZ_asdoc` adds `Z₁` back.

| cell | `auxZ_shipped` | `auxZ_asdoc` | paired difference |
|---|---|---|---|
| `zr_collider` | +0.49% | +0.19% | **−0.30 ± 0.10 pp** |
| `zrmar_collider` | +0.67% | +1.27% | **+0.60 ± 0.10 pp** |

**Detectable, not material.** Both are significant at 500 reps and both are under 1 pp, and
**the sign differs between the two mechanisms** — so there is no consistent direction to
claim, in either the pipeline's favour or against it.

**Recommendation: fix it, for consistency rather than accuracy.** Adding
`use_as_auxiliary` to `auto_preds` is a two-line change that makes the code do what the
documentation says. On these numbers it will not measurably change an answer, so it is not
urgent and it should not be sold as an accuracy improvement. The alternative — amending the
documentation to say auxiliaries reach the covariate block only — is worse, because the
censored draw is the one this strategy exists for.

---

## What this revises

| claim | status |
|---|---|
| `THEORY.md` congeniality condition, covariate half | **Gains its scope qualifier**: binds a covariate's draw when that covariate is adjusted for **and** on an open X–Y path. Confirmed on 4/4 gates |
| "…through the covariate's association with the exposure" (§8j, pre-run) | **Withdrawn** — the pipe refutes it |
| "MAR adds +4 to +6 pp on-path" | **Refuted** — measured +0.7, +0.7, −0.6, +2.3 |
| "The path structure dominates the mechanism" | **Confirmed, and more strongly than derived** |
| V13: "a Z-only fix moves bias away from truth" | **Scope-limited to an inert covariate.** V13 measured `noYz`-adjacent effects where that arm moves 1.2 pp; here it moves 25 |
| V12's Z-draw shape attribution (32–49% of the penalty) | **Scope-limited to the precision structure** for the same reason |
| Shipped default's accuracy under a causally active covariate | **New defect, unmeasured before**: +5% to +10.4%, invisible to coverage |

---

## Caveats (scope of this run)

- **One DGP per role, `n` = 800, `m` = 30, additive surface, 40% non-detects on the focal
  exposure, 40% covariate missingness.** Arrow strengths are fixed (`delta_zx` = `delta_xz` =
  `delta_yz` = 0.60); the +25 pp figures scale with them and no sweep was run.
- **The shipped-default bias is measured in four cells and is not attributed.** Whether it
  comes from the X block, the Z block or their alternation needs a V13-style decomposition
  under a fork — which the arms now exist for.
- **The `n` = 8,000 coverage projection is arithmetic on one cell, not a measurement.**
- **`mixed` presses the binary `Z₂` into service as the mediator.** A third continuous
  covariate would be a cleaner design; the fork half of that cell is the reliable part.
- **The MAR conjecture (BART recovers `Y` from proxies) is untested.** A `mice pmm` Z block
  under the same cells is the discriminating run.
- **The collider cell's analysis omits `Z₁` by construction** (`truth$z_in_model`). A user
  who leaves a collider in the model gets the −0.10 estimate from a true +0.40 that
  `test_v17_zrole.R` records — a specification error no imputation method can repair, and
  one the shipped `use_in_model = TRUE` convention invites.
