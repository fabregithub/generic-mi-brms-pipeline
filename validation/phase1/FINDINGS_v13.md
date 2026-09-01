# V13 findings — the exposure draw dominates, and a Z-only fix makes things worse

**Run:** 2026-09-01 12:32 → 15:06 JST (**153.3 min**, 8 scenarios × 5 arms + oracle × 1000
reps = 8,000 tasks, `n_ok` = 1000 in all 48 cells, zero task errors, zero non-finite
estimates). Results: `results/v13_latest.rds`. Criteria: `../PLAN_pipeline_validation.md` §8f.

**The question.** V12 confirmed the Z-draw *shape* mechanism but could attribute only 32–49%
of the ≈3 pp non-linear-outcome penalty to it, because its design varied the Z draw and the
exposure draw together. This run separates them.

**The answer overturns the plan.** The **exposure draw** is the dominant term, and fixing the
Z draw alone — which is exactly what roadmap item 08 proposed — **makes bias worse**.
Item 08 should not be built as scoped.

---

## Controls: one passed, one failed and bounded the run

**CONTROL 2 passed.** The `_xship` arms reproduce the shipped arm's between-imputation
variance, confirming they genuinely share its exposure draw:

| arm | `b` | `b_share` |
|---|---|---|
| `pipeline_bartMI` | 0.00097 | 0.341 |
| `smc_zgauss_xship` | 0.00096 | 0.291 |
| `smc_zexact_xship` | 0.00095 | 0.308 |
| `smc_zgauss` *(exact X)* | 0.00037 | 0.162 |
| `smc_zexact` *(exact X)* | 0.00037 | 0.167 |

**CONTROL 1 failed, in two of four linear cells** — and the failure localised the problem
rather than voiding the run:

| cell | missing outcome? | spread across the four SMC arms |
|---|---|---|
| `mcar_z40` | no | **0.10 pp** ✅ |
| `mar_z40` | no | **0.13 pp** ✅ |
| `combined` | 20% | 4.10 pp ❌ |
| `missing_y20` | 20% | **11.33 pp** ❌ |

**Cause, identified:** the SMC instrument does not impute `Y`. The pipeline does — its Z
block sets `impute_y = TRUE` so the outcome is complete when the X block runs, and MID
deletes those rows before the fit. Without that, the `leftcens` call receives a predictor
matrix containing `NA`, and the conditional it fits is wrong for *every* row, not only the
incomplete ones. The exact-X arms escape because a censored row with missing `Y` falls back
to the prior mean and `lm` drops it anyway.

**Consequence:** everything below is restricted to the two cells with a complete outcome,
where CONTROL 1 holds at 0.10–0.13 pp. Those are also the only cells where V12 found a
significant shape effect, so nothing of interest is lost — but the instrument's limitation is
real and is recorded in the caveats.

---

## Primary: the three-way decomposition

Signed contributions, paired per replication. **They sum to the measured total exactly.**

### `ynl_mcar_z40` — starting from `pipeline_bartMI` at −4.04%

| step | contribution | MC error | *p* |
|---|---|---|---|
| Z conditional **mean** (BART-estimated → exact) | **+1.85 pp** | 0.104 | <1e-4 |
| Z draw **shape** (Gaussian → exact) | **−1.43 pp** | 0.083 | <1e-4 |
| **Exposure draw** (`leftcens` → exact) | **+3.41 pp** | 0.177 | <1e-4 |
| total | +3.83 pp | 0.224 | |
| *sum of the three* | *+3.83 pp* | | |

### `ynl_mar_z40` — starting from −3.93%

| step | contribution | MC error | *p* |
|---|---|---|---|
| Z conditional **mean** | **+2.14 pp** | 0.145 | <1e-4 |
| Z draw **shape** | **−1.32 pp** | 0.090 | <1e-4 |
| **Exposure draw** | **+3.14 pp** | 0.176 | <1e-4 |
| total | +3.95 pp | 0.246 | |
| *sum of the three* | *+3.95 pp* | | |

**The exposure draw is the largest single term** — 3.1–3.4 pp of a 3.8–4.0 pp total, larger
than either Z component.

---

## The shape fix has a fixed direction, and that is the problem

The shape effect is a consistent **negative shift in signed bias, at both exposure settings**:

| cell | shape effect, X = shipped | shape effect, X = exact |
|---|---|---|
| `ynl_mcar_z40` | −1.43 pp | −1.89 pp |
| `ynl_mar_z40` | −1.32 pp | −1.28 pp |

Same direction, similar magnitude. What differs is **where you start from**:

| cell | Zgauss/Xship | Zexact/Xship | Zgauss/Xexact | Zexact/Xexact |
|---|---|---|---|---|
| `ynl_mcar_z40` | −2.19% | **−3.62%** | +1.69% | **−0.21%** |
| `ynl_mar_z40` | −1.79% | **−3.11%** | +1.30% | **+0.02%** |

With the **exact** exposure draw the bias sits *positive* (+1.69%, +1.30%), so a negative
shift carries it toward zero — the shape fix helps, which is what V12 saw. With the
**shipped** exposure draw the bias is already *negative* (−2.19%, −1.79%), so the same shift
carries it further away — the shape fix hurts.

**The Gaussian Z draw and the linear exposure draw were partly cancelling.** Removing one
error while leaving the other exposes it. This is the third time this pattern has appeared
in the study (V11's `forest`/`properBoot`, V12's reading of shape, and now this), and it is
the reason single-component fixes cannot be judged in isolation.

---

## Secondary: a Z-only fix is not sufficient — it is counterproductive

`smc_zexact_xship` applies **both** Z fixes (correct mean and correct shape) and leaves the
exposure draw as shipped — exactly what roadmap item 08 proposed to build:

| cell | shipped | **Z-only fix** | full fix | Z-only recovers |
|---|---|---|---|---|
| `ynl_mcar_z40` | −4.04% | **−3.62%** | −0.21% | **11%** of the gap |
| `ynl_mar_z40` | −3.93% | **−3.11%** | +0.02% | **21%** of the gap |

It recovers 11–21%, and only because the Z *mean* correction (+1.85, +2.14) outweighs the Z
*shape* correction working against it (−1.43, −1.32). **The shape component — the part V12
identified and item 08 was written around — is actively harmful in the shipped
configuration.**

---

## What this settles

- **Item 08 should not be built as scoped.** A shape-aware Z draw, on its own, moves bias
  away from truth in the configuration the pipeline actually ships.
- **The exposure draw is the priority**, at 3.1–3.4 pp of a ~3.9 pp penalty. That is item
  07's territory — the same `leftcens` linear conditional that V3/V9 found destroys mixture
  curvature. **The two defects share a cause after all**, which is the opposite of what V12's
  framing suggested.
- **V12's conclusion needs qualifying, not withdrawing.** Its shape measurement was correct;
  its implication — that a shape fix would buy 1.3–1.9 pp — held only because it measured
  with the exposure draw already exact. Stated as a standalone remedy it was wrong.

## Caveats (scope of this run)

1. **Two cells, not four.** CONTROL 1 failed wherever the outcome had missing values, because
   the SMC instrument does not impute `Y` as the pipeline does. Fixing that means adding a
   `Y` draw plus MID to the instrument — worth doing before any run that needs the
   missing-outcome cells.
2. **The Z "mean" step conflates two things**: BART's estimation error *and* any difference
   in what surface is assumed. `smc_zgauss_xship` is handed the true outcome model; BART
   estimates a mean without knowing it. So +1.85/+2.14 pp is an upper bound on what a better
   mean model could buy.
3. **The exposure-draw step is measured against an *exact* alternative** that a real
   implementation cannot have — it uses the generator's true surface. The 3.1–3.4 pp is
   therefore the ceiling for fixing that block, not a delivery estimate.
4. **One curvature form**, additive ERF, scalar estimand, `m` = 30, `n` = 800, complete
   outcome. All `smc_*` arms are harness instruments, not pipeline code.
5. **Signed, not absolute, decomposition.** The steps sum exactly because they are signed
   contributions to the bias. Reading them as improvements requires knowing the starting
   sign, which is precisely the trap this run exposed.
