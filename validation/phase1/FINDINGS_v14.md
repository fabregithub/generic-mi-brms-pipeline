# V14 findings — a shippable exposure draw recovers 83%, and calibrates better than the exact one

**Run:** 2026-09-01 16:52 → 21:03 JST (**250.3 min**, 8 scenarios × 3 arms + oracle × 1000
reps = 8,000 tasks, `n_ok` = 1000 in all 32 cells, zero task errors, zero non-finite
estimates). Results: `results/v14_latest.rds`. Criteria:
`../PLAN_pipeline_validation.md` §8g.

**The question.** V13 showed the exposure draw carries most of the non-linear-outcome
penalty; V9 that it carries the whole mixture failure. Both were measured against a sampler
that **knew the generator's surface**. Can a draw knowing only the analysis *formula* —
estimating everything else, as a shipped version must — match it?

**Answer: mostly, and with better intervals.** `smc_xgrid` removes **83%** of the penalty
and is the **best-calibrated arm of the three**. It misses the registered ±1 pp bar in one
of three cells, and it is compromised in a cell shape the instrument cannot yet handle.

*(Runtime note: 250 min against a predicted 3.2 h — 30% over, the largest cost miss so far.)*

---

## Control: passes, except where the instrument breaks

| cell | design | `bartMI` | `smc_zexact` | `smc_xgrid` | spread |
|---|---|---|---|---|---|
| `mcar_z40` | X1 censored | −1.67% | −0.29% | −0.47% | 1.38 pp |
| `missing_y20` | X1, 20% missing Y | −0.22% | +0.14% | +0.20% | 0.41 pp |
| `mar_z40` | X1 | −0.17% | +0.09% | +0.05% | 0.26 pp |
| **`combined`** | **all 3 censored** + missing Y | −0.16% | −0.13% | **−1.49%** | 1.36 pp |

**`combined` is the one that matters.** It is linear, so all three arms should agree, and
`smc_xgrid` is 1.36 pp adrift. In its non-linear twin the gap is worse: `ynl_combined` gives
`smc_zexact` −0.08% against `smc_xgrid` **−3.81%** — the candidate is worse there than the
*shipped* arm (−2.41%).

**Diagnosis.** `combined` is the only cell censoring **all three** exposures, and it also has
20% missing outcomes. The grid arm has no likelihood for a row whose `Y` is missing, so that
row falls back to a deterministic value instead of a draw; with three censored exposures the
fallback compounds through the FCS loop. `missing_y20` has the same missing-`Y` fraction but
censors one exposure and is fine (+0.20% vs +0.14%), so **missing `Y` alone is not
sufficient** to break it.

**This run cannot separate `censor_all` from `censor_all × missing Y`** — no cell has the
former without the latter. Both `combined` cells are excluded from what follows.

---

## Primary: within ±1 pp of the exact sampler in two of three cells

Paired, `smc_xgrid` − `smc_zexact`:

| cell | Δ (pp) | MC error | *p* | vs ±1 pp |
|---|---|---|---|---|
| `ynl_mcar_z40` | −0.556 | 0.167 | 0.0009 | **PASS** (CI inside) |
| `ynl_missing_y20` | +0.130 | 0.147 | 0.38 | **PASS** (CI inside) |
| `ynl_mar_z40` | **−1.363** | 0.182 | <1e-4 | **FAIL** |

The failure is real, not noise (7.5σ), and `ynl_mar_z40` is also where recovery is weakest.

## Secondary: 83% of the penalty removed

Degradation from each linear cell to its non-linear twin, paired within arm:

| arm | `mcar_z40` | `missing_y20` | `mar_z40` | **mean** |
|---|---|---|---|---|
| `pipeline_bartMI` | −2.37 | −3.17 | −3.76 | **−3.10 pp** |
| `smc_zexact` (knows the surface) | +0.08 | +0.01 | −0.07 | **+0.01 pp** |
| **`smc_xgrid`** (knows only the formula) | −0.30 | +0.09 | −1.39 | **−0.53 pp** |

Per cell, the share of the shipped arm's penalty removed: **86%** (`ynl_mcar_z40`),
**~100%** (`ynl_missing_y20`), **66%** (`ynl_mar_z40`).

## Secondary: the lower `b` is not a variance collapse

`smc_xgrid` has a smaller between-imputation variance than the shipped arm
(`b` 0.00068 vs 0.00102), which was pre-registered as a warning sign — a draw that buys bias
by collapsing `b` would be V4's improper-imputation defect wearing a new hat. **It is not:**

| arm | `b` | `b_share` | **width/SE** | **coverage** |
|---|---|---|---|---|
| `pipeline_bartMI` | 0.00102 | 0.344 | 1.049 | 0.954 |
| `smc_zexact` | 0.00037 | 0.175 | **1.247** | **0.982** |
| **`smc_xgrid`** | 0.00068 | 0.273 | **1.004** | **0.950** |

`smc_xgrid` is the **best-calibrated arm in the run** — intervals 0.4% wider than the true
sampling spread, coverage 0.950. The *exact* sampler over-covers badly (width/SE 1.247,
coverage 0.982): drawing from the correct tighter conditional shrinks `b` faster than it
shrinks the estimator's actual spread. A lower `b` is a defect only when the intervals end
up too narrow, and here they do not.

---

## What this revises in the theory

[`THEORY.md`](../THEORY.md) records V9's finding that **the functional *form* is what
matters, not the parameters** — its plug-in arm matched the oracle to 0.2 pp. **V14
qualifies that.** Here the same substitution costs **0.53 pp on average and 1.36 pp at
worst**, an order of magnitude more than V9 measured.

The two are not in contradiction — V9 measured BKMR mixture estimands, V14 a scalar
coefficient under a non-linear outcome — but the general claim "estimating the parameters is
free" is now known to be setting-dependent, and should not be carried forward unqualified.

## What this means for item 07

**The method works and is worth building**, with two conditions:

1. **A missing-`Y` path is required, not optional.** The `combined` failure is an
   implementation gap, and the fix is known: impute `Y` (as the pipeline's own Z block does
   with `impute_y = TRUE`) so the likelihood is defined for every row, then MID-delete before
   the fit. Until then the draw degrades exactly where multiple censored exposures meet
   missing outcomes — a realistic combination.
2. **83% is not 100%.** The residual is the cost of estimating the outcome model, and
   `ynl_mar_z40` shows it can reach 1.4 pp. Whether that matters depends on the alternative:
   against the shipped arm's −3.10 pp it is a large improvement; against the exact sampler it
   is a real gap.

## Caveats (scope of this run)

1. **Six cells, not eight.** Both `combined` cells excluded for the control failure above.
2. **`ynl_mar_z40` fails the registered bar** at −1.363 pp, and is the weakest cell on
   recovery (66%). Not explained.
3. **Additive path only.** Whether this draw recovers the **mixture** estimands needs the
   BKMR harness (V9's cells) — a separate ~10 h run, and the other half of item 07's claim.
4. **Harness instrument, not pipeline code.** Nothing shipped changes; wiring
   `.ce_smc_x_grid()` into `00_censored_exposure.R` behind a config flag, with its own tests,
   remains undone.
5. **One curvature form**, additive ERF, scalar estimand, `m` = 30, `n` = 800, `n_grid` = 512.
