# V12 findings — the Z-draw shape is *a* mechanism, worth about 40% of the penalty

**Run:** 2026-09-01 08:33 → 10:28 JST (**114.8 min**, 8 scenarios × 3 arms + oracle × 1000
reps = 8,000 tasks, `n_ok` = 1000 in all 32 cells, zero task errors, zero non-finite
estimates). Results: `results/v12_latest.rds`. Criteria: `../PLAN_pipeline_validation.md` §8e.

**The question.** V11 found a non-linear outcome costs ≈3 pp of bias for **every** imputer,
spread only 0.44 pp — strange if the conditional *mean* were the problem, since the arms
differ enormously in how they model it. The proposed explanation: they all draw a covariate
as *(fitted mean) + homoscedastic Gaussian noise*, which is exactly right under a linear
outcome and wrong under a non-linear one, where the true conditional is heteroscedastic and
skewed.

**Answer: confirmed, but it is roughly 40% of the story, not all of it.**

---

## Control: passes, twice over

In the linear cells the true conditional *is* Gaussian, so the two draws must agree.

| cell | `smc_zgauss` − `smc_zexact` | MC error | *p* |
|---|---|---|---|
| `mcar_z40` | −0.031 pp | 0.036 | 0.40 |
| `missing_y20` | **0.000 pp** | 0.000 | — |
| `combined` | +0.126 pp | 0.156 | 0.42 |
| `mar_z40` | +0.018 pp | 0.038 | 0.65 |

**A second, unplanned control fell out of the design.** `missing_y20` has a missing
*outcome* but complete *covariates* — so no Z draw happens at all and the two arms are
byte-identical, giving exactly 0.000. That the machinery produces an exact zero where it
should is a stronger check than the near-zeros elsewhere.

---

## Primary: the cost of getting the shape wrong

Paired; both arms share the same exact conditional **mean** and the same exposure draw, so
the difference is shape alone.

| cell | covariate missingness | Δ (pp) | MC error | *p* |
|---|---|---|---|---|
| `ynl_mcar_z40` | 40% MCAR | **+1.894** | 0.061 | <1e-4 |
| `ynl_mar_z40` | 40% MAR | **+1.278** | 0.073 | <1e-4 |
| `ynl_combined` | 20% MCAR | +0.088 | 0.152 | 0.56 |
| `ynl_missing_y20` | none | 0.000 | — | — |

**Shape is real where there is enough to impute**, and scales with how much: 1.9 pp at 40%
missingness, 1.3 pp at 40% MAR, indistinguishable from zero at 20%, exactly zero at none.
That monotone pattern is itself evidence the effect is the Z draw and not something else.

---

## Secondary: how much of V11's ≈3 pp does that explain?

Degradation from each linear cell to its non-linear twin, paired within arm:

| arm | `mcar_z40` | `missing_y20` | `combined` | `mar_z40` | **mean** |
|---|---|---|---|---|---|
| `pipeline_bartMI` | −2.37 | −3.17 | −2.26 | −3.76 | **−2.89 pp** |
| `smc_zgauss` | +2.00 | +0.01 | +0.01 | +1.19 | **+0.81 pp** |
| `smc_zexact` | +0.08 | +0.01 | +0.05 | −0.07 | **+0.02 pp** |

**The full remedy eliminates the penalty.** `smc_zexact` degrades by 0.02 pp where the
shipped path degrades by 2.89 — under a non-linear outcome it is essentially unbiased
(−0.21% and +0.02% in the two 40% cells).

**But shape alone is not the full remedy.** Decomposing the gap between the shipped arm and
the exact one:

| cell | `bartMI` | `smc_zgauss` | `smc_zexact` | shape's share |
|---|---|---|---|---|
| `ynl_mcar_z40` | −4.04% | +1.69% | −0.21% | 1.9 of 3.8 pp — **49%** |
| `ynl_mar_z40` | −3.93% | +1.30% | +0.02% | 1.3 of 4.0 pp — **32%** |

The remaining half to two-thirds is what the SMC arms *also* get exactly right: the
conditional **mean** (BART must estimate it; SMC is handed it) and the **exposure draw**
(`leftcens` linear conditional versus the exact one). **This run cannot separate those two**
— both differ between `bartMI` and `smc_zgauss`, and no arm isolates them.

---

## What this means for roadmap item 08

**The item is confirmed and should stay open**, with a realistic payoff:

- A shape-aware Z draw recovers **1.3–1.9 pp** in cells with substantial covariate
  missingness, and nothing where covariates are complete.
- It will **not** recover the full ≈3 pp. Half or more depends on the conditional mean and
  the exposure draw, which a shape fix does not touch.
- The gain is largest exactly where the pipeline is most used — heavy covariate missingness
  is the case multiple imputation exists for.

**The unresolved half is the obvious follow-up**: an arm with the exact Z draw but the
*shipped* exposure draw would separate the mean/X contributions. That is a small addition to
the existing harness.

## Caveats (scope of this run)

1. **`smc_zgauss` is handed the exact conditional mean**, which no real imputer has. So the
   1.3–1.9 pp is the shape cost *given a perfect mean* — an upper bound on what a shape fix
   buys if the mean is estimated well, and possibly optimistic if it is not.
2. **The remaining 51–68% is unattributed.** It is the conditional mean, the exposure draw,
   or both; this design conflates them.
3. **Both SMC arms over-cover**: coverage 0.982–0.984 with `width/SE` 1.22–1.26, against
   `bartMI`'s well-calibrated 0.955 at 1.03–1.07. The exact draws buy bias at the price of
   conservative intervals. A shipped version would need that examined — it is not obviously
   the right trade for every analysis.
4. **Shape needs missingness to bite**: null at 20% MCAR (0.09 pp, *p* = 0.56) and exactly
   zero with complete covariates. The 1.9 pp figure applies to heavy missingness.
5. **One curvature form** (`Z1^2`, coefficient 0.40), additive ERF, scalar estimand,
   `m` = 30, `n` = 800. Both SMC arms are harness instruments, not pipeline code.
