# V4 findings — the under-coverage attributed: improper imputation in the Z block

**Run:** 2026-08-25, 17:34–22:48 (**5.25 h**, 1500 tasks, 22 fork workers).
Track V4 of [`../PLAN_pipeline_validation.md`](../PLAN_pipeline_validation.md) §9.

**Design:** 5 scenarios × 4 arms × **300 replications**, `m`=30, `outer_sweeps`=3,
`margin`="shash", base seed **20260825** (identical to V2, so the reference arm sees the
same data). **Zero worker failures, zero task errors, 300/300 reps in every cell.**
Estimand: `b_logX1`, true value **0.40**.

**The question.** V2 found the shipped engine's point estimates sound everywhere
(|bias| ≤ 2.71%) but its 95% intervals ~17% too narrow when covariate *and* outcome
missingness are both substantial (coverage 0.893), and proved it is not a small-`m`
artefact. Two mechanisms were in play. This run changes one component at a time to
attribute the shortfall, and reports the Rubin variance decomposition V2 could not.

---

## Results

`width/SE` = (mean CI width ÷ 3.92) ÷ empirical SE. 1.0 = intervals match the true
sampling spread; below 1.0 = too narrow.

| scenario | arm | rel. bias | width/SE | `Ubar` | `B` | FMI | coverage |
|---|---|---|---|---|---|---|---|
| **base** | oracle | −0.62% | 1.029 | — | — | — | 0.953 |
| | block_fcs | −1.06% | 1.005 | 0.00164 | 0.00077 | 0.322 | 0.950 |
| | properZ | −1.06% | 1.005 | 0.00164 | 0.00077 | 0.322 | 0.950 |
| | noMID | −1.06% | 1.005 | 0.00164 | 0.00077 | 0.322 | 0.950 |
| | properZ_noMID | −1.06% | 1.005 | 0.00164 | 0.00077 | 0.322 | 0.950 |
| **mcar_z20** | block_fcs | −0.95% | 1.000 | 0.00156 | 0.00078 | 0.337 | 0.957 |
| | **properZ** | −1.37% | 1.043 | 0.00167 | 0.00079 | 0.325 | 0.963 |
| **mcar_z40** | block_fcs | −2.11% | **0.899** | 0.00147 | 0.00078 | 0.348 | 0.933 |
| | **properZ** | −3.18% | **1.085** | 0.00168 | 0.00099 | 0.374 | **0.960** |
| | noMID | −2.11% | 0.899 | 0.00147 | 0.00078 | 0.348 | 0.933 |
| **missing_y20** | block_fcs | +2.71% | 0.964 | 0.00208 | 0.00070 | 0.256 | 0.927 |
| | **properZ** | −0.95% | **1.003** | 0.00210 | 0.00083 | 0.286 | **0.947** |
| | **noMID** | +4.14% | **0.858** | 0.00146 | 0.00078 | 0.350 | **0.897** |
| | properZ_noMID | **−7.14%** | 1.051 | 0.00171 | 0.00120 | 0.413 | 0.933 |
| **combined** | block_fcs | +2.02% | **0.843** | 0.00188 | 0.00067 | 0.268 | 0.893 |
| | **properZ** | −1.34% | **0.973** | 0.00203 | 0.00099 | 0.329 | **0.933** |
| | **noMID** | +3.37% | **0.754** | 0.00132 | 0.00072 | 0.356 | **0.850** |
| | properZ_noMID | **−6.24%** | 1.003 | 0.00168 | 0.00133 | 0.443 | 0.920 |

**The negative control works.** In `base` — no covariate missingness, no outcome
missingness — all four arms return *identical* numbers to five decimals. Both
interventions are genuine no-ops there, which is exactly what makes the differences
elsewhere interpretable. The same holds within `mcar_z20`/`mcar_z40`, where `noMID` is
identical to the reference because there is no missing outcome for MID to act on.

---

## Mechanism (a) CONFIRMED, and localised to `B`

Swapping the Z-block imputer for a proper Bayesian draw restores calibration in every
cell where the Z block does real work:

| scenario | block_fcs → properZ (width/SE) | coverage |
|---|---|---|
| mcar_z40 | 0.899 → **1.085** | 0.933 → 0.960 |
| missing_y20 | 0.964 → **1.003** | 0.927 → 0.947 |
| combined | 0.843 → **0.973** | 0.893 → 0.933 |

Paired against identical data, the effect sits specifically in the **between**-imputation
term, not the within-imputation term:

| scenario | SE | **`B`** | `Ubar` | t(SE) |
|---|---|---|---|---|
| base | +0.0% | +0.0% | +0.0% | — |
| mcar_z20 | +2.9% | +8.0% | +7.5% | 7.9 |
| mcar_z40 | +9.4% | **+41.4%** | +14.4% | 20.0 |
| missing_y20 | +3.0% | **+29.2%** | +1.4% | 8.4 |
| combined | +9.0% | **+58.9%** | +8.4% | 25.7 |

`B` is understated by 8–59%, scaling with how much the Z block does, while `Ubar` moves
far less. That is the textbook signature of **improper multiple imputation**: miceRanger's
random-forest / predictive-mean-matching draw does not sample imputation-model parameters
from a posterior, so between-imputation variability is too small. The X block, which
*does* draw its margin and regression parameters properly, is calibrated wherever it works
alone (`base`: width/SE 1.005).

The `se_deficit` column tells the same story directly: the pooled SE understates the true
sampling SE by **16%** in `combined` and **10.6%** in `mcar_z40` under the shipped engine,
falling to **3.2%** and **−7.8%** with a proper draw.

---

## Mechanism (b) REFUTED — and the earlier reasoning was wrong

Disabling MID does not fix the intervals. It makes them **worse**, in both cells where it
is active:

| scenario | block_fcs → noMID (width/SE) | coverage | rel. bias |
|---|---|---|---|
| missing_y20 | 0.964 → **0.858** | 0.927 → 0.897 | +2.71% → +4.14% |
| combined | 0.843 → **0.754** | 0.893 → 0.850 | +2.02% → +3.37% |

And the both-off arm buys calibration at the cost of serious bias: **−7.14%** in
`missing_y20`, **−6.24%** in `combined`, against a 3% acceptance bar. Keeping imputed-`Y`
rows pulls the estimate toward the imputation model — which is precisely why MID exists.
**MID is doing its job and should stay.**

### Correcting the V2 reasoning

[`FINDINGS_v2.md`](FINDINGS_v2.md) argued a two-mechanism account on this basis:

> `missing_y20` has **complete covariates**, so the Z block barely works — yet width/SE
> is 0.96. The Z block cannot explain that cell.

**That premise was false.** `run_censored_exposure_block_fcs()` sets
`z_spec_as$imputation$impute_y <- TRUE` so the X block can condition on a complete `Y`.
In `missing_y20` the Z block is therefore busy imputing **`Y` itself** — improperly. That
is why `properZ` fixes that cell (`B` +29.2%) despite the covariates being complete, and
why `noMID`, which leaves the improper `Y` imputations in the analysis, makes it worse.

**One mechanism explains all five scenarios.** The two-mechanism account is withdrawn.

---

## What this does and does not license

**Established:** the interval under-coverage is caused by improper imputation in the
miceRanger Z block understating the between-imputation variance `B`. It is not `m`
(V2 refuted that on identical data), not MID, and not the censored-exposure X block.

**Scope beyond this feature.** The Z block is the pipeline's *general* row-level
imputation path. Nothing about this finding is specific to censored exposures — any
analysis run through this pipeline with substantial covariate or outcome missingness
should be expected to produce intervals that are too narrow, by a margin that grows with
the amount imputed.

**NOT established: that swapping the imputer is the fix.** Two reasons to be careful:

1. **`properZ` overshoots** — width/SE 1.085 in `mcar_z40` and 1.043 in `mcar_z20`, i.e.
   intervals ~4–9% *too wide* where the reference was calibrated or nearly so. It also
   shows `se_deficit` of −7.8% in `mcar_z40`.
2. **`properZ` bias is slightly worse** in `mcar_z40` (−3.18% against the reference's
   −2.11%), marginally outside the 3% bar.

`.v4_proper_row_imputation()` is a deliberately simple linear/logistic draw built as a
**validation instrument to isolate the mechanism**, not as a general-purpose imputer. It
has no non-linearity and no interactions, which is adequate for this DGP's linear-Gaussian
covariates and would be inadequate in general. The correct fix is to make the Z block
propagate parameter uncertainty properly — via a `mice`-style proper draw, a bootstrapped
random forest, or miceRanger's own settings if they can deliver it — not to adopt this
instrument.

---

## Caveats (scope of this run)

1. **MCAR, not MAR.** Covariate and outcome missingness are MCAR. MAR-on-`Z` is the more
   realistic mechanism and remains unswept.
2. **Additive ERF only.** All five scenarios are additive; mixture surfaces are Track V3.
3. **One estimand.** `b_logX1` only. Whether the `B` shortfall is uniform across
   parameters is untested.
4. **`properZ` is not validated as an imputer** — only as an instrument. Its overshoot
   (item 1 above) means it should not be read as a calibrated alternative.
5. **The proper draw is cheaper**, so the arms are not runtime-comparable: any real fix
   must be costed separately against the Z block's dominant runtime (V2 measured 90
   miceRanger fits per replication).
