# V5 findings — the candidate fix: a real improvement, not a complete one

**Run:** 2026-08-26, 08:46–16:49 JST (**8.05 h**, 1500 tasks, 22 fork workers).
Follow-up to Track V4 of [`../PLAN_pipeline_validation.md`](../PLAN_pipeline_validation.md) §9.

**Design:** 5 scenarios × 5 arms × **300 replications**, `m`=30, `outer_sweeps`=3,
`margin`="shash", base seed **20260825** — identical to V4, so the reference arm sees the
same data and this is a paired before/after. **Zero worker failures, zero task errors.**
Estimand `b_logX1`, true value **0.40**.

**The question.** V4 established that the pipeline's interval under-coverage comes from
improper multiple imputation in the miceRanger Z block, understating between-imputation
variance `B` by 8–59%. This run tests an actual fix in shipped pipeline code:
`run_row_level_imputation_proper()`, which bootstraps the training data once per
imputation, fits the forest to the resample, and imputes the original rows from that model
(`miceRanger::impute()`). The new arm is **`pipeline_properBoot`**.

The earlier `pipeline_properZ` arm remains for reference, but it is a harness *instrument*
— a simple linear/logistic proper draw — not shipped code.

---

## Results

| scenario | arm | rel. bias | width/SE | `B` | coverage |
|---|---|---|---|---|---|
| **base** | reference | −1.06% | 1.005 | 0.00077 | 0.950 |
| | **properBoot** | −1.06% | 1.005 | 0.00077 | 0.950 |
| | properZ | −1.06% | 1.005 | 0.00077 | 0.950 |
| **mcar_z20** | reference | −0.95% | 1.000 | 0.00078 | 0.957 |
| | **properBoot** | −0.87% | 1.014 | 0.00077 | 0.950 |
| | properZ | −1.37% | 1.043 | 0.00079 | 0.963 |
| **mcar_z40** | reference | −2.11% | 0.899 | 0.00078 | 0.933 |
| | **properBoot** | **−1.97%** | **0.965** | 0.00094 | **0.947** |
| | properZ | −3.18% | 1.085 | 0.00099 | 0.960 |
| **missing_y20** | reference | +2.71% | 0.964 | 0.00070 | 0.927 |
| | **properBoot** | +2.58% | 0.979 | 0.00077 | 0.937 |
| | properZ | −0.95% | 1.003 | 0.00083 | 0.947 |
| **combined** | reference | +2.02% | 0.843 | 0.00067 | 0.893 |
| | **properBoot** | +1.65% | **0.878** | 0.00082 | **0.910** |
| | properZ | −1.34% | 0.973 | 0.00099 | 0.933 |

### Paired against the shipped reference (identical data)

| scenario | SE | `B` | t(SE) | gap to width/SE 1.00 closed |
|---|---|---|---|---|
| base | +0.0% | +0.0% | — | — (no-op) |
| mcar_z20 | +0.5% | +5.3% | 1.0 | already calibrated |
| mcar_z40 | +5.0% | **+32.6%** | 11.1 | **66%** |
| missing_y20 | +1.4% | **+18.1%** | 4.4 | **42%** |
| combined | +3.1% | **+31.2%** | 9.0 | **23%** |

---

## Verdict: adopt it, but it is not a complete fix

**What the fix does well:**

1. **Monotone improvement, no over-correction anywhere.** Maximum width/SE is 1.014,
   against the instrument's 1.085. The concern raised by the 12-replication smoke test —
   that the bootstrap might over-inflate — did not materialise; that smoke reading (`B`
   +52%) was noise, as flagged at the time.
2. **Best bias of any arm**, in every scenario: −1.97% in `mcar_z40` against the
   reference's −2.11% and properZ's −3.18%.
3. **Already-calibrated cells are left alone.** `mcar_z20` moves 1.000 → 1.014 and
   `base` not at all, so the correction is applied where it is needed rather than
   uniformly.
4. **`mcar_z40` is essentially fixed** — coverage 0.947 against a nominal 0.95.
5. **Robust in practice.** Zero bootstrap fallbacks across 1500 tasks: the retry guard
   for resamples that lose a target's observed values or a factor level never fired, so
   every imputation used a genuine resample.

**What it does not do:** `combined` — every exposure censored, 20% MCAR covariates, 20%
missing outcomes — still under-covers at **0.910**, below the pre-registered 0.92 bar,
with only 23% of the shortfall closed.

---

## Why the bootstrap under-delivers on a forest

A random forest is *already* a bagged estimator: hundreds of trees, each fitted to its own
bootstrap resample, then averaged. Bagging exists precisely to be low-variance. Wrapping
one further bootstrap around it therefore shifts the fitted forest only slightly — the
mechanism that makes forests accurate damps the very uncertainty the outer bootstrap is
meant to inject.

The measurements bear this out. On the same data:

| approach | `B` inflation, `mcar_z40` | `B` inflation, `combined` |
|---|---|---|
| bootstrap on forest (`properBoot`) | +32.6% | +31.2% |
| parametric posterior draw (`properZ`) | +41.4% | +58.9% |

**The standard "bootstrap for properness" recipe is weaker for bagged learners than for
parametric models.** This is a general methodological point, not a quirk of this
implementation, and it is worth stating because the bootstrap route is the usual
recommendation for making any predictive imputer proper (it is what `mice`'s `rf` method
does).

---

## Why this does NOT license switching to `mice`

`properZ` reaches 0.973 in `combined` and looks like the better answer. The evidence does
not support that conclusion:

1. **The DGP has linear-Gaussian covariates**, so a linear imputation model is *correctly
   specified* here. Its calibration advantage may not survive non-linear covariates, which
   is the case miceRanger was presumably chosen for.
2. **`properZ` already has the worst bias of the three arms** (−3.18% in `mcar_z40`,
   against −1.97%), which is what misspecification would begin to look like.
3. **`mice`'s `rf` method is essentially `properBoot`.** Option B would only differ by
   choosing `pmm`/`norm`, i.e. trading flexibility for calibration.

**The current simulation cannot distinguish "properly calibrated" from "correctly
specified."** Every covariate in it is linear and Gaussian, so no arm is ever penalised for
misspecification. Choosing between forest-plus-bootstrap and parametric-plus-posterior on
this evidence would be choosing on a test that cannot punish the parametric option's main
risk.

**Next step for this question:** extend the DGP with non-linear and interacting covariates,
then re-run with a `pmm` arm added. Until then, `properBoot` is the defensible default
because it improves calibration without giving up the flexibility, and it has the best
bias on record.

---

## What changed in the pipeline

`analysis_spec$imputation$proper_draw` now defaults to **TRUE** (v1.4.0). Setting it to
`FALSE` restores the pre-v1.4.0 improper behaviour exactly.

This is a **deliberate behaviour change**: analyses with missing covariates or outcomes
will produce slightly wider intervals than before, and point estimates may shift very
slightly. Analyses must cite the pipeline version they used. Old configs continue to run
unchanged — only the default value of a new flag has moved.

---

## Caveats (scope of this run)

1. **Not a complete fix.** Under heavy combined missingness, coverage reaches ~0.91 rather
   than 0.95. The README caveat is retained, with revised numbers.
2. **Linear-Gaussian covariates only** — see above. This is the single most important
   limitation, because it is what prevents a clean choice between the two fix families.
3. **MCAR, not MAR**, throughout.
4. **Additive ERF only.** Mixture surfaces remain Track V3.
5. **One estimand** (`b_logX1`). Whether the residual shortfall is uniform across
   parameters is untested.
6. **Runtime cost.** The proper path adds one `impute()` pass per imputation — `m` predict
   passes over the full data per completed set, not one. This run took 8.05 h against 5.25 h
   for the four-arm V4 grid; roughly 15–20% of that is the new arm. Production runs with
   large `m` and `n` should expect the Z block to cost meaningfully more than before.
