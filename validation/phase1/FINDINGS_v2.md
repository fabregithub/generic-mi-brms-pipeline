# V2 findings — robustness sweep, and an under-coverage finding

**Run:** 2026-08-25, 09:10–12:07 (**2.95 h**, 63.7 CPU-hours). Track V2 of
[`../PLAN_pipeline_validation.md`](../PLAN_pipeline_validation.md) §7.

**Design:** 10 scenarios × 300 replications (`large_n` at 50), `m`=30, `p`=3 exposures,
`outer_sweeps`=3, `margin`="shash", base seed **20260825**, 22 fork workers.
**2750 tasks, zero worker failures, zero task errors, 98% worker efficiency.**
Estimand: `b_logX1` (true value **0.40**). Reproduce with `./run_v2_robustness.sh`.

**The question.** V1 validated the shipped engine in *one corner* of the design space:
only X1 censored, covariates complete, outcome complete. Two of the pipeline's code
paths were therefore never executed — MID (which needs a missing outcome) and the
miceRanger Z block (which needs missing covariates), so the block-FCS *alternation
itself* had never run under Monte Carlo. V2 sweeps the axes that corner left untouched.

---

## Both V1 dead spots are now genuinely exercised

From the procedure notes across all 2750 tasks:

```
cens=1                    1850     single censored exposure
cens=3                     600     every exposure censored -- X-block loop
cens=1; MID dropped 160    300     MID fires
cens=3; MID dropped 160    300     MID fires, all exposures censored
WARN notes                   0
```

MID fired in every scenario with a missing outcome, dropping exactly the 160 imputed-`Y`
rows (800 × 0.2) each time, and the row-count assertion passed everywhere. The
multi-exposure X-block loop ran 600 times without a failure.

**Calibration check.** `base` reproduces V1's additive 40%-ND cell: −1.06% / coverage
0.950 here against +0.48% / 0.957 in V1. These are independent samples (different task
seeds), and the difference is ~1.5 MC SE — V2's machinery did not change the answer.

---

## Results — pipeline_block_fcs across all ten scenarios

`width/SE` is (mean CI width ÷ 3.92) ÷ empirical SE: at 1.0 the intervals match the true
sampling spread, below 1.0 they are too narrow.

| scenario | Z-block work | rel. bias | rmse/oracle | width/SE | coverage |
|---|---|---|---|---|---|
| **combined** | Z + MID | +2.02% | 1.43 | **0.84** | **0.893** |
| mcar_z40 | Z only | −2.11% | 1.35 | 0.90 | 0.933 |
| missing_y20 | MID only | +2.71% | 1.41 | 0.96 | 0.927 |
| censor_all | neither | +0.96% | 1.21 | 1.00 | 0.957 |
| mcar_z20 | Z only | −0.95% | 1.25 | 1.00 | 0.957 |
| base | neither | −1.06% | 1.25 | 1.01 | 0.950 |
| rho00 | neither | −0.47% | 1.30 | 1.03 | 0.970 |
| rho80 | neither | −1.37% | 1.15 | 1.08 | 0.960 |
| large_n (n=80 000) | neither | −0.09% | 1.09 | 1.09 | 0.960 |
| skew075 | neither | +1.57% | 1.14 | 1.13 | 0.957 |

**Bias passes everywhere.** The largest relative bias across all ten scenarios is
**2.71%**, inside the pre-registered 3% bar — including with every exposure censored, at
n = 80 000, under right-skew, and at exposure correlation 0.8.

**Nine of ten scenarios pass outright.** `combined` — every exposure censored *and* 20%
MCAR covariates *and* 20% missing outcomes — has coverage **0.893**, below the 0.90
"investigate" line the plan set in advance.

---

## The finding: uncertainty is under-propagated, the point estimates are not

This is a variance problem, not a bias problem. `width/SE` correlates **0.837** with
coverage across the ten scenarios, and the ordering is not random:

| Z-block involvement | width/SE | coverage |
|---|---|---|
| none (6 scenarios) | 1.00 – 1.13 | 0.950 – 0.970 |
| MID only | 0.96 | 0.927 |
| MCAR covariates 40% | 0.90 | 0.933 |
| both (combined) | 0.84 | 0.893 |

Three things make this a real signal rather than an artefact:

1. **It is dose-dependent.** 20% MCAR covariates is perfectly calibrated (width/SE 1.00,
   coverage 0.957); 40% is not (0.90, 0.933).
2. **It compounds.** MID alone and MCAR alone each cost a little; together they cost
   more than either.
3. **The oracle stays calibrated throughout** — its width/SE spans 0.95–1.07 across every
   scenario. The harness, the metrics, and the analysis model are fine. This is the
   procedure.

### Leading hypothesis — CONFIRMED by Track V4 (see [`FINDINGS_v4.md`](FINDINGS_v4.md))

The X block draws its margin parameters and regression coefficients from their posterior
— that is *proper* multiple imputation in Rubin's sense. Every scenario where only the X
block does work is calibrated to within a percent.

miceRanger's random-forest / predictive-mean-matching imputation does **not** draw
imputation-model parameters from a posterior. Improper MI understates
between-imputation variance, which is exactly the observed signature: correct centre,
intervals too narrow, worsening as the Z block does more of the work.

If that is the cause, **the implication reaches beyond this feature**: it would affect
any analysis run through this pipeline with substantial covariate missingness, not just
censored exposures. That makes it worth resolving properly rather than noting and moving
on.

### The discriminating test

Improper MI is a *bias in the variance estimator*: raising `m` cannot fix it. Small-`m`
noise in Rubin's `B/m` correction, by contrast, shrinks with `m`. So re-running
`combined` at `m`=100 on **identical data** separates the two hypotheses cleanly:

- coverage stays ≈0.89 → improper MI (or another structural cause); `m` is not the lever
- coverage moves toward 0.95 → a small-`m` artefact, and the fix is a larger default `m`

### Result: the small-`m` hypothesis is REFUTED (2026-08-25, 71 min)

`combined` re-run at `m`=100 on **byte-identical data** — confirmed by the fact that the
`m`-independent procedures (`oracle`, `complete_case`) reproduce to `max |diff| = 0`:

| arm | est | rel. bias | emp SE | CI width | width/SE | coverage |
|---|---|---|---|---|---|---|
| m=30 (main run, reps 1–100) | 0.4079 | +1.98% | 0.0596 | 0.1997 | 0.855 | 0.88 |
| m=100 (diagnostic) | 0.4072 | +1.81% | 0.0591 | 0.1986 | 0.857 | 0.87 |

Tripling the imputations changed nothing: CI width moved 0.6%, width/SE moved 0.002,
coverage is flat. On the same 100 datasets the oracle sits at width/SE 1.10, coverage 0.97.

**The variance estimator is systematically too small, not merely noisy at small `m`.** To
reach nominal coverage the intervals would need to be ~**17% wider** (0.0591 × 3.92 =
0.2317 against the 0.1986 produced). Raising `m` is not the lever.

### Correction: there are TWO candidate mechanisms, not one — WITHDRAWN, see below

The single "improper MI in miceRanger" hypothesis stated above is **incomplete**. The
scenario pattern rules out either mechanism acting alone:

- `mcar_z40` has a **complete outcome**, so MID never fires — yet width/SE is 0.90. MID
  cannot explain that cell; something in the Z-block imputation understates variance.
- `missing_y20` has **complete covariates**, so the Z block barely works — yet width/SE
  is 0.96. The Z block cannot explain that cell. Relevant prior result: von Hippel (2007)
  derived a *corrected* variance estimator for multiple-imputation-then-deletion
  precisely because naive Rubin pooling on MID data understates variance. This pipeline
  uses naive Rubin pooling.
- `combined` has both and lands at 0.84 — consistent with the two effects compounding.

A two-mechanism account fits all ten scenarios; the single-mechanism one does not.

> **WITHDRAWN (2026-08-25, Track V4).** The premise of the `missing_y20` bullet above is
> **false**. `run_censored_exposure_block_fcs()` sets `impute_y <- TRUE` internally so the
> X block can condition on a complete `Y` — so in `missing_y20` the Z block is *not* idle,
> it is busy imputing **`Y` itself**. V4 confirms it: a proper Z-block draw fixes that cell
> (`B` +29.2%) even though the covariates are complete, while disabling MID makes it
> *worse* (width/SE 0.964 → 0.858) and adds bias. **One mechanism — improper imputation in
> the Z block understating `B` — explains all ten scenarios.** MID is doing its job and
> should stay. See [`FINDINGS_v4.md`](FINDINGS_v4.md).

### What cannot be concluded from the data as stored

`rubin_pool()` returns only the pooled SE, so it is impossible to see from these results
whether the shortfall sits in the between-imputation term `B` or the within-imputation
term `Ubar`. That decomposition is the entry point for Track V4, whose scope was widened
to cover it: report `B`/`Ubar`/FMI, isolate the Z block by swapping in a proper Bayesian
draw, and isolate MID by disabling it.

---

## A reproducibility bug found and fixed while setting up that test

The per-task seed was keyed on the scenario's index **in the filtered list**, so
`SCENARIOS=combined` generated different data than the same scenario inside a full run.
The `m`=30 vs `m`=100 comparison would then have confounded two changes at once.

Fixed to key on the scenario's position in the **canonical, unfiltered** grid. This is
strictly better than a name hash: subsets are now consistent with full runs, *and* the
completed run above stays exactly reproducible (in a full run, filtered position and
canonical position coincide). `combined` is canonical index 6 either way.

---

## Cost (measured, not estimated)

| scenario | s/rep (1 worker) | CPU-h |
|---|---|---|
| large_n (n = 80 000) | 367.3 | 5.10 |
| combined | 265.8 | 22.15 |
| mcar_z20 | 170.5 | 14.21 |
| mcar_z40 | 133.1 | 11.09 |
| missing_y20 | 98.0 | 8.17 |
| censor_all | 13.0 | 1.08 |
| base / rho00 / rho80 / skew075 | 5.6 – 6.1 | 1.95 total |

**Runtime is dominated by the miceRanger Z block, not the censored-exposure X block.**
Anything that gives the Z block work costs 15–40× the base case, because the block-FCS
runs `outer_sweeps` × `m` = 90 miceRanger fits per replication. The X-block scale test in
the design plan (§10, Phase 2) called itself "a rounding error next to `miceRanger` and
`brms`" — this quantifies that: at n = 80 000 with *no* covariate missingness a
replication costs 367 s, while at n = 800 *with* MAR covariates and a missing outcome it
costs 266 s. Covariate missingness, not sample size, is the cost driver on this path.

This is useful input for Track V4: raising `m` on the censored path is expensive
primarily because of the Z block.

---

## What this changes

**The V1 claim stands, with its scope now measured rather than assumed.** For additive
exposure–response functions the shipped engine is unbiased to within ~1–2% with nominal
coverage — and V2 shows that holds with every exposure censored, at n = 80 000, under
right-skew, and at high exposure correlation.

**A new caveat applies to interval coverage.** With *both* substantial covariate
missingness and missing outcomes, 95% intervals achieve ~89–93%. Point estimates remain
sound. Until the discriminating test resolves the cause, analyses with heavy covariate
and outcome missingness should treat the intervals as somewhat anti-conservative.

## Caveats (scope of this run)

1. **The under-coverage cause is unresolved** — the leading hypothesis is stated above,
   the test is running, and nothing here should be read as having demonstrated a
   mechanism.
2. **`cens_mi_y_shash` is absent from four scenarios** (`censor_all`, both MCAR cells,
   `missing_y20`, `combined`). It hardcodes X1 as the only censored variable and passes
   everything else as complete predictors, so it is undefined there. In those scenarios
   the pipeline has no prototype counterpart — which is the point of having built it.
3. **MCAR, not MAR.** Covariate and outcome missingness are MCAR here. MAR-on-`Z` is the
   more realistic mechanism and is not swept.
4. **Mixture estimands untouched** — all V2 scenarios are additive. Track V3.
5. **`large_n` is 50 replications**, so its MC error is ~2.4× the other cells'.
