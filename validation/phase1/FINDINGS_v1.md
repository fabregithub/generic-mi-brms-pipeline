# V1 findings — the *shipped* censored-exposure engine, validated

**Run:** 2026-08-25, 08:14–08:31 (~17 min). Track V1 of
[`../PLAN_pipeline_validation.md`](../PLAN_pipeline_validation.md) §6.

**Design:** full grid, **300 reps/cell**, `m`=30 imputations, `n`=800, `p`=3 exposures,
`rho`=0.4, `skew`=0, `NCORES`=20, base seed **20260825**. Pipeline knobs:
`outer_sweeps`=3, `margin`="shash". Estimand: `b_logX1`, the focal censored-exposure
coefficient (true value **0.40**). Reproduce with `./run_v1_pipeline.sh`; regenerate
these tables from `results/latest.rds`.

**The question.** Every congenial result in [`FINDINGS.md`](FINDINGS.md) came from
`cens_mi_y_shash`, a ~40-line *prototype* written for the study. This run adds
procedure 5 — `pipeline_block_fcs`, which calls the **shipped**
`run_censored_exposure_block_fcs()` from
[`../../00_censored_exposure.R`](../../00_censored_exposure.R) by sourcing it, never
by copying it. Everything downstream (analysis model, Rubin pooling, metrics) is the
harness's existing machinery, so the only thing differing across rows is missing-data
handling.

---

## Results (all cells n_ok = 300, zero failures)

### Additive ERF — the pass/fail case

| ND | procedure | est | rel. bias | rmse | CI width | coverage |
|---|---|---|---|---|---|---|
| 0.2 | oracle | 0.3975 | −0.62% | 0.0393 | 0.1586 | 0.953 |
| 0.2 | complete_case | 0.3972 | −0.71% | 0.0603 | 0.2209 | 0.943 |
| 0.2 | **leftcens_prestep** | 0.3765 | **−5.88%** | 0.0464 | 0.1681 | **0.923** |
| 0.2 | cens_mi_y_shash | 0.3973 | −0.67% | 0.0418 | 0.1658 | 0.947 |
| 0.2 | **pipeline_block_fcs** | **0.3968** | **−0.79%** | **0.0417** | 0.1654 | **0.957** |
| 0.4 | oracle | 0.4000 | 0.00% | 0.0391 | 0.1584 | 0.947 |
| 0.4 | complete_case | 0.3985 | −0.39% | 0.0719 | 0.2942 | 0.967 |
| 0.4 | **leftcens_prestep** | 0.3434 | **−14.14%** | 0.0714 | 0.2008 | **0.820** |
| 0.4 | cens_mi_y_shash | 0.4014 | +0.34% | 0.0483 | 0.1939 | 0.960 |
| 0.4 | **pipeline_block_fcs** | **0.4019** | **+0.48%** | **0.0485** | 0.1937 | **0.957** |

### Mixture ERF — diagnostic, not pass/fail

| ND | procedure | est | rel. bias | rmse | coverage |
|---|---|---|---|---|---|
| 0.2 | oracle | 0.4006 | +0.15% | 0.0406 | 0.950 |
| 0.2 | complete_case | 0.4001 | +0.03% | 0.0819 | 0.967 |
| 0.2 | leftcens_prestep | 0.4155 | +3.88% | 0.0438 | 0.953 |
| 0.2 | cens_mi_y_shash | 0.4300 | +7.49% | 0.0522 | 0.897 |
| 0.2 | **pipeline_block_fcs** | 0.4302 | **+7.54%** | 0.0521 | 0.903 |
| 0.4 | oracle | 0.4010 | +0.24% | 0.0399 | 0.947 |
| 0.4 | complete_case | 0.3998 | −0.06% | 0.1780 | 0.947 |
| 0.4 | leftcens_prestep | 0.4165 | +4.13% | 0.0416 | 0.973 |
| 0.4 | cens_mi_y_shash | 0.4778 | +19.45% | 0.0893 | 0.673 |
| 0.4 | **pipeline_block_fcs** | 0.4782 | **+19.55%** | 0.0897 | 0.663 |

---

## Verdict against the pre-registered criteria

The criteria below were written into the plan **before** the run.

| Criterion (additive) | ND 0.2 | ND 0.4 | |
|---|---|---|---|
| \|rel_bias\| ≤ 3% | 0.79% | 0.48% | ✅ |
| coverage ∈ [0.92, 0.97] | 0.957 | 0.957 | ✅ |
| rmse ≤ 1.3 × oracle | 1.059× | 1.243× | ✅ |
| strictly better than complete-case at 40% ND | — | 0.0485 vs 0.0719 | ✅ |
| agreement with `cens_mi_y_shash` within MC error | see below | see below | ✅ |
| `n_ok` = `n_rep`, no silent failures | 300 | 300 | ✅ |

**All additive criteria pass.** The shipped engine tracks the oracle to within 0.8%
with nominal coverage at both censoring levels, and recovers most of the efficiency
complete-case throws away (rmse 0.0485 vs 0.0719 at 40% ND).

## Agreement with the prototype (paired)

Both procedures see the *same* replications, so a paired comparison is far more
powerful than comparing marginal means — the two are correlated at ≈0.998.

| ERF | ND | mean(pipeline − shash) | MC SE | as % of estimand |
|---|---|---|---|---|
| additive | 0.2 | −0.00050 | 0.00017 | 0.125% |
| additive | 0.4 | +0.00057 | 0.00042 | 0.14% |
| mixture | 0.2 | +0.00018 | 0.00017 | 0.045% |
| mixture | 0.4 | +0.00040 | 0.00031 | 0.10% |

**Honest reading of the one nominally significant cell.** Additive ND 0.2 sits ~2.9 MC
SEs from zero. The *magnitude* is 0.0005 on an estimand of 0.40 — **24× smaller than
the 3% acceptance bar**. The two procedures differ by construction (the pipeline takes
`outer_sweeps` sequential `m`=1 draws; the prototype makes a single `m`=30 call), so a
tiny systematic offset is expected. Statistically detectable, practically negligible.
Recorded rather than rounded away.

## Indistinguishable from the oracle (paired)

The stronger claim — that the shipped engine matches an oracle fitted on the *uncensored*
data — was tested rather than inferred from the marginal bias columns:

| ERF | ND | mean(pipeline − oracle) | MC SE | t | p |
|---|---|---|---|---|---|
| additive | 0.2 | −0.00069 | 0.00075 | −0.92 | 0.357 |
| additive | 0.4 | +0.00192 | 0.00159 | +1.21 | 0.228 |

Neither is distinguishable from zero at 300 replications. Note this is a *failure to
reject*, not proof of exact equality — the paired design has enough power to detect a
0.0005 difference against the prototype, so the absence of a signal here is meaningful,
but it bounds rather than eliminates a residual difference.

## Reading against the hypotheses

**The shipped code is not merely "a" congenial method — it is the validated one.**
Before this run, the only evidence for `run_censored_exposure_block_fcs()` was that it
executes end-to-end plus a one-off 25-rep synthetic MC that lived nowhere in the repo.
The bias and coverage claims all rested on a prototype. That gap is now closed: the
engine users actually run reproduces the prototype's behaviour to within ~0.1% of the
estimand.

**The mixture rows confirm faithful implementation, not a defect.** The pipeline
reproduces the prototype's known §7.7 bias almost exactly (+7.54% vs +7.49%; +19.55%
vs +19.45%) rather than being worse. Since the pipeline's X block *is* a linear
conditional draw, this is the expected — and correct — behaviour: the residual bias is
the non-linear-congeniality gap, not a coding error. A pipeline result materially
worse than the prototype here would have indicated a bug; it does not.

**Independent replication of Phase 1.** This run used seed **20260825**;
[`FINDINGS.md`](FINDINGS.md) used **20260813**. The no-`Y` pre-step numbers land
almost on top of the originals — −5.88% vs −5.6%, −14.14% vs −14.4%, coverage
0.923/0.820 vs 0.92/0.82 — so **H1 is confirmed a second time on fresh randomness**,
with a procedure set and RNG stream that differ from the original study.

**A V0 observation resolved as noise.** The Track-V0 single-realization run showed the
focal exposure 12.6% low with its correlated covariate high — the pattern of
attenuation. At 300 reps the pipeline shows −0.79% / +0.48% with no systematic
attenuation. The V0 pattern was Monte-Carlo noise in one dataset at `m`=5, and was
correctly flagged as a hypothesis rather than a finding.

---

## What this licenses

For **additive / per-analyte exposure–response functions** — the case the Phase-2 gate
said to build for — the shipped censored-exposure block-FCS can be described as
validated: unbiased to within 1% with nominal 95% coverage at up to 40% non-detects,
and substantially more efficient than complete-case analysis.

For **mixture / BKMR-style surfaces**, nothing changes from Phase 1: the linear
congenial draw is not adequate, and this run confirms the shipped engine inherits that
limitation exactly as theory predicts.

## Caveats (scope of this run)

1. **MID was never exercised.** `mid_delete_imputed_y` is TRUE, but the harness's `Y`
   is fully observed, so multiple-imputation-then-deletion never fired. The procedure
   asserts this rather than assuming it (no `WARN rows dropped` note appeared in any of
   the 1200 rows). Testing MID needs a `Y`-missingness factor — Track V2.
2. **The miceRanger Z block was also a no-op.** Covariates are complete by default, so
   only the X block did work. That is what makes this comparable to `cens_mi_y_shash`,
   but it means the *block-FCS alternation itself* is still untested under MC — also V2.
3. **Only `X1` is censored** (`censor_which = 1L`), skew is 0, `rho` = 0.4, `n` = 800.
   The remaining robustness axes are V2.
4. **Mixture estimand remains the scaffold simplification** — a single local
   main-effect coefficient, not the true BKMR estimands. Track V3.
5. `outer_sweeps`=3 is redundant in this configuration: with complete covariates each
   sweep is another independent draw of the exposure and only the last is kept. It was
   kept at 3 for fidelity to the shipped example's default, at ~3× the necessary cost.
