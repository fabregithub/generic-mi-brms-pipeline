# Phase 1 harness — confirm the ERF bias and test the joint-model fix

Implements **Phase 1** of `../PLAN_leftcensored_exposure_integration.md` (§7, §10).
A standalone Monte-Carlo study with a **known exposure–response function (ERF)** that
measures the bias each imputation strategy incurs when a **censored exposure is the
focal predictor**. Because everything is simulated, the truth is known and bias is
exact.

It is deliberately built from **procedures that need no new package code**, so it can
run today. It is the **acceptance test** (decision gate, PLAN Phase 2) for whether the
bespoke §4 / block-FCS engine needs to be built at all.

## The question (H1)

> Does the *independent `leftcens` pre-step* (impute the exposure **without** the
> outcome `Y`) bias the focal exposure–response coefficient, and does the
> **congenial `brms` joint model** remove that bias?

The focal estimand is `b_logX1` — the main-effect coefficient of the censored focal
exposure `X1` on the log scale. Every procedure reports **the same estimand** with a
95% interval, so differences are attributable only to missing-data handling.

## Procedures compared (PLAN §7.3)

| # | Procedure | Hypothesis |
|---|---|---|
| 1 | `oracle` — complete data, no missingness | unbiased reference |
| 2 | `complete_case` — listwise-drop censored/missing rows | ~unbiased slope, loses low-dose region |
| 3 | `leftcens_prestep` — `leftcens::gsimp_mi` imputes `X` **without `Y`**, pool by Rubin | **biased** (current architecture) |
| 4 | `cens_mi_y` — congenial censored MI **including `Y`**: left-censored `survreg` draw of `X` conditioned on `Y,Z,other X`, truncated below the LOD, Rubin-pooled | **gold standard / low bias** (verified: tracks oracle to ~1% on the additive DGP) |

**Note on the gold standard.** The original plan named a `brms` joint model as procedure 4
(`bf(Y ~ mi(X)) + bf(X | mi() + cens())`). In the installed `brms`, that encoding creates
no latent per-observation parameters for censored rows, so `mi(X)` silently uses the
LOD-substituted values — i.e. it collapses to LOD-substitution and is biased. `proc_brms_joint`
is kept in the code for reference but is **deprecated in the scaffold and not in the default
set**; `cens_mi_y` is the working congenial reference and a direct preview of the §4 component.

## Layout

```
phase1/
├── R/
│   ├── dgp.R              # known-ERF generator: additive + mixture surface (PLAN §7.1)
│   ├── censoring.R        # left-censor X below a per-analyte LOD; optional MCAR on Z (§7.2)
│   ├── procedures.R       # the four procedures above (§7.3)
│   ├── procedures_pipeline.R  # arms that call the SHIPPED pipeline (V1, V4-V8)
│   ├── robustness.R       # V2 scenario grid + outcome/covariate missingness injectors
│   ├── proper_impute.R    # V4/V5 proper-draw Z-block instruments
│   ├── mice_impute.R      # V6 mice-pmm Z block
│   ├── bart_impute.R      # V7 BART Z block (harness instrument; V8 proved it exact)
│   ├── estimands_bkmr.R   # V3: the seven BKMR estimands + their ANALYTIC truth
│   ├── procedures_bkmr.R  # V3: the four BKMR arms
│   └── metrics.R          # bias / RMSE / coverage vs the known ERF (§7.4)
├── run_phase1.R       # runner: grid × reps, writes results/ (sourceable: run_phase1())
├── run_v2_robustness.R/.sh   # V2
├── run_v4_variance.R         # V4-V8 engine; arms/scenarios via ARMS/SCENARIOS
├── run_v4_variance.sh        # V4/V5
├── run_r8_overnight.sh       # V7
├── run_inner_iter_check.sh   # V8
├── run_v3_bkmr.R/.sh         # V3 (BKMR mixture estimands)
├── test_v3_harness.R         # V3 self-test: 16 checks, no MCMC, ~1 s
├── run_status.sh             # which runs are alive (PID-based; see its header)
└── results/           # generated (git-ignored) — regenerate by running
```

## Is a run still going?

```bash
bash run_status.sh          # every run this directory knows about
bash run_status.sh v3       # just the V3 runs
```

Every runner writes `logs/<tag>_<stamp>.pid`, and `run_status.sh` tests those with
`kill -0` — which sends no signal and only asks whether the process exists.

**Do not check liveness by grepping the process table.** The obvious version is subtly
broken — and broken in a way that survives testing:

```bash
until ! pgrep -f "run_v3_bkmr.R" >/dev/null; do sleep 10; done   # fine alone, hangs in pairs
```

`pgrep -f` matches **full command lines**, and this shell's command line contains the text
`run_v3_bkmr.R`. pgrep does not match its own ancestors, so *one* such loop works — which
is exactly why the bug hides. Two of them are **siblings**, not ancestors: each one's pgrep
finds the *other*, both conditions stay true forever, and both spin until killed. The
failure is silent and looks identical to a job that is still running.

If you genuinely need a pattern, bracket a character so the literal text does not appear in
the watchers' own command lines (`pgrep -f "[r]un_v3_bkmr.R"`). That protects against other
watchers, but not against unrelated processes that mention the same filename — so prefer
the PID file. PIDs are exact; patterns are approximate.

## Track V3 — BKMR mixture estimands

The one track that changes **what is measured** rather than how missingness is handled.
Seven contrasts on the exposure surface, each with truth derived analytically from the
generator rather than from an oracle fit — see `../PLAN_pipeline_validation.md` §8 for
the estimand table and the pre-registered criteria.

```bash
# a short pilot first (recommended)
N_REP=10 SCENARIOS=nd40 ./run_v3_bkmr.sh

# the gate on its own -- cheapest possible check that the harness is sound
ARMS=oracle_bkmr SCENARIOS=nd40 N_REP=44 ./run_v3_bkmr.sh

# the registered run: nd20, nd40, nd40_all x 200 reps  (~12.7 h)
./run_v3_bkmr.sh

# cheaper: drop the mildest cell               (~8.5 h)
SCENARIOS=nd40,nd40_all ./run_v3_bkmr.sh
```

**`nd40_all` censors all three exposures at 40%** and is in the default set. It is the
cell the track most needs: the mixture estimands are functions of all three exposures, and
the two sharp tests (`int_X1X2`, `curv_X1`) depend on X1 and X2 jointly — so this is the
only cell where both exposures entering them are imputed rather than one.

**The gate matters.** `oracle_bkmr` is not the yardstick here — the analytic truth is —
so the oracle's job is to test the harness. If it does not recover the analytic truth
with near-nominal coverage, nothing the other arms report is interpretable. Read it
first, every time.

**Check the machinery first — it takes a second and no MCMC:**

```bash
Rscript test_v3_harness.R
```

16 checks. The load-bearing one is the first: it rebuilds `eta` from
`simulate_complete()`'s own output, strips the covariate and intercept terms, and
requires what is left to equal `h_true()` exactly. If those two ever drift apart, every
bias V3 reports is wrong by the same silent amount and nothing in the output would show
it. The rest cover the contrast algebra, the closed-form cross-check at three
`(mu_x, sd_x)` settings, the additive generator zeroing the non-linear estimands, and
that `int_X1X2` and `curv_X1` respond to `b_int` and `b_quad` respectively and to
nothing else.

**A 50-knot Gaussian predictive process is the default** (`KNOTS=50`), because the full
GP costs 165 s per fit at n = 800 against 16-17 s with knots. The approximation is
*validated by the gate* rather than assumed; `KNOTS=0` fits the full GP if you want to
check it directly.

## Running

From this directory, with `leftcens` (≥ 0.8.0) installed:

```bash
# fast path — the three cheap procedures, both ERF forms, two ND fractions
Rscript run_phase1.R

# include the brms joint model (compiles Stan once per grid cell, then reuses it)
PROCS=oracle,complete_case,leftcens_prestep,brms_joint N_REP=20 Rscript run_phase1.R

# fuller grid
CONFIG=full Rscript run_phase1.R
```

**Environment variables:** `CONFIG` (`quick`|`full`), `N_REP`, `M` (imputations for
the pre-step), `N` (sample size), `PROCS` (comma list), `ERF` (`additive,mixture`),
`ND` (comma list of non-detect fractions), `SEED`.

`brms` is **excluded by default** because each fit compiles/samples Stan and dominates
runtime; add `brms_joint` to `PROCS` explicitly. The compiled Stan model is reused
across replications within a grid cell via `update(recompile = FALSE)`.

Outputs (`results/`): `latest.rds`, a timestamped `.rds`, `phase1_summary.csv`
(bias/coverage table), `phase1_raw.csv` (per-rep estimates).

### Reuse the functions

```r
source("run_phase1.R")                       # defines run_phase1(); does NOT auto-run
res <- run_phase1(config = "quick", n_rep = 100,
                  procs = c("oracle", "complete_case", "leftcens_prestep"))
print_phase1(res$summary)
```

## Reading the result

The summary table is grouped by `erf_form × nd_frac × procedure`, with `bias`,
`rel_bias`, `rmse`, `emp_se`, `mean_ci_w`, and **`coverage`** (95% interval coverage
of the true estimand; a calibrated procedure is ~0.95). The gate reads:

- If `leftcens_prestep` shows material bias / under-coverage that `brms_joint` removes
  → **H1 confirmed**; the joint model is the fix. Decide next (PLAN Phase 2) whether it
  scales, or whether the bespoke block-FCS is needed.
- Expect the pre-step bias to be **larger under the mixture surface** than the additive
  model (PLAN H2), and complete-case to be roughly unbiased for the additive slope but
  to lose the low-dose region (H3).

## Scaffold limitations (intended extensions, not bugs)

These are explicit simplifications to keep Phase 1 to *no new package code*; each is a
known follow-up, flagged in code:

1. ~~**Estimand = focal main-effect coefficient `b_logX1`** for both ERF forms.~~
   **Addressed 2026-08-28 by Track V3** (`R/estimands_bkmr.R`, `run_v3_bkmr.sh`): the
   richer mixture estimands — overall mixture effect, single-exposure effects, the
   X1:X2 interaction and the curvature in X1 — are now implemented with **analytic**
   truth, and BKMR is the analysis model rather than a matched linear one. The scalar
   tracks (V1–V2, V4–V8) still use `b_logX1`, which is what keeps them comparable.
2. **`cens_mi_y` uses a *linear* imputation model**, correctly specified for the additive
   (log-linear-Gaussian) DGP, where it tracks the oracle. For the **mixture surface** the
   linear model is misspecified, so `cens_mi_y` is itself biased there (empirically it is —
   the run surfaces this) — which is exactly the **§7.7 non-linear-congeniality** gap the
   substantive-model-compatible §4 approach is meant to close. It carries a
   `linear-approx (§7.7)` note in the mixture rows. The deprecated `brms_joint` (see above)
   is additive-oriented and not part of the default set.
3. **Only the focal exposure `X1` is censored** (`censor_which = 1L`), to isolate the
   mechanism. `inject_left_censoring(..., censor_which = seq_len(p))` censors all.
4. **MCAR on covariates defaults off** (`inject_mcar_covariates`, `mcar_frac = 0`). With
   MCAR on, the pre-step needs `miceRanger` for `Z` and the joint model needs `mi()` on
   `Z` — wired as a hook, not yet exercised.
5. Procedures 5 (block-FCS) and 6-as-primary are **Phase 3–4**, added only if the gate
   says build.

## Relationship to `leftcens/validation/`

Mirrors that study's conventions (env-var config, sourceable functions, git-ignored
`results/`). The difference: `leftcens/validation` validates *marginal* recovery of
censored analytes; this validates the *joint* exposure–response functional — the thing
the marginal pre-step is hypothesised to get wrong.
