# Validation plan: the *shipped* censored-exposure pipeline

**Status:** draft, 2026-08-25. Companion to
[`PLAN_leftcensored_exposure_integration.md`](PLAN_leftcensored_exposure_integration.md)
(the design/theory document) and [`phase1/FINDINGS.md`](phase1/FINDINGS.md) (the
Monte-Carlo results that motivated the build).

---

## 1. Why a separate plan

The existing plan is a *design* document: it argued the case for congenial
imputation of a censored exposure, gated the build behind a Monte-Carlo study, and
tracked the build itself through Phases 1–4. That plan has done its job — the
mechanism is demonstrated and the engine is written.

This plan covers what is left, which is a different question:

> The bias numbers we have validate a **prototype**. Do they hold for the **code we
> ship**?

Every headline result in `FINDINGS.md` comes from `cens_mi_y_shash`, a ~40-line
reference procedure inside the MC harness. The pipeline path that users will actually
run — `run_censored_exposure_block_fcs()` in
[`00_censored_exposure.R`](../00_censored_exposure.R) — has been shown to *run*
(end-to-end example, Steps 1→4) and was spot-checked by a one-off 25-replication
synthetic MC that lives nowhere in the repo. That is not the same evidence, and the
gap is the main thing this plan closes.

### Scope boundary

| In scope | Out of scope |
|---|---|
| Validating shipped pipeline code against known truth | Re-litigating the design (settled in the other plan) |
| Making the mixture verdict rest on real estimands | Building a substantive-model-compatible mixture engine |
| Documenting the required `m` | Production runs on real cohort data |
| Robustness axes not yet swept | `leftcens` package-internal validation (lives in that repo) |

---

## 2. Status snapshot

**Validated (prototype-level), from `FINDINGS.md` — 300 reps/cell, `n`=800, `M`=30:**

- H1 confirmed on the additive DGP: the no-`Y` pre-step attenuates the focal
  coefficient (−5.6% at 20% ND → −14.4% at 40% ND; coverage 0.92 → 0.82), and the
  congenial `Y`-aware method removes it (bias ≤0.4%, coverage ~0.96).
- H3 confirmed: complete-case ~unbiased for the slope, inefficient (rmse up to ~4.5×
  oracle).
- §7.7 gap material for mixtures: linear congenial imputation is itself biased on the
  mixture surface (+7.7% → +19.6%, coverage → 0.62).
- §7.5 skew: the reference must be skew-aware; a Gaussian/tobit `Y`-aware reference
  breaks under right-skew (coverage → 0.47), the shash margin restores it.
- Scale resolved: the X-block draw is ~`k²` in predictor width, not `n`; ~1.6–2.4 min
  at `k` ≤ 50 for 20 datasets × 5 sweeps at `n` = 80k.

**Not validated:**

| Gap | Consequence if left open |
|---|---|
| ~~Shipped `run_censored_exposure_block_fcs()` never MC-tested~~ | **CLOSED 2026-08-25 (V1)** — see `phase1/FINDINGS_v1.md` |
| ~~Steps 5–11 + Quarto never run on the censored path~~ | **CLOSED 2026-08-25 (V0)** |
| Mixture verdict rests on a single local coefficient | The §7.7 claim is not manuscript-final |
| No FMI / pilot-`m` figure for this path | Users have no principled `m` |
| ~~Robustness axes unswept~~ | **CLOSED 2026-08-25 (V2)** — 9/10 pass; one under-coverage finding open |
| ~~NEW: interval under-coverage under heavy covariate + outcome missingness~~ | **ATTRIBUTED 2026-08-25 (V4)** — improper MI in the miceRanger Z block understates `B`; affects the general imputation path, not just this feature |
| **OPEN: the Z block needs a proper draw** | Until fixed, intervals are too narrow wherever much is imputed. The V4 instrument is not the fix (overshoots) |

---

## 3. Working agreement

Per the arrangement for this project:

- **Claude writes the scripts.** Each track below produces files that are reviewed
  before anything is run.
- **The user runs them** on this machine, and pastes back the summary table or points
  at the results file.
- **We evaluate together** against the acceptance criteria stated per track — written
  *before* the run, so the read-out is not post-hoc.
- Findings land in `FINDINGS.md` (or a per-track findings file) with the same
  conventions as Phase 1: env-var config, sourceable functions, git-ignored `results/`.

Nothing in this plan should silently change pipeline behaviour. Where a track
uncovers a bug, the fix is a separate, reviewed change.

---

## 4. Environment (verified 2026-08-25)

```
leftcens    0.9.0   (tag v0.9.0 -> b3ebc8c, pushed)
brms        2.23.0
cmdstanr    0.9.0.9000  (CmdStan 2.38.0)
miceRanger  1.5.0
survival    3.8.11
future      1.75.0 / furrr 0.4.0
bkmr        NOT INSTALLED   <- required by track V3 only
cores       24
```

`NCORES=14` is the Phase-1 default and leaves headroom; with 24 cores, `NCORES=20`
is reasonable for the fork-parallel tracks. Keep BLAS single-threaded
(`run_phase1.sh` already exports the four `*_NUM_THREADS=1` variables) so it does not
fight the R-level loop.

---

## 5. Track V0 — Full-pipeline regression on the censored example

> **STATUS: PASSED — 2026-08-25**, ~78 s. All five criteria met: clean exit, populated
> `parameter_summary.csv`, Step 11 outputs, report rendered (HTML + DOCX), Step 12
> skipped cleanly. Steps 09/10 correctly omitted (no `mo()` terms), MID dropped 30
> imputed-`Y` rows per dataset, 0 divergences and 0 treedepth hits across 5 fits.
> Confirmed the `12_export_draws.R` harness fix — Step 12 now runs at all.

**Question:** do Steps 5–11 and the Quarto report work on the censored-exposure path?

Steps 1→4 are verified; 5–11 are generic but have never been exercised here. A
harness bug that aborted *every* `test/` run at Step 12 was fixed on 2026-08-25
(`12_export_draws.R` was missing from the copied file list), so this is the first run
that can complete.

**New code:** none. The script already exists.

**Run:**

```bash
bash test/test_censored_exposure_quick.sh
```

**Expected:** ~15–30 min at quick settings (`m`=5, chains=1). Output lands in
`test/runs/<timestamp>_censored_exposure/project/`.

**Acceptance criteria:**

1. Exit status 0 and no `pipeline_error.flag`.
2. `results/parameter_summary.csv` exists and contains `log(expo1)` / `log(expo2)`
   rows with finite estimates and CIs.
3. Step 11 stability outputs exist.
4. The combined Quarto report renders (HTML + DOCX) and its stability chapter is
   populated.
5. Step 12 is *skipped cleanly* (`export$cohort_id` is NULL) rather than erroring.

**What we evaluate together:** any step that errors or produces empty output; whether
the `mo()` chapter is correctly omitted; whether the pooled coefficients are sane
against the example's known generating values.

**Gate:** V0 must pass before V1 is worth running — V1 reuses the same imputation
engine, and a Steps 5–11 failure here is a pipeline bug, not a statistical one.

---

## 6. Track V1 — Procedure 5: the shipped engine in the MC harness

> **STATUS: PASSED — 2026-08-25.** 300 reps/cell, `m`=30, `n`=800, seed 20260825,
> ~17 min. Every pre-registered additive criterion met: rel. bias −0.79% (20% ND) and
> +0.48% (40% ND), coverage 0.957 at both, rmse 1.06×/1.24× oracle, and 0.0485 vs
> complete-case's 0.0719 at 40% ND. Paired agreement with the `cens_mi_y_shash`
> prototype is ≤0.14% of the estimand in every cell. Mixture cells reproduce the
> prototype's known §7.7 bias (+7.54% / +19.55%) rather than exceeding it —
> faithful implementation, not a defect. Also an independent replication of Phase 1's
> H1 on a fresh seed. Full table and reading:
> [`phase1/FINDINGS_v1.md`](phase1/FINDINGS_v1.md).
>
> *Correction applied during implementation:* the "force `n_cores = 1`" note below was
> based on a misreading. `run_phase1()`'s replication loop is **serial**, so there is
> no outer fork to nest inside; `n_cores` is passed through to the module's `mclapply`,
> which also exercises the shipped parallel path.

**The priority track.** Everything else is refinement; this one converts "the
prototype is unbiased" into "the shipped code is unbiased".

**Question:** does `run_censored_exposure_block_fcs()` recover the known
exposure–response coefficient with nominal coverage, on the same DGP, censoring, and
metrics as the existing four procedures?

### New files

```
validation/phase1/R/procedures_pipeline.R    # proc_pipeline_block_fcs() + adapters
validation/phase1/run_v1_pipeline.sh         # detached runner (mirrors run_phase1.sh)
```

plus small edits to `run_phase1.R` (registry) and `R/metrics.R` (see gotchas).

### Design

`proc_pipeline_block_fcs(bundle, m, seed)` must:

1. Adapt the harness bundle into the pipeline's input convention.
2. Call the **real** `run_censored_exposure_block_fcs()` — imported by sourcing
   `00_censored_exposure.R`, not copy-pasted.
3. Fit the matched ERF model to each of the `m` completed datasets with the harness's
   existing `fit_lm_estimand()`.
4. Pool with the harness's existing `rubin_pool()` and return `one_row(...)`.

Steps 3–4 deliberately reuse the harness helpers so the comparison against
`oracle` / `complete_case` / `leftcens_prestep` / `cens_mi_y_shash` is like-for-like:
the *only* thing that differs across procedures is the missing-data handling.

### Implementation gotchas (found while scoping — do not rediscover these)

1. **Scale mismatch.** The harness works in log space: `logX1` holds the true value
   and is `NA` on censored rows, with `logX1_lod` carrying the log LOD. The pipeline
   expects a **data-scale** exposure column plus `X_lo`/`X_hi` bounds. Adapter:

   | row | `X1` | `X1_lo` | `X1_hi` |
   |---|---|---|---|
   | observed | `exp(logX1)` | `exp(logX1)` | `exp(logX1)` |
   | censored | `NA` | `0` | `exp(logX1_lod)` |

   Set `censored_exposure$log_scale = TRUE` so the X block imputes on the log scale
   (where the shash margin lives), matching `cens_mi_y_shash`.

2. **`00_common_functions.R` is not standalone-sourceable.** It fails with
   `object 'paths' not found` — [line 1334](../00_common_functions.R:1334) reads
   `paths$objects` at top level for the runtime-override mechanism. Define stubs
   *before* sourcing:

   ```r
   paths <- list(objects = tempdir())
   analysis_spec <- list()
   source("00_common_functions.R")
   source("00_censored_exposure.R")
   ```

   (`log_msg()` works without `init_logging()`.)

3. **Nested forking.** `run_censored_exposure_block_fcs()` parallelises the `m`
   datasets with `mclapply`, and the harness already forks across replications via
   `NCORES`. Force `censored_exposure$n_cores = 1L` inside the procedure and
   parallelise only at the harness level — the same reasoning the module itself
   applies to its Z block.

4. **`metrics.R` will silently drop the new rows.** `summarise_phase1()` hardcodes
   `proc_order` as factor levels; a procedure name not in that vector becomes `NA`
   and sorts out of the table. Add `"pipeline_block_fcs"` to
   [`proc_order`](phase1/R/metrics.R:40).

5. **MID is a no-op as configured.** The harness's `Y` is fully observed, so
   `mid_delete_imputed_y` never fires. Exercising MID needs a `Y`-missingness factor —
   deferred to V4, and the procedure should assert MID did nothing rather than
   pretending it was tested.

6. **Synthetic `var_dict`.** The module reads `use_in_model`, `impute_target`, and
   `scale` from a dictionary data frame. Build a minimal one in the adapter with the
   exposures marked `impute_target = FALSE, use_in_model = TRUE`, matching what
   `01_validate_config.R` now enforces.

### Run

```bash
CONFIG=full NCORES=20 PROCS=oracle,complete_case,leftcens_prestep,cens_mi_y_shash,pipeline_block_fcs bash validation/phase1/run_v1_pipeline.sh
```

Start with a smoke pass (`CONFIG=quick N_REP=20`) to confirm the adapter before
committing to the full grid.

**Expected runtime:** the block-FCS runs `outer_sweeps` × (miceRanger + X-block) per
completed dataset, so it is the most expensive non-Stan procedure. Budget
significantly more than the ~5–8 min of the existing no-brms grid; the smoke pass
will give a per-rep timing to extrapolate from before the full run is launched.

### Acceptance criteria

Stated in advance. Monte-Carlo error at 300 reps with `emp_se` ≈ 0.04 on a true value
of 0.40 gives an MC standard error of ~0.6% on relative bias and ~±0.025 on coverage,
so these thresholds are ~5 MC SEs wide and are not measuring noise.

**Additive DGP (the case the gate said to build for) — this is a pass/fail:**

| Metric | Criterion |
|---|---|
| relative bias, ND 0.2 and 0.4 | \|rel_bias\| ≤ 3% |
| 95% interval coverage | 0.92 ≤ coverage ≤ 0.97 |
| rmse vs oracle | ≤ 1.3 × oracle rmse |
| rmse vs complete-case | strictly better at 40% ND |
| agreement with the prototype | within MC error of `cens_mi_y_shash` |
| `n_ok` | equals `n_rep` (no silent failures) |

**Mixture DGP — diagnostic, not pass/fail.** Phase 1 already established that *any*
linear congenial draw is biased here (+7.7% → +19.6%). The pipeline engine is a linear
draw, so it is expected to reproduce roughly that bias. The check is that it **tracks
`cens_mi_y_shash`** rather than being worse: agreement confirms the engine is a
faithful implementation, and the residual bias is the known §7.7 gap, not a coding
error. A pipeline result materially *worse* than the prototype on the mixture surface
would indicate an implementation bug.

**What we evaluate together:** the summary table row-by-row; any divergence from
`cens_mi_y_shash` beyond MC error; `n_ok` shortfalls and their `note` column; whether
the additive criteria are met at both ND levels.

**Outcome:** either a defensible bias/coverage claim for the shipped engine (which
goes into `README.md` and the manuscript), or a located bug.

---

## 7. Track V2 — Robustness sweep

> **STATUS: 9/10 PASSED, one finding — 2026-08-25.** 2750 tasks, 2.95 h, 63.7 CPU-h,
> zero failures, 98% worker efficiency. Bias passes everywhere (max 2.71%, bar 3%),
> including all exposures censored, n = 80 000, skew 0.75 and rho 0.8. **MID and the
> miceRanger Z block are now genuinely exercised** — both were dead code in V1.
> **Finding:** `combined` (all exposures censored + 20% MCAR covariates + 20% missing
> outcomes) has coverage **0.893**, below the 0.90 investigate line. It is a variance
> problem, not bias: interval width / empirical SE correlates 0.837 with coverage and
> falls to 0.84 exactly where the Z block does the most work, while the oracle stays
> calibrated (0.95–1.07) throughout. Leading hypothesis — miceRanger's RF/PMM draw is
> improper MI — is **not established**; a discriminating `m`=100 test is running.
> Full reading: [`phase1/FINDINGS_v2.md`](phase1/FINDINGS_v2.md).
>
> *Bug fixed during setup:* per-task seeds keyed on the scenario's index in the
> **filtered** list, so `SCENARIOS=` subsets silently generated different data than a
> full run. Now keyed on the canonical unfiltered position, which keeps subsets
> comparable to full runs and leaves the completed run reproducible.

**Question:** does the V1 result survive outside the single tested corner?

The Phase-1 caveats list four unswept axes. Each is a small change to the existing
generator/injector, and they can share one runner.

| Axis | Change | Why it matters |
|---|---|---|
| Censor **all** exposures | `inject_left_censoring(..., censor_which = seq_len(p))` | Real cohorts censor every analyte; the X block loops over exposures and that loop is untested under MC |
| MCAR covariates | `inject_mcar_covariates(..., mcar_frac = 0.2)` | Exercises the miceRanger Z block *and* the block-FCS alternation, currently a hook that has never been driven |
| Missing `Y` | new: MCAR a fraction of `Y` | The only way to exercise MID; asserts the imputed-`Y` rows are dropped and pooling is unaffected |
| Exposure correlation | vary `rho` (0.0 / 0.4 / 0.8) | High correlation is where the reduced predictor set could bite |
| Large `n` | one cell at `n` = 20k–80k, few reps | Confirms the scale test's timing at MC level |

**New files:** `validation/phase1/run_v2_robustness.R` + `.sh`, reusing V1's
procedure.

**Acceptance criteria:** the V1 additive criteria hold in every cell, with two
explicitly allowed degradations — coverage may drop toward 0.92 under heavy MCAR, and
rmse may rise with `rho`. Any cell where **bias** exceeds 3% or coverage falls below
0.90 is a finding to investigate, not a threshold to relax.

**Gate:** run only after V1 passes. A robustness sweep of a broken engine wastes core-hours.

---

## 8. Track V3 — Phase 2b: BKMR estimands

**Question:** does the §7.7 mixture conclusion hold for the estimands that actually
matter, rather than the scaffold's single local coefficient?

This is the one item the existing plan explicitly flags as needed before the mixture
verdict is manuscript-final. It changes *what is measured*, not *what is compared*.

**Prerequisite:** `install.packages("bkmr")` — not currently installed.

**New files:**

```
validation/phase1/R/estimands_mixture.R   # the four estimands below
validation/phase1/run_v3_estimands.R/.sh
```

**Estimands** (PLAN §7.4):

1. Overall mixture effect: all exposures moved q25 → q75 jointly.
2. Single-exposure effects with the others held at their medians.
3. Pairwise interactions.
4. `h` evaluated at representative exposure profiles.

**Design notes:**

- Each estimand needs a *known truth* computed analytically from `make_truth()`'s
  mixture surface, not estimated from an oracle fit — otherwise "bias" is contaminated
  by the reference's own error. This is the main piece of work in the track.
- `metrics.R` currently assumes one scalar estimand per row (`estimand_true`). Either
  add an `estimand` column and group by it, or keep a parallel summariser. Prefer the
  former so V1/V2 results stay comparable.
- BKMR is expensive. Expect a reduced grid: fewer reps, smaller `n`, and possibly only
  the 40% ND cell where the effect is largest.

**Acceptance criteria:** this track is *descriptive*, not pass/fail. It resolves a
research question. The deliverable is a table of bias/coverage per estimand ×
procedure, and an explicit written verdict on whether the linear congenial draw is
adequate for BKMR-style mixtures — with the answer allowed to be either yes or no.

**What we evaluate together:** whether the mixture bias seen on the scaffold estimand
persists, grows, or vanishes on the real estimands; and whether the ranking of
procedures changes.

---

## 9. Track V4 — Variance attribution (and pilot `m`)

> **STATUS: RESOLVED — 2026-08-25.** 5 scenarios × 4 arms × 300 reps, 1500 tasks, 5.25 h,
> zero failures. **The V2 under-coverage is caused by improper multiple imputation in the
> miceRanger Z block, understating the between-imputation variance `B` by 8–59%.**
> Swapping in a proper Bayesian draw restores calibration everywhere the Z block works
> (width/SE 0.899 → 1.085 in `mcar_z40`; 0.843 → 0.973 in `combined`). `Ubar` moves far
> less than `B`, which is the improper-MI signature. **MID is exonerated** — disabling it
> makes calibration *worse* and adds −6 to −7% bias, so it stays. The V2 two-mechanism
> account is **withdrawn**: one mechanism explains all scenarios. The fix is NOT to adopt
> the validation instrument (it overshoots to 1.085 and its bias is marginally worse) but
> to make the Z block propagate parameter uncertainty properly.
> Full reading: [`phase1/FINDINGS_v4.md`](phase1/FINDINGS_v4.md).
>
> *Scope note:* the Z block is the pipeline's **general** row-level imputation path, so
> this is not a censored-exposure finding — it applies to any analysis with substantial
> covariate or outcome missingness.
>
> *Operational note:* `resume` in `run_v4_variance.R` is **off by default and not
> trusted** — resuming from an accumulated checkpoint reproducibly killed the parent
> process after the first chunk, in both detached and foreground runs, while a fresh
> checkpoint ran fine. Cause unknown. Partial results are recoverable without resume via
> `summarise_v4(readRDS("results/v4_checkpoint.rds"))`.

**Question:** how many imputations does this path actually need?

**New file:** `validation/phase1/run_v4_fmi.R`

**Design:** `rubin_pool()` already returns the Barnard–Rubin `df`; FMI follows from
the between/within variance decomposition it computes. Run the pipeline engine at a
large `m` (say 50) on a few replications and report FMI for the focal coefficient,
then the `m` required for stable interval width per PLAN §6's two-`m` recipe.

Cross-check against the pipeline's own `mi_stability` auto-increment loop: if the loop
stops at an `m` well below the FMI-implied figure, that is a finding about the loop's
stopping rule.

**Acceptance criteria:** a documented FMI figure and a recommended default `m` for the
censored-exposure path, written into `README.md` and `00_config.R`'s commented block.

---

## 10. Sequencing

```
V0  regression ✅     ── gate ──▶  V1  procedure 5 ✅     ── gate ──▶  V2  robustness ✅
    PASSED 2026-08-25                  PASSED 2026-08-25                9/10 + 1 finding
                                           │                               │
                                           └──────────▶  V4 ✅  ◀──────────┘
                                                    RESOLVED 2026-08-25
                                                              
V3  BKMR estimands  ── independent; needs bkmr; can run any time after V0
```

**All four tracks are complete.** What remains is not validation but a fix and one open
research question:

1. **Make the Z block a proper draw** (engineering, in `00_common_functions.R`). This is
   now the highest-value change in the repo: it affects every analysis with substantial
   missingness, and V4 quantifies both the problem and the target. Needs its own
   validation run afterwards — the V4 harness measures it directly.
2. **Pilot `m` / FMI** — the original V4 scope, still undone. FMI is now reported per
   arm (0.26–0.44 across scenarios), so this is a short step from the existing results.
3. **V3 (BKMR estimands)** — independent, still needs `bkmr` installed. The only piece
   left before the mixture verdict is manuscript-final.

- **V0 → V1** is a hard gate: statistical validation of a pipeline whose reporting
  steps are broken is premature.
- **V1 → V2, V4** is a hard gate: both build on V1's procedure.
- **V3** is independent of V1 — it asks a different question (which estimand) and
  touches different files. It can be scheduled whenever `bkmr` is available.

## 10b. RNG conventions (learned the hard way)

Two properties of the harness's randomness that have already caused confusion:

1. **Task seeds are keyed on the scenario's position in the CANONICAL (unfiltered) grid**,
   not in whatever subset `SCENARIOS=` selected. Otherwise running one scenario alone
   generates different data than the same scenario in a full run, and the two cannot be
   compared. Verified by the deterministic arms reproducing to `max |diff| = 0`.
2. **An arm's numbers are only comparable across runs with an identical arm set.**
   Procedures consume the task's RNG stream *sequentially*, so adding or removing an arm
   changes the stream reaching every later arm. The *datasets* are identical; the
   imputation draws are not. `properBoot`'s `combined` coverage read 0.910 in V5 and 0.927
   in V6 for exactly this reason (~0.9 SE apart, i.e. noise). When comparing one arm across
   runs, hold the arm set fixed — or compare within a run, which is what every verdict here
   does.

## 11. Reporting conventions

Mirror Phase 1:

- Config via environment variables; runners are sourceable and do **not** auto-run
  when sourced (`sys.nframe() == 0L` guard).
- `results/` git-ignored; per-cell checkpointing so a late failure never loses
  finished work.
- Each track writes `results/<track>_summary.csv` + `<track>_raw.csv` + `latest.rds`.
- Findings go in a markdown file per track, in the `FINDINGS.md` house style: a table,
  a reading against the stated criteria, and an explicit caveats section naming what
  the run does *not* establish.
- Every acceptance criterion in this plan is written before its run. If a criterion
  turns out to be wrong, change it explicitly and say why — do not quietly re-read the
  result against a softer bar.
