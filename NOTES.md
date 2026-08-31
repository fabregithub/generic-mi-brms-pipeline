# Worklog

## Session summary

- Fixed incorrect multiple-imputation posterior pooling in `06_posterior_summary.R`: replaced naive
  draw-stacking with weighted pooling (`1/(m*K_i)` per draw) plus a conditional finite-`m` Rubin's-rule
  variance correction, applied only when a bimodality diagnostic confirms the pooled posterior is
  unimodal enough for the correction's symmetry assumption. Added audit columns
  (`m_imputations`, `between_var`, `within_var`, `variance_corrected`, `transform_used`,
  `bimodality_coef`) to `parameter_summary.csv`.
- Made `09_check_mo_parameter_columns.R` and `10_publication_mo_results.R` generic: they now discover
  `mo()` variables directly from the fitted model's formula instead of hardcoded names, with
  anchored regex matching that also handles brms's `idEQ<id>` naming infix when `mo()` is called with
  an explicit `id =` argument.
- Added an automatic m-increment stability loop (`analysis_spec$mi_stability$auto_increment`, now
  default `TRUE`) to `run_all.R`: fits imputations in `fit_workers`-aligned batches, stops increasing
  `m` once posterior summaries stabilise, and persists per-batch config overrides via
  `objects/mi_runtime_override.rds` (read by `00_common_functions.R`) since per-example configs
  replace `00_config.R` wholesale. The one-off smoke fit now runs only once per `run_all.R` invocation,
  not once per batch.
- Extended `03_impute.R`/`00_common_functions.R` to support safely extending an existing imputation
  run to a larger `m` (`allow_extend`), with deterministic batch seeding.
- Combined the previously separate `mi_stability` and `mo_effects` reports into the main Quarto
  report: `08_publication_results.R` now embeds Step 11's stability tables/figures and Step 9/10's
  `mo()` odds-ratio tables/figures as report chapters, and self-renders the report (HTML + DOCX)
  automatically. The `mo()` chapter is fully omitted (not just empty) for models without `mo()` terms.
  `run_all.R` was reordered so Steps 11 and 09/10 run before Step 8.
- Updated `00_config.R` with visible, documented placeholders for the new optional settings
  (`imputation$seed`/`allow_extend`, `mi_stability`, `mo_effects`) and flipped
  `mi_stability$auto_increment` to default `TRUE`.
- Rewrote large parts of `README.md`: new "Manuscript writing guide" section, documentation of the
  pooling math, the auto-increment loop, the combined report, and converted ~70 enumerated
  ```text``` blocks into proper Markdown bullets/numbered lists/blockquotes. Updated the four
  translated docs (`docs/README.*.md`) to match.
- Verified throughout via the `test/` example suite (quick + parallel, all three bundled examples)
  and the user's own real data (logistic model with repeated-Y measurements, both with and without
  `mo()` terms).
- Tagged `v0.9.5`; user has decided to hold `v1.0.0` until other users test the pipeline against their
  own (more varied) data and analysis patterns.
- Added `CLAUDE.md` documenting commands and architecture for future Claude Code sessions.

## Next steps (after v0.9.5, before v1.0.0)

- Recruit other users to test the pipeline on their own data/analysis patterns (different families,
  data structures, missingness mechanisms) before tagging `v1.0.0`.
- Fold any new analysis pattern that trips up an external tester into a fourth bundled example, so it
  becomes part of permanent regression coverage rather than a one-off fix.
- Revisit whether to package this as an installable R package later, once it's clear whether user
  friction is about workflow clarity (packaging won't help) or wanting a callable API without cloning
  the repo (packaging would help). Decided to keep it as a script-based template for now.

---

## Session: 2026-07-16 → 2026-07-17 (v1.0.0 → v1.1.0)

### Bug fixes

- **Cox PH prior detection** (`00_common_functions.R`): `make_default_priors()` was incorrectly
  detecting `time | cens(censored)` as a random-effects formula and adding an `sd` prior.
  Fixed with a regex that excludes brms response modifiers (`cens`, `trunc`, `mi`, etc.) from
  the `|` check.
- **Empty brmsprior crash** (`filter_priors_to_model()`): when all priors were filtered out, an
  empty 0-row `brmsprior` object was passed to `brm()`, causing `undefined columns selected`.
  Fixed by returning `NULL` instead of an empty object; wrapped `get_prior()` in `tryCatch`.
- **Cox censoring variable** (`examples/lung_cox/00_variable_dictionary_lung_cox.csv`): the
  `censored` column had `type = "binary"`, which triggers `as.factor()` in the data-prep step.
  brms `cens()` requires numeric 0/1. Fixed by changing to `type = "integer"`.
  General rule: survival time and censoring indicator columns must use `type = "integer"` or
  `type = "continuous"`, never `type = "binary"`.
- **pd display scale**: pd stored as 0–1 proportion; template sentences now multiply by 100 for display.

### New features (v1.0.0)

- **Cox proportional hazards model** (`family = "cox"`):
  - Formula syntax: `time | cens(censored) ~ ...`
  - brms convention: 0 = event (not censored), 1 = right-censored
  - Requires `survival` and `splines2` R packages
  - Effect scale: log hazard ratio; back-transformed to HR in publication outputs
  - Example 4: `examples/lung_cox/` using `survival::lung`
  - Test scripts: `test/test_lung_cox_quick.sh`, `test/test_lung_cox_parallel.sh`

- **Publication template sentences** (Step 8, `08_publication_results.R`):
  - bayestestR neutral format: "The effect of X has a probability of pd% of being [direction]
    (Median = ..., CI [...]) and can be considered as [label] (X% in ROPE)"
  - Covers fixed effects and `bsp_mo*` (monotonic main coefficients)
  - Output: `results/publication/tables/parameter_template_sentences.csv`
  - `analysis_spec$publication$template_sentences_scope`: `"all"` (default) or `"exposure_only"`
  - `"exposure_only"`: generates sentences only for `role == "exposure"` variables and their
    interaction terms; useful for causal inference analyses with a common DAG

- **Step 12 — federated draw export** (`12_export_draws.R`):
  - Long-format export: `cohort_id`, `parameter`, `imputation`, `draw_index`, `value`
  - Outputs: `results/export/cohort_draws.rds`, `cohort_draws.csv`, `cohort_metadata.json`
  - `analysis_spec$export$cohort_id` (NULL = skip step), `scope` (`"exposure_only"` default)
  - Companion repo: `generic-mi-brms-meta` (see its own NOTES.md)
  - Key design principle: cohorts share `00_variable_dictionary.csv` template from coordinating
    site → parameter names align automatically → no post-hoc renaming needed

- **`EXPORT_FORMAT.md`**: documents the export file schema, link-scale note, scope behaviour,
  and a pre-transfer checklist for cohorts

### New features (v1.1.0)

- **`launch.R`**: RStudio interactive menu launcher — open and click Source, no terminal needed.
  Designed for Windows users. Menu covers full pipeline, individual steps, clean all, clean fits.
- **`launch.sh`**: Bash equivalent for Mac/Linux (`bash launch.sh` or double-click after `chmod +x`)

### New features (v1.2.0)

- **`fit_single_imputation.R` added to both launchers** (menu option 6): prompts for an
  imputation number and calls `fit_single_imputation.R <n>`; useful for debugging a specific
  failed imputation without the terminal
- **README clarified**: launchers (`launch.R`, `launch.sh`) are the primary entry points;
  `run_all.R` described as the internal orchestrator that users should not call directly;
  `launch.sh` documented with `chmod +x` instructions alongside `launch.R`
- **Note on `launch.R` testing**: `readline()` blocks in non-interactive R sessions (e.g. piped
  input); `launch.R` must be sourced interactively in RStudio — only `launch.sh` is
  non-interactively testable

### New features (v1.3.0)

- **Launcher menu reordered**: validate config (option 1) now appears before run full pipeline
  (option 2) in both `launch.R` and `launch.sh`; encourages users to validate before committing
  to a full run
- **README updated**: `launch.sh` command-line section now shows `bash launch.sh` as the primary
  entry point (previously showed bare `Rscript` commands)

---

## Session: 2026-08-14 — Censored-exposure imputation (block-FCS)

### New feature: congenial imputation of censored exposures

- **`00_censored_exposure.R`** — new opt-in Step-3 strategy
  `strategy = "censored_exposure_block_fcs"` for the case where a focal *exposure*
  is left-/interval-censored below a reporting limit. Imputing such an exposure
  without the outcome biases the exposure–response coefficient (congeniality); this
  path imputes it **outcome-aware, skew-robust, and interval-respecting** via a
  two-engine block-FCS: `miceRanger` for MAR covariates (the reused
  `run_row_level_imputation()` helper) + `leftcens::impute_censored_conditional()`
  for the censored exposures. Handles multiple exposures and the three-tier
  (ND / DNQ / quantified) interval structure.
- **Dependency**: `leftcens` (>= 0.9.0), which was extended (in that package) to
  export `impute_censored_conditional()` plus the sinh-arcsinh margin API
  (`fit_shash_margin`, `draw_margin`, `x_to_z`/`z_to_x`).
- **Input convention**: each exposure `X` carries `X_lo`/`X_hi` interval columns
  (leftcens form). Dictionary: `role = exposure`, `impute_target = FALSE`,
  `use_in_model = TRUE`.
- **MID**: missing outcomes are imputed for the FCS predictor step, then the
  imputed-Y rows are deleted before the fit (von Hippel 2007). Asserted at runtime.
- **Dispatch**: one `else if` branch in `03_impute.R`; returns the standard
  `imputed_list`, so the file/manifest logic and Steps 4–12 are unchanged. Config
  documented in `00_config.R`.
- **Parallelism**: the `m` completed datasets are imputed in parallel
  (`censored_exposure$n_cores`, fork/`mclapply`); the miceRanger Z-block runs
  serially within each dataset to avoid nested PSOCK-in-fork failures.
- **Example**: `examples/censored_exposure/` (two exposures — one left-censored,
  one three-tier interval; logistic-agnostic Gaussian demo), registered in the
  shell test harness (`test/test_censored_exposure_quick.sh`). Verified end-to-end
  (Steps 1→4) and via a synthetic MC (recovers the exposure–response coefficient;
  coverage ~0.95).

### Validation study (`validation/`)

> **Status index:** [`validation/README.md`](validation/README.md) is the single place
> showing what is validated, what is open and what is next. Keep it current — it is
> the evidence base for the claims the README makes.

- `validation/PLAN_leftcensored_exposure_integration.md` — the design/decision
  document (congeniality theory, the two-engine block-FCS, the estimand/scale
  gates, the joint-model reference).
- `validation/phase1/` — a Monte-Carlo harness quantifying: (i) the no-Y pre-step
  bias (**H1** confirmed — attenuation growing with censoring), (ii) the
  §7.7 non-linear-congeniality gap for mixtures, and (iii) §7.5 skew robustness
  (a Gaussian/tobit reference breaks under right-skew; the shash margin fixes it).
  Also a scale test (X-block cost is ~`k^2`, cheap at n=80k with a reduced
  predictor set).

### Note

- The JECS PFAS→Kawasaki-disease **manuscript demonstration** (a synthetic cohort
  grounded in Lai 2025 / Iwata 2024 / Atagi 2024, the 4-method comparison, BKMR,
  and the manuscript skeleton) has been **moved out of this repo** to keep the
  pipeline general. It lives in the manuscript working directory. The generic
  censored-exposure feature and `examples/censored_exposure/` remain here.

---

## Session: 2026-08-25 — leftcens dependency + Step-1 validation for the censored path

### `leftcens` 0.9.0 was untagged

The pipeline documents a hard requirement on `leftcens >= 0.9.0` (for
`impute_censored_conditional()`), but the `leftcens` repo's newest tag was `v0.8.0`.
The 0.9.0 release was fully committed *and pushed* to `origin/main` — only the release
tag was missing, so there was no immutable ref to pin and no way to satisfy the
documented requirement from a tag. Tagged `v0.9.0` (at `b3ebc8c`) and pushed; the
locally installed 0.8.0 (built 2026-08-06, predating the new API) was upgraded.

### Fail-early validation for `censored_exposure_block_fcs` (`01_validate_config.R`)

Step 1 previously validated only the `subject_wide_with_repeated_y_auxiliary`
strategy, so a misconfigured censored-exposure run got through Step 1 and Step 2 and
failed inside Step 3. Added a matching validation branch that checks, before any real
work:

- `leftcens` is installed **and** exports `impute_censored_conditional()` — the error
  names the `install_github(...@v0.9.0)` command (this duplicates the runtime guard in
  `00_censored_exposure.R`, but fires at Step 1 instead of Step 3);
- `censored_exposure$exposure_vars` is set, and the exposures exist in the raw data;
- the `_lo`/`_hi` bound columns exist (honouring configured suffixes), are numeric, are
  not all-missing, and satisfy `lo <= hi`;
- `margin` is one of `"shash"` / `"gaussian"`;
- `predictors`, if given, exist — and **warns** if the outcome is not among them, since
  that silently reproduces the biased no-`Y` pre-step this whole strategy exists to
  avoid;
- dictionary flags: `use_in_model = TRUE` is required (**stop**); `impute_target = TRUE`
  **warns** (the X block owns those columns, not miceRanger).

It also logs the censored/interval-valued row count per exposure, which is a useful
sanity check on the input convention. Verified: all seven paths fire as intended, and
Step 1 still passes for all five bundled examples.

### Documentation gaps closed

- **`README.md`**: `leftcens` was required by the censored path but **absent from the
  dependency install block** entirely — it is a GitHub-only package, so there was no
  install instruction anywhere. Added a pinned
  `remotes::install_github("fabregithub/leftcens@v0.9.0")` block next to the `survival`
  one, plus a short install pointer in the feature section.
- **`00_config.R`** (and the censored example's config): the commented
  `censored_exposure` template was missing `n_cores`, even though the code reads it
  (`ce$n_cores %||% parallel$impute_workers`) and the README documents it. Added, with
  the unix-only/fork caveat. Also added the install command to the block's comment.

### Test-harness bug: Step 12 was never copied into the isolated run

`run_all.R` has sourced `12_export_draws.R` since 2026-07-16 (`8065903`), but
`test/test_example_common.sh`'s `root_files` list was never updated to copy it. Since
`run_step()` is a bare `source()` with no existence check, **every** `test/` run has
aborted at Step 12 with `cannot open the connection` since that date — including the
full-pipeline run the censored-exposure example still needs. Added
`"12_export_draws.R"` to `root_files`.

### Validation tracks V0 and V1 — the shipped censored-exposure engine is now validated

New plan document [`validation/PLAN_pipeline_validation.md`](validation/PLAN_pipeline_validation.md)
covering what the original integration plan left open. Its framing: every bias number in
`phase1/FINDINGS.md` validated a **prototype** (`cens_mi_y_shash`, ~40 lines written for
the study), not `run_censored_exposure_block_fcs()` — the code users actually run. Five
tracks (V0–V4), each with acceptance criteria written *before* the run.

**V0 — full-pipeline regression, PASSED (~78 s).** First complete 11-step run of the
censored-exposure example (`test/test_censored_exposure_quick.sh`), only possible after
the `12_export_draws.R` harness fix above. All five criteria met; Steps 09/10 correctly
omitted (no `mo()`), MID dropped 30 imputed-Y rows per dataset, 0 divergences.

**V1 — the shipped engine in the MC harness, PASSED (300 reps/cell, ~17 min).** New
procedure 5 (`pipeline_block_fcs`) *sources* `00_censored_exposure.R` and calls the real
`run_censored_exposure_block_fcs()`, then reuses the harness's own `fit_lm_estimand()` /
`rubin_pool()` / metrics — so the comparison against oracle / complete-case / no-Y
pre-step / prototype is like-for-like.

- **Additive DGP:** rel. bias −0.79% (20% ND), +0.48% (40% ND); coverage 0.957 at both;
  rmse 1.06×/1.24× oracle; 0.0485 vs complete-case's 0.0719 at 40% ND. Every
  pre-registered criterion met, 300/300 reps in every cell.
- **Paired agreement with the prototype** ≤0.14% of the estimand in all four cells. One
  cell (additive, 20% ND) is ~2.9 MC SEs from zero but only 0.0005 in magnitude — 24×
  under the acceptance bar; recorded honestly rather than rounded away.
- **Mixture DGP:** reproduces the prototype's known §7.7 bias (+7.54% / +19.55%) rather
  than exceeding it — evidence of faithful implementation, since the X block *is* a
  linear draw.
- **Bonus: independent replication of Phase 1.** Fresh seed (20260825 vs 20260813) puts
  the no-Y pre-step at −5.88% / −14.14%, coverage 0.923 / 0.820 — essentially on top of
  the original −5.6% / −14.4%, 0.92 / 0.82. H1 confirmed twice.
- **Bonus: a V0 observation resolved as noise.** V0's single realization showed the focal
  exposure 12.6% low with its correlated covariate high (the attenuation signature); at
  300 reps there is no systematic attenuation. Flagged as a hypothesis at the time, not
  a finding — correctly.

Findings: [`validation/phase1/FINDINGS_v1.md`](validation/phase1/FINDINGS_v1.md).

**Implementation notes worth keeping.** `00_common_functions.R` is *not*
standalone-sourceable — it reads `paths$objects` at top level for the runtime-override
mechanism, so the harness stubs `paths`/`analysis_spec` and sources into a private env
(not globalenv) to avoid shadowing harness functions. `metrics.R`'s `proc_order` is a
hardcoded factor-level vector: a procedure name missing from it becomes `NA` and silently
drops out of the summary. The harness's replication loop is serial, so the module's
`mclapply` over `m` datasets is safe and is used.

**Still not exercised (deferred to V2):** MID never fired (harness `Y` is complete) and
the miceRanger Z block was a no-op (covariates complete), so the block-FCS *alternation*
itself remains untested under Monte Carlo.

### Validation track V2 — robustness sweep, and an under-coverage finding

New harness pieces: [`validation/phase1/R/robustness.R`](validation/phase1/R/robustness.R)
(MCAR-outcome injector, the 10-scenario grid, per-scenario procedure gating) and
[`validation/phase1/run_v2_robustness.R`](validation/phase1/run_v2_robustness.R). The V1
adapter was generalised to any subset of censored exposures, MCAR covariates, and a
missing outcome.

**Parallelism rewrite.** `run_phase1.R` parallelises *inside* procedures and runs
replications serially, which idles most of a 24-core box when the inner work is small.
V2 inverts it: one fork per (scenario × rep) task, `mc.preschedule = FALSE` for load
balancing (task cost varies ~65×), every procedure told `n_cores = 1`, and the pipeline
sourced **once in the parent** so forks inherit it copy-on-write. Measured **98% worker
efficiency** across 2750 tasks.

**Result: 9/10 scenarios pass, 2.95 h / 63.7 CPU-h, zero failures.** Bias never exceeds
2.71% (bar: 3%) — including every exposure censored, n = 80 000, skew 0.75, rho 0.8.
MID and the miceRanger Z block are now genuinely exercised (both were dead code in V1):
600 multi-exposure tasks and 600 MID tasks, all with correct row counts, zero warnings.

**Finding — interval under-coverage, not bias.** `combined` (all exposures censored +
20% MCAR covariates + 20% missing outcomes) has coverage **0.893**, under the
pre-registered 0.90 investigate line. Diagnosis: intervals are too narrow, not
mis-centred. (CI width / 3.92) / empirical SE correlates **0.837** with coverage and
tracks Z-block involvement — 1.00–1.13 where the Z block is idle (coverage 0.950–0.970),
0.96 with MID only, 0.90 at 40% MCAR, 0.84 with both. Dose-dependent (20% MCAR is fine,
40% is not), compounding, and the **oracle stays calibrated (0.95–1.07) throughout**, so
it is the procedure and not the harness.

Leading hypothesis, **not established**: the X block draws its parameters from a
posterior (proper MI) and is calibrated wherever it works alone; miceRanger's RF/PMM draw
does not, and improper MI understates between-imputation variance — exactly the observed
signature. If so the implication reaches beyond this feature to **any** pipeline analysis
with substantial covariate missingness. Discriminating test running: `combined` at
`m`=100 on identical data — improper MI is a bias in the variance estimator and cannot be
fixed by raising `m`, whereas small-`m` Rubin noise would shrink.

**Cost is dominated by the Z block, not the censored X block.** 90 miceRanger fits per
replication (`outer_sweeps` × `m`). A replication at n = 80 000 with complete covariates
costs 367 s; at n = 800 *with* MAR covariates and a missing outcome it costs 266 s —
covariate missingness, not sample size, drives cost on this path.

**Reproducibility bug fixed.** Per-task seeds keyed on the scenario's index in the
*filtered* list, so `SCENARIOS=` subsets silently generated different data than the same
scenario in a full run — which would have confounded the `m`=30 vs `m`=100 comparison.
Now keyed on the canonical unfiltered position: subsets stay comparable to full runs and
the completed run remains reproducible.

Findings: [`validation/phase1/FINDINGS_v2.md`](validation/phase1/FINDINGS_v2.md).

### Validation track V4 — the under-coverage attributed: improper MI in the Z block

**Result (5 scenarios × 4 arms × 300 reps, 1500 tasks, 5.25 h, zero failures).** V2's
interval under-coverage is caused by **improper multiple imputation in the miceRanger Z
block**, which understates the between-imputation variance `B` by **8–59%**, scaling with
how much the Z block imputes. `Ubar` moves far less — the textbook improper-MI signature.
Swapping the Z block for a proper Bayesian draw restores calibration everywhere it works:
width/SE 0.899 → 1.085 (`mcar_z40`), 0.964 → 1.003 (`missing_y20`), 0.843 → 0.973
(`combined`). In `base`, where nothing is imputed by that block, all four arms return
identical numbers — the negative control that makes the rest interpretable.

**MID is exonerated and should stay.** Disabling it makes calibration *worse*
(width/SE 0.964 → 0.858; 0.843 → 0.754) and introduces **−6 to −7% bias**, against a 3%
bar. Keeping imputed-`Y` rows pulls the estimate toward the imputation model — exactly
what MID exists to prevent.

**The V2 two-mechanism account is WITHDRAWN.** It rested on "`missing_y20` has complete
covariates, so the Z block barely works." That premise was false:
`run_censored_exposure_block_fcs()` sets `impute_y <- TRUE` internally so the X block can
condition on a complete `Y`, so in that scenario the Z block is busy imputing **`Y`
itself** — improperly. One mechanism explains all ten V2 scenarios.

**Scope: this is not a censored-exposure finding.** The Z block is the pipeline's general
row-level imputation path, so any analysis with substantial covariate or outcome
missingness should be expected to produce intervals that are too narrow.

**Does this affect earlier work? No.** The defect predates the left-censored feature — it
lives in the general `run_row_level_imputation()` path — so the question was raised for the
companion study that used an earlier tagged release. Answer: **not affected**. That study
used a version tagged before the left-censored work began, and had very little
missingness. The error is strongly dose-dependent: V4 measured 20% MCAR covariates as
fully calibrated (width/SE 1.00, coverage 0.957), with degradation only appearing at 40%
and above, or when the outcome is also imputed. Recorded here so the question does not
have to be re-opened when the fix ships.

**Back-compatibility position.** Old configs must keep *running* (config compatibility);
identical *numbers* are not promised across versions — analyses cite the pipeline version
they used instead. Two Step-1 checks added during this session were over-strict and were
downgraded to warnings after review, because they would have rejected configs that
previously ran: an unknown name in `censored_exposure$predictors` (the X block has always
silently dropped these) and `use_in_model = FALSE` on a censored exposure (harmless when
`predictors` is set explicitly, since the flag only feeds the auto predictor set).

**The fix is NOT to adopt the V4 instrument.** `.v4_proper_row_imputation()` is a
deliberately simple linear/logistic draw built to isolate the mechanism. It overshoots
(width/SE 1.085 in `mcar_z40`, i.e. ~9% too wide) and its bias is marginally worse
(−3.18% vs −2.11%). The real fix is to make the Z block propagate parameter uncertainty
properly. The V4 harness measures any candidate directly.

Findings: [`validation/phase1/FINDINGS_v4.md`](validation/phase1/FINDINGS_v4.md).

**New harness pieces:** `R/proper_impute.R` (the proper-draw instrument),
`run_v4_variance.R` / `.sh`, three isolating variants in `procedures_pipeline.R`, and
`rubin_pool()` now returns the variance decomposition (`ubar`, `b`, `fmi`) — the piece V2
lacked, without which the shortfall could not have been localised.

**Operational note — `resume` is not trusted.** Resuming from an accumulated checkpoint
reproducibly killed the parent process right after the first chunk, in both detached and
foreground runs, while an identical configuration with a fresh checkpoint completed all 35
chunks. The checkpoint's structure looks correct (right columns, classes, procedures) and
the cause was not found; `resume` is therefore off by default with a warning. This costs
nothing operationally: the checkpoint still flushes after every chunk and partial results
are recoverable via `summarise_v4(readRDS("results/v4_checkpoint.rds"))`. Several hours
were lost to misdiagnosing this as an environment problem (nohup, screen, caffeinate, OOM,
sleep) before a user-launched run after a reboot ruled all of that out.

### v1.4.0 — proper multiple imputation is now the default (behaviour change)

**`analysis_spec$imputation$proper_draw` defaults to `TRUE`.** The Z block now bootstraps
its training data once per imputation (`run_row_level_imputation_proper()`), so the forest
varies across imputations and that uncertainty reaches the intervals. Set it to `FALSE` to
reproduce a pre-v1.4.0 analysis exactly.

**This changes results.** Analyses with missing covariates or outcomes get wider intervals,
and point estimates shift slightly. Config compatibility is preserved — old `00_config.R`
files run unchanged — but numbers are not promised across versions, so analyses must cite
the pipeline version used. The companion study is unaffected (pre-left-censored tag, very
little missingness).

**Evidence (V5, 5 scenarios × 5 arms × 300 reps, 8.05 h, zero failures):** the fix closes
**66%** of the calibration shortfall at 40% covariate missingness (coverage 0.933 → 0.947,
`B` +32.6%, t = 11.1), **42%** with a missing outcome, and **23%** under heavy combined
missingness (coverage 0.893 → **0.910**). It never over-corrects (max width/SE 1.014
against the linear instrument's 1.085), leaves already-calibrated cells alone, and has the
**best bias of any variant tested**. Zero bootstrap fallbacks across 1500 tasks.

**It is not a complete fix.** Under heavy combined missingness coverage reaches ~0.91, not
0.95. The reason is worth remembering: **a random forest is already bagged**, so an outer
bootstrap shifts it only slightly — bagging damps the very uncertainty the bootstrap is
meant to inject. Measured on the same data, the bootstrap-on-forest inflated `B` by
+31–33% where a parametric posterior draw inflated it by +41–59%. The standard
"bootstrap for properness" recipe is weaker for bagged learners than for parametric models.

**Why we did NOT switch to `mice`.** The linear instrument calibrates better (0.973 vs
0.878 in the worst cell) but the DGP has **linear-Gaussian covariates**, so a linear
imputation model is correctly specified there and is never penalised for misspecification —
and it already carries the worst bias of the three arms. `mice`'s `rf` method is
essentially what we implemented. The simulation cannot currently distinguish "properly
calibrated" from "correctly specified", so that choice is deferred until the DGP has
non-linear covariates. Findings:
[`validation/phase1/FINDINGS_v5.md`](validation/phase1/FINDINGS_v5.md).

**Cost.** The proper path adds one prediction pass *per imputation* (`m` passes per
completed set, not one). The V5 grid took 8.05 h against 5.25 h for the four-arm V4 grid.
Production runs with large `m` and `n` should expect the Z block to cost meaningfully more.

### v1.5.0 — BART is the Z-block imputer; R8 and R12 close (behaviour change)

> **No `v1.5.0` tag exists, and none can be added.** This release and v1.5.1 landed in the
> same commit (`5a66b7d`), so the tag on that commit is `v1.5.1`. To check out the tree in
> which BART became the default, use `v1.5.1`; `v1.4.0` is the last tag that predates it.
> Do not go looking for a `v1.5.0` tag.


**`analysis_spec$imputation$z_imputer` replaces the `proper_draw` flag**, defaulting to
`"bart"`. Values: `"bart"` (v1.5.0 default, needs `dbarts`), `"forest_boot"` (the v1.4.0
default), `"forest"` (plain miceRanger, improper, pre-v1.4.0). `proper_draw` is still
honoured — `TRUE` → `forest_boot`, `FALSE` → `forest` — so a config written against v1.4.0
keeps the imputer it was validated with; only configs specifying *neither* move. If `bart`
is selected and `dbarts` is missing, the pipeline warns loudly and falls back to
`forest_boot` rather than failing an upgraded install. The imputer used is logged each run.

**Why BART.** The Z block must be proper MI or `B` is understated and intervals are too
narrow. Two earlier attempts each got half of it: a bootstrapped forest is flexible but
under-disperses (a forest is *already* bagged, so an outer bootstrap barely moves it), and
`mice pmm` is properly dispersed but misspecified under non-linear covariates. BART is
both — a sum-of-trees model with priors on tree structure *and* leaf parameters, so a
posterior sample IS a proper draw.

**Evidence (V7, 5 chained phases, ~15 h, zero failures):**
- **P1** (7 scenarios × 4 arms × 300 reps): BART has the **best bias in all six diagnostic
  cells** (max 1.29% against `forest_boot`'s 2.67% and `mice pmm`'s 4.70%), coverage
  0.937–0.963, `base` control bit-identical across all four arms, `n_ok` 300 in all 35
  cells. Cost ~**4% of runtime** — cheaper than either forest arm.
- **P2** (`ntree` 50 → 200): bias stays ≤0.8%, width/SE moves ≤0.016, and `forest_boot` is
  identical at both settings (the control). **The verdict is not tuning-fragile**, and the
  overshoot persists at 200 trees so it is intrinsic, not a tuning artefact.
- **P3** (`m` ∈ {10,20,30,50}): **width stabilises by `m` ≈ 30** — 10 → 30 narrows ~2.7%,
  30 → 50 changes <0.2%. FMI flat at ≈0.30, and `m ≈ 100 × FMI` independently gives 30.
  **R12 / design-plan Phase 5 closes**; the shipped default was already right.

**The width/SE criterion was mis-designed — my error, not BART's.** BART failed the
pre-registered band (0.97–1.03) in three cells. But `width/SE` divides by the empirical SD
of the estimates, whose standard error is `sd/sqrt(2(n-1))` ≈ **4.1%** at 300 reps: a ±3%
band on a ±4%-noise quantity is unattainable by construction. Tested against that noise,
**none of BART's deviations is significant** (max z = +1.87); the only real deviation in the
whole table is `forest_boot` being *too narrow* in `combined` (z = −2.10) — the very
anti-conservatism R8 exists to remove. Two alternative explanations were tested and refuted
first: a `t`-vs-`z` bias in the 3.92 divisor (`t` = 1.965–1.974, moves the ratio ≤0.008) and
light-tailed estimates (`emp95/SD` 0.95–1.02, excess kurtosis −0.43 to +0.51).

**Lesson worth keeping:** `width/SE` is an excellent *mechanism* diagnostic — it is how V4
localised the shortfall to `B` rather than `Ubar` — but a poor *acceptance criterion*,
because its noise floor exceeds the effect size of interest. Coverage is the acceptance
criterion (MC SE 0.013). The claim is therefore that BART's width is consistent with
calibration to within ~±8%, not that it is calibrated.

**Follow-up opened:** a redundancy — block-FCS calls the Z block with `m`=1 per outer sweep
while the BART imputer runs its own 3 inner FCS iterations, so the censored path performs
~3× more BART fits than it needs (9 alternations where 3 were designed). Tracked in
`validation/INTEGRATION_SUMMARY.md` §2. **Mostly resolved in v1.5.1 below — the redundancy
was real but removing it was not free.** The interval-overshoot item raised alongside it was
**withdrawn** once the noise floor was computed — see above.

### v1.5.1 — the Z-block inner-iteration count is chosen per sweep (V8)

The v1.5.0 follow-up above assumed the outer block-FCS loop's alternation made the Z
block's own inner FCS iterations redundant, and defaulted them to 1. **V8 tested that
assumption at 1000 replications per cell and it is only half true.** The outer loop
alternates **Z ↔ X**; it does nothing to make the Z block's own targets condition on each
other. With three targets (two covariates plus an imputed `Y`), `inner = 1` cost **−0.72
percentage points** of bias against a pre-registered 0.5 pp bar — the entire confidence
interval outside it. With two non-outcome targets it was genuinely free (+0.04 pp, CI
±0.13). Evidence: [`validation/phase1/FINDINGS_v8.md`](validation/phase1/FINDINGS_v8.md).

**Coverage and width/SE moved by less than 0.02 in every cell while bias moved 0.72 pp.**
This was a bias-only defect, invisible to the calibration diagnostics — worth remembering
before trusting a coverage check to validate a change to the imputation chain.

`.ce_default_inner_iter(targets, y_var)` in `00_censored_exposure.R` now picks the count
per sweep from the Z block's actual target set: **1** when the block has ≤2 targets and
none is the outcome, **3** otherwise. `.ce_one_imputation()` calls it inside the sweep
loop, after `make_row_level_imputation_spec()` — the target set is not known before that.
An explicit `analysis_spec$imputation$bart_inner_iter` still wins.

`Y` is excluded from the cheap path beyond what the target count alone would require,
because the X block conditions on `Y` — an under-conditioned `Y` feeds straight into the
censored exposure draw. V8's cells confound "3 targets" with "`Y` is a target"; requiring
both conditions is the reading that stays safe under either explanation.

**`bart_inner_iter = 3` was removed from the shipped `00_config.R`** (left commented, with
the reason). An explicit value overrides the per-sweep choice, so leaving it set would have
made the new default dead code for anyone using the stock config. Off the censored path the
fallback in `run_row_level_imputation_bart()` is unchanged at 3.

**This is not a behaviour change for a stock censored-exposure analysis** — that path
imputes `Y` (MID), so it takes the `inner = 3` branch, which is the behaviour V7 validated.
Verified bit-for-bit against the stored V8 replications in both routing directions.

**Also settled in the same run:** V7's BART arm overrode `run_row_level_imputation` with the
validation harness's own instrument, so V7 had formally measured the instrument rather than
the shipped code. The two are **bit-identical** across all 4000 replications (max difference
exactly 0 on `estimate`, `se`, `ci_lo`, `ci_hi`, `ubar`, `b`). V7's conclusions apply to the
shipped pipeline without qualification.

### V3 / R11 — the mixture path is confirmed unusable (2026-08-29)

Not a code change: a **scope finding** that changes what the pipeline may be used for.

The censored-exposure strategy was always documented as additive-only, with the mixture
caveat resting on a scaffold estimand. V3 measured it on the estimands a BKMR analysis
actually reports, each with truth derived analytically from the generator, at 200
replications per cell. The result is worse than the caveat implied: the shipped engine is
the **worst of four arms** (mean 14.6% absolute paired excess over the oracle, against 6.4%
for LOD/√2 substitution and 10.5% for complete case), and it **destroys 57% of the
curvature** at 40% non-detects.

**Why**: the X block draws the censored exposure from a conditional that is *linear in the
predictors*, so it imposes linearity on precisely the rows where curvature would appear.
Imputing with the wrong functional form is worse than not imputing — on curvature, LOD/√2
substitution is statistically indistinguishable from the oracle (p = 0.42) while the
congenial engine is not.

**Interaction is more robust** than curvature: it survives censoring of a single exposure
(+2.6% to +5.7%) and breaks only when a second exposure is censored too (−12.5%). The
failure is *curvature*, not mixtures in general.

**Coverage was blind to it** — 0.945–0.960 alongside a −42% point bias, because the pooled
intervals ran 1.5–1.6× wider than the estimates' true sampling spread. This is the second
consecutive track where calibration diagnostics missed a real bias (see v1.5.1 above).
Coverage is not a substitute for a bias check against known truth.

The root README's mixture scope box has been rewritten from "remains biased" to an explicit
**do not use**. Evidence:
[`validation/phase1/FINDINGS_v3.md`](validation/phase1/FINDINGS_v3.md).

The fix — substantive-model-compatible imputation — is unbuilt and is now the only
substantive track left (`validation/ROADMAP.md` item 07).

### V9 — the mixture failure's cause is confirmed, and the fix is demonstrated (2026-08-29)

Again not a code change: a **mechanism result** that de-risks the remaining build.

V3 concluded that the censored-exposure X block destroys curvature *because* its
conditional is linear in the predictors — but that was an inference from a pattern, not a
demonstration. V9 replaced only that component, drawing the censored exposure from the
correct conditional under the same surface, and changed nothing else.

| `curv_X1`, paired excess over the oracle | `nd40` | `nd40_all` |
|---|---|---|
| shipped pipeline (linear conditional) | −55.6% | −59.2% |
| SMC, true surface | +5.9% | −0.3% |
| SMC, surface estimated each sweep | +6.2% | −0.1% |

**89–100% of the gap closes.** Across all fourteen cells the mean absolute excess over the
oracle falls from **17.2% to 0.9%**. The control arms were bit-identical to V3 on all 2,800
shared rows, so this is measured against exactly the behaviour V3 recorded.

**Estimating the surface rather than knowing it costs 0.2 pp.** The plug-in arm fits the
outcome model on the current completed data each sweep and draws the coefficients from
their posterior. The fix needs the right functional *form*, not the true parameters — which
is the useful news for anyone building it.

**What is still missing**: both arms were *handed* the generator's formula. A real BKMR
analysis does not know the surface. That gap is the whole content of
`validation/ROADMAP.md` item 07, now much better defined.

**Nothing shipped changes.** The SMC code is a validation harness instrument
(`validation/phase1/R/smc_impute.R`), deliberately not wired into
`00_censored_exposure.R`, and the root README's mixture restriction stands. Evidence:
[`validation/phase1/FINDINGS_v9.md`](validation/phase1/FINDINGS_v9.md).

### V10 — MAR covariate missingness tested; the additive claims hold (2026-08-29)

The pipeline's additive-path claims all rested on **MCAR** covariates — the easy,
unrealistic mechanism. V10 is the first test of MAR, the mechanism MI is actually built for.

| cell | rel. bias | coverage |
|---|---|---|
| `mcar_z40` | −1.668% | 0.957 |
| `mar_z40` | **−0.167%** | 0.967 |
| `mar_z40_strong` (doubled coefficients) | +0.252% | 0.966 |
| `nl_mar_combined` (MAR + non-linear covariates + 20% missing Y) | −0.945% | 0.958 |

**MAR does not degrade the engine — it is slightly easier than MCAR.** The strength=0
control reproduced the MCAR cell (−0.11 pp), so the injector changed only the mechanism,
and realised missing fractions matched to 0.1 pp, ruling out the obvious confound.

**The registered ±1 pp criterion is breached at 40% (+1.50 pp)** — in the favourable
direction. Recorded as a breach rather than rewritten after the fact; the conclusion does
not depend on the bar.

Consequence for the docs: the root README's scope note widens from MCAR to MAR, and the
MCAR-based figures are now known to be *conservative*. **MNAR remains untested** and is the
harder open question. Evidence:
[`validation/phase1/FINDINGS_v10.md`](validation/phase1/FINDINGS_v10.md).

### v1.5.2 — documentation restructure and a published site (no behaviour change)

**No pipeline code changed in this release.** It is versioned because for a template that
people clone, the documentation is much of the product.

> **The `v1.5.2` tag (`1c01b2e`) predates three commits that belong to this work**, because
> it was cut before they were written: the GitHub Actions bump to Node 24 majors
> (`c7796da`), the site URL added to the README (`f037541`), and this changelog entry. A
> checkout of `v1.5.2` therefore has a workflow that emits deprecation warnings and a
> README that does not link to the site the release introduces. The tag was left in place
> rather than force-moved, since it had already been pushed — moving a published tag
> silently leaves anyone who fetched it on different code. Use `main`, or the next tag, for
> the complete state.

**The README went from 3,078 lines to 615** (~66 minutes of reading to ~16). Ten focused
guides now live in [`docs/`](docs/), indexed from the README and from
[`docs/index.md`](docs/index.md):

| guide | covers |
|---|---|
| `setup.md` | R, CmdStan, Quarto and packages on macOS, Windows, Linux |
| `config-reference.md` | option-by-option `00_config.R` reference |
| `variable-dictionary.md` | dictionary fields, and when config may override them |
| `choosing-m.md` | how many imputations, the auto m-increment loop, the three run modes |
| `censored-exposures.md` | below-LOD exposures, their validation, and the mixture restriction |
| `repeated-measures.md` | subject-wide versus row-wise imputation |
| `performance.md` | sizing workers and chains; memory-versus-cores |
| `operations.md` | restarting, monitoring, CmdStan cache, debugging a fit |
| `reporting.md` | publication outputs and Methods/Results templates |
| `examples.md` | the four bundled examples and the automated tests |

**The mixture restriction is deliberately duplicated**, in full, in both the README and
`docs/censored-exposures.md` — the only content copied rather than linked. Acting on that
warning matters more than finding it, and a reader who never opens `docs/` must still not
miss it.

**Heading levels were repaired.** Sections 5, 8, 9 and 10 had been `#` rather than `##`,
giving the document five H1s and a broken outline (the links worked; the hierarchy did
not). Section numbers were then dropped entirely — the documentation index is the map now,
and numbering had become a holdover from when the README *was* the documentation.

**The four translated READMEs were deleted** (`docs/README.{de,es,fr,ja}.md`). They
predated the censored-exposure feature, never mentioned `leftcens` or the block-FCS
strategy, and — most seriously — omitted the restriction against mixture/BKMR analyses, so
non-English readers were missing a safety-critical caveat. Keeping four translations in
step with the pipeline was not sustainable; a translation that has silently fallen behind
is worse than none.

**A documentation site is published** at
<https://fabregithub.github.io/generic-mi-brms-pipeline/>, built by
[`_quarto.yml`](_quarto.yml) from the same Markdown the repository serves, so the two
cannot drift. Quarto rather than pkgdown: pkgdown builds from package metadata
(`DESCRIPTION`, `man/`, `NAMESPACE`), none of which this template has, and its central
feature — the function reference — would be empty. Quarto was already a pipeline
dependency. Deployment is via a GitHub Actions artifact, so there is **no `gh-pages`
branch** to maintain.

### V11 — the imputer default survives a non-linear outcome (2026-08-31)

A scope finding, not a code change.

`z_imputer = "bart"` was chosen in V6/V7 on a design where the outcome is **linear in
`(logX, Z)` by construction**, so a parametric imputer for `Z1` was correctly specified
against the outcome and could never be penalised for misspecification. `FINDINGS_v7.md`
named that confound and left it untested. V11 removed it by adding `b_zq (Z1^2 - 1)` to the
outcome, with the analysis model gaining the matching `I(Z1^2)` term so the focal estimand
stays exactly recoverable.

**The default survives.** `bartMI` is never significantly beaten in any non-linear cell —
where the summary table shows another arm ahead, the *paired* difference is
indistinguishable (*p* = 0.11–0.34) and favours `bartMI` in two of three. Those arms only
look better because they start at +1.6% to +1.8% bias and a uniform −3 pp shift carries them
through zero: cancellation, not robustness, and it would reverse under opposite curvature.

**The V7 hypothesis was wrong.** All four imputers degrade by 2.86–3.30 pp — a spread of
0.44 pp. The confound was real but inert, and the V6/V7 verdict stands on evidence that
could have overturned it.

**The finding with consequences is pipeline-wide**: a non-linear outcome costs ≈3 pp of bias
*whichever* imputer is used, taking the shipped default from 0.53% to 3.10% mean absolute
bias — outside the figures quoted for the additive path. The scope box in
[`docs/censored-exposures.md`](docs/censored-exposures.md) now says so, and notes the
mechanism is general: the Z block conditions on the outcome, so any analysis with missing
covariates is affected. **Coverage was unchanged** (0.953 → 0.953), since a 3% relative bias
is ≈0.28 SE at `n` = 800 — the third consecutive track where calibration diagnostics missed
a bias finding (V3, V8, V11).

Evidence: [`validation/phase1/FINDINGS_v11.md`](validation/phase1/FINDINGS_v11.md).

### Known gap (not addressed)

The four translated READMEs (`docs/README.{de,es,fr,ja}.md`) contain **zero** mentions
of `leftcens` or the censored-exposure strategy — they predate the whole feature. This
needs a proper translation pass, not a lone install line.

---
## Ideas for future development

### Pipeline
- [ ] **Sensitivity analysis runner**: re-run Step 4 with alternative priors/formula and compare
      posteriors against the primary run
- [ ] **Multi-outcome support**: structured config for fitting the same exposure to multiple outcomes
      in one `run_all.R` invocation
- [ ] **DAG-based confounder selection**: integrate `dagitty`/`ggdag` to auto-generate the
      adjustment set from a user-supplied DAG and populate `00_variable_dictionary.csv`
- [ ] **Windows clean helper**: `launch.R` handles cleaning from RStudio, but a pure R
      `clean_all()` in `00_common_functions.R` would make it callable from any context
- [ ] **Progress notifications**: email or Slack webhook on pipeline success/failure
      (hook into existing `pipeline_success.flag` / `pipeline_error.flag`)
- [ ] **R package** (low priority): thin package wrapping the scripts; would need refactoring
      to pass `paths` and `analysis_spec` as arguments rather than sourcing them

### Federated export (Step 12)
- [ ] **Draws subsampling**: option to export a random subsample of draws per imputation to
      reduce file size when `m` and `iter` are large
- [ ] **Differential privacy**: optional noise injection on exported draws before transfer
- [ ] **Encrypted export**: wrap `cohort_draws.rds` in an encrypted container for secure transfer


