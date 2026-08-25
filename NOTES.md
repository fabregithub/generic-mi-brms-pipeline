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


