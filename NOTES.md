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
