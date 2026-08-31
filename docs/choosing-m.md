# Choosing `m`, and running at scale

How many imputations you need, how the automatic m-increment loop decides for you, and the three run modes from a quick smoke test to a full production run.

> `m` is the setting people most often guess at. The pipeline can determine it for you by running batches until the pooled summaries stop moving.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

## Choosing the number of imputations adaptively

There is no universal value of `m` that is automatically sufficient for every analysis. The required number of imputations depends on the fraction of missing information, the amount and pattern of missingness, the target estimands, and how reproducible the posterior summaries need to be. A large dataset does not automatically require `m = 100`, and a small dataset can sometimes require more than `m = 100` if the fraction of missing information is high or the target estimates are unstable. With large datasets in particular, models tend to stabilise at a smaller `m` than with small datasets, while each fit also takes longer -- so paying for a fixed `m = 100` up front can waste a large amount of compute time.

There are two ways to use this adaptive idea: an **automatic loop** that `run_all.R` drives for you, and a **manual staged workflow** that you drive yourself, one batch at a time. The automatic loop is recommended for most analyses; the manual workflow remains useful for closer manual inspection between batches, or for resuming an already-started large run.

### Automatic m-increment loop (recommended)

`analysis_spec$mi_stability$auto_increment` is `TRUE` by default in `00_config.R`, so a typical setup looks like:

```r
analysis_spec$imputation$m <- 100  # acts as the ceiling, not the starting point

analysis_spec$mi_stability <- list(
  auto_increment = TRUE,

  # Defaults to analysis_spec$parallel$fit_workers, rounded up to a multiple
  # of fit_workers if you override it, so every batch fully occupies the
  # parallel workers with none left idle. For example, with fit_workers = 4,
  # the default batches are m = 4, 8, 12, ...
  increment_size = NULL,

  # Same tolerances used by the manual stability check below.
  parameter_regex = "^b_",
  exclude_intercept = TRUE,
  estimate_tolerance = 0.05,
  ci_endpoint_tolerance = 0.05,
  pd_tolerance = 0.02
)
```

Then run the pipeline as usual:

```bash
# Bash command block
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

`run_all.R` then:

1. Fits the first batch of imputations (size = `increment_size`).
2. Runs Steps 3, 4 and 6 for that batch.
3. Fits the next batch, extending (not replacing) the existing imputations and fits.
4. Compares the new cumulative results against the previous batch using the configured tolerances.
5. Stops increasing `m` as soon as every selected parameter is stable, or once the configured `analysis_spec$imputation$m` is reached, whichever comes first.
6. Runs Steps 5, 7, 8, (9-10 if `mo()` is present) and 11 once, on whatever `m` the loop settled on.

Each batch's imputations are seeded deterministically (`base_seed + n_existing`, where `base_seed` is `analysis_spec$imputation$seed`, defaulting to `analysis_spec$model$seed`), so repeated incremental runs reproduce the same imputed datasets. This is not bit-identical to a hypothetical single-shot `miceRanger(m = target_m)` call, because `miceRanger` generates all `m` chains of one call together with no native per-chain seed/resume support, but it is fully deterministic and reproducible across repeated incremental runs of this pipeline.

The loop only ever adds new imputations and fits; it never reduces or overwrites a finished batch. If you want to extend an `auto_increment = FALSE` run that already has some imputations on disk, you can also use the same extension mechanism manually -- see [If a larger run has already been started](#if-a-larger-run-has-already-been-started) below.

`analysis_spec$model$run_smoke_fit` does not need any special handling for the loop. The smoke fit's job is to catch formula/prior/data/CmdStan problems early, once -- so the loop runs it only for the first batch, then automatically disables it for every subsequent batch, regardless of how `run_smoke_fit` is set in `00_config.R`. You can safely leave `run_smoke_fit = TRUE` for the whole run; it will not be repeated before every batch's parallel fitting.

### Manual staged workflow

If you prefer to inspect the stability check yourself between batches rather than letting `run_all.R` decide automatically, set `analysis_spec$mi_stability$auto_increment <- FALSE` in `00_config.R` and use the manual workflow instead. The safest rule is:

> Start with a modest `m` and increase `m` only if needed. Do not reduce `m` within the same saved run.

Recommended staged workflow:

1. Start with a modest number of imputations, for example `m = 20` or `m = 24`, and set `analysis_spec$imputation$allow_extend <- TRUE` in `00_config.R`.
2. Run the pipeline through at least Step 6 so per-imputation posterior draws are available.
3. Run `11_check_imputation_stability.R`.
4. If the primary posterior summaries are stable, stop and use that `m`.
5. If the primary posterior summaries are not stable, increase `analysis_spec$imputation$m`, for example to 40.
6. Rerun the pipeline. With `allow_extend = TRUE`, existing valid imputations and fits are reused and only the new ones are added.
7. Repeat with larger values, for example 60, 80, 100 or beyond, only if needed.

Convenient sequences are `20 -> 40 -> 60 -> 80 -> 100`, or, if four models are fitted in parallel, `24 -> 40 -> 60 -> 80 -> 100`.

The batching values can reflect the available computing layout. For example, if four models are fitted at a time, values such as `m = 24`, `40`, `60`, `80` and `100` are convenient. The batching convenience is not itself the scientific justification; the scientific justification is the stability of the prespecified primary summaries.

For small samples, highly incomplete variables, weakly identified models or estimates near a decision boundary, more than `m = 100` may be needed. In those cases, either increase `m` further or report that the imputation-count stability check did not support a smaller value.

### Running the imputation-count stability script

After Step 6 has created per-imputation posterior draw files, run:

```bash
# Bash command block
Rscript 11_check_imputation_stability.R
```

This creates publication-ready stability outputs in `results/publication/mi_stability/`.

Typical outputs include:

- `tables/imputation_stability_all_batches.csv`
- `tables/imputation_stability_final_comparison_full.csv`
- `tables/imputation_stability_final_comparison_display.csv`
- `tables/imputation_stability_settings.csv`
- `tables/imputation_stability_stepwise_summary.csv`
- `tables/imputation_stability_stepwise_comparison_full.csv`
- `figures/imputation_stability_trajectories.png`
- `figures/imputation_stability_stepwise_change.png`
- `report/imputation_stability_report.qmd`

If Quarto is installed, the report can be rendered to HTML or DOCX.

The stepwise stability summary is especially useful for deciding whether results have already flattened at an early value of `m`. For example, it compares transitions such as `m = 8 to m = 12`, `m = 12 to m = 16`, `m = 16 to m = 20`, and `m = 20 to m = 24`.

Key quantities to inspect include:

- Stable, %
- Maximum absolute estimate change
- Median absolute estimate change
- Maximum absolute CrI-endpoint change
- Maximum relative odds-ratio change, %
- CrI exclusion changed, n

If these values are already negligible at an early transition and remain negligible afterwards, it is defensible to stop increasing `m`, provided the checked estimands are the prespecified primary estimands.

You can customise the stability check by adding an optional `mi_stability = list(...)` block in `00_config.R`:

```r
mi_stability = list(
  # If NULL, the script chooses batch sizes from available imputations.
  batch_sizes = c(12, 16, 20, 24),

  # Use exact posterior draw column names for primary estimands,
  # or leave NULL and select by parameter_regex.
  primary_parameters = NULL,

  # Default checks ordinary fixed-effect coefficients.
  # For all extracted parameters, use ".*".
  parameter_regex = "^b_",

  # Usually TRUE for publication summaries.
  exclude_intercept = TRUE,

  # Practical stability thresholds.
  # Replace these with thresholds appropriate for the scientific question.
  estimate_tolerance = 0.05,
  ci_endpoint_tolerance = 0.05,
  relative_transformed_tolerance_pct = 5,
  pd_tolerance = 0.02,

  max_plot_parameters = 12,

  # FALSE by default: this script's own tables/figures are also embedded in
  # the main report's "Imputation-count stability" chapter (see the reporting guide),
  # which run_all.R renders after this script runs. Set this TRUE only if
  # you also want this script's standalone, more detailed report rendered.
  render_quarto = FALSE
)
```

### If a larger run has already been started

If you already started a larger run, for example with `m = 100`, and then decide after 45 completed fits that the results are stable, do not simply reduce `m` in `00_config.R` and rerun `run_all.R` from the beginning. That may cause the saved imputation specification to be treated as incompatible with the current config.

Safer options are:

- **Option 1:** continue the larger run if the extra computation is acceptable.
- **Option 2:** run only downstream scripts on the completed fits, avoiding Step 3 imputation.
- **Option 3:** start a clean new run in a new output folder with the smaller final `m`.
- **Option 4:** extend the existing imputations to a larger `m` using `allow_extend`.

Option 4 only ever moves `m` upwards. If 20 imputations already exist on disk and you raise `analysis_spec$imputation$m` to 30, setting:

```r
analysis_spec$imputation$allow_extend <- TRUE
```

makes Step 3 generate only the 10 new imputations and append them to the existing manifest; the original 20 imputed-data files are never touched or regenerated. This is exactly the mechanism `run_all.R`'s automatic increment loop uses internally, and is also the right tool for manually extending a fixed-`m` run after inspecting Step 11's output yourself.

For Option 2, the downstream sequence is usually:

```bash
# Bash command block
Rscript 05_diagnostics.R
Rscript 06_posterior_summary.R
Rscript 07_posterior_prediction.R
Rscript 11_check_imputation_stability.R
Rscript 08_publication_results.R
```

Step 11 runs before Step 8 here for the same reason `run_all.R` orders them this way: Step 8's report embeds Step 11's tables/figures, so they need to exist on disk first.

This avoids re-running imputation. Use this only after confirming that the completed fit files and model-data files correspond to the same imputed datasets.

This stability check assesses numerical Monte Carlo stability as `m` increases. It does not validate the missing-data mechanism and does not make MNAR or censored data problems ignorable. For ready-to-adapt Methods text describing this adaptive-`m` justification, see [the manuscript writing guide](../docs/reporting.md#manuscript-writing-guide).

Recommended approach for a new analysis:

## 1. Quick test run

For the first test of a new dataset or model, use small settings. Edit these in `00_config.R` using RStudio or another text editor:

```r
imputation = list(
  ...
  m = 5,
  ...
)

model = list(
  ...
  chains = 1,
  iter = 500,
  warmup = 250,
  run_smoke_fit = TRUE,
  ...
)

parallel = list(
  ...
  fit_workers = 1,
  cores_per_fit = 1,
  summary_workers = 1,
  prediction_workers = 1,
  future_globals_maxsize_gb = 8
)

posterior_prediction = list(
  ...
  ndraws = 200,
  ...
)
```

This quick test is intended to check that:

- the data are read correctly
- the variable dictionary is valid
- imputation runs
- the `brms` formula is correct
- priors are compatible with the model
- one model can be fitted successfully
- posterior summaries and publication outputs are created

## 2. Parallel test run

After the quick test succeeds, test parallel fitting with modest settings:

```r
imputation = list(
  ...
  m = 10,
  ...
)

model = list(
  ...
  chains = 4,
  iter = 500,
  warmup = 250,
  run_smoke_fit = TRUE,
  ...
)

parallel = list(
  ...
  fit_workers = 2,
  cores_per_fit = 4,
  summary_workers = 2,
  prediction_workers = 2,
  future_globals_maxsize_gb = 20
)
```

This checks that multiple chains and multiple imputed datasets can run safely on the machine.

## 3. Full production run

For the final analysis, increase the imputation and MCMC settings. For a high-memory machine, a typical setting is:

```r
imputation = list(
  ...
  m = 100,
  ...
)

model = list(
  ...
  chains = 4,
  iter = 2000,
  warmup = 1000,
  run_smoke_fit = TRUE,
  ...
)

parallel = list(
  ...
  fit_workers = 4,
  cores_per_fit = 4,
  summary_workers = 4,
  prediction_workers = 4,
  future_globals_maxsize_gb = 80
)

posterior_prediction = list(
  ...
  ndraws = 1000,
  ...
)
```

For repeated runs of the same already-tested analysis, you may set:

```r
model = list(
  ...
  run_smoke_fit = FALSE,
  ...
)
```

but keep it as `TRUE` when changing the dataset, formula, priors, outcome family or imputation strategy.

For large analyses, avoid `brm_multiple()` unless you have a specific reason to use it. The one-fit-per-imputation design is usually safer and easier to restart.

---

---

*[← Back to the main README](../README.md)*
