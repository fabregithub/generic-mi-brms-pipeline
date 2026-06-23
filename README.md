# Generic MICE + brms Pipeline Template

This is a reusable R pipeline template for Bayesian regression analyses with optional multiple imputation.

It supports:

- data validation
- optional multiple imputation with `miceRanger`
- Bayesian regression with `brms` + `cmdstanr`
- one-fit-per-imputation checkpointing
- parallel model fitting across imputations
- diagnostics
- posterior summaries
- posterior prediction for rows with missing outcomes
- publication-ready tables, figures, methods/settings metadata and report templates

The default example uses the built-in public dataset `datasets::airquality`, so the template can be tested and demonstrated without private data.

---

## Contents

1. [Background and purpose](#1-background-and-purpose)
2. [Structure of the pipeline](#2-structure-of-the-pipeline)
3. [Quick start](#3-quick-start)
4. [Adapting the pipeline to private study data](#4-adapting-the-pipeline-to-private-study-data)
5. [Variable dictionary](#5-variable-dictionary)
6. [Parallelisation and performance tuning](#6-parallelisation-and-performance-tuning)
7. [Logging, monitoring, restarting, troubleshooting and debugging](#7-logging-monitoring-restarting-troubleshooting-and-debugging)
8. [Publication outputs and inference guidance](#8-publication-outputs-and-inference-guidance)
9. [Manuscript writing guide](#9-manuscript-writing-guide)
10. [Examples and tests](#10-examples-and-tests)
11. [Computing environment setup](#11-computing-environment-setup)

---

## 1. Background and purpose

This pipeline is intended for applied Bayesian regression analyses where the data may contain missing covariates, repeated outcomes, large models or models that need careful checkpointing. It combines multiple imputation, one-fit-per-imputation Bayesian modelling, diagnostics, posterior summaries, posterior prediction and publication-oriented outputs.

The design prioritises reproducibility and restartability over keeping all fitted models in memory. This is why the pipeline fits one `brms` model per imputed dataset, saves each fit immediately and reuses valid checkpoint files on rerun.

### Brief cautions and limitations

This repository is a workflow scaffold, not a substitute for statistical judgement. Before using it for a scientific analysis, check that the imputation strategy, model formula, priors, diagnostics and posterior summaries are appropriate for your study question.


The multiple-imputation steps are intended for variables where a standard MICE-style assumption is scientifically defensible, typically missing completely at random (MCAR) or missing at random (MAR) after conditioning on observed variables included in the imputation model. The pipeline does not automatically handle non-ignorable missingness, MNAR mechanisms, censoring, truncation, limit-of-detection problems or structural missingness.

If a variable has non-standard missingness, such as left-censored values below a detection limit, skip-pattern missingness or values missing for design reasons, process or model that missingness appropriately before using this pipeline. Do not simply code such values as ordinary `NA` and rely on the default MICE workflow unless that is justified for the study.

Some models can be computationally expensive. In particular, large mixed logistic models, spline terms, monotonic ordinal terms and many imputations can take substantial time. Always start with a small quick test before a production run.

---

### Important design notes

This template does **not** use `brm_multiple()` for model fitting.

Instead, it fits one `brms` model per imputed dataset and saves each fit immediately:

```text
fits/fit_imp_001.rds
fits/fit_imp_002.rds
...
fits/fit_imp_100.rds
```

This is intentionally safer for large datasets because:

- the main R session does not hold all fitted models in memory;
- completed fits are preserved if the run stops;
- failed or slow imputations can be re-run separately
- valid existing fits are skipped on re-run; and
- worker processes return only small status objects to the main session.

Parallelism happens across imputations using `future` / `furrr`, with dynamic scheduling. This is configured inside the R scripts:

```r
furrr::furrr_options(
  seed = TRUE,
  scheduling = Inf
)
```

This improves load balancing when some imputed datasets take longer than others.

---

### Supported analysis patterns

The template is designed to support:

- `single-time outcome, row-level covariates`
- `repeated outcome with subject-level covariates`
- `repeated outcome with time-varying covariates`
- `complete-case analysis without imputation`
- `row-level multiple imputation`
- `subject-level multiple imputation`
- `subject-wide imputation using repeated Y as auxiliary variables`

Supported model families include:

- `gaussian`
- `bernoulli`
- `poisson`
- `negbinomial`
- `beta`
- `ordinal`
- `categorical`

Model families and links are set in `00_config.R`.

---

## 2. Structure of the pipeline

The repository is organised around a small set of user-edited files and a sequence of numbered pipeline scripts. In most projects, users only need to edit:

- `00_config.R`
- `00_variable_dictionary.csv`

All other scripts should usually be treated as pipeline code.

---

### Pipeline scripts

The core pipeline is run by `run_all.R`, which calls these scripts in order:

1. `01_validate_config.R`
2. `02_prepare_data.R`
3. `03_impute.R`
4. `04_fit_models.R`
5. `05_diagnostics.R`
6. `06_posterior_summary.R`
7. `07_posterior_prediction.R`
8. `11_check_imputation_stability.R` -- always runs, producing the imputation-count stability tables/figures
9. `08_publication_results.R` -- writes and renders the combined main report

Step 11 runs *before* Step 8, not after, even though Step 11's checks are about whether `m` was large enough -- a question that only makes sense once fitting is done. The reason for the order is that Step 8's report template embeds Step 11's tables and figures directly, in an "Imputation-count stability" chapter (see [How Step 8's report embeds Step 11's results](#how-step-8s-report-embeds-step-11s-results) below), so Step 11's output files need to already exist on disk by the time Step 8 writes and renders the report.

`run_all.R` then automatically runs two more scripts:

- `09_check_mo_parameter_columns.R` -- runs automatically only if the fitted model's formula contains `mo()` terms
- `10_publication_mo_results.R` -- runs automatically only if the fitted model's formula contains `mo()` terms

`run_all.R` detects `mo()` terms directly from the fitted model's own formula (via `extract_special_term_vars()`), so steps 09-10 are skipped automatically for ordinary Gaussian, logistic, spline-only or factor-coded models, and run automatically whenever the formula contains `mo()`. You do not need to remember to run them manually, and you do not need to edit them per study; see [Section 8](#8-publication-outputs-and-inference-guidance) for how they discover monotonic-effect variables generically.

Step 11 runs on whichever final `m` the run settled on. It creates stepwise publication-ready stability summaries, including tables and figures that quantify how much estimates, credible intervals, posterior direction probabilities and transformed summaries change when `m` is increased.

If `analysis_spec$mi_stability$auto_increment = TRUE`, `run_all.R` also drives Steps 3/4/6 through an automatic imputation-count loop instead of fitting a single fixed `m` up front: it fits a small batch of imputations, checks stability against the previous batch, and stops increasing `m` as soon as the configured stability thresholds are met (or once the configured `m` is reached). See [Choosing the number of imputations adaptively](#choosing-the-number-of-imputations-adaptively) for details.

All of `01`-`07`, `11`, `08`, and `09`-`10` (when applicable) can also be run manually and standalone with `Rscript <script>.R`, which is useful for debugging or for resuming a partially completed run. If you run Step 8 manually before Step 11 has ever been run, its report's "Imputation-count stability" chapter is simply omitted, rather than failing.

---

### Main files to edit for a new project

Usually edit only:

- `00_config.R`
- `00_variable_dictionary.csv`

`00_config.R` defines the analysis structure, including:

- outcome variable
- model family and link
- data structure
- imputation strategy
- brms formula settings
- priors
- MCMC settings
- parallel settings
- memory guard settings

The `examples/` folder contains ready-to-run public-data configurations. These are useful for testing the pipeline and for learning how to structure a new analysis.

`00_variable_dictionary.csv` defines:

- variable labels
- variable roles
- variable types
- reference categories
- scaling
- imputation targets
- model inclusion

---

## 3. Quick start

This quick start uses the public Gaussian example based on `datasets::airquality`.

If R, CmdStan or Quarto are not yet installed, see [Section 11](#11-computing-environment-setup) first.

First, obtain a local copy of the template repository.

Using Git:

```bash
# Bash command block
git clone https://github.com/fabregithub/generic-mi-brms-pipeline.git
cd generic-mi-brms-pipeline
```

Alternatively, download the repository as a ZIP file from GitHub, unzip it, and open Terminal in the unzipped project folder.

Then run:

```bash
# Bash command block
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv
Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

Outputs are written to:

- `objects/`
- `fits/`
- `results/`
- `results/publication/`

### Common commands

Once a project is set up (its own `00_config.R` and `00_variable_dictionary.csv` in place), these are the commands you will use most often, whether for the quick-start example above or your own study data.

Validate config. Run in Terminal:

```bash
# Bash command block
Rscript 01_validate_config.R
```

Run all steps. Run in Terminal:

```bash
# Bash command block
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Fit one imputed dataset only. Run in Terminal:

```bash
# Bash command block
Rscript fit_single_imputation.R 1
```

For example, to fit the model only for imputed dataset 51, run in Terminal:

```bash
# Bash command block
Rscript fit_single_imputation.R 51
```

`run_all.R` (via `08_publication_results.R`) renders the main Quarto report automatically, so this step is not normally needed. If you want to re-render it manually, for example after editing the generated `.qmd` by hand, run in Terminal:

```bash
# Bash command block
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

For what to do if a run is interrupted, how to monitor a long run, or how to debug a problematic imputed dataset, see [Section 7](#7-logging-monitoring-restarting-troubleshooting-and-debugging).

---

## 4. Adapting the pipeline to private study data

For a new study, the recommended workflow is:

1. Prepare one clean analysis dataset and save it as an `.rds` file.
2. Edit `00_variable_dictionary.csv`.
3. Edit `00_config.R`.
4. Run validation.
5. Run a small quick test.
6. Run a modest parallel test.
7. Run the full production analysis.
8. Render and inspect publication outputs.

### Data preparation

Before running the pipeline, prepare one clean input dataset, for example `data/my_analysis_data.rds`.

The dataset should already contain consistent variable names, explicit ID variables if needed, explicit time or wave variables for repeated data, and any derived variables that are not created by the pipeline. Save the dataset in R with:

```r
saveRDS(my_data, "data/my_analysis_data.rds")
```

Then point `00_config.R` to it by editing the `data = list(...)` block:

```r
data = list(
  raw_data_file = "data/my_analysis_data.rds",
  ...
)
```

### Check the missing-data mechanism before imputation

Before using the imputation step, review why each variable is missing.

The pipeline uses MICE-style multiple imputation through `miceRanger`. This is suitable when treating the missing values as MCAR or MAR is reasonable after conditioning on observed variables in the imputation model.

It is not a general solution for non-ignorable missingness or censoring. In particular, do not pass the following directly to the pipeline as ordinary `NA` values without careful pre-processing:

- left-censored measurements
- values below a detection limit
- right-censored or interval-censored measurements
- structural missingness
- skip-pattern missingness
- not-applicable responses
- missingness caused by study design
- known MNAR variables

For these cases, first create an analysis-ready representation outside the pipeline. Depending on the scientific context, this might involve:

- using a censored-data model outside this pipeline
- creating a below-detection-limit indicator
- using an appropriate substitution or interval representation
- creating explicit "not applicable" categories
- separating structural missingness from true missingness
- conducting sensitivity analyses for MNAR assumptions

Only variables that can reasonably be treated as ordinary missing values under the chosen imputation model should be set as imputation targets in `00_variable_dictionary.csv`.

### Decision tree: where to start

```text
If your data are one row per person or analytic unit
  -> use data_structure = "single_time"

If your data are one row per subject-time observation
  -> use a repeated-data structure and set id_var and time_var

If your covariates are mostly measured once per subject
  -> consider subject_wide_with_repeated_y_auxiliary imputation

If your predictors change over time
  -> use a repeated/time-varying pattern and check that timing is set correctly

If there are no missing covariates to impute, or if you want a complete-data analysis
  -> edit the imputation = list(...) block so enabled = FALSE,
     strategy = "none", and m = 1

If you need ordinary linear/logistic effects
  -> use fixed_effects = "auto" and control variables through the dictionary

If you need nonlinear continuous effects
  -> use custom_formula with s()

If you need ordinal monotonic effects
  -> mark the variable as type = ordinal, use custom_formula with mo(),
     and run optional scripts 09 and 10 for derived odds-ratio summaries
```

### Important `00_config.R` options

The example configs are the easiest reference, but the following settings are worth checking whenever you adapt the pipeline.

#### Outcome family and link

Set the outcome family and link in `analysis_spec$outcome`:

```r
outcome = list(
  y_var = "low",
  family = "bernoulli",
  link = "logit",
  y_prefix = NULL,
  y_wide_regex = NULL,
  predict_missing_y = TRUE
)
```

Common link examples:

| Family | Common links |
|---|---|
| `gaussian` | `identity`, `log` |
| `bernoulli` | `logit`, `probit`, `cloglog` |
| `poisson` | `log`, `identity` |
| `negbinomial` | `log`, `identity` |
| `beta` | `logit`, `probit`, `cloglog`, `log` |
| `ordinal` | often `logit` or `probit`, depending on the ordinal family |
| `categorical` | often `logit` |

Choose a link that is valid for the selected family and appropriate for the scientific question.

#### `y_prefix` and `y_wide_regex`

These are only needed for repeated-outcome wide imputation. For ordinary single-time analyses, keep both as `NULL`.

For example, if the long outcome is `ps` and `time` is `1, 2, ..., 6`, the subject-wide auxiliary outcome columns may be `ps_1`, `ps_2`, `ps_3`, `ps_4`, `ps_5`, `ps_6`.

Then use:

```r
y_prefix = "ps_"
y_wide_regex = "^ps_"
```

#### Imputation targets and missingness assumptions

The imputation step is intended for ordinary missing data where a MICE-style imputation model is defensible. In the variable dictionary, set `impute_target = TRUE` only for variables to be imputed under that assumption.

If a variable is left-censored, below a detection limit, structurally missing, not applicable by design or likely not missing at random (MNAR), handle that issue before running the pipeline. Such values should not be treated as ordinary `NA` values unless that choice is explicitly justified.

#### Skipping imputation when there is no missing data

If the analysis dataset has no missing covariate data, or if you deliberately want to run a complete-data analysis without multiple imputation, edit the `imputation = list(...)` block in `00_config.R` directly:

```r
imputation = list(
  enabled = FALSE,
  strategy = "none",

  # These are ignored when enabled = FALSE, but keeping m = 1
  # makes the intended one-dataset workflow explicit.
  m = 1,

  maxiter = 0,
  mean_match_k = 5,
  verbose = FALSE,
  impute_y = FALSE,
  extra_exclude_targets = character(0)
)
```


With this setting, the pipeline prepares one model dataset and fits one `brms` model, rather than creating multiple imputed datasets.

#### Imputation seeding and extending an existing run

Two further `imputation = list(...)` settings control reproducibility and safe extension:

```r
imputation = list(
  ...
  # Base seed for miceRanger. Defaults to analysis_spec$model$seed if not set.
  # miceRanger has no native per-chain seed/resume support, so each batch of
  # new imputations (the initial m, or a later extension batch) is seeded
  # deterministically as seed + (number of imputations that already existed
  # before that batch). This reproduces the same imputed datasets across
  # repeated runs, though it is not bit-identical to a hypothetical
  # single-shot miceRanger(m = target_m) call.
  seed = 12345,

  # FALSE by default. Step 3 normally refuses to fit more imputations than
  # already exist on disk, to avoid silently treating a config edit as a
  # request to overwrite existing imputed data. Set this to TRUE to let
  # Step 3 generate only the additional imputations needed to reach the
  # current m, appending them to the existing manifest without touching
  # any existing imputed-data file. run_all.R's automatic m-increment loop
  # (analysis_spec$mi_stability$auto_increment) sets this internally.
  allow_extend = FALSE,
  ...
)
```

For a no-missing-data analysis, the cleanest dictionary setup is usually `impute_target = FALSE` for all variables.

Before disabling imputation, confirm that there are no missing predictor values that the model needs:

```r
source("00_config.R")
d <- readRDS(paths$raw_data)
colSums(is.na(d))
```

If missing predictor values remain and imputation is disabled, those rows may be dropped or model-data checks may fail, depending on the variables and model formula.

#### Variable groups

The recommended default is:

```r
variables = NULL
```

This means variable roles, types, timing, scaling, reference categories, imputation targets, model inclusion and auxiliary-variable status are read from `00_variable_dictionary.csv`.

Advanced users can override selected derived groups in `00_config.R`:

```r
variables = list(
  exposure_vars = c("exposure"),
  covariate_vars = c("age", "sex", "income"),
  auxiliary_vars = c("baseline_score"),

  continuous_vars = c("age", "baseline_score"),
  categorical_vars = c("sex"),
  ordinal_vars = c("income"),

  subject_level_vars = c("age", "sex", "income", "baseline_score"),
  time_varying_vars = c("current_exposure"),

  scale_vars = c("age", "baseline_score")
)
```

For most analyses, it is simpler and safer to keep `variables = NULL` and control the analysis through the dictionary.

#### Initial values

The default example setting is:

```r
init = "random"
```

This uses ordinary Stan-style random initialisation.

A useful alternative for some difficult logistic models is:

```r
init = 0
```

This starts parameters at zero and can sometimes avoid early warm-up overflow warnings. It is not required for ordinary analyses.

Advanced users may also supply a custom initialisation function:

```r
init = function() {
  list()
}
```

Use `init = "random"` unless there is a clear reason to use `init = 0`.

#### Posterior draw extraction regex

Step 6 extracts posterior draws using `analysis_spec$model$parameter_draw_regex`.

For ordinary fixed, random-effect SD and residual parameters only:

```r
parameter_draw_regex = "^(b_|sd_|sigma)"
```

The recommended general default is:

```r
parameter_draw_regex = "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)"
```

This includes:

| Prefix | Meaning |
|---|---|
| `b_` | ordinary fixed effects |
| `bsp_` | monotonic-effect coefficients used by `brms::mo()` |
| `sd_` | group-level standard deviations |
| `sigma` | residual SD for Gaussian models |
| `sds_` | smooth-term standard deviations |
| `bs_` | smooth basis coefficients |
| `simo_` | simplex parameters used by `brms::mo()` |

Use the general default if the model contains `s()` or `mo()` terms.

#### How Step 6 pools draws across imputations

Step 6 does not simply concatenate every imputation's posterior draws and summarise the pooled sample directly. Two corrections are applied per parameter:

1. Each draw is weighted by `1 / (m * K_i)`, where `K_i` is the number of finite draws from imputation `i`. This makes every imputation contribute exactly `1/m` to every summary statistic, regardless of how many draws a particular imputation happened to produce after filtering non-finite values (e.g. for monotonic simplex parameters).
2. A finite-`m` variance correction (the `B/m` term from Rubin's combining rule) is applied on top of the weighted pooled sample, but only when a per-parameter shape diagnostic (a bimodality coefficient) indicates the pooled posterior is unimodal enough for the correction's symmetry assumption to be safe, and only on a support-respecting transform (log for strictly positive parameters, logit for (0, 1)-bounded parameters). For visibly multimodal parameters, the correction is skipped and the uncorrected weighted pooled sample is reported instead, since a symmetric rescale would distort genuine between-imputation structure rather than represent it.

`results/parameter_summary.csv`/`.rds` record this per parameter via the `m_imputations`, `between_var`, `within_var`, `variance_corrected`, `transform_used` and `bimodality_coef` columns, so you can audit whether and how each parameter was corrected.

#### Posterior summaries and ROPE

The `summary` block controls posterior summaries:

```r
summary = list(
  effects = "fixed",
  component = "conditional",
  centrality = "median",
  ci = 0.95,
  ci_method = "HDI",
  test = c("p_direction", "rope"),

  rope = list(
    method = "fixed",
    fixed_range = c(-1, 1),
    width_probability = 0.05
  )
)
```

`p_direction` is the posterior probability of the dominant sign. ROPE means region of practical equivalence and should ideally be defined using scientific or clinical relevance.

Common ROPE options:

```r
# No ROPE summary
rope = list(
  method = "none"
)
```

```r
# Fixed ROPE on the model coefficient scale
rope = list(
  method = "fixed",
  fixed_range = c(-1, 1)
)
```

```r
# Logistic/logit example: odds ratio between 0.95 and 1.05
rope = list(
  method = "fixed",
  fixed_range = log(c(0.95, 1.05))
)
```

```r
# Convenience option for workflows that implement family-specific defaults
rope = list(
  method = "auto",
  width_probability = 0.05
)
```

For Gaussian models, the ROPE is on the outcome/model scale. For Bernoulli/logit models, the ROPE is on the log-odds scale; using `log(c(0.95, 1.05))` corresponds to an odds-ratio range from 0.95 to 1.05.

If no meaningful ROPE has been defined, use:

```r
test = c("p_direction")
rope = list(method = "none")
```

---

### Notes for adapting to private study data

For a new project:

1. Replace or edit `00_variable_dictionary.csv`.
2. Edit `00_config.R`.
3. Run in Terminal:

```bash
# Bash command block
Rscript 01_validate_config.R
```

4. If validation passes, run in Terminal:

```bash
# Bash command block
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```


### Choosing the number of imputations adaptively

There is no universal value of `m` that is automatically sufficient for every analysis. The required number of imputations depends on the fraction of missing information, the amount and pattern of missingness, the target estimands, and how reproducible the posterior summaries need to be. A large dataset does not automatically require `m = 100`, and a small dataset can sometimes require more than `m = 100` if the fraction of missing information is high or the target estimates are unstable. With large datasets in particular, models tend to stabilise at a smaller `m` than with small datasets, while each fit also takes longer -- so paying for a fixed `m = 100` up front can waste a large amount of compute time.

There are two ways to use this adaptive idea: an **automatic loop** that `run_all.R` drives for you, and a **manual staged workflow** that you drive yourself, one batch at a time. The automatic loop is recommended for most analyses; the manual workflow remains useful for closer manual inspection between batches, or for resuming an already-started large run.

#### Automatic m-increment loop (recommended)

Set `analysis_spec$mi_stability$auto_increment <- TRUE` in `00_config.R`:

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

#### Manual staged workflow

If you prefer to inspect the stability check yourself between batches rather than letting `run_all.R` decide automatically, use the manual workflow instead. The safest rule is:

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

#### Running the imputation-count stability script

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
  # the main report's "Imputation-count stability" chapter (see Section 8),
  # which run_all.R renders after this script runs. Set this TRUE only if
  # you also want this script's standalone, more detailed report rendered.
  render_quarto = FALSE
)
```

#### If a larger run has already been started

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

This stability check assesses numerical Monte Carlo stability as `m` increases. It does not validate the missing-data mechanism and does not make MNAR or censored data problems ignorable. For ready-to-adapt Methods text describing this adaptive-`m` justification, see [Section 9](#9-manuscript-writing-guide).

Recommended approach for a new analysis:

### 1. Quick test run

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

### 2. Parallel test run

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

### 3. Full production run

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

### Repeated-outcome data and subject-wide imputation

The generic pipeline can also handle repeated-outcome data where the outcome is measured multiple times per subject and most covariates are subject-level or baseline variables.

Use this imputation strategy in `00_config.R`:

```r
imputation = list(
  ...
  strategy = "subject_wide_with_repeated_y_auxiliary",
  ...
)
```

This strategy does the following:

1. Starts from long data, one row per subject-time observation.
2. Creates a subject-wide imputation dataset, one row per subject.
3. Converts the repeated outcome into wide auxiliary variables, for example `ps_1, ps_2, ps_3, ps_4, ps_5, ps_6`.
4. Imputes subject-level covariates once per subject.
5. Joins the imputed subject-level variables back to the original long outcome data.
6. Fits the `brms` model on the long data.

Required `00_config.R` fields:

- `analysis_spec$data$id_var`
- `analysis_spec$data$time_var`
- `analysis_spec$outcome$y_var`
- `analysis_spec$outcome$y_prefix`
- `analysis_spec$outcome$y_wide_regex`

Example:

```r
data = list(
  raw_data_file = "data/my_repeated_data.rds",
  id_var = "ID",
  row_id_var = "row_id",
  data_structure = "repeated_y_subject_covariates",
  time_var = "time"
)

outcome = list(
  y_var = "ps",
  family = "bernoulli",
  link = "logit",
  y_prefix = "ps_",
  y_wide_regex = "^ps_",
  predict_missing_y = TRUE
)

imputation = list(
  enabled = TRUE,
  strategy = "subject_wide_with_repeated_y_auxiliary",
  m = 24,
  maxiter = 5,
  mean_match_k = 5,
  verbose = FALSE,
  impute_y = FALSE,
  extra_exclude_targets = character(0)
)
```

In `00_variable_dictionary.csv`, subject-level covariates should have:

```r
timing = single
```

or:

```r
timing = baseline
```

Repeated or time-varying variables should have:

```r
timing = repeated
```

or:

```r
timing = time_varying
```

If `impute_y = FALSE`, repeated outcome values are not imputed as targets. They are used as wide auxiliary predictors when possible, and missing outcome rows can later be summarised using posterior prediction.

### Ordinal predictors and `mo()`

Ordinal variables should be marked in the dictionary:

```r
type = ordinal
```

The pipeline converts these variables to ordered factors before modelling. They can then be used in a custom `brms` formula:

```r
custom_formula = brms::bf(
  y ~ mo(education) + mo(income) + age_z + sex + (1 | id)
)
```

For large mixed logistic models, `mo()` can be much slower than ordinary factor coding. A practical approach is to use the factor-coded model as the main analysis and use `mo()` as a sensitivity analysis with fewer imputations.

You do not need to hardcode `mo()` variable names anywhere. `09_check_mo_parameter_columns.R` and `10_publication_mo_results.R` discover them directly from the fitted model's own formula, and `run_all.R` runs both scripts automatically only when `mo()` terms are present. See [Section 8](#8-publication-outputs-and-inference-guidance) for how to optionally supply nicer labels and category levels.

---

### Variable roles: dictionary by default, config as optional override

The pipeline uses `00_variable_dictionary.csv` as the default source of truth for variable metadata.

This means the following information should normally be specified only in the dictionary:

- `role`
- `type`
- `timing`
- `scale`
- `reference`
- `impute_target`
- `use_in_model`
- `use_as_auxiliary`

Therefore, the standard setting in `00_config.R` is:

```r
variables = NULL
```

The pipeline automatically derives internal variable groups from the dictionary, including:

- `exposure_vars`
- `covariate_vars`
- `auxiliary_vars`
- `continuous_vars`
- `categorical_vars`
- `ordinal_vars`
- `subject_level_vars`
- `time_varying_vars`
- `scale_vars`

For example, if `00_variable_dictionary.csv` contains:

```text
var,role,type,timing,scale,use_in_model
Solar.R,covariate,continuous,single,z,TRUE
Wind,covariate,continuous,single,z,TRUE
Temp,covariate,continuous,single,z,TRUE
Month,covariate,categorical,single,no,TRUE
```

and `analysis_spec$model$fixed_effects` is:

```r
fixed_effects = "auto"
```

the model uses `Solar.R_z + Wind_z + Temp_z + Month`.

### Optional overrides

Advanced users can override selected derived groups in `00_config.R`.

Example:

```r
variables = list(
  scale_vars = c("age", "income"),
  auxiliary_vars = c("baseline_score")
)
```

If `variables = NULL`, no override is applied.

If `exposure_vars` or `covariate_vars` are explicitly supplied in `analysis_spec$variables`, they are used by `fixed_effects = "auto"` instead of the dictionary-derived `use_in_model` list. For most analyses, it is simpler and safer to leave `variables = NULL` and control the model through the dictionary.

---

# 5. Variable dictionary

The file `00_variable_dictionary.csv` is the main machine-readable description of the analysis variables.

Expected columns are:

- `var`
- `label`
- `role`
- `type`
- `timing`
- `scale`
- `reference`
- `impute_target`
- `use_in_model`
- `use_as_auxiliary`

Some fields are straightforward:

- `var`: exact variable name in the dataset
- `label`: human-readable label for tables and reports
- `impute_target`: whether this variable should be imputed when missing
- `use_in_model`: whether this variable should appear in the final brms model

Set `impute_target = TRUE` only for variables whose missing values should be treated as ordinary missing data in the MICE-style imputation model. Do not set `impute_target = TRUE` for left-censored values, below-detection-limit values, structural missingness, skip-pattern missingness or known MNAR variables unless those issues have already been handled appropriately before the pipeline.

The other fields need more explanation.

### `role`

`role` describes the analytical purpose of the variable.

Recommended values:

| Value | Meaning | Example |
|---|---|---|
| `outcome` | Continuous or general outcome | `Ozone` |
| `binary_outcome` | Binary outcome for Bernoulli/logistic models | `low` |
| `exposure` | Main exposure or predictor of scientific interest | treatment group, air pollution |
| `covariate` | Adjustment variable / confounder / predictor | age, sex, income |
| `auxiliary` | Used for imputation only, not included in final model | extra baseline score |
| `id` | Subject, cluster or row identifier | `ID`, `row_id` |
| `time` | Measurement occasion, wave, visit or follow-up time | `time`, `wave` |
| `cluster` | Grouping variable for random effects or clustering | school ID, hospital ID |
| `strata` | Stratification variable, if relevant | site, cohort |

For most standard regression analyses, the most common roles are:

- `outcome`
- `binary_outcome`
- `exposure`
- `covariate`
- `id`
- `time`
- `auxiliary`

Example:

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

### `type`

`type` describes the statistical data type.

Recommended values:

| Value | Meaning | Typical R class |
|---|---|---|
| `continuous` | Numeric variable on a continuous scale | numeric |
| `integer` | Count-like integer variable | integer/numeric |
| `binary` | Two-level variable, often 0/1 | factor or integer |
| `categorical` | Unordered categorical variable with 3+ levels | factor |
| `ordinal` | Ordered categorical variable | ordered factor |
| `date` | Calendar date | Date |
| `id` | Identifier, not modelled as numeric | character/factor/integer |

Examples:

```text
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
```

Important notes:

- Binary and categorical variables should usually be factors before modelling.
- Continuous variables with `scale = z` get a new `_z` version for modelling.
- ID variables should not be scaled or imputed unless there is a specific reason.

### `timing`

`timing` describes when the variable was measured.

Recommended values:

| Value | Meaning | Example |
|---|---|---|
| `single` | Measured once per analytic unit | sex, baseline income |
| `baseline` | Measured at baseline only | baseline age |
| `repeated` | Measured repeatedly across time/visits | repeated outcome |
| `time_varying` | Predictor changes over time | current exposure, current medication |
| `derived` | Created from other variables | standardised score |
| `id` | Identifier variable | subject ID |

Examples:

```text
ps,Psychological stress,binary_outcome,binary,repeated,no,0,FALSE,TRUE,TRUE
time,Measurement wave,time,integer,repeated,no,,FALSE,TRUE,FALSE
c_sex,Child sex,covariate,binary,baseline,no,0,TRUE,TRUE,FALSE
```

How this matters:

- For `single_time` datasets, most variables usually have `timing = single`.
- For repeated-outcome data with subject-level covariates, the outcome and time variable are usually `repeated`, while most covariates are `baseline` or `single`.
- For repeated data with time-varying predictors, mark those predictors as `time_varying`.

### `scale`

`scale` tells the pipeline whether and how to transform a variable before modelling.

Recommended values:

| Value | Meaning | Result |
|---|---|---|
| `no` | No scaling or transformation | original variable used |
| `z` | Standardise to mean 0 and SD 1 | creates `var_z` |
| `centre` | Mean-centre only | creates centred version if supported |
| `log` | Log-transform | creates log version if supported |
| `custom` | User-defined transform outside dictionary | handled in config/functions |

Currently, the most important supported value is `z`.

Example:

```text
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
```

This creates `age_z`, and the model should use `age_z`, not `age`.

Use `z` for continuous predictors when coefficients should represent the effect per 1 SD increase. This is usually helpful for Bayesian models because weakly informative priors like `normal(0, 1)` or `normal(0, 1.5)` are easier to interpret on standardised predictors.

Do not use `z` for:

- outcomes
- ID variables
- binary/categorical variables
- already standardised variables

### `reference`

`reference` defines the reference category for binary, categorical or ordinal variables.

Examples:

```text
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
```

For a binary variable coded `0/1`, setting:

```r
reference = 0
```

means the coefficient compares level `1` against level `0`.

For categorical variables, choose the scientifically meaningful or most common category as the reference.

Important notes:

- The value in `reference` must match an actual value/level in the data.
- Leave `reference` blank for continuous variables.
- Reference levels affect coefficient interpretation but not overall model fit.

### `use_as_auxiliary`

`use_as_auxiliary` controls whether a variable is used in imputation but excluded from the final analysis model.

This is useful for variables that help predict missingness or missing values, but are not part of the scientific model.

Examples:

```text
extra_score,Extra baseline score,auxiliary,continuous,baseline,z,,TRUE,FALSE,TRUE
hospital_id,Hospital ID,auxiliary,categorical,single,no,,FALSE,FALSE,TRUE
```

Interpretation:

| impute_target | use_in_model | use_as_auxiliary | Meaning |
|---|---|---|---|
| `TRUE` | `TRUE` | `FALSE` | Impute if missing and include in model |
| `TRUE` | `FALSE` | `TRUE` | Impute/use for imputation but exclude from model |
| `FALSE` | `FALSE` | `TRUE` | Use as imputation predictor only |
| `FALSE` | `TRUE` | `FALSE` | Include in model but do not impute |
| `FALSE` | `FALSE` | `FALSE` | Keep as metadata/ID or ignore analytically |

Common use cases:

1. Auxiliary predictor only:

```text
baseline_score,Baseline questionnaire score,auxiliary,continuous,baseline,z,,FALSE,FALSE,TRUE
```

2. Variable to impute and include:

```text
income,Household income,covariate,categorical,baseline,no,1,TRUE,TRUE,FALSE
```

3. Outcome not imputed:

```text
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
```

4. Repeated outcome as imputation auxiliary:

```text
ps,Psychological stress,binary_outcome,binary,repeated,no,0,FALSE,TRUE,TRUE
```

Whether repeated outcomes are used as auxiliary predictors depends on the imputation strategy in `00_config.R`.

### Example: Gaussian demo

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
Ozone,Ozone concentration,outcome,continuous,single,no,,FALSE,TRUE,FALSE
Solar.R,Solar radiation,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Wind,Wind speed,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Temp,Temperature,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Month,Month,covariate,categorical,single,no,5,FALSE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

### Example: Logistic demo

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
lwt,Maternal weight at last menstrual period,covariate,continuous,single,z,,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ptl,Previous premature labours,covariate,continuous,single,z,,TRUE,TRUE,FALSE
ht,History of hypertension,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ui,Uterine irritability,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ftv,Physician visits during first trimester,covariate,continuous,single,z,,TRUE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

### Example: Spline + monotonic demo

The `birthwt_spline_monotonic` example is used to test custom `brms` formula terms: `s(age_z, k = 5)` and `mo(lwt_q)`.

The corresponding dictionary demonstrates two important ideas:

- `scale = z` creates `age_z`, `ptl_z` and `ftv_z`
- `type = ordinal` creates an ordered factor suitable for `mo()`

Example:

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
lwt_q,Maternal weight quintile,covariate,ordinal,single,no,1,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ptl,Previous premature labours,covariate,continuous,single,z,,TRUE,TRUE,FALSE
ht,History of hypertension,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ui,Uterine irritability,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ftv,Physician visits during first trimester,covariate,continuous,single,z,,TRUE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

In `00_config.R`, this dictionary is paired with a custom formula:

```r
custom_formula = brms::bf(
  low ~ s(age_z, k = 5) + mo(lwt_q) + race + smoke + ptl_z + ht + ui + ftv_z
)
```

Notes:

- `lwt_q` is created by the example-data script before the pipeline runs.
- `lwt_q` is marked as `type = ordinal`, so the pipeline converts it to an ordered factor.
- `mo(lwt_q)` then models the ordinal effect as monotonic but not necessarily equally spaced.
- The posterior draw regex for `mo()` models should include both `bsp_` and `simo_` parameters:

```r
model = list(
  ...
  parameter_draw_regex = "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)",
  ...
)
```

---

## 6. Parallelisation and performance tuning

The pipeline has several levels of parallelisation.

| Setting | Used in | Meaning |
|---|---|---|
| `impute_workers` | Step 3 | Number of parallel `miceRanger` workers |
| `num_impute_threads_per_worker` | Step 3 | Threads used by each imputation worker |
| `fit_workers` | Step 4 | Number of imputed datasets fitted in parallel |
| `cores_per_fit` | Step 4 | Number of chains/cores per `brms` fit |
| `summary_workers` | Step 6 | Number of workers used for posterior-draw extraction |
| `prediction_workers` | Step 7 | Number of workers used for posterior prediction |
| `future_globals_maxsize_gb` | Steps 4, 6, 7 | Maximum future globals size |

Recommended order inside `analysis_spec$parallel`:

```r
parallel = list(
  # miceRanger imputation workers and threads
  impute_workers = 2,
  num_impute_threads_per_worker = 2,

  # Backward-compatible fallback used by older scripts.
  # This can be removed later if no scripts reference it.
  num_impute_threads = 2,

  # Model fitting workers and cores
  fit_workers = 4,
  cores_per_fit = 4,

  # Step 6 and Step 7 workers
  summary_workers = 2,
  prediction_workers = 2,

  future_globals_maxsize_gb = 80
)
```

For large `brmsfit` objects, start conservatively with `summary_workers = 2` and `prediction_workers = 2`. Increase only if memory is comfortable.

If you use the automatic m-increment loop (`analysis_spec$mi_stability$auto_increment = TRUE`, see [Choosing the number of imputations adaptively](#choosing-the-number-of-imputations-adaptively)), `fit_workers` also determines the default batch size: each batch defaults to `fit_workers` imputations, so every batch fully occupies the parallel workers with none left idle. With `fit_workers = 4`, the default batches are `m = 4, 8, 12, ...`; with `fit_workers = 6`, they are `m = 6, 12, 18, ...`. A custom `increment_size` that is not already a multiple of `fit_workers` is rounded up automatically.

`num_impute_threads` is a backward-compatible fallback for older scripts. If all scripts in the repository use `impute_workers` and `num_impute_threads_per_worker`, it can eventually be removed.

### Suggested settings by computing environment

#### Laptop or low-memory desktop

```r
parallel = list(
  ...
  impute_workers = 1,
  num_impute_threads_per_worker = 1,
  fit_workers = 1,
  cores_per_fit = 1,
  summary_workers = 1,
  prediction_workers = 1,
  future_globals_maxsize_gb = 8
)
```

#### Standard desktop or small workstation

```r
parallel = list(
  ...
  impute_workers = 2,
  num_impute_threads_per_worker = 2,
  fit_workers = 2,
  cores_per_fit = 4,
  summary_workers = 2,
  prediction_workers = 2,
  future_globals_maxsize_gb = 20
)
```

#### High-memory workstation

```r
parallel = list(
  ...
  impute_workers = 4,
  num_impute_threads_per_worker = 4,
  fit_workers = 4,
  cores_per_fit = 4,
  summary_workers = 4,
  prediction_workers = 4,
  future_globals_maxsize_gb = 80
)
```

#### Shared server or high-performance computing environment

Use conservative per-job settings unless your scheduler explicitly allocates more resources. Avoid requesting more threads than the scheduler has allocated to the job.

---

### Parallel miceRanger imputation

The pipeline can parallelise `miceRanger` imputation using `doParallel` and `foreach`.

In `00_config.R`, the relevant settings are:

```r
parallel = list(
  ...
  impute_workers = 4,
  num_impute_threads_per_worker = 4,
  ...
)
```

The approximate CPU demand during imputation is `impute_workers * num_impute_threads_per_worker`. For example, 4 workers * 4 threads per worker = about 16 active threads.

Start conservatively, especially on laptops or when the imputation data are large, because parallel `miceRanger` can copy data to worker processes and increase memory use.

Recommended starting values:

```r
# Public examples or ordinary laptops
parallel = list(
  ...
  impute_workers = 1,
  num_impute_threads_per_worker = 1,
  ...
)

# High-memory workstation
parallel = list(
  ...
  impute_workers = 4,
  num_impute_threads_per_worker = 4,
  ...
)
```

To test parallel imputation from scratch, remove old imputation outputs first:

```bash
# Bash command block
rm -f objects/imputation_manifest.rds
rm -f objects/imputation_spec.rds
rm -rf objects/imputed_data
rm -rf objects/imputed_wide
rm -rf objects/model_data
rm -rf fits results
rm -f pipeline_error.flag pipeline_success.flag
```

Then re-run the pipeline.

---

### Recommended high-performance settings

The following settings were tested successfully on a high-memory machine. Edit these in `00_config.R` using RStudio or another text editor:

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
```

This runs:

- 100 imputations
- 4 chains per model
- 4 imputed datasets fitted in parallel
- 16 active chains total

For a new analysis, keep the following setting in `00_config.R`:

```r
model = list(
  ...
  run_smoke_fit = TRUE,
  ...
)
```

The smoke fit runs one sequential model first, before launching parallel workers. This catches formula, prior, data or CmdStan problems early.

If you use the automatic m-increment loop (`analysis_spec$mi_stability$auto_increment = TRUE`), you do not need to manually toggle this setting between batches: the loop runs the smoke fit only once, for the first batch, and disables it automatically for every later batch in the same `run_all.R` invocation. See [Automatic m-increment loop](#automatic-m-increment-loop-recommended).

After a configuration has been tested successfully, you may set the following in `00_config.R`:

```r
model = list(
  ...
  run_smoke_fit = FALSE,
  ...
)
```

to save a little time.

---

## 7. Logging, monitoring, restarting, troubleshooting and debugging

This section covers what to do once a run is already underway or has already produced some output: restarting after an interruption, monitoring progress, recovering from CmdStan cache problems, and debugging a specific imputed dataset. For the commands used to start a run in the first place, see [Quick start](#3-quick-start).

---

### Restarting after interruption

The pipeline is checkpointed. If a run is interrupted, run the following in Terminal:

```bash
# Bash command block
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Existing valid fit files are skipped.

To clean only fitting and downstream results while keeping prepared/imputed data, run in Terminal:

```bash
# Bash command block
rm -f fits/fit_imp_*.rds
rm -f objects/fit_manifest.rds
rm -f objects/fit_status.rds
rm -f objects/fit_smoke_status.rds
rm -f results/fit_status.csv
rm -f results/fit_smoke_status.csv
rm -f results/worker_logs/fit_worker_imp_*.log

rm -f results/parameter_draws.rds
rm -f results/parameter_summary.rds
rm -f results/parameter_summary.csv
rm -f results/parameter_draws_imp_*.rds
rm -f objects/parameter_manifest.rds

rm -f results/missing_y_draws.rds
rm -f results/missing_y_summary.rds
rm -f results/missing_y_summary.csv
rm -f results/missing_y_draws_imp_*.rds

rm -f pipeline_error.flag
rm -f pipeline_success.flag
rm -f run_all_stdout.log
```
or run the Bash command `bash 99_clean_fitting_results.sh`.

To clean imputation and all downstream outputs, run in Terminal:

```bash
# Bash command block
rm -f objects/imputation_manifest.rds
rm -f objects/imputed_data/*.rds
rm -f objects/imputed_wide/*.rds

rm -f objects/model_data_manifest.rds
rm -f objects/model_data/*.rds

rm -f fits/fit_imp_*.rds
rm -f objects/fit_manifest.rds
rm -f objects/fit_status.rds
rm -f objects/fit_smoke_status.rds
rm -f results/fit_status.csv
rm -f results/fit_smoke_status.csv
rm -f results/worker_logs/fit_worker_imp_*.log

rm -f results/parameter_draws.rds
rm -f results/parameter_summary.rds
rm -f results/parameter_summary.csv
rm -f results/parameter_draws_imp_*.rds
rm -f objects/parameter_manifest.rds

rm -f results/missing_y_draws.rds
rm -f results/missing_y_summary.rds
rm -f results/missing_y_summary.csv
rm -f results/missing_y_draws_imp_*.rds

rm -f pipeline_error.flag
rm -f pipeline_success.flag
rm -f run_all_stdout.log
```
or run the Bash command `bash 99_cleanall.sh`.

---

### Logging and monitoring

The pipeline writes logs and status files:

- `pipeline_progress.log`
- `pipeline_stdout.log`
- `pipeline_heartbeat.txt`
- `pipeline_success.flag`
- `pipeline_error.flag`
- `results/worker_logs/`

Monitor progress by running the following in Terminal:

```bash
# Bash command block
tail -f pipeline_progress.log
```

Inspect worker-level fitting logs by running the following in Terminal:

```bash
# Bash command block
ls -lh results/worker_logs
cat results/worker_logs/fit_worker_imp_001.log
```

If the pipeline completes successfully, the file `pipeline_success.flag` is created. If an R-level error is caught instead, `pipeline_error.flag` is created.

If the machine crashes or restarts unexpectedly, there may be no error flag. In that case, check the heartbeat and logs, then re-run. The checkpoint system should skip completed valid fits.

---

### Troubleshooting CmdStan cache issues

If a simple model fails with an error like:

```text
Fitting failed. Unable to retrieve the metadata.
No chains finished successfully. Unable to retrieve the fit.
```

and the data look valid, try clearing the CmdStanR cache. Run in Terminal:

```bash
# Bash command block
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

Then re-run the direct debug fit or the pipeline.

This can resolve stale or corrupted compiled-model cache issues.

---

### Debugging model fitting

A good sequence is:

1. Validate the config. Run in Terminal:

```bash
# Bash command block
Rscript 01_validate_config.R
```

2. Run one direct or smoke fit.

3. If successful, run the full pipeline.

4. If a specific imputed dataset is slow or problematic during model fitting, fit it alone. Run in Terminal:

```bash
# Bash command block
Rscript fit_single_imputation.R 51
```

You can temporarily skip or restrict **model fitting for selected imputed datasets** in `00_config.R`.

These options do **not** skip the imputation step itself. They only control which already-created imputed datasets are passed to Step 4, where `brms` models are fitted.

To fit all imputations except imputations 45 and 51:

```r
model = list(
  ...
  skip_imputations = c(45, 51),
  only_imputations = integer(0),
  ...
)
```

To fit only imputation 51, for example as a diagnostic re-run:

```r
model = list(
  ...
  skip_imputations = integer(0),
  only_imputations = c(51),
  ...
)
```

Interpretation:

- `skip_imputations`: do not fit `brms` models for these imputed datasets
- `only_imputations`: fit `brms` models only for these imputed datasets

These settings are useful when one imputed dataset is unusually slow or problematic. Completed valid fit files are still preserved and skipped on re-run.

---

# 8. Publication outputs and inference guidance

After successful completion, publication outputs are written to `results/publication/`.

Typical outputs include:

- `results/publication/tables/main_effect_table_display.csv`
- `results/publication/tables/main_effect_table_full.csv`
- `results/publication/tables/diagnostics_summary.csv`
- `results/publication/tables/analysis_metadata.csv`
- `results/publication/tables/analysis_metadata.rds`
- `results/publication/figures/forest_plot_odds_ratios.png`
- `results/publication/report/bayesian_mi_report_template.qmd`

The generated Quarto report includes posterior results, diagnostics, figures, an imputation-count stability chapter (see below), and a methods/settings table based on `analysis_metadata.csv`. This table records key analysis settings such as the imputation strategy, target number of imputations, fitted imputations used in posterior summaries, model formula, family/link, priors, MCMC settings, parallel settings, posterior-summary settings and predictive-draw settings.

`08_publication_results.R` renders this report itself, to both HTML and DOCX, as the last thing it does. Running `run_all.R` therefore always produces a finished, rendered report with no separate manual step. If you need to re-render it by hand (for example after manually editing the `.qmd`), run:

```bash
# Bash command block
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

Self-rendering can be turned off with `analysis_spec$publication$render_quarto <- FALSE` in `00_config.R`, if you only want the `.qmd` written and prefer to render it yourself.

---

### How Step 8's report embeds Step 11's results

The main report's "Imputation-count stability" chapter is not a separate report -- it is built from the same tables and figures that `11_check_imputation_stability.R` writes to `results/publication/mi_stability/`, referenced from the main report's `.qmd` by relative path. There is one report, one `.qmd`, and one rendered HTML/DOCX pair; `results/publication/mi_stability/` still exists alongside it as the underlying data (and as a `.qmd` of its own, for anyone who wants the full, unabridged stepwise detail across every evaluated batch), but it is not rendered separately.

This is why `run_all.R` runs Step 11 before Step 8 (see [Pipeline scripts](#pipeline-scripts)): the embedded chunks are evaluated when the report is rendered, so Step 11's files need to already be on disk at that point. `11_check_imputation_stability.R`'s own `render_quarto` option now defaults to `FALSE` for this reason -- its standalone report would just duplicate what is already in the main report's chapter. Set `analysis_spec$mi_stability$render_quarto <- TRUE` if you specifically want that standalone, more detailed report rendered as well.

---

### Optional monotonic-effect post-processing

For models that use `brms::mo()`, the standard posterior parameter table is not always the most interpretable summary. Monotonic effects are parameterised using an overall monotonic coefficient and simplex parameters, so category-specific odds ratios should be derived from posterior draws.

Two scripts are provided for this purpose:

- `09_check_mo_parameter_columns.R`
- `10_publication_mo_results.R`

`run_all.R` runs both of these automatically, but only when the fitted model's formula actually contains `mo()` terms (detected directly from `model_spec$formula`, not from any hardcoded variable list). For ordinary Gaussian, logistic, spline-only or factor-coded models, they are skipped automatically. You can also run them manually after the main pipeline has completed and `results/parameter_draws.rds` has been created:

```bash
# Bash command block
Rscript 09_check_mo_parameter_columns.R
Rscript 10_publication_mo_results.R
```

`09_check_mo_parameter_columns.R` discovers every `mo()` variable from the model formula, reports the matching `bsp_`/`simo_` columns and the implied number of ordered categories, and prints a ready-to-paste `analysis_spec$mo_effects` config block.

By default, `10_publication_mo_results.R` runs with **zero configuration**: it discovers the same `mo()` variables from the formula and labels categories generically as `"Level 1"`, `"Level 2"`, and so on. For publication-ready labels and category names, paste the block 09 printed into `00_config.R` and edit it, for example:

```r
analysis_spec$mo_effects <- list(
  vars = list(
    lwt_q = list(
      label = "Maternal weight quintile",
      levels = c("Q1", "Q2", "Q3", "Q4", "Q5")
    )
  ),

  # Only needed if your formula contains a time * mo(variable) interaction.
  time_var = NULL,
  time_values = NULL
)
```

Any `mo()` variable not listed in `analysis_spec$mo_effects$vars` still gets generic "Level N" labels automatically, so this config is optional and additive -- you only need to add entries for the variables where you want nicer labels.

Optional Quarto rendering:

```bash
# Bash command block
quarto render results/publication/mo_effects/report/mo_effects_report.qmd
```

The main output table is `results/publication/mo_effects/tables/mo_cumulative_or_table.csv`.

Additional outputs include:

- `results/publication/mo_effects/tables/mo_adjacent_or_table.csv`
- `results/publication/mo_effects/tables/mo_average_or_table.csv`
- `results/publication/mo_effects/tables/mo_simplex_table.csv`
- `results/publication/mo_effects/figures/mo_cumulative_or_plot.png`
- `results/publication/mo_effects/figures/mo_adjacent_or_plot.png`
- `results/publication/mo_effects/report/mo_effects_report.qmd`

For models using `mo()`, make sure `00_config.R` includes `bsp_` and `simo_` in the posterior draw extraction regex:

```r
model = list(
  ...
  parameter_draw_regex = "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)",
  ...
)
```

The `10_publication_mo_results.R` script can handle monotonic effects with interactions such as `time * mo(ordinal_variable)`.

In that case, set `analysis_spec$mo_effects$time_var`/`time_values`, and it calculates category-specific monotonic-effect odds ratios at the configured time values.

These scripts are not required for ordinary Gaussian, logistic, spline-only or factor-coded models, and `run_all.R` skips them automatically in that case. They only run when derived monotonic-effect odds ratios are needed.


### Optional imputation-count stability outputs

If `11_check_imputation_stability.R` is run, publication-ready imputation-count stability outputs are written to `results/publication/mi_stability/`.

These outputs summarise how selected posterior summaries change as the number of completed imputations increases. They are intended to support transparent reporting when an adaptive number of imputations is used, especially for expensive Bayesian models where blindly running `m = 100` or more may be impractical.

Main outputs include:

- `results/publication/mi_stability/tables/imputation_stability_all_batches.csv`
- `results/publication/mi_stability/tables/imputation_stability_final_comparison_display.csv`
- `results/publication/mi_stability/tables/imputation_stability_settings.csv`
- `results/publication/mi_stability/tables/imputation_stability_stepwise_summary.csv`
- `results/publication/mi_stability/tables/imputation_stability_stepwise_comparison_full.csv`
- `results/publication/mi_stability/figures/imputation_stability_trajectories.png`
- `results/publication/mi_stability/figures/imputation_stability_stepwise_change.png`
- `results/publication/mi_stability/report/imputation_stability_report.qmd`

The stability check should focus on prespecified primary estimands rather than every nuisance parameter.

---

# 9. Manuscript writing guide

The publication outputs in [Section 8](#8-publication-outputs-and-inference-guidance) are designed to be quoted directly rather than re-derived by hand. This section maps each manuscript paragraph to the file that supplies its numbers, and gives fill-in-the-blank templates for Methods and Results text.

These templates are starting points, not finished prose. They still need scientific interpretation, comparison to prior literature, and journal-specific formatting, all of which require human judgement that this pipeline does not provide.

### Where each paragraph's numbers come from

| Manuscript paragraph | Source file | What it gives you |
|---|---|---|
| Missing data and imputation | `results/publication/tables/analysis_metadata.csv` | Imputation strategy, target/used `m`, imputation iterations, mean-matching candidates |
| Model specification | `results/publication/tables/analysis_metadata.csv` | Model formula, family/link, priors |
| MCMC settings | `results/publication/tables/analysis_metadata.csv` | Chains, iterations, warm-up, seed, `adapt_delta`, `max_treedepth` |
| Convergence/diagnostics | `results/publication/tables/diagnostics_summary.csv` | Divergences, treedepth hits, Rhat/ESS-type checks per fit |
| Adaptive imputation-count justification | `results/publication/mi_stability/tables/imputation_stability_stepwise_summary_display.csv` | Whether/when posterior summaries stopped changing as `m` increased |
| Main fixed-effect results | `results/publication/tables/main_effect_table_display.csv` | Estimate, 95% CrI, `pd`, ROPE %, and the transformed effect (e.g. odds ratio) per variable |
| Smooth/monotonic supplementary results | `results/publication/tables/special_parameter_table.csv` | Auxiliary `s()`/`mo()` parameters not shown in the main table |
| Monotonic-effect (`mo()`) results | `results/publication/mo_effects/tables/mo_cumulative_or_table.csv` | Category-specific odds ratios for monotonic ordinal predictors |
| Posterior prediction for missing outcomes | `results/missing_y_summary.csv` | Summaries of predicted values for rows with missing outcome data |

Open `analysis_metadata.csv` directly to read off exact values for the templates below, for example:

```r
metadata <- readr::read_csv("results/publication/tables/analysis_metadata.csv")
print(metadata, n = Inf)
```

### Methods text templates

**Missing data and imputation.** Adapt to the actual `Imputation strategy` and `m` values from `analysis_metadata.csv`:

> Missing covariate data were handled using multiple imputation by chained
equations implemented with random forests (miceRanger). Variables with
missing values were imputed under a missing-at-random assumption,
conditional on the variables included in the imputation model. We
generated m = XX imputed datasets, fitted a separate Bayesian model to
each, and combined posterior draws across imputations.

If you used the subject-wide repeated-outcome strategy, adapt instead:

> For the repeated outcome, subject-level covariates were imputed once per
subject using a subject-wide imputation dataset, then merged back onto
the long-format repeated-measures data before model fitting.

**Adaptive number of imputations.** If you used the automatic m-increment loop or the manual staged workflow described in [Section 4](#4-adapting-the-pipeline-to-private-study-data):

> We used an adaptive multiple-imputation strategy because each imputed-data
Bayesian model was computationally expensive. We first fitted an initial
set of imputed datasets and assessed the stability of prespecified
primary posterior summaries. We increased the number of imputations until
the pooled posterior medians, credible intervals, posterior direction
probabilities and substantive conclusions changed negligibly with
additional imputations. The final analysis used m = XX imputations.

If the stepwise stability table (`imputation_stability_stepwise_summary_display.csv`) shows early flattening, a more specific statement can be used instead:

> The stability assessment showed that posterior summaries changed
negligibly after m = XX. The largest subsequent changes in posterior
medians, credible-interval endpoints and odds-ratio summaries were below
prespecified practical thresholds. The final analysis therefore used
m = XX imputations.

**Model specification.** Read `Model formula`, `Model family/link` and `Priors` from `analysis_metadata.csv`:

> We fitted a Bayesian [family] regression model with a [link] link function
using brms (formula: [paste Model formula here]) and the cmdstanr backend
for Stan. Priors were [paste Priors here]. One model was fitted separately
to each imputed dataset rather than using brm_multiple(), so that fitted
models could be checkpointed and restarted independently.

**MCMC settings.** Read `MCMC chains`, `Total/Warm-up/Post-warm-up iterations`, `Seed`, `adapt_delta` and `max_treedepth` from `analysis_metadata.csv`:

>Each model was fitted using XX chains of XX total iterations, including
XX warm-up iterations, yielding XX post-warm-up draws per chain per
imputed dataset (XX total post-warm-up draws pooled across XX imputations).

**Posterior pooling across imputations.** Worth including given the corrected pooling method (see [Section 4](#4-adapting-the-pipeline-to-private-study-data), "How Step 6 pools draws across imputations"):

> Posterior draws were pooled across imputed datasets with equal weight per
imputation, and an additional finite-imputation variance correction
(following Rubin's combining rules) was applied where the pooled
posterior distribution was sufficiently unimodal for this correction to
be appropriate.

**Software and reproducibility.** Fill in package versions from your R session (`sessionInfo()`) and the `Seed` value from `analysis_metadata.csv`:

> Analyses were conducted in R version XX using brms version XX with the
cmdstanr backend (CmdStan version XX), and miceRanger version XX for
multiple imputation. A fixed random seed (XX) was used for model fitting;
imputation batches were seeded deterministically for reproducibility.

### Results text templates

**Reporting a main effect**, using a row from `main_effect_table_display.csv`:

> [Variable label] was associated with [outcome]: estimate = XX (95% credible interval [CrI]: XX to XX; posterior probability of direction = XX). [If a transformed effect column is present, e.g. odds ratio:] On the odds-ratio scale, this corresponds to an odds ratio of XX (95% CrI: XX to XX).


If a ROPE was defined and is reported in the `ROPE %` column:

> The posterior probability that this effect fell within the prespecified region of practical equivalence (ROPE: XX to XX) was XX%, [supporting / not supporting] practical equivalence to no effect.


**Reporting a monotonic (`mo()`) effect**, using a row from `mo_cumulative_or_table.csv` (see [Section 8](#8-publication-outputs-and-inference-guidance) for how these are derived):

> [Variable label] showed a monotonic ordinal association with [outcome]. Compared with the lowest category, the odds ratio for [highest category] was XX (95% CrI: XX to XX; Pr(OR > 1) = XX), with intermediate categories shown in [Table/Figure XX].


### What this guide does not write for you

- Scientific interpretation of effect sizes
- Clinical or substantive significance
- Comparison with prior literature
- Study limitations
- Sample size or power justification
- Ethics, consent and data-availability statements
- Citations for software and methods (cite R, brms, Stan, miceRanger and the relevant statistical methods papers directly)

---

# 10. Examples and tests

The repository includes three public example analyses. These are intended to help users test the full pipeline before applying it to private study data.

The example folders are:

```text
examples
├── airquality_gaussian
│   ├── 00_config_airquality_gaussian.R
│   ├── 00_create_airquality_example_data.R
│   └── 00_variable_dictionary_airquality_gaussian.csv
├── birthwt_logistic
│   ├── 00_config_birthwt_logistic.R
│   ├── 00_create_birthwt_logistic_example_data.R
│   └── 00_variable_dictionary_birthwt_logistic.csv
└── birthwt_spline_monotonic
    ├── 00_config_birthwt_spline_monotonic.R
    ├── 00_create_birthwt_spline_monotonic_example_data.R
    ├── 00_variable_dictionary_birthwt_spline_monotonic.csv
    └── README_birthwt_spline_monotonic.md
```

Each example contains:

- a config file
- a variable dictionary
- a data-creation script

To use an example, copy its config file and variable dictionary to the project root as `00_config.R` and `00_variable_dictionary.csv`.

Then run the example data-creation script.

### Example 1: Gaussian model using `datasets::airquality`

This is the default Gaussian example.

It tests:

- continuous outcome
- `gaussian(identity)` model
- row-level imputation
- ordinary fixed-effect reporting
- posterior prediction for rows with missing outcome values

Model outline:

- **Dataset:** `datasets::airquality`
- **Outcome:** `Ozone`
- **Model family:** `gaussian(identity)`
- **Model:** `Ozone ~ Solar.R_z + Wind_z + Temp_z + Month`

Run in Terminal from the project root:

```bash
# Bash command block
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv

Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

`run_all.R` already rendered `results/publication/report/bayesian_mi_report_template.html`/`.docx`; no separate render step is needed.

### Example 2: Logistic model using `MASS::birthwt`

This example tests the Bernoulli/logit workflow with a binary outcome.

It tests:

- binary outcome
- `bernoulli(logit)` model
- row-level imputation
- odds-ratio reporting
- ordinary fixed-effect reporting

Model outline:

- **Dataset:** `MASS::birthwt`
- **Outcome:** `low`
- **Model family:** `bernoulli(logit)`
- **Model:** `low ~ age_z + lwt_z + race + smoke + ptl_z + ht + ui + ftv_z`

Run in Terminal from the project root:

```bash
# Bash command block
cp examples/birthwt_logistic/00_config_birthwt_logistic.R 00_config.R
cp examples/birthwt_logistic/00_variable_dictionary_birthwt_logistic.csv 00_variable_dictionary.csv

Rscript examples/birthwt_logistic/00_create_birthwt_logistic_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_birthwt_logistic_stdout.log
```

`run_all.R` already rendered `results/publication/report/bayesian_mi_report_template.html`/`.docx`; no separate render step is needed.

### Example 3: Spline and monotonic effects using `MASS::birthwt`

This example tests custom `brms` formulae with `s()` and `mo()` terms.

It tests:

- `custom_formula`
- `s()` smooth terms
- `mo()` monotonic effects
- ordered categorical predictors
- special `brms` parameter summaries
- conditional-effect plots
- odds-ratio reporting for ordinary fixed effects

Model outline:

- **Dataset:** `MASS::birthwt`
- **Outcome:** `low`
- **Model family:** `bernoulli(logit)`
- **Model:** `low ~ s(age_z, k = 5) + mo(lwt_q) + race + smoke + ptl_z + ht + ui + ftv_z`

Here, `lwt_q` is an ordered quintile version of maternal weight, created by the example data script.

Run in Terminal from the project root:

```bash
# Bash command block
cp examples/birthwt_spline_monotonic/00_config_birthwt_spline_monotonic.R 00_config.R
cp examples/birthwt_spline_monotonic/00_variable_dictionary_birthwt_spline_monotonic.csv 00_variable_dictionary.csv

Rscript examples/birthwt_spline_monotonic/00_create_birthwt_spline_monotonic_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_birthwt_spline_monotonic_stdout.log
```

`run_all.R` already rendered `results/publication/report/bayesian_mi_report_template.html`/`.docx`; no separate render step is needed.

For this example, the model includes mo(lwt_q). To create derived monotonic-effect odds-ratio summaries, also run:

```bash
# Bash command block
Rscript 09_check_mo_parameter_columns.R
Rscript 10_publication_mo_results.R
quarto render results/publication/mo_effects/report/mo_effects_report.qmd
```

The main monotonic-effect table is `results/publication/mo_effects/tables/mo_cumulative_or_table.csv`.

Supplementary monotonic-effect outputs include:

- `results/publication/mo_effects/tables/mo_adjacent_or_table.csv`
- `results/publication/mo_effects/tables/mo_simplex_table.csv`
- `results/publication/mo_effects/report/mo_effects_report.html`
- `results/publication/mo_effects/report/mo_effects_report.docx`

See also `examples/birthwt_spline_monotonic/README_birthwt_spline_monotonic.md`.

### Cleaning outputs before switching examples

When switching from one example to another, clean the previous outputs first. Run in Terminal:

```bash
# Bash command block
rm -rf objects fits results
rm -f pipeline_error.flag
rm -f pipeline_success.flag
rm -f pipeline_progress.log
rm -f pipeline_heartbeat.txt
rm -f pipeline_stdout.log
rm -f run_all_stdout.log
```

If the repository includes a cleaning script, you can instead run:

```bash
# Bash command block
bash 99_cleanall.sh
```

### Automated example tests

The repository may include optional bash scripts to test the examples automatically.

Quick tests use small settings such as:

- `m = 5`
- `chains = 1`
- `iter = 500`
- `warmup = 250`
- `fit_workers = 1`
- `cores_per_fit = 1`
- `summary_workers = 1`
- `prediction_workers = 1`

Run in Terminal:

```bash
# Bash command block
bash test/test_all_examples_quick.sh
```

Parallel tests use modest parallel settings such as:

- `m = 10`
- `chains = 4`
- `iter = 500`
- `warmup = 250`
- `fit_workers = 2`
- `cores_per_fit = 4`
- `summary_workers = 2`
- `prediction_workers = 2`

Run in Terminal:

```bash
# Bash command block
bash test/test_all_examples_parallel.sh
```

A successful example test should create:

- `results/diagnostics.rds`
- `results/parameter_summary.rds`
- `results/publication/tables/main_effect_table_display.csv`
- `results/publication/tables/analysis_metadata.csv`
- `results/publication/report/bayesian_mi_report_template.qmd`
- `results/publication/report/bayesian_mi_report_template.html`
- `results/publication/report/bayesian_mi_report_template.docx`

---

### Test outputs

The bash test scripts write isolated preserved runs under `test/runs/`.

This means example tests do not delete or overwrite root-level `objects/`, `fits/` or `results/`.

Run quick tests:

```bash
# Bash command block
bash test/test_all_examples_quick.sh
```

Run parallel tests:

```bash
# Bash command block
bash test/test_all_examples_parallel.sh
```

List preserved test runs:

```bash
# Bash command block
bash test/list_test_runs.sh
```

To clean all preserved test runs manually:

```bash
# Bash command block
rm -rf test/runs
```

Add `test/runs/` to `.gitignore` if it is not already present.

The all-example test scripts cover:

- `airquality_gaussian`
- `birthwt_logistic`
- `birthwt_spline_monotonic`

The `birthwt_spline_monotonic` example exercises custom `brms` formula support for `s()` and `mo()` terms.

---

## 11. Computing environment setup

Before running the pipeline, prepare the R environment and CmdStan toolchain.

The exact setup depends on the operating system and computing environment. The sections below cover macOS, Windows and Linux. If you are using a managed workstation, shared server or high-performance computing cluster, some system tools may need to be installed by an administrator or loaded through environment modules.

---

### 11.1 macOS setup

Before running the pipeline, prepare the R environment and CmdStan toolchain. These instructions are for macOS. If running Windows, set up the Windows toolchain and R environment appropriately.

### 1. Install system tools on macOS

Run in Terminal to install the Apple command line tools:

```bash
# Bash command block
xcode-select --install
```

Run in Terminal to confirm that `make` and a C++ compiler are available:

```bash
# Bash command block
make --version
clang++ --version
```

If these commands fail, restart the Terminal and try again.

### 2. Install required R packages

Run in R or RStudio to install the packages required for this pipeline:

```r
install.packages(c(
  "tidyverse",
  "miceRanger",
  "brms",
  "posterior",
  "bayestestR",
  "future",
  "furrr",
  "doParallel",
  "foreach",
  "gt",
  "flextable",
  "officer",
  "forcats",
  "glue",
  "readr",
  "tibble",
  "dplyr",
  "stringr",
  "purrr",
  "rlang"
))
```

Run in R or RStudio to install `cmdstanr` from the Stan R-universe repository:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
```

### 3. Install CmdStan

Run in R or RStudio:

```r
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

cmdstanr::install_cmdstan()
```

This may take several minutes because CmdStan is compiled locally.

Run in R or RStudio to check the CmdStan path:

```r
cmdstanr::cmdstan_path()
```

You should see something like:

```text
/Users/yourname/.cmdstan/cmdstan-2.xx.x
```

### 4. Verify CmdStan works

Run this small CmdStanR test in R or RStudio:

```r
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

stan_file <- file.path(
  cmdstanr::cmdstan_path(),
  "examples",
  "bernoulli",
  "bernoulli.stan"
)

mod <- cmdstanr::cmdstan_model(stan_file)

fit <- mod$sample(
  data = list(
    N = 10,
    y = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0)
  ),
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 10,
  iter_sampling = 10,
  refresh = 1
)

fit$summary()
```

If this runs successfully, CmdStan is ready.

### 5. Verify brms + cmdstanr works

Run this minimal `brms` test in R or RStudio:

```r
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")

fit_test <- brms::brm(
  mpg ~ wt,
  data = mtcars,
  family = gaussian(),
  chains = 1,
  iter = 500,
  warmup = 250,
  cores = 1,
  backend = "cmdstanr",
  refresh = 10,
  silent = 0
)

summary(fit_test)
```

If this succeeds, the local Bayesian modelling environment is ready for the pipeline.

### 6. Install Quarto for report rendering

The pipeline can generate a Quarto report template in `results/publication/report/bayesian_mi_report_template.qmd`.

To render this report to HTML or DOCX, install Quarto.

Run in Terminal on macOS with Homebrew:

```bash
# Bash command block
brew install --cask quarto
```

Alternatively, download and install Quarto from the Quarto website.

Run in Terminal to verify the installation:

```bash
# Bash command block
quarto --version
```

If this command prints a version number, Quarto is ready.

### 7. Optional: clear CmdStanR cache if fitting fails unexpectedly

If a simple model fails with an error such as:

```text
Fitting failed. Unable to retrieve the metadata.
No chains finished successfully. Unable to retrieve the fit.
```

and the data look valid, run the following in Terminal to clear the CmdStanR cache:

```bash
# Bash command block
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

Then re-run the validation or pipeline.

### 8. Recommended reproducibility option: renv

For long-term reproducibility, consider using `renv`. Run in R or RStudio:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

This creates a project-specific package lockfile so the same package versions can be restored later. Run in R or RStudio:

```r
renv::restore()
```

---

### 11.2 Windows setup

Install:

- R
- RStudio
- Rtools
- Quarto

Use the version of Rtools that matches your R version. After installing Rtools, open R or RStudio and run:

```r
Sys.which("make")
Sys.which("g++")
```

Both should return valid paths.

Install the required R packages and `cmdstanr` using the same R commands shown in the macOS section. Then run:

```r
library(cmdstanr)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_path()
```

If toolchain problems persist, restart RStudio after installing Rtools. Then run the same CmdStan and `brms` verification examples shown in the macOS section.

To verify Quarto from Command Prompt or PowerShell:

```bash
# Command Prompt / PowerShell
quarto --version
```

---

### 11.3 Linux setup

Linux distributions differ, but the required system tools are generally:

- `make`
- `g++`
- `tar`
- `gzip`

For Ubuntu/Debian:

```bash
# Bash command block
sudo apt update
sudo apt install -y build-essential gfortran make
```

For Fedora:

```bash
# Bash command block
sudo dnf groupinstall -y "Development Tools"
sudo dnf install -y gcc-c++ gcc-gfortran make
```

For Red Hat Enterprise Linux or compatible systems, the exact commands may depend on the system configuration and repositories.

On a shared server or cluster, system tools may be provided through environment modules, for example:

```bash
# Bash command block
module avail gcc
module load gcc
```

Check with your system administrator or cluster documentation.

Install the required R packages and `cmdstanr` using the same R commands shown in the macOS section. Then run:

```r
library(cmdstanr)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_path()
```

Verify `brms` + `cmdstanr` using the same test shown in the macOS section.

Install Quarto using the installer appropriate for your distribution. On Ubuntu/Debian, for example:

```bash
# Bash command block
sudo dpkg -i quarto-*-linux-amd64.deb
```

Then check:

```bash
# Bash command block
quarto --version
```

On a shared server, Quarto may be available as a module or may need to be installed in a user directory.

---

### Dependencies

The pipeline requires both R packages and the Quarto command-line tool for report rendering.

Core R packages:

```r
install.packages(c(
  "tidyverse",
  "miceRanger",
  "brms",
  "posterior",
  "bayestestR",
  "future",
  "furrr",
  "doParallel",
  "foreach",
  "gt",
  "flextable",
  "officer",
  "forcats",
  "glue"
))
```

Quarto is also required if you want to render the generated `.qmd` report. Run in Terminal:

```bash
# Bash command block
# macOS with Homebrew
brew install --cask quarto

# Check installation
quarto --version
```

Install `cmdstanr` and CmdStan. Run in R or RStudio:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

cmdstanr::install_cmdstan()
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
```

If CmdStan is already installed, confirm the path in R or RStudio:

```r
cmdstanr::cmdstan_path()
```

---