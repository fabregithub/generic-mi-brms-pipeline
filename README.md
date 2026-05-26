# Generic MI + brms Pipeline Template

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
- publication-ready tables, figures, methods/settings metadata, and report templates

The default example uses the built-in public dataset `datasets::airquality`, so the template can be tested and demonstrated without private data.

---

## Initial R environment setup

Before running the pipeline, prepare the R environment and CmdStan toolchain. These instructions are for macOS. If running Windows, set up the Windows toolchain and R environment appropriately.

### 1. Install system tools on macOS

Run in Terminal to install the Apple command line tools:

```
xcode-select --install
```

Run in Terminal to confirm that `make` and a C++ compiler are available:

```
make --version
clang++ --version
```

If these commands fail, restart the Terminal and try again.

### 2. Install required R packages

Run in R or RStudio to install the packages required for this pipeline:

```
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

```
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
```

### 3. Install CmdStan

Run in R or RStudio:

```
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

cmdstanr::install_cmdstan()
```

This may take several minutes because CmdStan is compiled locally.

Run in R or RStudio to check the CmdStan path:

```
cmdstanr::cmdstan_path()
```

You should see something like:

```
/Users/yourname/.cmdstan/cmdstan-2.xx.x
```

### 4. Verify CmdStan works

Run this small CmdStanR test in R or RStudio:

```
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

```
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

The pipeline can generate a Quarto report template in:

```
results/publication/report/bayesian_mi_report_template.qmd
```

To render this report to HTML or DOCX, install Quarto.

Run in Terminal on macOS with Homebrew:

```
brew install --cask quarto
```

Alternatively, download and install Quarto from the Quarto website.

Run in Terminal to verify the installation:

```
quarto --version
```

If this command prints a version number, Quarto is ready.

### 7. Optional: clear CmdStanR cache if fitting fails unexpectedly

If a simple model fails with an error such as:

```
Fitting failed. Unable to retrieve the metadata.
No chains finished successfully. Unable to retrieve the fit.
```

and the data look valid, run the following in Terminal to clear the CmdStanR cache:

```
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

Then re-run the validation or pipeline.

### 8. Recommended reproducibility option: renv

For long-term reproducibility, consider using `renv`. Run in R or RStudio:

```
install.packages("renv")
renv::init()
renv::snapshot()
```

This creates a project-specific package lockfile so the same package versions can be restored later. Run in R or RStudio:

```
renv::restore()
```


## Quick start

From a fresh template folder, run in Terminal:

```
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv
Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

Outputs are written to:

```
objects/
fits/
results/
results/publication/
```

---

## Main files to edit for a new project

Usually edit only:

```
00_config.R
00_variable_dictionary.csv
```

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

---

---

## Repeated-outcome data and subject-wide imputation

The generic pipeline can also handle repeated-outcome data where the outcome is measured multiple times per subject and most covariates are subject-level or baseline variables.

Use this imputation strategy in `00_config.R`:

```
analysis_spec$imputation$strategy <- "subject_wide_with_repeated_y_auxiliary"
```

This strategy does the following:

```
1. Starts from long data, one row per subject-time observation.
2. Creates a subject-wide imputation dataset, one row per subject.
3. Converts the repeated outcome into wide auxiliary variables, for example:
   ps_1, ps_2, ps_3, ps_4, ps_5, ps_6
4. Imputes subject-level covariates once per subject.
5. Joins the imputed subject-level variables back to the original long outcome data.
6. Fits the brms model on the long data.
```

Required `00_config.R` fields:

```
analysis_spec$data$id_var
analysis_spec$data$time_var
analysis_spec$outcome$y_var
analysis_spec$outcome$y_prefix
analysis_spec$outcome$y_wide_regex
```

Example:

```
data = list(
  raw_data_file = "data/my_repeated_data.rds",
  id_var = "NO",
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

```
timing = single
```

or:

```
timing = baseline
```

Repeated or time-varying variables should have:

```
timing = repeated
```

or:

```
timing = time_varying
```

If `impute_y = FALSE`, repeated outcome values are not imputed as targets. They are used as wide auxiliary predictors when possible, and missing outcome rows can later be summarised using posterior prediction.

### Ordinal predictors and `mo()`

Ordinal variables should be marked in the dictionary:

```
type = ordinal
```

The pipeline converts these variables to ordered factors before modelling. They can then be used in a custom `brms` formula:

```
custom_formula = brms::bf(
  y ~ mo(education) + mo(income) + age_z + sex + (1 | id)
)
```

For large mixed logistic models, `mo()` can be much slower than ordinary factor coding. A practical approach is to use the factor-coded model as the main analysis and use `mo()` as a sensitivity analysis with fewer imputations.



## Variable roles: dictionary by default, config as optional override

The pipeline uses `00_variable_dictionary.csv` as the default source of truth for variable metadata.

This means the following information should normally be specified only in the dictionary:

```
role
type
timing
scale
reference
impute_target
use_in_model
use_as_auxiliary
```

Therefore, the standard setting in `00_config.R` is:

```
variables = NULL
```

The pipeline automatically derives internal variable groups from the dictionary, including:

```
exposure_vars
covariate_vars
auxiliary_vars
continuous_vars
categorical_vars
ordinal_vars
subject_level_vars
time_varying_vars
scale_vars
```

For example, if `00_variable_dictionary.csv` contains:

```
var,role,type,timing,scale,use_in_model
Solar.R,covariate,continuous,single,z,TRUE
Wind,covariate,continuous,single,z,TRUE
Temp,covariate,continuous,single,z,TRUE
Month,covariate,categorical,single,no,TRUE
```

and `analysis_spec$model$fixed_effects` is:

```
fixed_effects = "auto"
```

the model uses:

```
Solar.R_z + Wind_z + Temp_z + Month
```

### Optional overrides

Advanced users can override selected derived groups in `00_config.R`.

Example:

```
variables = list(
  scale_vars = c("age", "income"),
  auxiliary_vars = c("baseline_score")
)
```

If `variables = NULL`, no override is applied.

If `exposure_vars` or `covariate_vars` are explicitly supplied in `analysis_spec$variables`, they are used by `fixed_effects = "auto"` instead of the dictionary-derived `use_in_model` list. For most analyses, it is simpler and safer to leave `variables = NULL` and control the model through the dictionary.



## Variable dictionary guide

The file `00_variable_dictionary.csv` is the main machine-readable description of the analysis variables.

Expected columns are:

```
var
label
role
type
timing
scale
reference
impute_target
use_in_model
use_as_auxiliary
```

Some fields are straightforward:

- `var`: exact variable name in the dataset
- `label`: human-readable label for tables and reports
- `impute_target`: whether this variable should be imputed when missing
- `use_in_model`: whether this variable should appear in the final brms model

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
| `id` | Subject, cluster, or row identifier | `ID`, `row_id` |
| `time` | Measurement occasion, wave, visit, or follow-up time | `time`, `wave` |
| `cluster` | Grouping variable for random effects or clustering | school ID, hospital ID |
| `strata` | Stratification variable, if relevant | site, cohort |

For most standard regression analyses, the most common roles are:

```
outcome
binary_outcome
exposure
covariate
id
time
auxiliary
```

Example:

```
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

```
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

```
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
| `z` | Standardize to mean 0 and SD 1 | creates `var_z` |
| `center` | Mean-center only | creates centered version if supported |
| `log` | Log-transform | creates log version if supported |
| `custom` | User-defined transform outside dictionary | handled in config/functions |

Currently, the most important supported value is:

```
z
```

Example:

```
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
```

This creates:

```
age_z
```

and the model should use `age_z`, not `age`.

Use `z` for continuous predictors when coefficients should represent the effect per 1 SD increase. This is usually helpful for Bayesian models because weakly informative priors like `normal(0, 1)` or `normal(0, 1.5)` are easier to interpret on standardised predictors.

Do not use `z` for:

```
outcomes
ID variables
binary/categorical variables
already standardised variables
```

### `reference`

`reference` defines the reference category for binary, categorical, or ordinal variables.

Examples:

```
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
```

For a binary variable coded `0/1`, setting:

```
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

```
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

```
baseline_score,Baseline questionnaire score,auxiliary,continuous,baseline,z,,FALSE,FALSE,TRUE
```

2. Variable to impute and include:

```
income,Household income,covariate,categorical,baseline,no,1,TRUE,TRUE,FALSE
```

3. Outcome not imputed:

```
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
```

4. Repeated outcome as imputation auxiliary:

```
ps,Psychological stress,binary_outcome,binary,repeated,no,0,FALSE,TRUE,TRUE
```

Whether repeated outcomes are used as auxiliary predictors depends on the imputation strategy in `00_config.R`.

### Example: Gaussian demo

```
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
Ozone,Ozone concentration,outcome,continuous,single,no,,FALSE,TRUE,FALSE
Solar.R,Solar radiation,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Wind,Wind speed,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Temp,Temperature,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Month,Month,covariate,categorical,single,no,5,FALSE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

### Example: Logistic demo

```
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


---

---

## Example datasets

The repository includes three public example analyses. These are intended to help users test the full pipeline before applying it to private study data.

The example folders are:

```
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

```
a config file
a variable dictionary
a data-creation script
```

To use an example, copy its config file and variable dictionary to the project root as:

```
00_config.R
00_variable_dictionary.csv
```

Then run the example data-creation script.

### Example 1: Gaussian model using `datasets::airquality`

This is the default Gaussian example.

It tests:

```
continuous outcome
gaussian(identity) model
row-level imputation
ordinary fixed-effect reporting
posterior prediction for rows with missing outcome values
```

Model outline:

```
Dataset: datasets::airquality
Outcome: Ozone
Model family: gaussian(identity)
Model: Ozone ~ Solar.R_z + Wind_z + Temp_z + Month
```

Run in Terminal from the project root:

```
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv

Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

Render the report:

```
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

### Example 2: Logistic model using `MASS::birthwt`

This example tests the Bernoulli/logit workflow with a binary outcome.

It tests:

```
binary outcome
bernoulli(logit) model
row-level imputation
odds-ratio reporting
ordinary fixed-effect reporting
```

Model outline:

```
Dataset: MASS::birthwt
Outcome: low
Model family: bernoulli(logit)
Model: low ~ age_z + lwt_z + race + smoke + ptl_z + ht + ui + ftv_z
```

Run in Terminal from the project root:

```
cp examples/birthwt_logistic/00_config_birthwt_logistic.R 00_config.R
cp examples/birthwt_logistic/00_variable_dictionary_birthwt_logistic.csv 00_variable_dictionary.csv

Rscript examples/birthwt_logistic/00_create_birthwt_logistic_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_birthwt_logistic_stdout.log
```

Render the report:

```
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

### Example 3: Spline and monotonic effects using `MASS::birthwt`

This example tests custom `brms` formulae with `s()` and `mo()` terms.

It tests:

```
custom_formula
s() smooth terms
mo() monotonic effects
ordered categorical predictors
special brms parameter summaries
conditional-effect plots
odds-ratio reporting for ordinary fixed effects
```

Model outline:

```
Dataset: MASS::birthwt
Outcome: low
Model family: bernoulli(logit)
Model: low ~ s(age_z, k = 5) + mo(lwt_q) + race + smoke + ptl_z + ht + ui + ftv_z
```

Here, `lwt_q` is an ordered quintile version of maternal weight, created by the example data script.

Run in Terminal from the project root:

```
cp examples/birthwt_spline_monotonic/00_config_birthwt_spline_monotonic.R 00_config.R
cp examples/birthwt_spline_monotonic/00_variable_dictionary_birthwt_spline_monotonic.csv 00_variable_dictionary.csv

Rscript examples/birthwt_spline_monotonic/00_create_birthwt_spline_monotonic_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_birthwt_spline_monotonic_stdout.log
```

Render the report:

```
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

See also:

```
examples/birthwt_spline_monotonic/README_birthwt_spline_monotonic.md
```

### Cleaning outputs before switching examples

When switching from one example to another, clean the previous outputs first. Run in Terminal:

```
rm -rf objects fits results
rm -f pipeline_error.flag
rm -f pipeline_success.flag
rm -f pipeline_progress.log
rm -f pipeline_heartbeat.txt
rm -f pipeline_stdout.log
rm -f run_all_stdout.log
```

If the repository includes a cleaning script, you can instead run:

```
bash 99_cleanall.sh
```

### Automated example tests

The repository may include optional bash scripts to test the examples automatically.

Quick tests use small settings such as:

```
m = 5
chains = 1
iter = 500
warmup = 250
fit_workers = 1
cores_per_fit = 1
```

Run in Terminal:

```
bash test/test_all_examples_quick.sh
```

Parallel tests use modest parallel settings such as:

```
m = 10
chains = 4
iter = 500
warmup = 250
fit_workers = 2
cores_per_fit = 4
```

Run in Terminal:

```
bash test/test_all_examples_parallel.sh
```

A successful example test should create:

```
results/diagnostics.rds
results/parameter_summary.rds
results/publication/tables/main_effect_table_display.csv
results/publication/tables/analysis_metadata.csv
results/publication/report/bayesian_mi_report_template.qmd
results/publication/report/bayesian_mi_report_template.html
results/publication/report/bayesian_mi_report_template.docx
```


---

## Parallel miceRanger imputation

The pipeline can parallelise `miceRanger` imputation using `doParallel` and `foreach`.

In `00_config.R`, the relevant settings are:

```
analysis_spec$parallel$impute_workers <- 4
analysis_spec$parallel$num_impute_threads_per_worker <- 4
```

The approximate CPU demand during imputation is:

```
impute_workers * num_impute_threads_per_worker
```

For example:

```
4 workers * 4 threads per worker = about 16 active threads
```

Start conservatively, especially on laptops or when the imputation data are large, because parallel `miceRanger` can copy data to worker processes and increase memory use.

Recommended starting values:

```
# Public examples or ordinary laptops
analysis_spec$parallel$impute_workers <- 1
analysis_spec$parallel$num_impute_threads_per_worker <- 1

# High-memory workstation
analysis_spec$parallel$impute_workers <- 4
analysis_spec$parallel$num_impute_threads_per_worker <- 4
```

To test parallel imputation from scratch, remove old imputation outputs first:

```
rm -f objects/imputation_manifest.rds
rm -f objects/imputation_spec.rds
rm -rf objects/imputed_data
rm -rf objects/imputed_wide
rm -rf objects/model_data
rm -rf fits results
rm -f pipeline_error.flag pipeline_success.flag
```

Then re-run the pipeline.


## Recommended high-performance settings

The following settings were tested successfully on a high-memory machine. Edit these in `00_config.R` using RStudio or another text editor:

```
analysis_spec$imputation$m <- 100

analysis_spec$model$chains <- 4
analysis_spec$model$iter <- 2000
analysis_spec$model$warmup <- 1000
analysis_spec$model$run_smoke_fit <- TRUE

analysis_spec$parallel$fit_workers <- 4
analysis_spec$parallel$cores_per_fit <- 4
analysis_spec$parallel$future_globals_maxsize_gb <- 80
```

This runs:

```
100 imputations
4 chains per model
4 imputed datasets fitted in parallel
16 active chains total
```

For a new analysis, keep the following setting in `00_config.R`:

```
analysis_spec$model$run_smoke_fit <- TRUE
```

The smoke fit runs one sequential model first, before launching parallel workers. This catches formula, prior, data, or CmdStan problems early.

After a configuration has been tested successfully, you may set the following in `00_config.R`:

```
analysis_spec$model$run_smoke_fit <- FALSE
```

to save a little time.

---

## Important design notes

This template does **not** use `brm_multiple()` for model fitting.

Instead, it fits one `brms` model per imputed dataset and saves each fit immediately:

```
fits/fit_imp_001.rds
fits/fit_imp_002.rds
...
fits/fit_imp_100.rds
```

This is intentionally safer for large datasets because:

- the main R session does not hold all fitted models in memory
- completed fits are preserved if the run stops
- failed or slow imputations can be re-run separately
- valid existing fits are skipped on re-run
- worker processes return only small status objects to the main session

Parallelism happens across imputations using `future` / `furrr`, with dynamic scheduling. This is configured inside the R scripts:

```
furrr::furrr_options(
  seed = TRUE,
  scheduling = Inf
)
```

This improves load balancing when some imputed datasets take longer than others.

---

## Common commands

Validate config. Run in Terminal:

```
Rscript 01_validate_config.R
```

Run all steps. Run in Terminal:

```
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Fit one imputed dataset only. Run in Terminal:

```
Rscript fit_single_imputation.R 1
```

For example, to fit the model only for imputed dataset 51, run in Terminal:

```
Rscript fit_single_imputation.R 51
```

Render the generated Quarto report after installing Quarto. Run in Terminal:

```
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

---

## Restarting after interruption

The pipeline is checkpointed. If a run is interrupted, run the following in Terminal:

```
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Existing valid fit files are skipped.

To clean only fitting and downstream results while keeping prepared/imputed data, run in Terminal:

```
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
or run ```bash 99_clean_fitting_results.sh```.

To clean imputation and all downstream outputs, run in Terminal:

```
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
or run ```bash 99_cleanall.sh```.

---

## Logging and monitoring

The pipeline writes logs and status files:

```
pipeline_progress.log
pipeline_stdout.log
pipeline_heartbeat.txt
pipeline_success.flag
pipeline_error.flag
results/worker_logs/
```

Monitor progress by running the following in Terminal:

```
tail -f pipeline_progress.log
```

Inspect worker-level fitting logs by running the following in Terminal:

```
ls -lh results/worker_logs
cat results/worker_logs/fit_worker_imp_001.log
```

If the pipeline completes successfully, this file is created:

```
pipeline_success.flag
```

If an R-level error is caught, this file is created:

```
pipeline_error.flag
```

If the machine crashes or restarts unexpectedly, there may be no error flag. In that case, check the heartbeat and logs, then re-run. The checkpoint system should skip completed valid fits.

---

## Troubleshooting CmdStan cache issues

If a simple model fails with an error like:

```
Fitting failed. Unable to retrieve the metadata.
No chains finished successfully. Unable to retrieve the fit.
```

and the data look valid, try clearing the CmdStanR cache. Run in Terminal:

```
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

Then re-run the direct debug fit or the pipeline.

This can resolve stale or corrupted compiled-model cache issues.

---

## Debugging model fitting

A good sequence is:

1. Validate the config. Run in Terminal:

```
Rscript 01_validate_config.R
```

2. Run one direct or smoke fit.

3. If successful, run the full pipeline.

4. If a specific imputed dataset is slow or problematic during model fitting, fit it alone. Run in Terminal:

```
Rscript fit_single_imputation.R 51
```

You can temporarily skip or restrict **model fitting for selected imputed datasets** in `00_config.R`.

These options do **not** skip the imputation step itself. They only control which already-created imputed datasets are passed to Step 4, where `brms` models are fitted.

To fit all imputations except imputations 45 and 51:

```
analysis_spec$model$skip_imputations <- c(45, 51)
analysis_spec$model$only_imputations <- integer(0)
```

To fit only imputation 51, for example as a diagnostic re-run:

```
analysis_spec$model$skip_imputations <- integer(0)
analysis_spec$model$only_imputations <- c(51)
```

Interpretation:

```
skip_imputations = do not fit brms models for these imputed datasets
only_imputations = fit brms models only for these imputed datasets
```

These settings are useful when one imputed dataset is unusually slow or problematic. Completed valid fit files are still preserved and skipped on re-run.

---

## Dependencies

The pipeline requires both R packages and the Quarto command-line tool for report rendering.

Core R packages:

```
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

```
# macOS with Homebrew
brew install --cask quarto

# Check installation
quarto --version
```

Install `cmdstanr` and CmdStan. Run in R or RStudio:

```
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

cmdstanr::install_cmdstan()
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
```

If CmdStan is already installed, confirm the path in R or RStudio:

```
cmdstanr::cmdstan_path()
```

---

## Supported analysis patterns

The template is designed to support:

```
single-time outcome, row-level covariates
repeated outcome with subject-level covariates
repeated outcome with time-varying covariates
complete-case analysis without imputation
row-level multiple imputation
subject-level multiple imputation
subject-wide imputation using repeated Y as auxiliary variables
```

Supported model families include:

```
gaussian
bernoulli
poisson
negbinomial
beta
ordinal
categorical
```

Model families and links are set in `00_config.R`.

---

## Publication outputs

After successful completion, publication outputs are written to:

```
results/publication/
```

Typical outputs include:

```
results/publication/tables/main_effect_table_display.csv
results/publication/tables/main_effect_table_full.csv
results/publication/tables/diagnostics_summary.csv
results/publication/tables/analysis_metadata.csv
results/publication/tables/analysis_metadata.rds
results/publication/figures/forest_plot_odds_ratios.png
results/publication/report/bayesian_mi_report_template.qmd
```

The generated Quarto report includes posterior results, diagnostics, figures, and a methods/settings table based on `analysis_metadata.csv`. This table records key analysis settings such as the imputation strategy, target number of imputations, fitted imputations used in posterior summaries, model formula, family/link, priors, MCMC settings, parallel settings, posterior-summary settings, and predictive-draw settings.

Render the report with:

```
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

---

## Notes for adapting to private study data

For a new project:

1. Replace or edit `00_variable_dictionary.csv`.
2. Edit `00_config.R`.
3. Run in Terminal:

```
Rscript 01_validate_config.R
```

4. If validation passes, run in Terminal:

```
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Recommended approach for a new analysis:

### 1. Quick test run

For the first test of a new dataset or model, use small settings. Edit these in `00_config.R` using RStudio or another text editor:

```
analysis_spec$imputation$m <- 5

analysis_spec$model$chains <- 1
analysis_spec$model$iter <- 500
analysis_spec$model$warmup <- 250
analysis_spec$model$run_smoke_fit <- TRUE

analysis_spec$parallel$fit_workers <- 1
analysis_spec$parallel$cores_per_fit <- 1
analysis_spec$parallel$future_globals_maxsize_gb <- 8

analysis_spec$posterior_prediction$ndraws <- 200
```

This quick test is intended to check that:

```
the data are read correctly
the variable dictionary is valid
imputation runs
the brms formula is correct
priors are compatible with the model
one model can be fitted successfully
posterior summaries and publication outputs are created
```

### 2. Parallel test run

After the quick test succeeds, test parallel fitting with modest settings:

```
analysis_spec$imputation$m <- 10

analysis_spec$model$chains <- 4
analysis_spec$model$iter <- 500
analysis_spec$model$warmup <- 250
analysis_spec$model$run_smoke_fit <- TRUE

analysis_spec$parallel$fit_workers <- 2
analysis_spec$parallel$cores_per_fit <- 4
analysis_spec$parallel$future_globals_maxsize_gb <- 20
```

This checks that multiple chains and multiple imputed datasets can run safely on the machine.

### 3. Full production run

For the final analysis, increase the imputation and MCMC settings. For a high-memory machine, a typical setting is:

```
analysis_spec$imputation$m <- 100

analysis_spec$model$chains <- 4
analysis_spec$model$iter <- 2000
analysis_spec$model$warmup <- 1000
analysis_spec$model$run_smoke_fit <- TRUE

analysis_spec$parallel$fit_workers <- 4
analysis_spec$parallel$cores_per_fit <- 4
analysis_spec$parallel$future_globals_maxsize_gb <- 80

analysis_spec$posterior_prediction$ndraws <- 1000
```

For repeated runs of the same already-tested analysis, you may set:

```
analysis_spec$model$run_smoke_fit <- FALSE
```

but keep it as `TRUE` when changing the dataset, formula, priors, outcome family, or imputation strategy.

For large analyses, avoid `brm_multiple()` unless you have a specific reason to use it. The one-fit-per-imputation design is usually safer and easier to restart.
