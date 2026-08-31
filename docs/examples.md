# Examples and tests

Four worked analyses on bundled public datasets, and the automated test scripts that run them end to end.

> The fastest way to see the pipeline work before pointing it at your own data. Each example is a complete, runnable configuration.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

The repository includes four public example analyses. These are intended to help users test the full pipeline before applying it to private study data.

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
├── birthwt_spline_monotonic
│   ├── 00_config_birthwt_spline_monotonic.R
│   ├── 00_create_birthwt_spline_monotonic_example_data.R
│   ├── 00_variable_dictionary_birthwt_spline_monotonic.csv
│   └── README_birthwt_spline_monotonic.md
└── lung_cox
    ├── 00_config_lung_cox.R
    ├── 00_create_lung_cox_example_data.R
    └── 00_variable_dictionary_lung_cox.csv
```

Each example contains:

- a config file
- a variable dictionary
- a data-creation script

To use an example, copy its config file and variable dictionary to the project root as `00_config.R` and `00_variable_dictionary.csv`.

Then run the example data-creation script.

## Example 1: Gaussian model using `datasets::airquality`

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

## Example 2: Logistic model using `MASS::birthwt`

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

## Example 3: Spline and monotonic effects using `MASS::birthwt`

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

## Example 4: Cox proportional hazards model using `survival::lung`

This example tests the Cox PH workflow for time-to-event data with right censoring and missing predictor values.

It tests:

- survival outcome (time to death, right-censored)
- `cox(log)` family with `time | cens(censored)` response
- `custom_formula` to specify the brms survival response
- hazard-ratio reporting
- multiple imputation of missing predictor values

Model outline:

- **Dataset:** `survival::lung` (Loprinzi et al., 1994)
- **Outcome:** days to death, right-censored (`time | cens(censored)`)
- **Model family:** `cox(log)`
- **Model:** `time | cens(censored) ~ age_z + sex + ph_ecog + ph_karno_z + wt_loss_z + meal_cal_z`

The censoring indicator is recoded so that `0` = event (death) and `1` = right-censored, which is the convention brms expects for `cens()`.

> **Package requirement:** Cox PH models in brms require the `survival` package at model-fitting time. `survival` ships with R as a recommended package and is normally already available; if not, install it with `install.packages("survival")`.

Run in Terminal from the project root:

```bash
cp examples/lung_cox/00_config_lung_cox.R 00_config.R
cp examples/lung_cox/00_variable_dictionary_lung_cox.csv 00_variable_dictionary.csv

Rscript examples/lung_cox/00_create_lung_cox_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_lung_cox_stdout.log
```

`run_all.R` already rendered `results/publication/report/bayesian_mi_report_template.html`/`.docx`; no separate render step is needed.

## Cleaning outputs before switching examples

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

## Automated example tests

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

## Test outputs

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

---

*[← Back to the main README](../README.md)*
