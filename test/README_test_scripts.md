# Test scripts

The test scripts live in the `test/` folder to keep the project root clean.

## Important behaviour

Each test runs in an isolated copy of the project under:

```text
test/runs/
```

This means:

```text
root-level objects/, fits/, and results/ are not deleted
root-level outputs from previous analyses are not overwritten
each example test result is preserved
each quick/parallel run can be inspected later
```

Example run folder:

```text
test/runs/20260526_101530_birthwt_spline_monotonic_quick/project/
```

Inside each run's `project/` folder, you will find:

```text
objects/
fits/
results/
pipeline_progress.log
```

The run-level stdout logs are kept outside the runtime project root:

```text
test/runs/<run_id>/logs/run_all_<example>_<mode>_stdout.log
```

The parent run folder also contains:

```text
RUN_INFO.txt
```

## Examples covered

The test suite now covers all bundled examples:

```text
airquality_gaussian
  gaussian(identity), simple public-data example

birthwt_logistic
  bernoulli(logit), logistic model example

birthwt_spline_monotonic
  bernoulli(logit), custom brms formula with s() and mo()
```

The `birthwt_spline_monotonic` tests specifically check that the pipeline can:

```text
validate custom brms formulae
fit a model containing s()
fit a model containing mo()
extract special brms parameters such as simo_ and smooth terms
create publication outputs for special parameters
render the Quarto report
```

## Quick tests

Quick tests are intended to check that the pipeline runs end-to-end.

They use conservative settings:

```text
m = 5
chains = 1
iter = 500
warmup = 250
impute_workers = 1
num_impute_threads_per_worker = 1
fit_workers = 1
cores_per_fit = 1
```

Run all quick tests in Terminal from the project root:

```text
bash test/test_all_examples_quick.sh
```

Or run examples individually:

```text
bash test/test_airquality_quick.sh
bash test/test_birthwt_logistic_quick.sh
bash test/test_birthwt_spline_monotonic_quick.sh
```

## Parallel tests

Parallel tests are intended to check both:

```text
parallel miceRanger imputation
parallel brms fitting
```

They use modest settings:

```text
m = 10
chains = 4
iter = 500
warmup = 250
impute_workers = 2
num_impute_threads_per_worker = 2
fit_workers = 2
cores_per_fit = 2
```

Run all parallel tests in Terminal from the project root:

```text
bash test/test_all_examples_parallel.sh
```

Or run examples individually:

```text
bash test/test_airquality_parallel.sh
bash test/test_birthwt_logistic_parallel.sh
bash test/test_birthwt_spline_monotonic_parallel.sh
```

## List previous test runs

Run:

```text
bash test/list_test_runs.sh
```

## Clean test runs manually

The scripts do not delete previous test runs.

To remove all preserved test runs manually:

```text
rm -rf test/runs
```

To remove one run:

```text
rm -rf test/runs/<run_folder_name>
```

## Successful outputs

A successful test should create these inside the run's `project/` folder:

```text
results/diagnostics.rds
results/parameter_summary.rds
results/publication/tables/main_effect_table_display.csv
results/publication/tables/analysis_metadata.csv
results/publication/report/bayesian_mi_report_template.qmd
results/publication/report/bayesian_mi_report_template.html
results/publication/report/bayesian_mi_report_template.docx
```

For the `birthwt_spline_monotonic` example, the tests also expect:

```text
results/special_parameter_summary.rds
```

or:

```text
results/special_parameter_summary.csv
```

## Notes

The scripts do not require `pipeline_success.flag`, because some versions of the pipeline print `Pipeline completed successfully` without writing that flag. Instead, the scripts check for key output files and fail if `pipeline_error.flag` exists.

## Stdout logs

The `run_all` stdout log is written to the run-level `logs/` folder, not to the runtime project root:

```text
test/runs/<run_id>/logs/
```

This keeps each runtime project root cleaner while preserving the full command output for debugging.
