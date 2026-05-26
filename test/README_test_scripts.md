# Test scripts

The test scripts live in the `test/` folder to keep the project root clean.

## Important behaviour

Each test now runs in an isolated copy of the project under:

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
test/runs/20260526_101530_airquality_quick/project/
```

Inside each run folder, you will find:

```text
objects/
fits/
results/
run_all_<example>_<mode>_stdout.log
pipeline_progress.log
```

The parent folder also contains:

```text
RUN_INFO.txt
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

Run in Terminal from the project root:

```text
bash test/test_all_examples_quick.sh
```

Or run examples individually:

```text
bash test/test_airquality_quick.sh
bash test/test_birthwt_logistic_quick.sh
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

Run in Terminal from the project root:

```text
bash test/test_all_examples_parallel.sh
```

Or run examples individually:

```text
bash test/test_airquality_parallel.sh
bash test/test_birthwt_logistic_parallel.sh
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

## What the scripts do

Each test script:

```text
creates an isolated runtime copy of the project under test/runs/
copies the example config to runtime 00_config.R
copies the example variable dictionary to runtime 00_variable_dictionary.csv
creates the public example data inside the runtime copy
appends temporary test overrides to runtime 00_config.R
runs validation
runs the full pipeline
renders the Quarto report
checks key outputs exist
leaves all outputs in the run folder
```

The temporary override block is clearly marked inside the runtime copy:

```text
# ---- BEGIN automated test overrides ----
...
# ---- END automated test overrides ----
```

The root project files are not edited.

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

## Notes

The scripts do not require `pipeline_success.flag`, because some versions of the pipeline print `Pipeline completed successfully` without writing that flag. Instead, the scripts check for key output files and fail if `pipeline_error.flag` exists.

## Implementation note

The helper function `log()` writes to stderr. This is intentional: some helper
functions return paths through stdout and are called by command substitution.
Keeping logs on stderr prevents paths from being polluted by timestamped log
messages.
