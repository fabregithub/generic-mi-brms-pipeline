# Test scripts

The test scripts live in the `test/` folder to keep the project root clean.

Project layout:

```
generic_mi_brms_pipeline
├── run_all.R
├── 00_config.R
├── examples
└── test
    ├── test_airquality_quick.sh
    ├── test_birthwt_logistic_quick.sh
    ├── test_all_examples_quick.sh
    ├── test_airquality_parallel.sh
    ├── test_birthwt_logistic_parallel.sh
    ├── test_all_examples_parallel.sh
    └── test_example_common.sh
```

The scripts automatically detect the project root as the parent of `test/`, so they can be run from the project root or from another working directory.

## Quick tests

Quick tests are intended to check that the pipeline runs end-to-end.

They use conservative settings:

```
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

```
bash test/test_all_examples_quick.sh
```

Or run examples individually:

```
bash test/test_airquality_quick.sh
bash test/test_birthwt_logistic_quick.sh
```

## Parallel tests

Parallel tests are intended to check both:

```
parallel miceRanger imputation
parallel brms fitting
```

They use modest settings:

```
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

```
bash test/test_all_examples_parallel.sh
```

Or run examples individually:

```
bash test/test_airquality_parallel.sh
bash test/test_birthwt_logistic_parallel.sh
```

The parallel scripts check the stdout log for evidence that the parallel imputation path was used, such as:

```
miceRanger parallel: TRUE
impute_workers: 2
```

If the log-check warning appears but the pipeline succeeds, inspect the stdout log manually. It may simply mean the log format has changed.

## What the scripts do

Each test script:

```
copies the example config to 00_config.R
copies the example variable dictionary to 00_variable_dictionary.csv
creates the public example data
cleans old outputs
appends temporary test overrides to 00_config.R
runs validation
runs the full pipeline
renders the Quarto report
checks key outputs exist
```

The temporary override block is clearly marked:

```
# ---- BEGIN automated test overrides ----
...
# ---- END automated test overrides ----
```

The next test run removes the previous override block before adding a new one.

## Successful outputs

A successful test should create:

```
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
