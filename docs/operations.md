# Running, monitoring and recovering

Restarting an interrupted run, following progress, clearing a stale CmdStan cache, and debugging a model that will not fit.

> Reach for this when a run stops, stalls, or produces something you did not expect.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

This section covers what to do once a run is already underway or has already produced some output: restarting after an interruption, monitoring progress, recovering from CmdStan cache problems, and debugging a specific imputed dataset. For the commands used to start a run in the first place, see [Quick start](../README.md#quick-start).

---

## Restarting after interruption

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

## Logging and monitoring

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

## Troubleshooting CmdStan cache issues

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

## Debugging model fitting

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

---

*[← Back to the main README](../README.md)*
