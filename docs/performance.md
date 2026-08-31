# Parallelisation and performance tuning

How the pipeline uses cores and memory, and how to size `impute_workers`, `fit_workers` and `chains` for the machine you actually have.

> Start here if a run is slower than expected, or if it exhausted memory. The defaults are deliberately conservative.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

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

If you use the automatic m-increment loop (`analysis_spec$mi_stability$auto_increment = TRUE`, see [Choosing the number of imputations adaptively](choosing-m.md)), `fit_workers` also determines the default batch size: each batch defaults to `fit_workers` imputations, so every batch fully occupies the parallel workers with none left idle. With `fit_workers = 4`, the default batches are `m = 4, 8, 12, ...`; with `fit_workers = 6`, they are `m = 6, 12, 18, ...`. A custom `increment_size` that is not already a multiple of `fit_workers` is rounded up automatically.

`num_impute_threads` is a backward-compatible fallback for older scripts. If all scripts in the repository use `impute_workers` and `num_impute_threads_per_worker`, it can eventually be removed.

## Suggested settings by computing environment

### Laptop or low-memory desktop

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

### Standard desktop or small workstation

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

### High-memory workstation

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

### Shared server or high-performance computing environment

Use conservative per-job settings unless your scheduler explicitly allocates more resources. Avoid requesting more threads than the scheduler has allocated to the job.

---

## Parallel miceRanger imputation

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

## Recommended high-performance settings

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

If you use the automatic m-increment loop (`analysis_spec$mi_stability$auto_increment = TRUE`), you do not need to manually toggle this setting between batches: the loop runs the smoke fit only once, for the first batch, and disables it automatically for every later batch in the same `run_all.R` invocation. See [Automatic m-increment loop](choosing-m.md#automatic-m-increment-loop-recommended).

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

---

*[← Back to the main README](../README.md)*
