# For the public airquality example, run:
# Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
# once before running the full pipeline.

run_step <- function(s) {
  cat("\n\n================ RUNNING", s, "================\n\n")
  source(s)
}

set_runtime_override <- function(paths, override) {
  saveRDS(override, file.path(paths$objects, "mi_runtime_override.rds"))
}

clear_runtime_override <- function(paths) {
  f <- file.path(paths$objects, "mi_runtime_override.rds")
  if (file.exists(f)) file.remove(f)
}

run_step("01_validate_config.R")

# Defensive: wipe any override left behind by a previous run that was
# interrupted mid-loop, so this run always starts from the static config.
clear_runtime_override(paths)

run_step("02_prepare_data.R")

# ------------------------------------------------------------
# Steps 3/4/6: imputation, fitting, posterior-draw extraction
#
# If analysis_spec$mi_stability$auto_increment is TRUE, m is increased in
# batches (sized by analysis_spec$mi_stability$increment_size, defaulting to
# fit_workers so each round exactly fills the parallel workers) rather than
# fitting analysis_spec$imputation$m models up front. After each new batch,
# a lightweight two-batch stability check compares the new cumulative
# results against the previous batch; the loop stops as soon as every
# selected parameter is stable by the configured thresholds, or once the
# configured m is reached, whichever comes first. This avoids paying for
# m = 100 model fits when m = 40 would already have given the same answer.
#
# If auto_increment is FALSE (default), all m imputations/fits/draws are
# produced in a single pass, exactly as before this loop was introduced.
# ------------------------------------------------------------

mi_loop_enabled <- isTRUE(analysis_spec$mi_stability$auto_increment)

if (mi_loop_enabled) {
  max_m <- as.integer(analysis_spec$imputation$m %||% 1)

  fit_workers <- as.integer(analysis_spec$parallel$fit_workers %||% 1)
  fit_workers <- max(1L, fit_workers)

  # Default increment_size to fit_workers (already an exact multiple). If
  # the user supplies a custom increment_size that isn't a multiple of
  # fit_workers, round it up to the nearest multiple -- purely for
  # efficiency, so every batch fully occupies the parallel workers with no
  # idle capacity on the final imputation of a batch.
  increment_size <- analysis_spec$mi_stability$increment_size %||% fit_workers
  increment_size <- max(1L, as.integer(increment_size))

  if (increment_size %% fit_workers != 0) {
    rounded_increment_size <- ceiling(increment_size / fit_workers) * fit_workers

    cat(
      "\nRounding mi_stability$increment_size from", increment_size,
      "up to", rounded_increment_size,
      "so each batch is a multiple of fit_workers =", fit_workers, "\n"
    )

    increment_size <- rounded_increment_size
  }

  cat(
    "\n\nAutomatic imputation-count stability loop enabled.\n",
    "Increment size:", increment_size, "| fit_workers:", fit_workers,
    "| Maximum m:", max_m, "\n\n"
  )

  m_current <- min(increment_size, max_m)

  # Resume support: a previous run may already have grown m past the first
  # batch (and then stopped, or failed in a later step). Step 3's overwrite
  # guard rightly refuses to shrink m below the existing imputation count,
  # so start the loop from that count instead of the first batch size.
  imputation_manifest_file <- file.path(paths$objects, "imputation_manifest.rds")
  n_existing_imputations <- 0L

  if (file.exists(imputation_manifest_file)) {
    existing_manifest <- readRDS(imputation_manifest_file)

    if (is.data.frame(existing_manifest) && "imputed_file" %in% names(existing_manifest)) {
      existing_files <- as.character(existing_manifest$imputed_file)
      n_existing_imputations <- sum(
        !is.na(existing_files) & nzchar(existing_files) & file.exists(existing_files)
      )
    }
  }

  resumed_past_first_batch <- n_existing_imputations > m_current

  if (resumed_past_first_batch) {
    m_current <- min(n_existing_imputations, max_m)

    cat(
      "\nFound", n_existing_imputations,
      "existing imputation(s); resuming the stability loop at m =", m_current, "\n"
    )
  }

  # The smoke fit's job is to catch formula/prior/data/CmdStan problems
  # early, once. Re-running it before every batch's parallel fitting would
  # just re-fit imputation 1's model again and again for no benefit, so it
  # runs only for this first batch; the override below disables it for all
  # subsequent batches regardless of analysis_spec$model$run_smoke_fit.
  set_runtime_override(
    paths,
    list(imputation = list(m = m_current, allow_extend = TRUE))
  )
  run_step("03_impute.R")
  run_step("04_fit_models.R")
  run_step("06_posterior_summary.R")

  # On resume, the previous run may have stopped exactly here because
  # stability was already reached; that decision is not persisted anywhere,
  # so re-evaluate the last batch transition before fitting new batches.
  stability_reached <- FALSE

  if (resumed_past_first_batch && m_current > increment_size && m_current < max_m) {
    m_previous_check <- m_current - increment_size
    stability <- evaluate_mi_stability_batches(paths, analysis_spec, m_previous_check, m_current)

    cat(
      "\nResume stability check m =", m_previous_check, "-> m =", m_current, ": ",
      stability$n_stable, "/", stability$n_parameters, "parameter(s) stable.\n"
    )

    if (stability$all_stable) {
      stability_reached <- TRUE
      cat("Stability already reached at m =", m_current, ". No further batches needed.\n")
    }
  }

  while (!stability_reached && m_current < max_m) {
    m_previous <- m_current
    m_current <- min(m_current + increment_size, max_m)

    set_runtime_override(
      paths,
      list(
        imputation = list(m = m_current, allow_extend = TRUE),
        model = list(run_smoke_fit = FALSE)
      )
    )
    run_step("03_impute.R")
    run_step("04_fit_models.R")
    run_step("06_posterior_summary.R")

    stability <- evaluate_mi_stability_batches(paths, analysis_spec, m_previous, m_current)

    cat(
      "\nStability check m =", m_previous, "-> m =", m_current, ": ",
      stability$n_stable, "/", stability$n_parameters, "parameter(s) stable.\n"
    )

    if (stability$all_stable) {
      cat("Stability reached at m =", m_current, ". Stopping the increment loop early.\n")
      break
    }

    if (m_current >= max_m) {
      cat(
        "Reached the configured maximum m =", max_m,
        "without every parameter meeting the stability thresholds.\n",
        "Consider raising analysis_spec$imputation$m if this matters for your primary parameters.\n"
      )
    }
  }

  clear_runtime_override(paths)
} else {
  run_step("03_impute.R")
  run_step("04_fit_models.R")
  run_step("06_posterior_summary.R")
}

run_step("05_diagnostics.R")
run_step("07_posterior_prediction.R")

# Step 11 produces the full multi-batch imputation-count stability tables
# and trajectory plots for whatever m the run above settled on, whether or
# not the auto-increment loop was used. It runs before Step 8 (not after,
# as in the original manual-incrementation design) because Step 8's report
# template embeds these tables/figures directly in its "Imputation-count
# stability" chapter, and needs them to already exist on disk -- both at
# write time (so the chapter content reflects this run) and at the
# self-render Step 8 performs at the end of its own script.
run_step("11_check_imputation_stability.R")

# Steps 09-10 are mo()-specific publication helpers. Run them only if the
# fitted model's own formula actually contains mo() terms, so a model
# without monotonic effects doesn't hit 10's "no main coefficient column
# found" error. Like Step 11, they run before Step 8 because Step 8's
# report template embeds their tables/figures in a "Monotonic-effect
# (mo())" chapter, and needs them to already exist on disk before that
# chapter is rendered.
model_spec <- readRDS(file.path(paths$objects, "model_spec.rds"))
mo_vars_detected <- extract_special_term_vars(model_spec$formula, fun = "mo")

if (length(mo_vars_detected) > 0) {
  cat(
    "\n\nDetected mo() term(s) in the model formula:",
    paste(mo_vars_detected, collapse = ", "), "\n"
  )

  run_step("09_check_mo_parameter_columns.R")
  run_step("10_publication_mo_results.R")
} else {
  cat("\nNo mo() terms detected in the model formula; skipping steps 09-10.\n")
}

run_step("08_publication_results.R")

# Step 12: export per-imputation draws for federated meta-analysis.
# Skipped automatically when analysis_spec$export$cohort_id is NULL.
run_step("12_export_draws.R")

cat("\nPipeline completed successfully.\n")
