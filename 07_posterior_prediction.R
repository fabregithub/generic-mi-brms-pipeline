# ============================================================
# 07_posterior_prediction.R
# Posterior prediction for missing outcome rows
#
# Parallel version:
# - predicts per imputation in parallel
# - writes one RDS file per imputation
# - combines prediction files after workers finish
# ============================================================

source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 7: Posterior prediction for missing outcome rows", {

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  get_parallel_workers <- function(analysis_spec, field, fallback_field = "fit_workers") {
    workers <- analysis_spec$parallel[[field]] %||%
      analysis_spec$parallel[[fallback_field]] %||%
      1L

    workers <- as.integer(workers)

    if (is.na(workers) || workers < 1L) {
      workers <- 1L
    }

    workers
  }

  if (!isTRUE(analysis_spec$posterior_prediction$enabled) ||
      !isTRUE(analysis_spec$outcome$predict_missing_y)) {
    log_msg("Posterior prediction disabled; skipping.")
  } else {

    fit_manifest <- readRDS(file.path(paths$objects, "fit_manifest.rds")) %>%
      dplyr::filter(purrr::map_lgl(fit_file, rds_ok))

    model_spec <- readRDS(file.path(paths$objects, "model_spec.rds"))

    if (nrow(fit_manifest) == 0) {
      stop("No valid fit files found.")
    }

    if (!"pred_file" %in% names(fit_manifest)) {
      stop("fit_manifest does not contain pred_file. Run Step 4/prepare model data with prediction manifest support.")
    }

    fit_manifest <- fit_manifest %>%
      dplyr::filter(purrr::map_lgl(pred_file, rds_ok))

    if (nrow(fit_manifest) == 0) {
      stop("No valid prediction data files found.")
    }

    missing_y_draw_files <- file.path(
      paths$results,
      sprintf("missing_y_draws_imp_%03d.rds", fit_manifest$imputation)
    )

    prediction_workers <- get_parallel_workers(
      analysis_spec = analysis_spec,
      field = "prediction_workers",
      fallback_field = "fit_workers"
    )

    prediction_workers <- min(prediction_workers, nrow(fit_manifest))

    log_msg("STEP 7 prediction_workers:", prediction_workers)

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    old_maxsize <- getOption("future.globals.maxSize")
    on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)

    options(
      future.globals.maxSize =
        (analysis_spec$parallel$future_globals_maxsize_gb %||% 8) * 1024^3
    )

    future::plan(
      future::multisession,
      workers = prediction_workers
    )

    prediction_status <- furrr::future_pmap_dfr(
      list(
        imp_i = fit_manifest$imputation,
        fit_file_i = fit_manifest$fit_file,
        pred_file_i = fit_manifest$pred_file,
        draw_file_i = missing_y_draw_files
      ),
      function(imp_i, fit_file_i, pred_file_i, draw_file_i) {
        if (rds_ok(draw_file_i)) {
          return(
            tibble::tibble(
              imputation = imp_i,
              status = "skipped_existing_valid_draws",
              missing_y_draw_file = draw_file_i
            )
          )
        }

        fit_i <- readRDS(fit_file_i)
        pred_i <- readRDS(pred_file_i)

        draws_i <- predict_missing_y_draws_one(
          fit_i,
          pred_i,
          analysis_spec,
          model_spec,
          ndraws = analysis_spec$posterior_prediction$ndraws
        ) %>%
          dplyr::mutate(imputation = imp_i, .before = 1)

        saveRDS(
          draws_i,
          draw_file_i,
          compress = FALSE
        )

        rm(fit_i, pred_i, draws_i)
        gc()

        tibble::tibble(
          imputation = imp_i,
          status = "completed",
          missing_y_draw_file = draw_file_i
        )
      },
      .options = furrr::furrr_options(
        seed = TRUE,
        packages = c(
          "brms",
          "posterior",
          "tibble",
          "dplyr",
          "purrr"
        )
      )
    )

    readr::write_csv(
      prediction_status,
      file.path(paths$results, "missing_y_prediction_status.csv")
    )

    prediction_manifest <- fit_manifest %>%
      dplyr::mutate(missing_y_draw_file = missing_y_draw_files)

    saveRDS(
      prediction_manifest,
      file.path(paths$objects, "prediction_manifest.rds"),
      compress = FALSE
    )

    valid_missing_y_draw_files <- missing_y_draw_files[
      purrr::map_lgl(missing_y_draw_files, rds_ok)
    ]

    if (length(valid_missing_y_draw_files) == 0) {
      stop("No valid missing-y posterior draw files found.")
    }

    missing_y_draws <- purrr::map_dfr(
      valid_missing_y_draw_files,
      readRDS
    )

    saveRDS(
      missing_y_draws,
      file.path(paths$results, "missing_y_draws.rds"),
      compress = FALSE
    )

    missing_y_summary <- summarise_missing_y_draws(
      missing_y_draws,
      analysis_spec,
      analysis_spec$summary
    )

    saveRDS(
      missing_y_summary,
      file.path(paths$results, "missing_y_summary.rds"),
      compress = FALSE
    )

    readr::write_csv(
      missing_y_summary,
      file.path(paths$results, "missing_y_summary.csv")
    )

    log_msg("Saved posterior predictions for missing outcome rows.")

    rm(
      fit_manifest,
      model_spec,
      missing_y_draw_files,
      prediction_manifest,
      valid_missing_y_draw_files,
      missing_y_draws,
      missing_y_summary
    )

    gc()

    guard_memory("after STEP 7 cleanup")
  }
}, analysis_spec)
