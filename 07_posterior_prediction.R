source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 7: Posterior prediction for missing outcome rows", {
  if (!isTRUE(analysis_spec$posterior_prediction$enabled) || !isTRUE(analysis_spec$outcome$predict_missing_y)) {
    log_msg("Posterior prediction disabled; skipping.")
  } else {
    fit_manifest <- readRDS(file.path(paths$objects, "fit_manifest.rds")) %>%
      filter(purrr::map_lgl(fit_file, rds_ok))
    model_spec <- readRDS(file.path(paths$objects, "model_spec.rds"))

    if (nrow(fit_manifest) == 0) stop("No valid fit files found.")

    missing_y_draw_files <- file.path(paths$results, sprintf("missing_y_draws_imp_%03d.rds", fit_manifest$imputation))

    for (i in seq_len(nrow(fit_manifest))) {
      if (rds_ok(missing_y_draw_files[i])) next
      fit_i <- readRDS(fit_manifest$fit_file[i])
      pred_i <- readRDS(fit_manifest$pred_file[i])
      draws_i <- predict_missing_y_draws_one(
        fit_i,
        pred_i,
        analysis_spec,
        model_spec,
        ndraws = analysis_spec$posterior_prediction$ndraws
      ) %>% mutate(imputation = fit_manifest$imputation[i], .before = 1)
      saveRDS(draws_i, missing_y_draw_files[i], compress = FALSE)
      rm(fit_i, pred_i, draws_i); gc()
    }

    prediction_manifest <- fit_manifest %>% mutate(missing_y_draw_file = missing_y_draw_files)
    saveRDS(prediction_manifest, file.path(paths$objects, "prediction_manifest.rds"), compress = FALSE)

    missing_y_draws <- purrr::map_dfr(missing_y_draw_files[file.exists(missing_y_draw_files)], readRDS)
    saveRDS(missing_y_draws, file.path(paths$results, "missing_y_draws.rds"), compress = FALSE)

    missing_y_summary <- summarise_missing_y_draws(missing_y_draws, analysis_spec, analysis_spec$summary)
    saveRDS(missing_y_summary, file.path(paths$results, "missing_y_summary.rds"), compress = FALSE)
    readr::write_csv(missing_y_summary, file.path(paths$results, "missing_y_summary.csv"))

    log_msg("Saved posterior predictions for missing outcome rows.")
  }
}, analysis_spec)
