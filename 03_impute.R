source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 3: Imputation", {
  dat <- readRDS(file.path(paths$objects, "dat_prepared.rds"))
  var_dict <- readRDS(file.path(paths$objects, "variable_dictionary.rds"))

  if (!isTRUE(analysis_spec$imputation$enabled) || identical(analysis_spec$imputation$strategy, "none")) {
    imputed_list <- list(dat)
  } else if (identical(analysis_spec$imputation$strategy, "row_level")) {
    if (!requireNamespace("miceRanger", quietly = TRUE)) {
      stop("Package 'miceRanger' is required for imputation.")
    }
    imputation_spec <- make_row_level_imputation_spec(dat, analysis_spec, var_dict)
    saveRDS(imputation_spec, file.path(paths$objects, "imputation_spec.rds"), compress = FALSE)

    log_msg(
      "Imputation parallel settings | impute_workers:",
      analysis_spec$parallel$impute_workers %||% 1,
      "| num_impute_threads_per_worker:",
      analysis_spec$parallel$num_impute_threads_per_worker %||%
        analysis_spec$parallel$num_impute_threads %||%
        1
    )

    imputed_list <- run_row_level_imputation(dat, imputation_spec, analysis_spec)
  } else if (identical(analysis_spec$imputation$strategy, "subject_wide_with_repeated_y_auxiliary")) {
    if (!requireNamespace("miceRanger", quietly = TRUE)) {
      stop("Package 'miceRanger' is required for imputation.")
    }

    log_msg("Building subject-wide imputation data")

    wide_data <- make_subject_wide_imputation_data(
      data = dat,
      analysis_spec = analysis_spec,
      var_dict = var_dict
    )

    saveRDS(
      wide_data$subject_wide,
      file.path(paths$objects, "subject_wide.rds"),
      compress = FALSE
    )

    saveRDS(
      wide_data$long_base,
      file.path(paths$objects, "long_base.rds"),
      compress = FALSE
    )

    log_msg("Subject-wide rows:", nrow(wide_data$subject_wide))
    log_msg("Long rows:", nrow(wide_data$long_base))

    imputation_spec <- make_subject_wide_imputation_spec(
      subject_wide = wide_data$subject_wide,
      analysis_spec = analysis_spec,
      var_dict = var_dict
    )

    saveRDS(imputation_spec, file.path(paths$objects, "imputation_spec.rds"), compress = FALSE)

    log_msg(
      "Imputation parallel settings | impute_workers:",
      analysis_spec$parallel$impute_workers %||% 1,
      "| num_impute_threads_per_worker:",
      analysis_spec$parallel$num_impute_threads_per_worker %||%
        analysis_spec$parallel$num_impute_threads %||%
        1
    )

    imputed_wide_list <- run_row_level_imputation(
      data = wide_data$subject_wide,
      imputation_spec = imputation_spec,
      analysis_spec = analysis_spec
    )

    imputed_list <- prepare_long_imputed_from_subject_wide(
      imputed_wide_list = imputed_wide_list,
      long_base = wide_data$long_base,
      analysis_spec = analysis_spec
    )

    rm(wide_data, imputed_wide_list)
    gc()
  } else {
    stop(
      "This template currently implements strategy = 'row_level', ",
      "'subject_wide_with_repeated_y_auxiliary', and 'none'. Requested: ",
      analysis_spec$imputation$strategy
    )
  }

  imputed_files <- file.path(paths$imputed_data, sprintf("imputed_%03d.rds", seq_along(imputed_list)))
  for (i in seq_along(imputed_list)) {
    saveRDS(imputed_list[[i]], imputed_files[i], compress = FALSE)
  }

  imputation_manifest <- tibble(imputation = seq_along(imputed_files), imputed_file = imputed_files)
  saveRDS(imputation_manifest, file.path(paths$objects, "imputation_manifest.rds"), compress = FALSE)
  log_msg("Saved imputation manifest with", nrow(imputation_manifest), "imputation(s).")
}, analysis_spec)
