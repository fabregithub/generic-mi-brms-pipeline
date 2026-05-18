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
    imputed_list <- run_row_level_imputation(dat, imputation_spec, analysis_spec)
  } else {
    stop("This template currently implements strategy = 'row_level' and 'none'. Requested: ", analysis_spec$imputation$strategy)
  }

  imputed_files <- file.path(paths$imputed_data, sprintf("imputed_%03d.rds", seq_along(imputed_list)))
  for (i in seq_along(imputed_list)) {
    saveRDS(imputed_list[[i]], imputed_files[i], compress = FALSE)
  }

  imputation_manifest <- tibble(imputation = seq_along(imputed_files), imputed_file = imputed_files)
  saveRDS(imputation_manifest, file.path(paths$objects, "imputation_manifest.rds"), compress = FALSE)
  log_msg("Saved imputation manifest with", nrow(imputation_manifest), "imputation(s).")
}, analysis_spec)
