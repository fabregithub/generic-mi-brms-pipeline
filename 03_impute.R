source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

# -------------------------------------------------------------------------
# Step 3 overwrite/reuse guard
# -------------------------------------------------------------------------
#
# This step deliberately avoids silently overwriting existing imputed datasets.
#
# Important behaviour:
# - If the existing imputation manifest already contains exactly the requested
#   number of valid imputed-data files, Step 3 reuses them and exits.
# - If the existing manifest contains more imputations than requested, Step 3
#   stops rather than reducing m and overwriting files.
# - If the existing manifest contains fewer imputations than requested, Step 3
#   stops rather than pretending to extend the run. Extension should be handled
#   deliberately in a separate extension workflow.
#
# This prevents accidental damage when a user changes, for example, m = 100
# to m = 45 in 00_config.R and reruns run_all.R.

valid_existing_imputation_files <- function(paths) {
  manifest_file <- file.path(paths$objects, "imputation_manifest.rds")

  if (!file.exists(manifest_file)) {
    return(character(0))
  }

  manifest <- readRDS(manifest_file)

  if (!is.data.frame(manifest) || !"imputed_file" %in% names(manifest)) {
    return(character(0))
  }

  files <- as.character(manifest$imputed_file)
  files <- files[!is.na(files) & nzchar(files)]
  files[file.exists(files)]
}

safe_step("STEP 3: Imputation", {
  dat <- readRDS(file.path(paths$objects, "dat_prepared.rds"))
  var_dict <- readRDS(file.path(paths$objects, "variable_dictionary.rds"))

  target_m <- analysis_spec$imputation$m %||% 1
  target_m <- as.integer(target_m)

  if (!is.finite(target_m) || target_m < 1) {
    stop("analysis_spec$imputation$m must be a positive integer.")
  }

  existing_imputed_files <- valid_existing_imputation_files(paths)
  n_existing <- length(existing_imputed_files)

  if (n_existing > 0) {
    log_msg("Existing valid imputed dataset files found:", n_existing)
    log_msg("Current requested m:", target_m)

    if (n_existing == target_m) {
      log_msg(
        "Existing imputed datasets match requested m. ",
        "Skipping imputation and reusing existing files."
      )

      imputation_manifest <- tibble::tibble(
        imputation = seq_len(n_existing),
        imputed_file = existing_imputed_files
      )

      saveRDS(
        imputation_manifest,
        file.path(paths$objects, "imputation_manifest.rds"),
        compress = FALSE
      )

      log_msg("Re-saved imputation manifest with", nrow(imputation_manifest), "existing imputation(s).")
      return(invisible(imputation_manifest))
    }

    if (n_existing > target_m) {
      stop(
        "Existing imputed datasets were created with a larger m than the current config. ",
        "Refusing to reduce m and overwrite existing imputed data. ",
        "Existing valid imputations: ", n_existing, "; requested m: ", target_m, ". ",
        "Safer options: keep the larger m, run only downstream scripts on completed fits, ",
        "or start a clean new run in a separate output folder."
      )
    }

    if (n_existing < target_m) {
      stop(
        "Existing imputed datasets were found, but fewer than the current requested m. ",
        "Refusing to overwrite them automatically. ",
        "Existing valid imputations: ", n_existing, "; requested m: ", target_m, ". ",
        "To extend m safely, use a deliberate extension workflow or start a clean new run. ",
        "Do not rely on Step 3 to overwrite an existing imputation set."
      )
    }
  }

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

  # Final defensive check before writing.
  existing_target_files <- imputed_files[file.exists(imputed_files)]
  if (length(existing_target_files) > 0) {
    stop(
      "Refusing to overwrite existing imputed dataset file(s). ",
      "First existing file: ", existing_target_files[[1]], ". ",
      "Move or delete old imputation files only if you intentionally want a clean re-imputation."
    )
  }

  for (i in seq_along(imputed_list)) {
    saveRDS(imputed_list[[i]], imputed_files[i], compress = FALSE)
  }

  imputation_manifest <- tibble::tibble(imputation = seq_along(imputed_files), imputed_file = imputed_files)
  saveRDS(imputation_manifest, file.path(paths$objects, "imputation_manifest.rds"), compress = FALSE)
  log_msg("Saved imputation manifest with", nrow(imputation_manifest), "imputation(s).")
}, analysis_spec)
