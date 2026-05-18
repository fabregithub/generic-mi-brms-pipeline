source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 5: Diagnostics", {
  fit_manifest <- readRDS(file.path(paths$objects, "fit_manifest.rds"))
  model_spec <- readRDS(file.path(paths$objects, "model_spec.rds"))

  fit_manifest <- fit_manifest %>% filter(purrr::map_lgl(fit_file, rds_ok))
  if (nrow(fit_manifest) == 0) stop("No valid fit files found.")

  diagnostics <- purrr::map2_dfr(
    fit_manifest$fit_file,
    fit_manifest$imputation,
    ~ diagnose_one_fit(.x, .y, model_spec$max_treedepth)
  )

  saveRDS(diagnostics, file.path(paths$results, "diagnostics.rds"), compress = FALSE)
  readr::write_csv(diagnostics, file.path(paths$results, "diagnostics.csv"))
  log_msg("Saved diagnostics for", nrow(diagnostics), "fit(s).")
}, analysis_spec)
