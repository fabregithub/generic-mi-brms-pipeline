# Fit one imputation only, e.g.:
# Rscript fit_single_imputation.R 3

source("00_config.R")
source("00_common_functions.R")

init_logging("single_fit")
setup_project_dirs(paths)
setup_info <- setup_brms_cmdstan(paths$cache)

args <- commandArgs(trailingOnly = TRUE)
target_imp <- if (length(args) >= 1) as.integer(args[1]) else 1

overwrite_existing <- FALSE

safe_step(paste("Fit single imputation", target_imp), {
  dat <- readRDS(file.path(paths$objects, "dat_prepared.rds"))
  var_dict <- readRDS(file.path(paths$objects, "variable_dictionary.rds"))
  imputation_manifest <- readRDS(file.path(paths$objects, "imputation_manifest.rds"))
  model_spec <- build_model_spec(analysis_spec, var_dict, dat)
  model_data_manifest <- prepare_model_data_files(imputation_manifest, analysis_spec, model_spec, paths$model_data)

  row <- model_data_manifest %>% filter(imputation == target_imp)
  if (nrow(row) != 1) stop("Could not find imputation ", target_imp)

  fit_file <- file.path(paths$fits, sprintf("fit_imp_%03d.rds", target_imp))
  fit_manifest <- row %>% mutate(fit_file = fit_file)

  if (rds_ok(fit_file) && !overwrite_existing) {
    log_msg("Existing valid fit found:", fit_file)
  } else {
    status <- fit_one_brm_file(1, fit_manifest, model_spec, analysis_spec, setup_info, log_file = paste0("fit_imp_", target_imp, "_only.log"))
    print(status)
  }
}, analysis_spec)
