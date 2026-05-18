source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 2: Prepare data", {
  raw_data <- readRDS(paths$raw_data)
  var_dict <- read_var_dict(paths$variable_dictionary)
  dat <- prepare_raw_data(raw_data, analysis_spec, var_dict)

  saveRDS(dat, file.path(paths$objects, "dat_prepared.rds"), compress = FALSE)
  saveRDS(analysis_spec, file.path(paths$objects, "analysis_spec.rds"), compress = FALSE)
  saveRDS(var_dict, file.path(paths$objects, "variable_dictionary.rds"), compress = FALSE)

  log_msg("Saved dat_prepared.rds with", nrow(dat), "rows and", ncol(dat), "columns.")
}, analysis_spec)
