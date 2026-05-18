source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)
setup_info <- setup_brms_cmdstan(paths$cache)

safe_step("STEP 1: Validate configuration", {
  if (!file.exists(paths$raw_data)) {
    stop("Raw data file not found: ", paths$raw_data, ". Run 00_create_airquality_example_data.R first.")
  }
  if (!file.exists(paths$variable_dictionary)) {
    stop("Variable dictionary not found: ", paths$variable_dictionary)
  }

  log_msg("Reading raw data:", paths$raw_data)
  raw_data <- readRDS(paths$raw_data)
  var_dict <- read_var_dict(paths$variable_dictionary)

  required <- unique(c(var_dict$var, analysis_spec$outcome$y_var, analysis_spec$data$id_var, analysis_spec$data$time_var))
  check_required_vars(raw_data, required, "variables listed in config/dictionary")

  dat <- prepare_raw_data(raw_data, analysis_spec, var_dict)
  model_spec <- build_model_spec(analysis_spec, var_dict, dat)

  dat2 <- apply_z_stats(dat, model_spec$z_stats)
  formula_vars <- all.vars(model_spec$formula)
  check_required_vars(dat2, formula_vars, "formula variables after transformation")

  log_msg("Formula:", paste(deparse(model_spec$formula), collapse = " "))
  log_msg("Family:", analysis_spec$outcome$family, "link:", analysis_spec$outcome$link)
  log_msg("Validation completed.")
}, analysis_spec)
