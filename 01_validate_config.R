source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)
setup_info <- setup_brms_cmdstan(paths$cache)

safe_step("STEP 1: Validate configuration", {
  if (!file.exists(paths$raw_data)) {
    stop(
      "Raw data file not found: ",
      paths$raw_data,
      ". Create the input data file specified in 00_config.R first. ",
      "For the airquality example, run: ",
      "Rscript examples/airquality_gaussian/00_create_airquality_example_data.R"
    )
  }
  if (!file.exists(paths$variable_dictionary)) {
    stop("Variable dictionary not found: ", paths$variable_dictionary)
  }

  log_msg("Reading raw data:", paths$raw_data)
  raw_data <- readRDS(paths$raw_data)
  var_dict <- read_var_dict(paths$variable_dictionary)

  required <- unique(c(var_dict$var, analysis_spec$outcome$y_var, analysis_spec$data$id_var, analysis_spec$data$time_var))
  check_required_vars(raw_data, required, "variables listed in config/dictionary")

  if (identical(analysis_spec$imputation$strategy, "subject_wide_with_repeated_y_auxiliary")) {
    if (is.null(analysis_spec$data$id_var) ||
        is.null(analysis_spec$data$time_var) ||
        is.na(analysis_spec$data$id_var) ||
        is.na(analysis_spec$data$time_var) ||
        !nzchar(analysis_spec$data$id_var) ||
        !nzchar(analysis_spec$data$time_var)) {
      stop(
        "strategy = 'subject_wide_with_repeated_y_auxiliary' requires ",
        "analysis_spec$data$id_var and analysis_spec$data$time_var."
      )
    }

    check_required_vars(
      raw_data,
      c(
        analysis_spec$data$id_var,
        analysis_spec$data$time_var,
        analysis_spec$outcome$y_var
      ),
      "subject-wide repeated outcome variables"
    )

    if (!"timing" %in% names(var_dict)) {
      stop(
        "strategy = 'subject_wide_with_repeated_y_auxiliary' requires ",
        "a timing column in 00_variable_dictionary.csv."
      )
    }

    n_subject_level <- var_dict %>%
      dplyr::filter(.data$timing %in% c("single", "baseline")) %>%
      nrow()

    if (n_subject_level == 0) {
      stop(
        "No subject-level variables found in 00_variable_dictionary.csv. ",
        "For subject-wide imputation, mark baseline/single-measure covariates ",
        "with timing = 'single' or timing = 'baseline'."
      )
    }

    log_msg(
      "Validated subject-wide repeated-outcome imputation strategy with",
      n_subject_level,
      "subject-level variable row(s) in dictionary"
    )
  }


  dat <- prepare_raw_data(raw_data, analysis_spec, var_dict)
  model_spec <- build_model_spec(analysis_spec, var_dict, dat)

  dat2 <- apply_z_stats(dat, model_spec$z_stats)
  formula_vars <- model_spec$formula_vars %||% get_brms_formula_vars(model_spec$formula)
  check_required_vars(dat2, formula_vars, "formula variables after transformation")

  # ------------------------------------------------------------
  # Additional validation for brms custom formula terms
  # ------------------------------------------------------------

  smooth_vars <- extract_special_term_vars(model_spec$formula, "s")
  monotonic_vars <- extract_special_term_vars(model_spec$formula, "mo")

  if (length(smooth_vars) > 0) {
    check_required_vars(dat2, smooth_vars, "variables used inside s()")

    non_numeric_smooths <- smooth_vars[
      !vapply(dat2[smooth_vars], is.numeric, logical(1))
    ]

    if (length(non_numeric_smooths) > 0) {
      stop(
        "Variables used inside s() should be numeric after transformation: ",
        paste(non_numeric_smooths, collapse = ", ")
      )
    }

    log_msg("Detected smooth term(s) s() for:", paste(smooth_vars, collapse = ", "))
  }

  if (length(monotonic_vars) > 0) {
    check_required_vars(dat2, monotonic_vars, "variables used inside mo()")

    bad_mo_vars <- monotonic_vars[
      !vapply(
        dat2[monotonic_vars],
        function(x) {
          is.ordered(x) ||
            is.integer(x) ||
            (is.numeric(x) && all(is.na(x) | abs(x - round(x)) < .Machine$double.eps^0.5))
        },
        logical(1)
      )
    ]

    if (length(bad_mo_vars) > 0) {
      stop(
        "Variables used inside mo() should be ordered factors or integer-like numeric variables: ",
        paste(bad_mo_vars, collapse = ", "),
        ". In the variable dictionary, these are usually type = ordinal or type = integer."
      )
    }

    mo_dict <- var_dict %>%
      dplyr::filter(.data$var %in% monotonic_vars)

    non_ordinal_mo <- mo_dict %>%
      dplyr::filter(!.data$type %in% c("ordinal", "integer")) %>%
      dplyr::pull(.data$var)

    if (length(non_ordinal_mo) > 0) {
      warning(
        "The following variables are used inside mo(), but their dictionary type is not ordinal/integer: ",
        paste(non_ordinal_mo, collapse = ", "),
        ". The model may still fit if the variables are ordered factors or integer-like."
      )
    }

    log_msg("Detected monotonic term(s) mo() for:", paste(monotonic_vars, collapse = ", "))
  }

  log_msg("Formula:", paste(deparse(model_spec$formula), collapse = " "))
  log_msg("Family:", analysis_spec$outcome$family, "link:", analysis_spec$outcome$link)
  log_msg("Validation completed.")
}, analysis_spec)
