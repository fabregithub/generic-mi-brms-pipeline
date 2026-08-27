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

  if (identical(analysis_spec$imputation$strategy, "censored_exposure_block_fcs")) {
    if (!requireNamespace("leftcens", quietly = TRUE) ||
        !"impute_censored_conditional" %in% getNamespaceExports("leftcens")) {
      stop(
        "strategy = 'censored_exposure_block_fcs' requires leftcens >= 0.9.0 ",
        "(exporting impute_censored_conditional()). Install it with: ",
        "remotes::install_github(\"fabregithub/leftcens@v0.9.0\")"
      )
    }

    ce <- analysis_spec$imputation$censored_exposure

    if (is.null(ce) || length(ce$exposure_vars) == 0) {
      stop(
        "strategy = 'censored_exposure_block_fcs' requires ",
        "analysis_spec$imputation$censored_exposure$exposure_vars."
      )
    }

    lo_suffix <- ce$lo_suffix %||% "_lo"
    hi_suffix <- ce$hi_suffix %||% "_hi"

    check_required_vars(raw_data, ce$exposure_vars, "censored exposure variables")
    check_required_vars(
      raw_data,
      c(paste0(ce$exposure_vars, lo_suffix), paste0(ce$exposure_vars, hi_suffix)),
      paste0("censored exposure bound columns (", lo_suffix, "/", hi_suffix, ")")
    )

    margin <- ce$margin %||% "shash"
    if (!margin %in% c("shash", "gaussian")) {
      stop(
        "censored_exposure$margin must be 'shash' or 'gaussian', not '",
        margin, "'."
      )
    }

    # The X block conditions on this predictor set. Including the outcome is the
    # whole point of this strategy (congeniality); omitting it reproduces the
    # biased no-Y pre-step. Auto (NULL) always includes the outcome.
    if (!is.null(ce$predictors)) {
      # BACK-COMPATIBILITY: the X block does `intersect(predictors, names(data))`,
      # so a predictor that is not in the data has always been silently dropped
      # rather than being an error. Warn (making the drop visible) instead of
      # stopping, so a config that ran before still runs.
      unknown_preds <- setdiff(ce$predictors, names(raw_data))
      if (length(unknown_preds) > 0) {
        warning(
          "censored_exposure$predictors names not present in the data: ",
          paste(unknown_preds, collapse = ", "),
          ". These are silently dropped by the X block; remove them from the ",
          "config or correct the spelling."
        )
      }

      if (!analysis_spec$outcome$y_var %in% ce$predictors) {
        warning(
          "censored_exposure$predictors does not include the outcome '",
          analysis_spec$outcome$y_var,
          "'. The X block will then impute the exposure WITHOUT the outcome, ",
          "which is the uncongenial pre-step this strategy exists to avoid ",
          "(attenuates the exposure-response coefficient)."
        )
      }
    }

    for (x in ce$exposure_vars) {
      lo_col <- paste0(x, lo_suffix)
      hi_col <- paste0(x, hi_suffix)
      lo <- raw_data[[lo_col]]
      hi <- raw_data[[hi_col]]

      if (!is.numeric(lo) || !is.numeric(hi)) {
        stop(
          "Censored exposure bound columns must be numeric: ",
          lo_col, " / ", hi_col, "."
        )
      }

      if (all(is.na(lo)) || all(is.na(hi))) {
        stop(
          "Censored exposure '", x, "' has an all-missing bound column (",
          lo_col, " / ", hi_col, ")."
        )
      }

      inverted <- is.finite(lo) & is.finite(hi) & lo > hi
      if (any(inverted)) {
        stop(
          "Censored exposure '", x, "' has ", sum(inverted),
          " row(s) where ", lo_col, " > ", hi_col,
          ". Bounds must satisfy lo <= hi."
        )
      }

      n_cens <- sum(!(is.finite(lo) & is.finite(hi) & abs(hi - lo) < 1e-9))
      log_msg(
        "Censored exposure", x, ":", n_cens, "of", nrow(raw_data),
        "row(s) censored/interval-valued"
      )
    }

    # Dictionary flags: the exposure must reach the model, and must NOT be a
    # miceRanger target (the X block owns it).
    ce_dict <- var_dict %>%
      dplyr::filter(.data$var %in% ce$exposure_vars)

    not_in_model <- ce_dict %>%
      dplyr::filter(!.data$use_in_model %in% TRUE) %>%
      dplyr::pull(.data$var)

    if (length(not_in_model) > 0) {
      # BACK-COMPATIBILITY: `use_in_model` feeds only the AUTO predictor set
      # (00_censored_exposure.R builds `auto_preds` from it). When `predictors`
      # is given explicitly the flag is never consulted, so a config with
      # use_in_model = FALSE plus an explicit predictor set and custom_formula
      # has always worked. Only fail in the auto case, where it genuinely breaks.
      if (is.null(ce$predictors)) {
        stop(
          "Censored exposure(s) must be marked use_in_model = TRUE in ",
          "00_variable_dictionary.csv when censored_exposure$predictors is not ",
          "set explicitly (the X block derives its predictors from that flag): ",
          paste(not_in_model, collapse = ", ")
        )
      }
      warning(
        "Censored exposure(s) marked use_in_model = FALSE: ",
        paste(not_in_model, collapse = ", "),
        ". Harmless here because censored_exposure$predictors is set explicitly, ",
        "but check the exposure still reaches the model formula."
      )
    }

    as_target <- ce_dict %>%
      dplyr::filter(.data$impute_target %in% TRUE) %>%
      dplyr::pull(.data$var)

    if (length(as_target) > 0) {
      warning(
        "Censored exposure(s) marked impute_target = TRUE in the dictionary: ",
        paste(as_target, collapse = ", "),
        ". The censored X block imputes these, not miceRanger; set ",
        "impute_target = FALSE to make that explicit."
      )
    }

    log_msg(
      "Validated censored-exposure block-FCS strategy |",
      "exposures:", paste(ce$exposure_vars, collapse = ", "),
      "| margin:", margin,
      "| outer_sweeps:", ce$outer_sweeps %||% 5L,
      "| MID:", isTRUE(ce$mid_delete_imputed_y %||% TRUE)
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
