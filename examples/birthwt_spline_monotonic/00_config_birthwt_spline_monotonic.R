# ============================================================
# When running, copy or rename this file to 00_config.R
# Generic MI + brms pipeline configuration
#
# Public custom-formula example:
# MASS::birthwt
#
# Outcome:
# low, binary indicator of low birth weight
#
# Purpose:
# Test custom brms formula terms:
#   - s(age_z): spline smooth term
#   - mo(lwt_q): monotonic ordered quintile effect
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(future)
  library(furrr)
  library(posterior)
  library(bayestestR)
})

analysis_spec <- list(

  analysis_id = "birthwt_spline_monotonic_demo",

  project_label = "Public birthwt spline + monotonic MI + brms demo",

  data = list(
    raw_data_file = "data/birthwt_spline_monotonic_example.rds",
    id_var = NULL,
    row_id_var = "row_id",
    data_structure = "single_time",
    time_var = NULL
  ),

  outcome = list(
    y_var = "low",
    family = "bernoulli",
    link = "logit",
    y_prefix = NULL,
    y_wide_regex = NULL,
    predict_missing_y = FALSE
  ),

  variables = list(
    exposure_vars = c(
      "age",
      "lwt_q",
      "race",
      "smoke",
      "ptl",
      "ht",
      "ui",
      "ftv"
    ),

    covariate_vars = character(0),

    auxiliary_vars = c("lwt"),

    continuous_vars = c(
      "age",
      "lwt",
      "ptl",
      "ftv"
    ),

    categorical_vars = c(
      "race",
      "smoke",
      "ht",
      "ui"
    ),

    ordinal_vars = c(
      "lwt_q"
    ),

    subject_level_vars = c(
      "age",
      "lwt",
      "lwt_q",
      "race",
      "smoke",
      "ptl",
      "ht",
      "ui",
      "ftv"
    ),

    time_varying_vars = character(0),

    scale_vars = c(
      "age",
      "lwt",
      "ptl",
      "ftv"
    )
  ),

  imputation = list(
    enabled = TRUE,
    strategy = "row_level",

    # Quick-test default. Increase for production testing.
    m = 5,

    maxiter = 5,
    mean_match_k = 5,
    verbose = FALSE,

    # Do not impute the outcome.
    impute_y = FALSE,

    extra_exclude_targets = character(0)
  ),

  model = list(
    fixed_effects = "auto",

    random_effects = list(
      subject_intercept = FALSE,
      subject_slope_vars = character(0)
    ),

    # Custom brms formula used to test the modified reporting scripts.
    # s(age_z) tests spline reporting.
    # mo(lwt_q) tests monotonic-effect reporting.
    custom_formula = brms::bf(
      low ~ s(age_z, k = 5) + mo(lwt_q) + race + smoke + ptl_z + ht + ui + ftv_z
    ),

    priors = "default_weakly_regularizing",

    # Quick-test defaults.
    chains = 1,
    iter = 500,
    warmup = 250,

    seed = 54321,

    adapt_delta = 0.95,
    max_treedepth = 12,

    # Include ordinary coefficients plus smooth and monotonic parameters.
    parameter_draw_regex = "^(b_|sd_|sigma|sds_|bs_|simo_|bsp_)",

    init = 0,
    silent = 0,

    run_smoke_fit = TRUE,

    skip_imputations = integer(0),
    only_imputations = integer(0)
  ),

  posterior_prediction = list(
    enabled = FALSE,
    ndraws = 200,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  ),

  summary = list(
    effects = "fixed",
    component = "conditional",
    centrality = "median",
    ci = 0.95,
    ci_method = "HDI",
    test = c("p_direction", "rope"),

    rope = list(
      method = "fixed",
      fixed_range = c(-0.1, 0.1),
      width_probability = 0.05
    )
  ),

  reporting = list(
    conditional_effects = list(
      enabled = TRUE,

      # "auto" detects variables used inside s() and mo() in custom_formula.
      # Expected here: age_z and lwt_q.
      effects = "auto",

      re_formula = NA,
      resolution = 100,
      representative_fit = "first_valid"
    )
  ),

  parallel = list(
    # ----------------------------------------------------------
    # miceRanger imputation parallelisation
    # ----------------------------------------------------------
    # impute_workers controls how many parallel miceRanger workers
    # are used.  num_impute_threads_per_worker controls the number
    # of ranger threads used inside each worker.
    #
    # Approximate CPU demand during imputation:
    # impute_workers * num_impute_threads_per_worker
    #
    # Public example default is deliberately conservative.
    # For large analyses on a high-memory machine, try:
    # impute_workers = 4
    # num_impute_threads_per_worker = 4
    impute_workers = 1,
    num_impute_threads_per_worker = 1,

    # Backward-compatible fallback used by older scripts.
    num_impute_threads = 1,

    # ----------------------------------------------------------
    # brms model fitting parallelisation
    # ----------------------------------------------------------
    # Public example default is deliberately conservative.
    # For real analyses, increase after the smoke fit succeeds.
    fit_workers = 1,
    cores_per_fit = 1,

    future_globals_maxsize_gb = 8
  ),

  memory_guard = list(
    enabled = TRUE,
    max_r_memory_gb = 20,
    min_mac_available_gb = 5,
    min_mac_available_before_brm_gb = 5,
    gc_before_check = TRUE
  )
)

paths <- list(
  root = ".",

  data = "data",

  raw_data = analysis_spec$data$raw_data_file,
  variable_dictionary = "00_variable_dictionary.csv",

  objects = "objects",

  imputed_data = "objects/imputed_data",
  imputed = "objects/imputed_data",
  imputed_wide = "objects/imputed_wide",

  model_data = "objects/model_data",

  fits = "fits",
  results = "results",

  publication = "results/publication",

  cache = file.path(path.expand("~"), ".cmdstanr-cache")
)

options(brms.backend = "cmdstanr")
options(scipen = 999)
