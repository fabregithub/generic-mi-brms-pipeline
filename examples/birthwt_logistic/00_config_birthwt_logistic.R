# ============================================================
# When running, rename this file to 00_config.R
# Generic MI + brms pipeline configuration
#
# Public logistic example:
# MASS::birthwt
#
# Outcome:
# low, binary indicator of low birth weight
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

  analysis_id = "birthwt_logistic_demo",

  project_label = "Public birthwt logistic MI + brms demo",

  data = list(
    raw_data_file = "data/birthwt_logistic_example.rds",
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

  variables = NULL,

  imputation = list(
    enabled = TRUE,
    strategy = "row_level",

    # Use small m for demo speed.
    # For stress testing, set m = 100.
    m = 100,

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

    custom_formula = NULL,

    priors = "default_weakly_regularizing",

    chains = 4,
    iter = 2000,
    warmup = 1000,

    seed = 54321,

    adapt_delta = 0.95,
    max_treedepth = 12,

    parameter_draw_regex = "^(b_|sd_|sigma)",

    init = 0,

    # Set to 0 for visible Stan output during testing.
    # Set to 2 for quieter production runs.
    silent = 0,

    run_smoke_fit = TRUE,

    skip_imputations = integer(0),
    only_imputations = integer(0)
  ),

  posterior_prediction = list(
    enabled = FALSE,
    ndraws = 500,
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

      # Log-odds ROPE placeholder for demo.
      # For real analyses, define this based on practical relevance.
      fixed_range = c(-0.1, 0.1),

      width_probability = 0.05
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
