# ============================================================
# When running, rename this file to 00_config.R
# Generic MI + brms pipeline configuration
#
# Public Cox proportional hazards example:
# survival::lung
#
# Outcome:
# time to death (days), right-censored
# Predictors: age, sex, ECOG score, Karnofsky score,
#             weight loss, meal calories
# Reference: Loprinzi et al. (1994) J Clin Oncol
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

  analysis_id = "lung_cox_demo",

  project_label = "Public lung survival Cox MI + brms demo",

  data = list(
    raw_data_file = "data/lung_cox_example.rds",
    id_var = NULL,
    row_id_var = "row_id",
    data_structure = "single_time",
    time_var = NULL
  ),

  outcome = list(
    # For Cox PH the response is specified via custom_formula below as
    # time | cens(censored) ~ ..., so y_var points to the time column.
    y_var = "time",

    # Options supported by the template:
    # "gaussian", "bernoulli", "poisson", "negbinomial",
    # "beta", "ordinal", "categorical", "cox"
    family = "cox",

    # Cox uses a log link (parameters are log hazard ratios).
    link = "log",

    y_prefix = NULL,
    y_wide_regex = NULL,

    # Survival time is always observed; there are no missing-outcome rows
    # to predict for Cox models.
    predict_missing_y = FALSE
  ),

  variables = NULL,

  imputation = list(
    enabled = TRUE,
    strategy = "row_level",

    # Use small m for demo speed; increase to 100 for production.
    m = 100,

    maxiter = 5,
    mean_match_k = 5,
    verbose = FALSE,

    # Do not impute the survival time or censoring indicator.
    impute_y = FALSE,

    extra_exclude_targets = character(0)
  ),

  model = list(
    fixed_effects = "auto",

    random_effects = list(
      subject_intercept = FALSE,
      subject_slope_vars = character(0)
    ),

    # Cox PH response must be specified as time | cens(censored).
    # The pipeline passes this formula directly to brms.
    custom_formula = brmsformula(
      time | cens(censored) ~ age_z + sex + ph_ecog + ph_karno_z + wt_loss_z + meal_cal_z
    ),

    priors = "default_weakly_regularizing",

    chains = 4,
    iter = 2000,
    warmup = 1000,

    seed = 54321,

    adapt_delta = 0.95,
    max_treedepth = 12,

    # Cox PH parameters: b_ (log hazard ratios)
    parameter_draw_regex = "^b_",

    init = "random",
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
      # log(c(0.90, 1/0.90)) corresponds to a +/-10% hazard ratio band,
      # a reasonable default for a Cox survival model.
      method = "fixed",
      fixed_range = log(c(0.90, 1 / 0.90)),
      width_probability = 0.10
    )
  ),

  parallel = list(
    impute_workers = 1,
    num_impute_threads_per_worker = 1,
    num_impute_threads = 1,
    fit_workers = 1,
    cores_per_fit = 1,
    summary_workers = 2,
    prediction_workers = 2,
    future_globals_maxsize_gb = 8
  ),

  mi_stability = list(
    auto_increment = TRUE,
    increment_size = NULL,
    parameter_regex = "^b_",
    exclude_intercept = TRUE,
    estimate_tolerance = 0.05,
    ci_endpoint_tolerance = 0.05,
    relative_transformed_tolerance_pct = 5,
    pd_tolerance = 0.02
  ),

  mo_effects = NULL,

  publication = list(
    template_sentences_scope = "all"
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
