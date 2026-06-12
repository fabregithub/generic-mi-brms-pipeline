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

    # Options supported by the template:
    # "gaussian", "bernoulli", "poisson", "negbinomial",
    # "beta", "ordinal", "categorical"
    family = "bernoulli",

    # Link function passed to brms family construction.
    #
    # Common examples:
    # gaussian    -> "identity", "log"
    # bernoulli   -> "logit", "probit", "cloglog"
    # poisson     -> "log", "identity"
    # negbinomial -> "log", "identity"
    # beta        -> "logit", "probit", "cloglog", "log"
    # ordinal     -> commonly "logit" or "probit", depending on the family
    # categorical -> commonly "logit"
    #
    # Choose a link that is valid for the selected family and appropriate
    # for the scientific question.
    link = "logit",

    # Only used for repeated-Y wide imputation.
    #
    # If the long outcome is ps and time is 1, 2, ..., 6,
    # the subject-wide auxiliary outcome columns may be:
    # ps_1, ps_2, ps_3, ps_4, ps_5, ps_6
    #
    # Example:
    # y_prefix = "ps_"
    # y_wide_regex = "^ps_"
    #
    # For single-time analyses, keep both as NULL.
    y_prefix = NULL,
    y_wide_regex = NULL,

    # If TRUE, rows with missing outcome are excluded from model fitting
    # but posterior predictions are generated for them later.
    predict_missing_y = FALSE
  ),

  # ------------------------------------------------------------
  # Variable roles
  # ------------------------------------------------------------
  # Variable roles, types, timing, scaling, imputation targets,
  # model inclusion, reference categories, and auxiliary-variable status are
  # read from 00_variable_dictionary.csv.
  #
  # Leave this as NULL for standard use. This keeps the dictionary as the
  # default source of truth.
  #
  # Advanced users may provide selected overrides here. If supplied,
  # these values override the groups derived from the dictionary.
  #
  # Complete override template:
  #
  # variables = list(
  #   exposure_vars = c("exposure"),
  #   covariate_vars = c("age", "sex", "income"),
  #   auxiliary_vars = c("baseline_score"),
  #
  #   continuous_vars = c("age", "baseline_score"),
  #   categorical_vars = c("sex"),
  #   ordinal_vars = c("income"),
  #
  #   subject_level_vars = c("age", "sex", "income", "baseline_score"),
  #   time_varying_vars = c("current_exposure"),
  #
  #   # Variables to z-scale for modelling. The model uses age_z, etc.
  #   scale_vars = c("age", "baseline_score")
  # )

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
    # Posterior draw extraction regex used in Step 6.
    #
    # Ordinary fixed/random/residual parameters only:
    # parameter_draw_regex = "^(b_|sd_|sigma)"
    #
    # Recommended general default, including spline and mo() parameters:
    # - b_      ordinary fixed effects
    # - bsp_    monotonic effect coefficients used by brms::mo()
    # - sd_     group-level standard deviations
    # - sigma   residual SD for Gaussian models
    # - sds_    smooth-term standard deviations
    # - bs_     smooth basis coefficients
    # - simo_   simplex parameters used by brms::mo()
    parameter_draw_regex = "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)",

    # Initial values passed to brms/Stan.
    #
    # Options:
    # init = "random"
    #   Default Stan-style random initialisation.
    #
    # init = 0
    #   Starts parameters at zero. This can help some difficult logistic
    #   models, but is not required for ordinary analyses.
    #
    # init = function() ...
    #   Advanced custom initialisation function.
    #
    # For the public examples, use "random" unless there is a clear reason
    # to use init = 0.
    init = "random",

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
    # Effects and component passed to posterior summary helpers.
    # Common choices:
    # effects = "fixed"
    # component = "conditional"
    effects = "fixed",
    component = "conditional",

    # Point summary:
    # "median" is robust and usually preferred for Bayesian summaries.
    # "mean" can also be used when appropriate.
    centrality = "median",

    # Credible interval probability and method.
    ci = 0.95,
    ci_method = "HDI",

    # Common bayestestR-style tests.
    #
    # "p_direction" gives the posterior probability of the dominant sign.
    # "rope" gives the percentage of the posterior inside a region of
    # practical equivalence, if a ROPE is defined.
    #
    # If you do not want ROPE summaries, use:
    # test = c("p_direction")
    test = c("p_direction", "rope"),

    rope = list(
      # Options:
      #
      # "none"
      #   Do not calculate or report ROPE summaries.
      #
      # "fixed"
      #   Use fixed_range exactly as supplied below.
      #   The range is on the model coefficient scale.
      #
      # "auto"
      #   Convenience option for workflows that implement family-specific
      #   ROPE defaults. For Bernoulli/logit models, width_probability = 0.05
      #   is commonly interpreted as approximately OR 0.95 to 1.05:
      #   log(c(0.95, 1.05)).
      #
      # Scientific analyses should define the ROPE from practical relevance
      # whenever possible.
      method = "fixed",

      # Used if method = "fixed".
      #
      # Logistic/logit example:
      # fixed_range = log(c(0.95, 1.05))
      #   Interpreted as an odds-ratio range from 0.95 to 1.05.
      #
      # Wider logistic/logit example:
      # fixed_range = log(c(0.90, 1.10))
      #
      # Gaussian example:
      # fixed_range = c(-1, 1)
      #   Interpreted on the outcome/model scale.
      fixed_range = log(c(0.95, 1.05)),

      # Used by method = "auto" in compatible binary/logit workflows.
      # width_probability = 0.05 means a +/- 5% odds-ratio band when
      # translated as log(c(1 - 0.05, 1 + 0.05)).
      width_probability = 0.05
    )
  ),

  parallel = list(
    # ----------------------------------------------------------
    # miceRanger imputation parallelisation
    # ----------------------------------------------------------
    # impute_workers controls how many parallel miceRanger workers
    # are used. num_impute_threads_per_worker controls the number
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
    # Keep this while supporting older pipeline versions. If all scripts
    # in your repository use impute_workers and num_impute_threads_per_worker,
    # this field can eventually be removed.
    num_impute_threads = 1,

    # ----------------------------------------------------------
    # brms model fitting parallelisation
    # ----------------------------------------------------------
    # fit_workers controls how many imputed datasets/models are fitted
    # at the same time. cores_per_fit controls how many chains/cores
    # each model can use.
    #
    # Approximate CPU demand during model fitting:
    # fit_workers * cores_per_fit
    #
    # Public example default is deliberately conservative.
    # For real analyses, increase after the smoke fit succeeds.
    fit_workers = 1,
    cores_per_fit = 1,

    # Parallel workers for Step 6 posterior draw extraction.
    # Start conservatively if brmsfit objects are large.
    summary_workers = 2,

    # Parallel workers for Step 7 posterior prediction.
    # Start conservatively if posterior prediction data are large.
    prediction_workers = 2,

    # Maximum future global size, in GB.
    # Increase for large models/data if future/furrr complains about globals.
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
