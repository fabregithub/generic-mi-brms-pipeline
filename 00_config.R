# ============================================================
# 00_config.R
# Generic MI + brms pipeline configuration
#
# Public example dataset:
# datasets::airquality
#
# Outcome:
# Ozone, Gaussian identity model
#
# Main model:
# Ozone ~ Solar.R_z + Wind_z + Temp_z + Month
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

# ============================================================
# Analysis specification
# ============================================================

analysis_spec <- list(

  # ------------------------------------------------------------
  # Basic analysis identity
  # ------------------------------------------------------------

  analysis_id = "airquality_gaussian_demo",

  project_label = "Public airquality Gaussian MI + brms demo",

  # ------------------------------------------------------------
  # Data specification
  # ------------------------------------------------------------

  data = list(
    raw_data_file = "data/airquality_example.rds",

    # airquality has no natural subject ID.
    # The data-creation script adds row_id.
    id_var = NULL,

    row_id_var = "row_id",

    # Options:
    # "single_time"
    #   One row per analytic unit.
    #
    # "repeated_y_subject_covariates"
    #   Long repeated-outcome data where most covariates are measured once
    #   per subject or at baseline.
    #
    # "repeated_y_timevarying_covariates"
    #   Long repeated-outcome data with time-varying predictors.
    data_structure = "single_time",

    # Repeated-measure time variable.
    # Use NULL for single-time analyses.
    # Example for repeated data:
    # time_var = "time"
    time_var = NULL
  ),

  # ------------------------------------------------------------
  # Outcome specification
  # ------------------------------------------------------------

  outcome = list(
    y_var = "Ozone",

    # Options supported by the template:
    # "gaussian", "bernoulli", "poisson", "negbinomial",
    # "beta", "ordinal", "categorical"
    family = "gaussian",

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
    link = "identity",

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
    predict_missing_y = TRUE
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
  #
  # Public airquality example, if written explicitly:
  #
  # variables = list(
  #   exposure_vars = c("Solar.R", "Wind", "Temp", "Month"),
  #   covariate_vars = character(0),
  #   auxiliary_vars = character(0),
  #   continuous_vars = c("Solar.R", "Wind", "Temp"),
  #   categorical_vars = c("Month"),
  #   ordinal_vars = character(0),
  #   subject_level_vars = c("Solar.R", "Wind", "Temp", "Month"),
  #   time_varying_vars = character(0),
  #   scale_vars = c("Solar.R", "Wind", "Temp")
  # )

  variables = NULL,

  # ------------------------------------------------------------
  # Imputation specification
  # ------------------------------------------------------------

  imputation = list(
    enabled = TRUE,

    # Options:
    # "none"
    # "row_level"
    # "subject_level"
    # "subject_wide_with_repeated_y_auxiliary"
    # "long_row_level"
    strategy = "row_level",

    # Public demo uses small m for speed.
    # For real analyses, use something like 50-100.
    m = 100,

    maxiter = 5,

    mean_match_k = 5,

    verbose = FALSE,

    # For this demo, do not impute the outcome as a target.
    # Missing Ozone rows will be predicted later using posterior_predict().
    impute_y = FALSE,

    extra_exclude_targets = character(0),

    # ----------------------------------------------------------
    # Optional: reproducibility and safe extension
    # ----------------------------------------------------------
    # Base seed for miceRanger. Defaults to analysis_spec$model$seed if left
    # NULL. Each batch of new imputations (the initial m, or a later
    # extension batch) is seeded deterministically as seed + (number of
    # imputations that already existed before that batch).
    seed = NULL,

    # FALSE by default. Step 3 normally refuses to fit more imputations than
    # already exist on disk. Set to TRUE to let Step 3 generate only the
    # additional imputations needed to reach the current m, without
    # touching any existing imputed-data file. Used automatically by the
    # mi_stability$auto_increment loop below.
    allow_extend = FALSE
  ),

  # ------------------------------------------------------------
  # Model specification
  # ------------------------------------------------------------

  model = list(
    # If "auto", fixed effects are built from variable dictionary/model roles.
    # In this demo, the final fixed effects should be:
    # Solar.R_z + Wind_z + Temp_z + Month
    fixed_effects = "auto",

    random_effects = list(
      subject_intercept = FALSE,
      subject_slope_vars = character(0)
    ),

    # Optional manual formula override.
    #
    # Examples:
    #
    # Standard Gaussian/logistic formula:
    # custom_formula = brms::bf(y ~ age_z + sex + exposure)
    #
    # Random intercept:
    # custom_formula = brms::bf(y ~ age_z + sex + exposure + (1 | id))
    #
    # Spline:
    # custom_formula = brms::bf(y ~ s(age_z, k = 5) + sex + exposure)
    #
    # Monotonic ordinal term:
    # custom_formula = brms::bf(y ~ mo(income) + age_z + sex)
    custom_formula = NULL,

    # Options:
    # "default_weakly_regularizing"
    # or a brms prior vector
    priors = "default_weakly_regularizing",

    # Public demo settings.
    # For real analyses, use e.g. chains = 4, iter = 2000, warmup = 1000.
    chains = 4,
    iter = 2000,
    warmup = 1000,

    seed = 12345,

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
    # For the public template, use "random" unless there is a clear reason
    # to use init = 0.
    init = "random",

    # 0 = show Stan output.
    # 2 = quieter.
    # For debugging the template, keep this at 0.
    silent = 0,

    # Run one sequential fit before parallel fitting.
    # This helps catch brms/CmdStan/config problems early.
    run_smoke_fit = TRUE,

    # For production analyses, you can temporarily skip problematic imputations.
    # Example: skip_imputations = c(45, 51)
    skip_imputations = integer(0),

    # To run only specific imputations.
    # Example: only_imputations = c(51)
    only_imputations = integer(0)
  ),

  # ------------------------------------------------------------
  # Posterior prediction specification
  # ------------------------------------------------------------

  posterior_prediction = list(
    enabled = TRUE,

    # Small for demo speed.
    # For real analyses, use 1000 or more.
    ndraws = 200,

    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  ),

  # ------------------------------------------------------------
  # Posterior summary specification
  # ------------------------------------------------------------

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
      # Gaussian example:
      # fixed_range = c(-1, 1)
      #   Interpreted on the outcome/model scale.
      #
      # Logistic/logit example:
      # fixed_range = log(c(0.95, 1.05))
      #   Interpreted as an odds-ratio range from 0.95 to 1.05.
      #
      # Wider logistic/logit example:
      # fixed_range = log(c(0.90, 1.10))
      fixed_range = c(-1, 1),

      # Used by method = "auto" in compatible binary/logit workflows.
      # width_probability = 0.05 means a +/- 5% odds-ratio band when
      # translated as log(c(1 - 0.05, 1 + 0.05)).
      width_probability = 0.05
    )
  ),

  # ------------------------------------------------------------
  # Parallel and memory settings
  # ------------------------------------------------------------

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

  # ------------------------------------------------------------
  # Memory guard settings
  # ------------------------------------------------------------

  memory_guard = list(
    enabled = TRUE,

    # Public demo values are conservative.
    # For a high-memory workstation, you may use larger values.
    max_r_memory_gb = 20,

    # macOS-specific available-memory checks. On non-macOS systems these
    # checks are usually skipped by the helper functions.
    min_mac_available_gb = 5,

    min_mac_available_before_brm_gb = 5,

    gc_before_check = TRUE
  ),

  # ------------------------------------------------------------
  # Optional: automatic imputation-count stability loop
  # ------------------------------------------------------------
  # Placeholder, inactive by default (auto_increment = FALSE). With this
  # left FALSE, run_all.R fits a single fixed m, exactly as if this block
  # were absent. Set auto_increment <- TRUE to let run_all.R fit imputations
  # in batches and stop increasing m automatically once posterior summaries
  # are stable, instead of fitting analysis_spec$imputation$m up front.
  # analysis_spec$imputation$m is then used as the ceiling, not the
  # starting point. See README.md, "Choosing the number of imputations
  # adaptively".
  mi_stability = list(
    auto_increment = FALSE,

    # NULL defaults to analysis_spec$parallel$fit_workers, rounded up to a
    # multiple of fit_workers if you override it here.
    increment_size = NULL,

    # The remaining settings are shared with the manual
    # 11_check_imputation_stability.R script and only matter if
    # auto_increment = TRUE or you run that script directly.
    parameter_regex = "^b_",
    exclude_intercept = TRUE,
    estimate_tolerance = 0.05,
    ci_endpoint_tolerance = 0.05,
    relative_transformed_tolerance_pct = 5,
    pd_tolerance = 0.02
  ),

  # ------------------------------------------------------------
  # Optional: monotonic-effect (mo()) labels
  # ------------------------------------------------------------
  # Placeholder, inactive unless your model formula contains mo(). If your
  # formula has no mo() term, run_all.R skips 09/10 automatically and this
  # block is simply ignored. If it does contain mo() and this is left NULL,
  # 10_publication_mo_results.R still runs, using generic "Level 1",
  # "Level 2", ... category labels. Run 09_check_mo_parameter_columns.R for
  # a ready-to-paste vars = list(...) block with the right number of levels
  # for each detected mo() variable, then add real labels/levels below.
  mo_effects = NULL

  # Example, once you have run 09_check_mo_parameter_columns.R:
  # mo_effects = list(
  #   vars = list(
  #     your_ordinal_var = list(
  #       label = "Human-readable label",
  #       levels = c("Level 1", "Level 2", "Level 3")
  #     )
  #   ),
  #   time_var = NULL,    # set only if your formula has time * mo(variable)
  #   time_values = NULL  # e.g. 1:6
  # )
)

# ============================================================
# Project paths
# ============================================================

paths <- list(
  root = ".",

  # Folders
  data = "data",
  objects = "objects",

  # Input files
  raw_data = analysis_spec$data$raw_data_file,
  variable_dictionary = "00_variable_dictionary.csv",

  # Imputation outputs
  imputed_data = "objects/imputed_data",
  imputed = "objects/imputed_data",          # backward-compatible alias
  imputed_wide = "objects/imputed_wide",

  # Model data and fits
  model_data = "objects/model_data",
  fits = "fits",

  # Results
  results = "results",
  publication = "results/publication",

  # CmdStan cache
  cache = file.path(path.expand("~"), ".cmdstanr-cache")
)

# ============================================================
# Convenience settings
# ============================================================

options(
  brms.backend = "cmdstanr"
)

# Avoid accidental scientific notation in output tables.
options(
  scipen = 999
)
