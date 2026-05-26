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
    # "repeated_y_subject_covariates"
    # "repeated_y_timevarying_covariates"
    data_structure = "single_time",
    
    # No repeated-measure time variable in this demo.
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
    
    link = "identity",
    
    # Only used for repeated-Y wide imputation.
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
  # This block is kept only for optional project-specific overrides.
  # Leave as NULL for standard use.
  #
  # Example override:
  # variables = list(
  #   scale_vars = c("age", "income"),
  #   auxiliary_vars = c("baseline_score")
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
    
    extra_exclude_targets = character(0)
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
    # Example:
    # custom_formula = brms::bf(Ozone ~ Solar.R_z + Wind_z + Temp_z + Month)
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
    
    parameter_draw_regex = "^(b_|sd_|sigma)",
    
    # Usually helpful for avoiding odd random initializations in demos.
    init = 0,
    
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
    effects = "fixed",
    component = "conditional",
    centrality = "median",
    ci = 0.95,
    ci_method = "HDI",
    
    # Common bayestestR tests.
    # For Gaussian demo this is okay; for real analyses choose deliberately.
    test = c("p_direction", "rope"),
    
    rope = list(
      # Options:
      # "none"
      # "fixed"
      # "auto"
      method = "fixed",
      
      # Used if method = "fixed".
      # For Gaussian Ozone scale this is just a demo placeholder.
      fixed_range = c(-1, 1),
      
      # Used by some custom binary/logit workflows.
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
  
  # ------------------------------------------------------------
  # Memory guard settings
  # ------------------------------------------------------------
  
  memory_guard = list(
    enabled = TRUE,
    
    # Public demo values are conservative.
    # For a 192 GB machine, you may use larger values.
    max_r_memory_gb = 20,
    
    min_mac_available_gb = 5,
    
    min_mac_available_before_brm_gb = 5,
    
    gc_before_check = TRUE
  )
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