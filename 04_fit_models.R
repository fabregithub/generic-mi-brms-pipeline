# ============================================================
# 04_fit_models.R
# Prepare model data and fit brms models in parallel
#
# Design:
# - no brm_multiple()
# - one imputation per worker task
# - each worker saves its own brmsfit RDS
# - main R session receives only small status rows
# - scheduling = Inf for load balancing
# ============================================================

source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

setup_info <- setup_brms_cmdstan(paths$cache)

safe_step("STEP 4: Prepare model data and fit brms models", {
  
  # ------------------------------------------------------------
  # Load prepared objects
  # ------------------------------------------------------------
  
  dat <- readRDS(file.path(paths$objects, "dat_prepared.rds"))
  var_dict <- readRDS(file.path(paths$objects, "variable_dictionary.rds"))
  imputation_manifest <- readRDS(file.path(paths$objects, "imputation_manifest.rds"))
  
  # ------------------------------------------------------------
  # Build model spec
  # IMPORTANT:
  # build_model_spec() uses argument name reference_data, not data.
  # ------------------------------------------------------------
  
  model_spec <- build_model_spec(
    analysis_spec = analysis_spec,
    var_dict = var_dict,
    reference_data = dat
  )
  
  saveRDS(
    model_spec,
    file.path(paths$objects, "model_spec.rds"),
    compress = FALSE
  )
  
  saveRDS(
    model_spec$z_stats,
    file.path(paths$objects, "z_stats.rds"),
    compress = FALSE
  )
  
  log_msg("Model formula:", paste(deparse(model_spec$formula), collapse = " "))
  log_msg("Model family:", analysis_spec$outcome$family)
  log_msg("Model link:", analysis_spec$outcome$link)
  
  # ------------------------------------------------------------
  # Prepare model data files
  # ------------------------------------------------------------
  
  model_data_manifest <- prepare_model_data_files(
    imputation_manifest = imputation_manifest,
    analysis_spec = analysis_spec,
    model_spec = model_spec,
    out_dir = paths$model_data
  )
  
  saveRDS(
    model_data_manifest,
    file.path(paths$objects, "model_data_manifest.rds"),
    compress = FALSE
  )
  
  log_msg("Saved model_data_manifest.rds.")
  
  # ------------------------------------------------------------
  # Optional skip / only filters
  # ------------------------------------------------------------
  
  skip_imputations <- analysis_spec$model$skip_imputations %||% integer(0)
  only_imputations <- analysis_spec$model$only_imputations %||% integer(0)
  
  if (length(skip_imputations) > 0) {
    log_msg(
      "Skipping imputations:",
      paste(skip_imputations, collapse = ", ")
    )
    
    model_data_manifest <- model_data_manifest %>%
      dplyr::filter(!imputation %in% skip_imputations)
  }
  
  if (length(only_imputations) > 0) {
    log_msg(
      "Running only imputations:",
      paste(only_imputations, collapse = ", ")
    )
    
    model_data_manifest <- model_data_manifest %>%
      dplyr::filter(imputation %in% only_imputations)
  }
  
  if (nrow(model_data_manifest) == 0) {
    stop("No imputations remain after skip/only filtering.")
  }
  
  # ------------------------------------------------------------
  # Fit manifest
  # ------------------------------------------------------------
  
  fit_files <- file.path(
    paths$fits,
    sprintf("fit_imp_%03d.rds", model_data_manifest$imputation)
  )
  
  fit_manifest <- model_data_manifest %>%
    dplyr::mutate(fit_file = fit_files)
  
  saveRDS(
    fit_manifest,
    file.path(paths$objects, "fit_manifest.rds"),
    compress = FALSE
  )
  
  log_msg("Saved fit_manifest.rds with", nrow(fit_manifest), "row(s).")
  
  # ------------------------------------------------------------
  # Worker logging
  # ------------------------------------------------------------
  
  worker_log_dir <- file.path(paths$results, "worker_logs")
  dir.create(worker_log_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------
  # Parallel settings
  # ------------------------------------------------------------
  
  fit_workers <- analysis_spec$parallel$fit_workers %||% 1
  cores_per_fit <- analysis_spec$parallel$cores_per_fit %||% 1
  
  log_msg("fit_workers:", fit_workers)
  log_msg("cores_per_fit:", cores_per_fit)
  
  options(
    future.globals.maxSize =
      (analysis_spec$parallel$future_globals_maxsize_gb %||% 8) * 1024^3
  )
  
  cmdstanr::set_cmdstan_path(setup_info$cmdstan_path)
  
  options(brms.backend = "cmdstanr")
  options(cmdstanr_write_stan_file_dir = setup_info$cache_dir)
  
  # ------------------------------------------------------------
  # Optional sequential smoke fit
  # ------------------------------------------------------------
  
  run_smoke_fit <- analysis_spec$model$run_smoke_fit %||% TRUE
  
  if (isTRUE(run_smoke_fit)) {
    first_missing_idx <- which(!purrr::map_lgl(fit_manifest$fit_file, rds_ok))[1]
    
    if (!is.na(first_missing_idx)) {
      imp_smoke <- fit_manifest$imputation[first_missing_idx]
      
      log_msg(
        "Running sequential smoke fit for imputation",
        imp_smoke,
        "before parallel fitting."
      )
      
      smoke_log_file <- file.path(
        worker_log_dir,
        sprintf("fit_worker_imp_%03d_smoke.log", imp_smoke)
      )
      
      smoke_status <- fit_one_brm_file(
        ii = first_missing_idx,
        fit_manifest = fit_manifest,
        model_spec = model_spec,
        analysis_spec = analysis_spec,
        setup_info = setup_info,
        log_file = smoke_log_file
      )
      
      saveRDS(
        smoke_status,
        file.path(paths$objects, "fit_smoke_status.rds"),
        compress = FALSE
      )
      
      readr::write_csv(
        smoke_status,
        file.path(paths$results, "fit_smoke_status.csv")
      )
      
      if (!smoke_status$status[[1]] %in% c("completed", "skipped_existing_valid_fit")) {
        stop(
          "Smoke fit failed for imputation ",
          imp_smoke,
          ". Check ",
          smoke_log_file
        )
      }
      
      log_msg("Smoke fit completed successfully.")
    } else {
      log_msg("All fits already exist; smoke fit skipped.")
    }
  }
  
  # ------------------------------------------------------------
  # Parallel fitting for remaining missing fits
  # ------------------------------------------------------------
  
  fit_manifest_to_run <- fit_manifest %>%
    dplyr::mutate(
      fit_valid = purrr::map_lgl(fit_file, rds_ok)
    ) %>%
    dplyr::filter(!fit_valid) %>%
    dplyr::select(-fit_valid)
  
  if (nrow(fit_manifest_to_run) == 0) {
    log_msg("All fits already exist. No parallel fitting needed.")
    
    fit_status <- fit_manifest %>%
      dplyr::transmute(
        imputation = imputation,
        status = "skipped_existing_valid_fit",
        fit_file = fit_file,
        error = NA_character_
      )
    
  } else {
    log_msg(
      "Starting parallel fitting for",
      nrow(fit_manifest_to_run),
      "remaining imputation(s)."
    )
    
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    if (fit_workers > 1) {
      future::plan(
        future::multisession,
        workers = fit_workers
      )
    } else {
      future::plan(future::sequential)
    }
    
    fit_status_new <- furrr::future_map_dfr(
      seq_len(nrow(fit_manifest_to_run)),
      function(ii) {
        imp_i <- fit_manifest_to_run$imputation[ii]
        
        worker_log_file <- file.path(
          worker_log_dir,
          sprintf("fit_worker_imp_%03d.log", imp_i)
        )
        
        tryCatch(
          {
            status_i <- fit_one_brm_file(
              ii = ii,
              fit_manifest = fit_manifest_to_run,
              model_spec = model_spec,
              analysis_spec = analysis_spec,
              setup_info = setup_info,
              log_file = worker_log_file
            )
            
            if (!"error" %in% names(status_i)) {
              status_i$error <- NA_character_
            }
            
            status_i
          },
          error = function(e) {
            msg <- conditionMessage(e)
            
            cat(
              paste0(
                "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                "ERROR for imputation ", imp_i, ": ", msg, "\n"
              ),
              file = worker_log_file,
              append = TRUE
            )
            
            tibble::tibble(
              imputation = imp_i,
              status = "failed",
              fit_file = fit_manifest_to_run$fit_file[ii],
              error = msg
            )
          }
        )
      },
      .options = furrr::furrr_options(
        seed = TRUE,
        scheduling = Inf,
        packages = c(
          "brms",
          "cmdstanr",
          "posterior",
          "tibble",
          "dplyr",
          "readr"
        )
      )
    )
    
    existing_status <- fit_manifest %>%
      dplyr::filter(purrr::map_lgl(fit_file, rds_ok)) %>%
      dplyr::transmute(
        imputation = imputation,
        status = "skipped_existing_valid_fit",
        fit_file = fit_file,
        error = NA_character_
      )
    
    fit_status <- dplyr::bind_rows(
      existing_status,
      fit_status_new
    ) %>%
      dplyr::arrange(imputation)
  }
  
  # ------------------------------------------------------------
  # Save status
  # ------------------------------------------------------------
  
  saveRDS(
    fit_status,
    file.path(paths$objects, "fit_status.rds"),
    compress = FALSE
  )
  
  readr::write_csv(
    fit_status,
    file.path(paths$results, "fit_status.csv")
  )
  
  log_msg("Saved fit_status for", nrow(fit_status), "imputation(s).")
  
  n_failed <- sum(fit_status$status == "failed", na.rm = TRUE)
  
  if (n_failed > 0) {
    failed_imps <- fit_status %>%
      dplyr::filter(status == "failed") %>%
      dplyr::pull(imputation)
    
    log_msg(
      "WARNING:",
      n_failed,
      "fit(s) failed:",
      paste(failed_imps, collapse = ", ")
    )
    
    stop(
      "Some fits failed. Check results/fit_status.csv and results/worker_logs."
    )
  }
  
  log_msg("All requested fits completed or were already available.")
  
  rm(
    dat,
    var_dict,
    imputation_manifest,
    model_spec,
    model_data_manifest,
    fit_manifest,
    fit_manifest_to_run,
    fit_status
  )
  
  gc()
}, analysis_spec)