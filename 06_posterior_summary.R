# ============================================================
# 06_posterior_summary.R
# Extract posterior parameter draws and create summaries
#
# Robust version:
# - reads each saved brms fit one by one
# - extracts parameter draws
# - saves per-imputation draw files
# - combines draw files
# - creates posterior summaries manually
# - avoids fragile describe_posterior(data.frame) method dispatch
# ============================================================

source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 6: Posterior parameter summaries and draws", {
  
  # ------------------------------------------------------------
  # Load manifests/specs
  # ------------------------------------------------------------
  
  fit_manifest <- readRDS(file.path(paths$objects, "fit_manifest.rds"))
  model_spec <- readRDS(file.path(paths$objects, "model_spec.rds"))
  
  # Keep only valid completed fits
  fit_manifest <- fit_manifest %>%
    dplyr::mutate(
      fit_valid = purrr::map_lgl(fit_file, rds_ok)
    ) %>%
    dplyr::filter(fit_valid) %>%
    dplyr::select(-fit_valid)
  
  if (nrow(fit_manifest) == 0) {
    stop("No valid brms fit files found.")
  }
  
  log_msg("Using", nrow(fit_manifest), "valid fit(s) for posterior summaries.")
  
  guard_memory("after reading inputs for STEP 6")
  
  # ------------------------------------------------------------
  # Extract parameter draws fit-by-fit
  # ------------------------------------------------------------
  
  parameter_draw_files <- file.path(
    paths$results,
    sprintf("parameter_draws_imp_%03d.rds", fit_manifest$imputation)
  )
  
  draw_regex <- model_spec$parameter_draw_regex %||% "^(b_|sd_|sigma|sds_|bs_|simo_|bsp_)"
  
  for (ii in seq_len(nrow(fit_manifest))) {
    imp_i <- fit_manifest$imputation[ii]
    
    update_heartbeat(
      paste(
        "STEP 6 extracting parameter draws imputation",
        imp_i,
        "of",
        nrow(fit_manifest)
      )
    )
    
    if (rds_ok(parameter_draw_files[ii])) {
      log_msg("Existing parameter draws found; skipping imputation", imp_i)
      next
    }
    
    log_msg("Reading fit for parameter draws, imputation", imp_i)
    
    fit_i <- readRDS(fit_manifest$fit_file[ii])
    
    draws_i <- posterior::as_draws_df(
      fit_i,
      variable = draw_regex,
      regex = TRUE
    ) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(imputation = imp_i, .before = 1)
    
    saveRDS(
      draws_i,
      parameter_draw_files[ii],
      compress = FALSE
    )
    
    log_msg("Saved parameter draws for imputation", imp_i)
    
    rm(fit_i, draws_i)
    gc()
    
    guard_memory(paste("after extracting parameter draws for imputation", imp_i))
  }
  
  parameter_manifest <- fit_manifest %>%
    dplyr::mutate(parameter_draw_file = parameter_draw_files)
  
  saveRDS(
    parameter_manifest,
    file.path(paths$objects, "parameter_manifest.rds"),
    compress = FALSE
  )
  
  # ------------------------------------------------------------
  # Combine parameter draw files
  # ------------------------------------------------------------
  
  valid_draw_files <- parameter_draw_files[
    purrr::map_lgl(parameter_draw_files, rds_ok)
  ]
  
  if (length(valid_draw_files) == 0) {
    stop("No valid parameter draw files found.")
  }
  
  log_msg("Combining", length(valid_draw_files), "parameter draw file(s).")
  
  parameter_draws <- purrr::map_dfr(
    valid_draw_files,
    readRDS
  )
  
  saveRDS(
    parameter_draws,
    file.path(paths$results, "parameter_draws.rds"),
    compress = FALSE
  )
  
  log_msg("Saved combined parameter_draws.rds")
  
  # ------------------------------------------------------------
  # Manual posterior summary
  # ------------------------------------------------------------
  
  summary_spec <- analysis_spec$summary
  
  ci <- summary_spec$ci %||% 0.95
  alpha <- (1 - ci) / 2
  
  centrality <- summary_spec$centrality %||% "median"
  
  meta_cols <- c(
    "imputation",
    ".chain",
    ".iteration",
    ".draw"
  )
  
  parameter_cols <- setdiff(names(parameter_draws), meta_cols)
  
  # Keep only numeric parameter columns
  parameter_cols <- parameter_cols[
    purrr::map_lgl(parameter_draws[parameter_cols], is.numeric)
  ]
  
  if (length(parameter_cols) == 0) {
    stop("No numeric parameter columns found in parameter_draws.")
  }
  
  log_msg("Summarising", length(parameter_cols), "parameter(s).")
  
  # ------------------------------------------------------------
  # ROPE handling
  # ------------------------------------------------------------
  # Support both old and new config styles:
  #
  # Old/list style:
  #   summary$rope <- list(method = "fixed", fixed_range = c(-0.1, 0.1))
  #
  # Current/simple style:
  #   summary$test <- c("p_direction", "rope")
  #   summary$rope_range <- "auto_logit_5pct"
  #
  # Also tolerate summary$rope being an atomic vector/string.

  rope_method <- "none"
  rope_range <- NULL

  rope_cfg <- summary_spec$rope %||% NULL

  if (is.list(rope_cfg)) {
    rope_method <- rope_cfg$method %||% "none"

    if (identical(rope_method, "fixed")) {
      rope_range <- rope_cfg$fixed_range %||% NULL
    }
  } else if (is.character(rope_cfg) && length(rope_cfg) >= 1) {
    rope_method <- rope_cfg[[1]]
  }

  rope_range_cfg <- summary_spec$rope_range %||% NULL

  if (!is.null(rope_range_cfg)) {
    if (is.numeric(rope_range_cfg) && length(rope_range_cfg) == 2) {
      rope_method <- "fixed"
      rope_range <- as.numeric(rope_range_cfg)
    } else if (is.character(rope_range_cfg) && length(rope_range_cfg) >= 1) {
      if (identical(rope_range_cfg[[1]], "auto_logit_5pct")) {
        # On the log-odds scale, this corresponds approximately to an
        # odds-ratio equivalence region of 0.95 to 1.05.
        rope_method <- "fixed"
        rope_range <- log(c(0.95, 1.05))
      } else if (identical(rope_range_cfg[[1]], "none")) {
        rope_method <- "none"
        rope_range <- NULL
      } else {
        warning(
          "Unknown summary$rope_range value: ",
          rope_range_cfg[[1]],
          ". ROPE percentages will be set to NA."
        )
      }
    }
  }

  if (!identical(rope_method, "fixed")) {
    rope_range <- NULL
  }

  classify_parameter <- function(param) {
    dplyr::case_when(
      param == "b_Intercept" ~ "intercept",
      grepl("^b_", param) ~ "fixed",
      grepl("^sd_", param) ~ "group_sd",
      grepl("^sigma", param) ~ "residual_sigma",
      grepl("^sds_", param) ~ "smooth_sd",
      grepl("^bs_", param) ~ "smooth_basis",
      grepl("^simo_", param) ~ "monotonic_simplex",
      grepl("^bsp_", param) ~ "special_brms",
      TRUE ~ "other"
    )
  }
  
  parameter_summary <- purrr::map_dfr(
    parameter_cols,
    function(param) {
      x <- parameter_draws[[param]]
      x <- x[is.finite(x)]
      
      if (length(x) == 0) {
        return(
          tibble::tibble(
            Parameter = param,
            Median = NA_real_,
            Mean = NA_real_,
            SD = NA_real_,
            CI = ci,
            CI_low = NA_real_,
            CI_high = NA_real_,
            pd = NA_real_,
            ROPE_Percentage = NA_real_,
            n_draws = 0L
          )
        )
      }
      
      central_value <- if (identical(centrality, "mean")) {
        mean(x)
      } else {
        stats::median(x)
      }
      
      ci_low <- as.numeric(stats::quantile(x, probs = alpha, names = FALSE))
      ci_high <- as.numeric(stats::quantile(x, probs = 1 - alpha, names = FALSE))
      
      pd <- max(
        mean(x > 0, na.rm = TRUE),
        mean(x < 0, na.rm = TRUE)
      )
      
      rope_percentage <- NA_real_
      
      if (!is.null(rope_range) && length(rope_range) == 2) {
        rope_percentage <- mean(
          x >= rope_range[1] & x <= rope_range[2],
          na.rm = TRUE
        ) * 100
      }
      
      tibble::tibble(
        Parameter = param,
        Median = stats::median(x),
        Mean = mean(x),
        SD = stats::sd(x),
        CI = ci,
        CI_low = ci_low,
        CI_high = ci_high,
        pd = pd,
        ROPE_Percentage = rope_percentage,
        n_draws = length(x)
      )
    }
  )
  
  parameter_summary <- parameter_summary %>%
    dplyr::mutate(
      Parameter_Class = classify_parameter(Parameter),
      .after = "Parameter"
    )
  
  saveRDS(
    parameter_summary,
    file.path(paths$results, "parameter_summary.rds"),
    compress = FALSE
  )
  
  readr::write_csv(
    parameter_summary,
    file.path(paths$results, "parameter_summary.csv")
  )
  
  special_parameter_summary <- parameter_summary %>%
    dplyr::filter(Parameter_Class %in% c(
      "smooth_sd",
      "smooth_basis",
      "monotonic_simplex",
      "special_brms"
    ))
  
  if (nrow(special_parameter_summary) > 0) {
    saveRDS(
      special_parameter_summary,
      file.path(paths$results, "special_parameter_summary.rds"),
      compress = FALSE
    )
    
    readr::write_csv(
      special_parameter_summary,
      file.path(paths$results, "special_parameter_summary.csv")
    )
    
    log_msg("Saved special_parameter_summary.rds and special_parameter_summary.csv")
  }
  
  log_msg("Saved parameter_summary.rds and parameter_summary.csv")
  
  rm(
    fit_manifest,
    model_spec,
    parameter_draw_files,
    parameter_manifest,
    valid_draw_files,
    parameter_draws,
    parameter_summary,
    special_parameter_summary
  )
  
  gc()
  
  guard_memory("after STEP 6 cleanup")
}, analysis_spec)