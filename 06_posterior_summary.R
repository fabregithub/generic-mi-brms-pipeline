# ============================================================
# 06_posterior_summary.R
# Extract posterior parameter draws and create summaries
#
# Parallel version:
# - extracts per-imputation parameter draws in parallel
# - writes one RDS file per imputation
# - combines draw files after extraction
# - creates posterior summaries without fragile describe_posterior()
# - supports brms special parameters used by s() and mo()
# ============================================================

source("00_config.R")
source("00_common_functions.R")

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 6: Posterior parameter summaries and draws", {

  # ------------------------------------------------------------
  # Small local helpers
  # ------------------------------------------------------------

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  get_parallel_workers <- function(analysis_spec, field, fallback_field = "fit_workers") {
    workers <- analysis_spec$parallel[[field]] %||%
      analysis_spec$parallel[[fallback_field]] %||%
      1L

    workers <- as.integer(workers)

    if (is.na(workers) || workers < 1L) {
      workers <- 1L
    }

    workers
  }

  resolve_rope <- function(summary_spec) {
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
          rope_method <- "fixed"
          rope_range <- log(c(0.95, 1.05))
        } else if (identical(rope_range_cfg[[1]], "none")) {
          rope_method <- "none"
          rope_range <- NULL
        }
      }
    }

    if (!identical(rope_method, "fixed")) {
      rope_range <- NULL
    }

    rope_range
  }

  classify_parameter <- function(param) {
    dplyr::case_when(
      stringr::str_detect(param, "^b_") ~ "fixed",
      stringr::str_detect(param, "^bsp_") ~ "special_brms",
      stringr::str_detect(param, "^simo_") ~ "monotonic_simplex",
      stringr::str_detect(param, "^sds_") ~ "smooth_sd",
      stringr::str_detect(param, "^bs_") ~ "smooth_basis",
      stringr::str_detect(param, "^sd_") ~ "random_sd",
      stringr::str_detect(param, "^sigma") ~ "residual_sigma",
      TRUE ~ "other"
    )
  }

  summarise_one_parameter <- function(parameter_draws, param, ci, alpha, centrality, rope_range) {
    x <- parameter_draws[[param]]
    x <- x[is.finite(x)]

    if (length(x) == 0) {
      return(
        tibble::tibble(
          Parameter = param,
          Parameter_Class = classify_parameter(param),
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
      Parameter_Class = classify_parameter(param),
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
  # Extract parameter draws fit-by-fit, in parallel
  # ------------------------------------------------------------

  parameter_draw_files <- file.path(
    paths$results,
    sprintf("parameter_draws_imp_%03d.rds", fit_manifest$imputation)
  )

  # Include brms special parameters by default so s() and mo() models work
  # without needing users to remember bsp_/simo_/sds_/bs_.
  draw_regex <- model_spec$parameter_draw_regex %||%
    "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)"

  summary_workers <- get_parallel_workers(
    analysis_spec = analysis_spec,
    field = "summary_workers",
    fallback_field = "fit_workers"
  )

  summary_workers <- min(summary_workers, nrow(fit_manifest))

  log_msg("STEP 6 summary_workers:", summary_workers)
  log_msg("STEP 6 parameter draw regex:", draw_regex)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  old_maxsize <- getOption("future.globals.maxSize")
  on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)

  options(
    future.globals.maxSize =
      (analysis_spec$parallel$future_globals_maxsize_gb %||% 8) * 1024^3
  )

  future::plan(
    future::multisession,
    workers = summary_workers
  )

  extraction_status <- furrr::future_pmap_dfr(
    list(
      imp_i = fit_manifest$imputation,
      fit_file_i = fit_manifest$fit_file,
      draw_file_i = parameter_draw_files
    ),
    function(imp_i, fit_file_i, draw_file_i) {
      if (rds_ok(draw_file_i)) {
        return(
          tibble::tibble(
            imputation = imp_i,
            status = "skipped_existing_valid_draws",
            parameter_draw_file = draw_file_i
          )
        )
      }

      fit_i <- readRDS(fit_file_i)

      draws_i <- posterior::as_draws_df(
        fit_i,
        variable = draw_regex,
        regex = TRUE
      ) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(imputation = imp_i, .before = 1)

      saveRDS(
        draws_i,
        draw_file_i,
        compress = FALSE
      )

      rm(fit_i, draws_i)
      gc()

      tibble::tibble(
        imputation = imp_i,
        status = "completed",
        parameter_draw_file = draw_file_i
      )
    },
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c(
        "posterior",
        "tibble",
        "dplyr",
        "purrr"
      )
    )
  )

  readr::write_csv(
    extraction_status,
    file.path(paths$results, "parameter_extraction_status.csv")
  )

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

  rope_range <- resolve_rope(summary_spec)

  # This is usually much cheaper than extracting draws, but use map_dfr rather
  # than a for-loop to keep the code vectorised and readable.
  parameter_summary <- purrr::map_dfr(
    parameter_cols,
    ~ summarise_one_parameter(
      parameter_draws = parameter_draws,
      param = .x,
      ci = ci,
      alpha = alpha,
      centrality = centrality,
      rope_range = rope_range
    )
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

  log_msg("Saved parameter_summary.rds and parameter_summary.csv")

  special_parameter_summary <- parameter_summary %>%
    dplyr::filter(.data$Parameter_Class != "fixed")

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

  rm(
    fit_manifest,
    model_spec,
    parameter_draw_files,
    parameter_manifest,
    valid_draw_files,
    parameter_draws,
    parameter_summary
  )

  gc()

  guard_memory("after STEP 6 cleanup")
}, analysis_spec)
