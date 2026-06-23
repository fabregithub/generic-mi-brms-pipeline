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

  # ------------------------------------------------------------
  # Weighted-statistics helpers for MI pooling.
  #
  # Each draw from imputation i carries weight 1/(m * K_i), so every
  # imputation contributes exactly 1/m regardless of how many finite
  # draws it happened to produce (fixes unequal-K_i weighting).
  # ------------------------------------------------------------

  weighted_mean_ <- function(x, w) sum(x * w) / sum(w)

  weighted_var_ <- function(x, w) {
    mu <- weighted_mean_(x, w)
    sum(w * (x - mu)^2) / sum(w)
  }

  weighted_quantile_ <- function(x, w, probs) {
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cw <- (cumsum(w) - 0.5 * w) / sum(w)
    as.numeric(stats::approx(cw, x, xout = probs, rule = 2)$y)
  }

  # Sarle's bimodality coefficient (no extra package dependency).
  # BC > ~0.555 (uniform-distribution threshold) flags likely
  # bimodal/multimodal shape, in which case we must NOT apply a
  # symmetric variance-inflation correction (it would smear distinct
  # modes instead of preserving genuine between-imputation structure).
  bimodality_coefficient_ <- function(x, w) {
    mu <- weighted_mean_(x, w)
    s2 <- weighted_var_(x, w)
    if (!is.finite(s2) || s2 <= 0) {
      return(NA_real_)
    }
    skew <- weighted_mean_((x - mu)^3, w) / s2^1.5
    kurt <- weighted_mean_((x - mu)^4, w) / s2^2
    (skew^2 + 1) / kurt
  }

  # Choose a support-respecting transform so the correction (which
  # assumes rough symmetry) is applied on a scale where that
  # assumption is defensible.
  choose_transform_ <- function(x) {
    if (all(x > 0)) {
      list(
        name = "log",
        fwd = log,
        inv = exp
      )
    } else if (all(x > 0 & x < 1)) {
      list(
        name = "logit",
        fwd = function(v) log(v / (1 - v)),
        inv = function(v) 1 / (1 + exp(-v))
      )
    } else {
      list(
        name = "identity",
        fwd = identity,
        inv = identity
      )
    }
  }

  pool_one_parameter <- function(parameter_draws, param, ci, alpha, rope_range, m_total) {
    df <- tibble::tibble(
      imputation = parameter_draws$imputation,
      value = parameter_draws[[param]]
    ) %>%
      dplyr::filter(is.finite(.data$value))

    if (nrow(df) == 0) {
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
          n_draws = 0L,
          m_imputations = 0L,
          between_var = NA_real_,
          within_var = NA_real_,
          variance_corrected = FALSE,
          transform_used = NA_character_,
          bimodality_coef = NA_real_
        )
      )
    }

    # Per-imputation statistics (Q_i, U_i, K_i) drive both the
    # weighting and the Rubin's-rule finite-m correction.
    per_imp <- df %>%
      dplyr::group_by(.data$imputation) %>%
      dplyr::summarise(
        Q_i = mean(.data$value),
        U_i = stats::var(.data$value),
        K_i = dplyr::n(),
        .groups = "drop"
      )

    m <- nrow(per_imp)

    if (m < m_total) {
      log_msg(
        "Note:", param, "has draws from only", m, "of", m_total,
        "imputations (e.g. a special-term parameter not present in every fit)."
      )
    }

    df <- df %>%
      dplyr::left_join(per_imp %>% dplyr::select(imputation, K_i), by = "imputation") %>%
      dplyr::mutate(weight = 1 / (m * .data$K_i))

    Qbar <- sum(per_imp$Q_i) / m
    Ubar <- mean(per_imp$U_i, na.rm = TRUE)
    B <- if (m > 1) stats::var(per_imp$Q_i) else 0
    V_mix <- Ubar + B
    T_var <- V_mix + B / m

    bc <- bimodality_coefficient_(df$value, df$weight)
    unimodal <- is.finite(bc) && bc <= 0.555

    scale_factor <- if (is.finite(V_mix) && V_mix > 0) sqrt(T_var / V_mix) else NA_real_

    apply_correction <- unimodal && is.finite(scale_factor) && m > 1

    transform_used <- "none"
    corrected_value <- df$value

    if (apply_correction) {
      tr <- choose_transform_(df$value)
      y <- tr$fwd(df$value)

      if (all(is.finite(y))) {
        ybar <- weighted_mean_(y, df$weight)
        y_corrected <- ybar + scale_factor * (y - ybar)
        corrected_value <- tr$inv(y_corrected)
        transform_used <- tr$name
      } else {
        apply_correction <- FALSE
      }
    }

    w <- df$weight
    x <- corrected_value

    ci_low <- weighted_quantile_(x, w, alpha)
    ci_high <- weighted_quantile_(x, w, 1 - alpha)
    med <- weighted_quantile_(x, w, 0.5)

    pd <- max(
      sum(w[x > 0]) / sum(w),
      sum(w[x < 0]) / sum(w)
    )

    rope_percentage <- NA_real_

    if (!is.null(rope_range) && length(rope_range) == 2) {
      in_rope <- x >= rope_range[1] & x <= rope_range[2]
      rope_percentage <- sum(w[in_rope]) / sum(w) * 100
    }

    tibble::tibble(
      Parameter = param,
      Parameter_Class = classify_parameter(param),
      Median = med,
      Mean = weighted_mean_(x, w),
      SD = sqrt(weighted_var_(x, w)),
      CI = ci,
      CI_low = ci_low,
      CI_high = ci_high,
      pd = pd,
      ROPE_Percentage = rope_percentage,
      n_draws = nrow(df),
      m_imputations = m,
      between_var = B,
      within_var = Ubar,
      variance_corrected = apply_correction,
      transform_used = transform_used,
      bimodality_coef = bc
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

  m_total <- dplyr::n_distinct(parameter_draws$imputation)

  # Proper MI pooling per parameter: weight draws by 1/(m*K_i) so every
  # imputation contributes equally regardless of how many finite draws
  # it produced, then apply the Rubin's-rule finite-m variance
  # correction (B/m) only where the pooled shape is unimodal enough for
  # the correction's symmetry assumption to be safe.
  parameter_summary <- purrr::map_dfr(
    parameter_cols,
    ~ pool_one_parameter(
      parameter_draws = parameter_draws,
      param = .x,
      ci = ci,
      alpha = alpha,
      rope_range = rope_range,
      m_total = m_total
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
