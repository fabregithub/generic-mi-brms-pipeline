# ============================================================
# 11_check_imputation_stability.R
# Imputation-count stability checks for MI + brms analyses
#
# Run after Step 6 has created per-imputation parameter draw files:
#   Rscript 11_check_imputation_stability.R
# ============================================================

source("00_config.R")
source("00_common_functions.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(gt)
  library(flextable)
  library(officer)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!exists("rds_ok", mode = "function")) {
  rds_ok <- function(path) {
    if (length(path) != 1 || is.na(path) || !file.exists(path)) return(FALSE)
    tryCatch({ readRDS(path); TRUE }, error = function(e) FALSE)
  }
}
if (!exists("init_logging", mode = "function")) init_logging <- function(...) invisible(NULL)
if (!exists("log_msg", mode = "function")) log_msg <- function(...) message(paste(..., collapse = " "))
if (!exists("setup_project_dirs", mode = "function")) {
  setup_project_dirs <- function(paths) {
    dir.create(paths$objects, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results, recursive = TRUE, showWarnings = FALSE)
  }
}
if (!exists("safe_step", mode = "function")) {
  safe_step <- function(step_name, expr, analysis_spec = NULL) {
    log_msg(step_name)
    force(expr)
  }
}

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 11: Imputation-count stability checks", {

  cfg <- analysis_spec$mi_stability %||% list()
  set.seed(cfg$seed %||% 12345)

  parameter_regex <- cfg$parameter_regex %||% "^b_"
  primary_parameters <- cfg$primary_parameters %||% NULL
  exclude_intercept <- cfg$exclude_intercept %||% TRUE

  ci <- analysis_spec$summary$ci %||% 0.95
  alpha <- (1 - ci) / 2

  estimate_tolerance <- cfg$estimate_tolerance %||% 0.05
  ci_endpoint_tolerance <- cfg$ci_endpoint_tolerance %||% 0.05
  relative_transformed_tolerance_pct <- cfg$relative_transformed_tolerance_pct %||% 5
  pd_tolerance <- cfg$pd_tolerance %||% 0.02
  max_plot_parameters <- cfg$max_plot_parameters %||% 12
  # Defaults to FALSE: the tables/figures this script produces are also
  # embedded directly into the main report's "Imputation-count stability"
  # chapter (see 08_publication_results.R), which run_all.R renders once
  # everything has finished. Set this TRUE if you want this script's own
  # standalone report rendered in addition to the combined main report,
  # e.g. to inspect every evaluated batch in more detail.
  render_quarto <- cfg$render_quarto %||% FALSE

  stability_dir <- file.path(paths$results, "publication", "mi_stability")
  table_dir <- file.path(stability_dir, "tables")
  figure_dir <- file.path(stability_dir, "figures")
  report_dir <- file.path(stability_dir, "report")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  fmt_num <- function(x, digits = 3) {
    ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
  }
  fmt_ci <- function(low, high, digits = 3) {
    paste0("[", fmt_num(low, digits), ", ", fmt_num(high, digits), "]")
  }
  clean_parameter_name <- function(x) {
    x %>%
      stringr::str_replace("^b_", "") %>%
      stringr::str_replace("^sd_", "SD: ") %>%
      stringr::str_replace("^sds_", "Smooth SD: ") %>%
      stringr::str_replace("^bs_", "Smooth basis: ") %>%
      stringr::str_replace("^bsp_", "Monotonic effect: ") %>%
      stringr::str_replace("^simo_", "Monotonic simplex: ") %>%
      stringr::str_replace_all(":", " x ") %>%
      stringr::str_replace_all("_z$", " per SD") %>%
      stringr::str_replace_all("_", " ")
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
  is_intercept <- function(x) x %in% c("b_Intercept", "Intercept", "(Intercept)")
  effect_transform_type <- function(analysis_spec) {
    fam <- analysis_spec$outcome$family
    link <- analysis_spec$outcome$link
    if (fam %in% c("bernoulli", "ordinal", "categorical") && link == "logit") return("odds_ratio")
    if (fam %in% c("poisson", "negbinomial") && link == "log") return("rate_ratio")
    if (fam == "gaussian" && link == "log") return("multiplicative_ratio")
    "none"
  }
  transform_label <- function(transform_type) {
    switch(transform_type, odds_ratio = "Odds ratio", rate_ratio = "Rate ratio", multiplicative_ratio = "Ratio", "Transformed estimate")
  }
  resolve_rope <- function(summary_spec) {
    rope_method <- "none"
    rope_range <- NULL
    rope_cfg <- summary_spec$rope %||% NULL
    if (is.list(rope_cfg)) {
      rope_method <- rope_cfg$method %||% "none"
      if (identical(rope_method, "fixed")) rope_range <- rope_cfg$fixed_range %||% NULL
      if (identical(rope_method, "auto")) {
        width <- rope_cfg$width_probability %||% 0.05
        fam <- analysis_spec$outcome$family
        link <- analysis_spec$outcome$link
        if (fam %in% c("bernoulli", "ordinal", "categorical") && link == "logit") {
          rope_method <- "fixed"
          rope_range <- log(c(1 - width, 1 + width))
        } else {
          rope_method <- "none"
          rope_range <- NULL
        }
      }
    } else if (is.character(rope_cfg) && length(rope_cfg) >= 1) {
      rope_method <- rope_cfg[[1]]
    }
    rope_range_cfg <- summary_spec$rope_range %||% NULL
    if (!is.null(rope_range_cfg)) {
      if (is.numeric(rope_range_cfg) && length(rope_range_cfg) == 2) {
        rope_method <- "fixed"; rope_range <- as.numeric(rope_range_cfg)
      } else if (is.character(rope_range_cfg) && length(rope_range_cfg) >= 1) {
        if (identical(rope_range_cfg[[1]], "auto_logit_5pct")) {
          rope_method <- "fixed"; rope_range <- log(c(0.95, 1.05))
        } else if (identical(rope_range_cfg[[1]], "none")) {
          rope_method <- "none"; rope_range <- NULL
        }
      }
    }
    if (!identical(rope_method, "fixed")) rope_range <- NULL
    rope_range
  }
  find_parameter_draw_files <- function(paths) {
    parameter_manifest_file <- file.path(paths$objects, "parameter_manifest.rds")
    if (rds_ok(parameter_manifest_file)) {
      pm <- readRDS(parameter_manifest_file)
      file_col <- intersect(c("parameter_draw_file", "draw_file", "file", "parameter_file"), names(pm))[1]
      if (!is.na(file_col)) {
        return(pm %>% dplyr::transmute(imputation = as.integer(.data$imputation), parameter_draw_file = as.character(.data[[file_col]])) %>% dplyr::filter(!is.na(.data$imputation), !is.na(.data$parameter_draw_file)))
      }
    }
    files <- list.files(paths$results, pattern = "^parameter_draws_imp_[0-9]+\\.rds$", full.names = TRUE)
    if (length(files) == 0) return(tibble::tibble(imputation = integer(), parameter_draw_file = character()))
    tibble::tibble(
      imputation = as.integer(stringr::str_match(basename(files), "^parameter_draws_imp_([0-9]+)\\.rds$")[, 2]),
      parameter_draw_file = files
    )
  }
  choose_batch_sizes <- function(n_imputations, cfg) {
    user_batches <- cfg$batch_sizes %||% NULL
    if (!is.null(user_batches)) {
      batches <- as.integer(user_batches)
      batches <- batches[!is.na(batches)]
      batches <- batches[batches >= 2 & batches <= n_imputations]
      return(sort(unique(c(batches, n_imputations))))
    }
    if (n_imputations <= 8) {
      batches <- seq(2, n_imputations)
    } else if (n_imputations <= 24) {
      batches <- sort(unique(c(seq(4, n_imputations, by = 4), n_imputations)))
    } else {
      batches <- c(8, 12, 16, 20, 24, 32, 40, 48, 60, 80, 100, 150, n_imputations)
      batches <- sort(unique(batches[batches >= 2 & batches <= n_imputations]))
    }
    batches
  }
  select_parameter_columns <- function(draws, parameter_regex, primary_parameters, exclude_intercept) {
    meta_cols <- c("imputation", ".chain", ".iteration", ".draw")
    parameter_cols <- setdiff(names(draws), meta_cols)
    parameter_cols <- parameter_cols[vapply(draws[parameter_cols], is.numeric, logical(1))]
    if (!is.null(primary_parameters) && length(primary_parameters) > 0) {
      parameter_cols <- intersect(parameter_cols, primary_parameters)
    } else if (!is.null(parameter_regex) && nzchar(parameter_regex)) {
      parameter_cols <- parameter_cols[stringr::str_detect(parameter_cols, parameter_regex)]
    }
    if (isTRUE(exclude_intercept)) parameter_cols <- parameter_cols[!is_intercept(parameter_cols)]
    parameter_cols
  }
  summarise_one_parameter <- function(x, param, batch_n, batch_imputations, ci, alpha, rope_range) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(tibble::tibble(n_imputations = batch_n, imputation_min = min(batch_imputations), imputation_max = max(batch_imputations), Parameter = param, Parameter_clean = clean_parameter_name(param), Parameter_Class = classify_parameter(param), Median = NA_real_, Mean = NA_real_, SD = NA_real_, CI = ci, CI_low = NA_real_, CI_high = NA_real_, pd = NA_real_, ROPE_Percentage = NA_real_, n_draws = 0L))
    }
    ci_low <- as.numeric(stats::quantile(x, probs = alpha, names = FALSE))
    ci_high <- as.numeric(stats::quantile(x, probs = 1 - alpha, names = FALSE))
    pd <- max(mean(x > 0, na.rm = TRUE), mean(x < 0, na.rm = TRUE))
    rope_percentage <- NA_real_
    if (!is.null(rope_range) && length(rope_range) == 2) {
      rope_percentage <- mean(x >= rope_range[1] & x <= rope_range[2], na.rm = TRUE) * 100
    }
    tibble::tibble(n_imputations = batch_n, imputation_min = min(batch_imputations), imputation_max = max(batch_imputations), Parameter = param, Parameter_clean = clean_parameter_name(param), Parameter_Class = classify_parameter(param), Median = stats::median(x), Mean = mean(x), SD = stats::sd(x), CI = ci, CI_low = ci_low, CI_high = ci_high, pd = pd, ROPE_Percentage = rope_percentage, n_draws = length(x))
  }
  summarise_batch <- function(batch_n, draw_file_manifest, parameter_cols, ci, alpha, rope_range) {
    batch_manifest <- draw_file_manifest %>% dplyr::arrange(.data$imputation) %>% dplyr::slice_head(n = batch_n)
    batch_imputations <- batch_manifest$imputation
    log_msg("Summarising stability batch m =", batch_n)
    batch_draws <- purrr::map_dfr(batch_manifest$parameter_draw_file, function(path_i) {
      draws_i <- readRDS(path_i)
      keep_cols <- intersect(c("imputation", ".chain", ".iteration", ".draw", parameter_cols), names(draws_i))
      draws_i[, keep_cols, drop = FALSE]
    })
    purrr::map_dfr(parameter_cols, ~ summarise_one_parameter(batch_draws[[.x]], .x, batch_n, batch_imputations, ci, alpha, rope_range))
  }
  add_transform_columns <- function(df, transform_type) {
    if (transform_type == "none") {
      df %>% dplyr::mutate(Transformed = NA_real_, Transformed_low = NA_real_, Transformed_high = NA_real_)
    } else {
      df %>% dplyr::mutate(Transformed = exp(.data$Median), Transformed_low = exp(.data$CI_low), Transformed_high = exp(.data$CI_high))
    }
  }
  make_comparison <- function(stability_summary, transform_type) {
    batches <- sort(unique(stability_summary$n_imputations))
    if (length(batches) < 2) stop("At least two batch sizes are required for a stability comparison.")
    m_previous <- batches[length(batches) - 1]
    m_final <- batches[length(batches)]
    metric_cols <- c("Median", "Mean", "SD", "CI_low", "CI_high", "pd", "ROPE_Percentage", "Transformed", "Transformed_low", "Transformed_high", "n_draws")
    previous_df <- stability_summary %>% dplyr::filter(.data$n_imputations == m_previous) %>% dplyr::select(Parameter, Parameter_clean, Parameter_Class, dplyr::all_of(metric_cols)) %>% dplyr::rename_with(~ paste0(.x, "_previous"), dplyr::all_of(metric_cols))
    final_df <- stability_summary %>% dplyr::filter(.data$n_imputations == m_final) %>% dplyr::select(Parameter, Parameter_clean, Parameter_Class, dplyr::all_of(metric_cols)) %>% dplyr::rename_with(~ paste0(.x, "_final"), dplyr::all_of(metric_cols))
    comparison <- final_df %>% dplyr::left_join(previous_df, by = c("Parameter", "Parameter_clean", "Parameter_Class")) %>%
      dplyr::mutate(
        m_previous = m_previous, m_final = m_final,
        Median_change = .data$Median_final - .data$Median_previous,
        abs_Median_change = abs(.data$Median_change),
        CI_low_change = .data$CI_low_final - .data$CI_low_previous,
        CI_high_change = .data$CI_high_final - .data$CI_high_previous,
        max_abs_CI_endpoint_change = pmax(abs(.data$CI_low_change), abs(.data$CI_high_change), na.rm = TRUE),
        pd_change = .data$pd_final - .data$pd_previous,
        abs_pd_change = abs(.data$pd_change),
        CI_excludes_null_previous = .data$CI_low_previous > 0 | .data$CI_high_previous < 0,
        CI_excludes_null_final = .data$CI_low_final > 0 | .data$CI_high_final < 0,
        CI_exclusion_changed = .data$CI_excludes_null_previous != .data$CI_excludes_null_final
      )
    if (transform_type != "none") {
      comparison <- comparison %>% dplyr::mutate(Transformed_relative_change_pct = 100 * ((.data$Transformed_final / .data$Transformed_previous) - 1), abs_Transformed_relative_change_pct = abs(.data$Transformed_relative_change_pct))
    } else {
      comparison <- comparison %>% dplyr::mutate(Transformed_relative_change_pct = NA_real_, abs_Transformed_relative_change_pct = NA_real_)
    }
    comparison %>% dplyr::mutate(
      stable_estimate = .data$abs_Median_change <= estimate_tolerance,
      stable_ci = .data$max_abs_CI_endpoint_change <= ci_endpoint_tolerance,
      stable_transformed = if (transform_type == "none") TRUE else .data$abs_Transformed_relative_change_pct <= relative_transformed_tolerance_pct,
      stable_pd = .data$abs_pd_change <= pd_tolerance,
      stable_interpretation = !.data$CI_exclusion_changed,
      stable_by_thresholds = .data$stable_estimate & .data$stable_ci & .data$stable_transformed & .data$stable_pd & .data$stable_interpretation
    )
  }

  make_stepwise_comparisons <- function(stability_summary, transform_type) {
    batches <- sort(unique(stability_summary$n_imputations))

    if (length(batches) < 2) {
      stop("At least two batch sizes are required for stepwise stability comparisons.")
    }

    purrr::map_dfr(
      seq_len(length(batches) - 1),
      function(i) {
        m_previous <- batches[i]
        m_final <- batches[i + 1]

        previous_df <- stability_summary %>%
          dplyr::filter(.data$n_imputations == m_previous) %>%
          dplyr::select(
            Parameter, Parameter_clean, Parameter_Class,
            Median, Mean, SD, CI_low, CI_high, pd, ROPE_Percentage,
            Transformed, Transformed_low, Transformed_high, n_draws
          ) %>%
          dplyr::rename_with(
            ~ paste0(.x, "_previous"),
            -c(Parameter, Parameter_clean, Parameter_Class)
          )

        final_df <- stability_summary %>%
          dplyr::filter(.data$n_imputations == m_final) %>%
          dplyr::select(
            Parameter, Parameter_clean, Parameter_Class,
            Median, Mean, SD, CI_low, CI_high, pd, ROPE_Percentage,
            Transformed, Transformed_low, Transformed_high, n_draws
          ) %>%
          dplyr::rename_with(
            ~ paste0(.x, "_final"),
            -c(Parameter, Parameter_clean, Parameter_Class)
          )

        comparison_i <- final_df %>%
          dplyr::left_join(
            previous_df,
            by = c("Parameter", "Parameter_clean", "Parameter_Class")
          ) %>%
          dplyr::mutate(
            m_previous = m_previous,
            m_final = m_final,
            batch_transition = paste0("m=", m_previous, " to m=", m_final),
            Median_change = .data$Median_final - .data$Median_previous,
            abs_Median_change = abs(.data$Median_change),
            CI_low_change = .data$CI_low_final - .data$CI_low_previous,
            CI_high_change = .data$CI_high_final - .data$CI_high_previous,
            max_abs_CI_endpoint_change = pmax(
              abs(.data$CI_low_change),
              abs(.data$CI_high_change),
              na.rm = TRUE
            ),
            pd_change = .data$pd_final - .data$pd_previous,
            abs_pd_change = abs(.data$pd_change),
            CI_excludes_null_previous = .data$CI_low_previous > 0 | .data$CI_high_previous < 0,
            CI_excludes_null_final = .data$CI_low_final > 0 | .data$CI_high_final < 0,
            CI_exclusion_changed = .data$CI_excludes_null_previous != .data$CI_excludes_null_final
          )

        if (transform_type != "none") {
          comparison_i <- comparison_i %>%
            dplyr::mutate(
              Transformed_relative_change_pct =
                100 * ((.data$Transformed_final / .data$Transformed_previous) - 1),
              abs_Transformed_relative_change_pct =
                abs(.data$Transformed_relative_change_pct)
            )
        } else {
          comparison_i <- comparison_i %>%
            dplyr::mutate(
              Transformed_relative_change_pct = NA_real_,
              abs_Transformed_relative_change_pct = NA_real_
            )
        }

        comparison_i %>%
          dplyr::mutate(
            stable_estimate = .data$abs_Median_change <= estimate_tolerance,
            stable_ci = .data$max_abs_CI_endpoint_change <= ci_endpoint_tolerance,
            stable_transformed = if (transform_type == "none") {
              TRUE
            } else {
              .data$abs_Transformed_relative_change_pct <= relative_transformed_tolerance_pct
            },
            stable_pd = .data$abs_pd_change <= pd_tolerance,
            stable_interpretation = !.data$CI_exclusion_changed,
            stable_by_thresholds =
              .data$stable_estimate &
              .data$stable_ci &
              .data$stable_transformed &
              .data$stable_pd &
              .data$stable_interpretation
          )
      }
    )
  }

  summarise_stepwise_comparisons <- function(stepwise_comparison) {
    stepwise_comparison %>%
      dplyr::group_by(.data$m_previous, .data$m_final, .data$batch_transition) %>%
      dplyr::summarise(
        n_parameters = dplyr::n(),
        n_stable = sum(.data$stable_by_thresholds, na.rm = TRUE),
        n_unstable = sum(!.data$stable_by_thresholds, na.rm = TRUE),
        percent_stable = 100 * .data$n_stable / .data$n_parameters,
        max_abs_Median_change = max(.data$abs_Median_change, na.rm = TRUE),
        median_abs_Median_change = stats::median(.data$abs_Median_change, na.rm = TRUE),
        max_abs_CI_endpoint_change = max(.data$max_abs_CI_endpoint_change, na.rm = TRUE),
        median_abs_CI_endpoint_change = stats::median(.data$max_abs_CI_endpoint_change, na.rm = TRUE),
        max_abs_pd_change = max(.data$abs_pd_change, na.rm = TRUE),
        median_abs_pd_change = stats::median(.data$abs_pd_change, na.rm = TRUE),
        max_abs_Transformed_relative_change_pct =
          max(.data$abs_Transformed_relative_change_pct, na.rm = TRUE),
        median_abs_Transformed_relative_change_pct =
          stats::median(.data$abs_Transformed_relative_change_pct, na.rm = TRUE),
        n_CI_exclusion_changed = sum(.data$CI_exclusion_changed, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$m_final)
  }

  save_gt_table <- function(data, filename_base, title, subtitle = NULL) {
    out_file <- file.path(table_dir, paste0(filename_base, ".html"))
    gt_obj <- data %>% gt::gt() %>% gt::tab_header(title = title, subtitle = subtitle)
    gt::gtsave(gt_obj, out_file)
    out_file
  }
  save_flextable_docx <- function(data, filename_base, title) {
    out_file <- file.path(table_dir, paste0(filename_base, ".docx"))
    ft <- data %>% flextable::flextable() %>% flextable::autofit() %>% flextable::align(align = "center", part = "all") %>% flextable::align(j = 1, align = "left", part = "all") %>% flextable::bold(part = "header") %>% flextable::set_caption(title)
    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, title, style = "heading 1")
    doc <- flextable::body_add_flextable(doc, ft)
    print(doc, target = out_file)
    out_file
  }

  draw_file_manifest <- find_parameter_draw_files(paths) %>%
    dplyr::mutate(draw_file_valid = purrr::map_lgl(.data$parameter_draw_file, rds_ok)) %>%
    dplyr::filter(.data$draw_file_valid) %>%
    dplyr::arrange(.data$imputation)
  if (nrow(draw_file_manifest) < 2) stop("At least two valid per-imputation parameter draw files are needed. Run Step 6 first.")
  log_msg("Found", nrow(draw_file_manifest), "valid per-imputation draw file(s).")

  first_draws <- readRDS(draw_file_manifest$parameter_draw_file[1])
  parameter_cols <- select_parameter_columns(first_draws, parameter_regex, primary_parameters, exclude_intercept)
  rm(first_draws); gc()
  if (length(parameter_cols) == 0) stop("No parameter columns selected for stability checking.")

  n_available <- nrow(draw_file_manifest)
  batch_sizes <- choose_batch_sizes(n_available, cfg)
  if (length(batch_sizes) < 2) stop("Need at least two batch sizes to evaluate imputation-count stability.")

  rope_range <- resolve_rope(analysis_spec$summary)
  transform_type <- cfg$transform %||% "auto"
  if (identical(transform_type, "auto")) transform_type <- effect_transform_type(analysis_spec)
  transformed_label <- transform_label(transform_type)

  stability_summary <- purrr::map_dfr(batch_sizes, ~ summarise_batch(.x, draw_file_manifest, parameter_cols, ci, alpha, rope_range)) %>%
    add_transform_columns(transform_type)
  saveRDS(stability_summary, file.path(table_dir, "imputation_stability_all_batches.rds"), compress = FALSE)
  readr::write_csv(stability_summary, file.path(table_dir, "imputation_stability_all_batches.csv"))

  stability_comparison <- make_comparison(stability_summary, transform_type)
  saveRDS(stability_comparison, file.path(table_dir, "imputation_stability_final_comparison_full.rds"), compress = FALSE)
  readr::write_csv(stability_comparison, file.path(table_dir, "imputation_stability_final_comparison_full.csv"))

  stability_comparison_display <- stability_comparison %>%
    dplyr::transmute(
      Parameter = .data$Parameter_clean,
      `Previous m` = .data$m_previous,
      `Final m` = .data$m_final,
      `Estimate, previous` = fmt_num(.data$Median_previous, 3),
      `Estimate, final` = fmt_num(.data$Median_final, 3),
      `Estimate change` = fmt_num(.data$Median_change, 3),
      `95% CrI, previous` = fmt_ci(.data$CI_low_previous, .data$CI_high_previous, 3),
      `95% CrI, final` = fmt_ci(.data$CI_low_final, .data$CI_high_final, 3),
      pd_previous = fmt_num(.data$pd_previous, 3),
      pd_final = fmt_num(.data$pd_final, 3),
      `ROPE %, previous` = fmt_num(.data$ROPE_Percentage_previous, 1),
      `ROPE %, final` = fmt_num(.data$ROPE_Percentage_final, 1),
      `Stable by thresholds` = ifelse(.data$stable_by_thresholds, "Yes", "No")
    )
  if (transform_type != "none") {
    transformed_display <- stability_comparison %>% dplyr::transmute(
      Parameter = .data$Parameter_clean,
      !!paste0(transformed_label, ", previous") := fmt_num(.data$Transformed_previous, 3),
      !!paste0(transformed_label, ", final") := fmt_num(.data$Transformed_final, 3),
      `Relative change, %` = fmt_num(.data$Transformed_relative_change_pct, 2)
    )
    stability_comparison_display <- stability_comparison_display %>% dplyr::left_join(transformed_display, by = "Parameter")
  }
  readr::write_csv(stability_comparison_display, file.path(table_dir, "imputation_stability_final_comparison_display.csv"))
  save_gt_table(stability_comparison_display, "imputation_stability_final_comparison", "Imputation-count Stability: Final Comparison", glue("Comparison of m = {max(batch_sizes)} with the previous evaluated batch"))
  save_flextable_docx(stability_comparison_display, "imputation_stability_final_comparison", "Imputation-count Stability: Final Comparison")

  # ------------------------------------------------------------
  # Stepwise comparisons across every increase in m
  # ------------------------------------------------------------
  stepwise_comparison <- make_stepwise_comparisons(
    stability_summary = stability_summary,
    transform_type = transform_type
  )

  stepwise_summary <- summarise_stepwise_comparisons(stepwise_comparison)

  saveRDS(stepwise_comparison, file.path(table_dir, "imputation_stability_stepwise_comparison_full.rds"), compress = FALSE)
  readr::write_csv(stepwise_comparison, file.path(table_dir, "imputation_stability_stepwise_comparison_full.csv"))
  readr::write_csv(stepwise_summary, file.path(table_dir, "imputation_stability_stepwise_summary.csv"))

  stepwise_summary_display <- stepwise_summary %>%
    dplyr::transmute(
      `Batch transition` = .data$batch_transition,
      `Parameters checked` = .data$n_parameters,
      `Stable, n` = .data$n_stable,
      `Unstable, n` = .data$n_unstable,
      `Stable, %` = fmt_num(.data$percent_stable, 1),
      `Maximum absolute estimate change` = fmt_num(.data$max_abs_Median_change, 3),
      `Median absolute estimate change` = fmt_num(.data$median_abs_Median_change, 3),
      `Maximum absolute CrI-endpoint change` = fmt_num(.data$max_abs_CI_endpoint_change, 3),
      `Median absolute CrI-endpoint change` = fmt_num(.data$median_abs_CI_endpoint_change, 3),
      `Maximum absolute pd change` = fmt_num(.data$max_abs_pd_change, 3),
      `Maximum relative OR change, %` = fmt_num(.data$max_abs_Transformed_relative_change_pct, 2),
      `Median relative OR change, %` = fmt_num(.data$median_abs_Transformed_relative_change_pct, 2),
      `CrI exclusion changed, n` = .data$n_CI_exclusion_changed
    )

  readr::write_csv(stepwise_summary_display, file.path(table_dir, "imputation_stability_stepwise_summary_display.csv"))

  save_gt_table(
    stepwise_summary_display,
    "imputation_stability_stepwise_summary",
    "Imputation-count Stability: Stepwise Summary",
    "Quantitative changes at each increase in the number of imputations"
  )

  save_flextable_docx(
    stepwise_summary_display,
    "imputation_stability_stepwise_summary",
    "Imputation-count Stability: Stepwise Summary"
  )

  stepwise_change_plot_data <- stepwise_summary %>%
    dplyr::select(
      .data$m_final,
      .data$batch_transition,
      .data$max_abs_Median_change,
      .data$median_abs_Median_change,
      .data$max_abs_CI_endpoint_change,
      .data$median_abs_CI_endpoint_change
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        .data$max_abs_Median_change,
        .data$median_abs_Median_change,
        .data$max_abs_CI_endpoint_change,
        .data$median_abs_CI_endpoint_change
      ),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(
        .data$metric,
        max_abs_Median_change = "Maximum absolute estimate change",
        median_abs_Median_change = "Median absolute estimate change",
        max_abs_CI_endpoint_change = "Maximum absolute CrI-endpoint change",
        median_abs_CI_endpoint_change = "Median absolute CrI-endpoint change"
      )
    )

  stepwise_change_plot <- ggplot2::ggplot(
    stepwise_change_plot_data,
    ggplot2::aes(x = .data$m_final, y = .data$value, linetype = .data$metric, shape = .data$metric)
  ) +
    ggplot2::geom_hline(yintercept = estimate_tolerance, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(stepwise_summary$m_final))) +
    ggplot2::labs(
      title = "Stepwise Imputation-count Stability",
      subtitle = "Change when increasing from the previous m to the current m",
      x = "Current number of imputations",
      y = "Absolute change on coefficient scale",
      linetype = NULL,
      shape = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    filename = file.path(figure_dir, "imputation_stability_stepwise_change.png"),
    plot = stepwise_change_plot,
    width = 10,
    height = 6,
    dpi = 300
  )

  ggplot2::ggsave(
    filename = file.path(figure_dir, "imputation_stability_stepwise_change.pdf"),
    plot = stepwise_change_plot,
    width = 10,
    height = 6
  )

  stability_settings <- tibble::tibble(
    setting = c("available_valid_imputations", "batch_sizes", "selected_parameters", "parameter_regex", "exclude_intercept", "transform_type", "estimate_tolerance", "ci_endpoint_tolerance", "relative_transformed_tolerance_pct", "pd_tolerance", "ci_probability", "rope_range"),
    value = c(as.character(n_available), paste(batch_sizes, collapse = ", "), paste(parameter_cols, collapse = ", "), parameter_regex, as.character(exclude_intercept), transform_type, as.character(estimate_tolerance), as.character(ci_endpoint_tolerance), as.character(relative_transformed_tolerance_pct), as.character(pd_tolerance), as.character(ci), ifelse(is.null(rope_range), "not specified", paste(round(rope_range, 5), collapse = " to ")))
  )
  readr::write_csv(stability_settings, file.path(table_dir, "imputation_stability_settings.csv"))
  save_gt_table(stability_settings, "imputation_stability_settings", "Imputation-count Stability Settings")
  save_flextable_docx(stability_settings, "imputation_stability_settings", "Imputation-count Stability Settings")

  plot_parameter_order <- stability_comparison %>% dplyr::arrange(dplyr::desc(.data$stable_by_thresholds == FALSE), dplyr::desc(.data$abs_Median_change)) %>% dplyr::slice_head(n = max_plot_parameters) %>% dplyr::pull(.data$Parameter)
  plot_data <- stability_summary %>% dplyr::filter(.data$Parameter %in% plot_parameter_order)
  if (nrow(plot_data) > 0) {
    if (transform_type != "none") {
      plot_data <- plot_data %>% dplyr::mutate(plot_estimate = .data$Transformed, plot_low = .data$Transformed_low, plot_high = .data$Transformed_high, null_value = 1)
      y_label <- transformed_label
    } else {
      plot_data <- plot_data %>% dplyr::mutate(plot_estimate = .data$Median, plot_low = .data$CI_low, plot_high = .data$CI_high, null_value = 0)
      y_label <- "Estimate"
    }
    trajectory_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$n_imputations, y = .data$plot_estimate)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$null_value), linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$plot_low, ymax = .data$plot_high), alpha = 0.6) +
      ggplot2::geom_line() + ggplot2::geom_point(size = 1.7) +
      ggplot2::facet_wrap(~ Parameter_clean, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = batch_sizes) +
      ggplot2::labs(title = "Imputation-count Stability of Posterior Summaries", subtitle = "Nested subsets of completed imputations", x = "Number of imputations included", y = y_label) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), strip.text = ggplot2::element_text(size = 8))
    plot_height <- max(5, 2.2 * ceiling(length(unique(plot_data$Parameter)) / 2))
    ggplot2::ggsave(file.path(figure_dir, "imputation_stability_trajectories.png"), plot = trajectory_plot, width = 10, height = plot_height, dpi = 300)
    ggplot2::ggsave(file.path(figure_dir, "imputation_stability_trajectories.pdf"), plot = trajectory_plot, width = 10, height = plot_height)
  }

  report_qmd <- file.path(report_dir, "imputation_stability_report.qmd")

  # Copy the trajectory figure into the report folder before rendering.
  # DOCX often embeds images from outside the report folder successfully, but HTML
  # usually links to image files. Keeping a copy beside the .qmd, and using
  # embed-resources: true, makes the HTML output portable and avoids broken images.
  trajectory_png <- file.path(figure_dir, "imputation_stability_trajectories.png")
  trajectory_png_report <- file.path(report_dir, "imputation_stability_trajectories.png")

  if (file.exists(trajectory_png)) {
    file.copy(
      from = trajectory_png,
      to = trajectory_png_report,
      overwrite = TRUE
    )
  }

  stepwise_png <- file.path(figure_dir, "imputation_stability_stepwise_change.png")
  stepwise_png_report <- file.path(report_dir, "imputation_stability_stepwise_change.png")

  if (file.exists(stepwise_png)) {
    file.copy(
      from = stepwise_png,
      to = stepwise_png_report,
      overwrite = TRUE
    )
  }

  report_text <- c(
    "---",
    "title: \"Imputation-count Stability Check\"",
    "format:",
    "  html:",
    "    toc: true",
    "    embed-resources: true",
    "  docx: default",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "```{r setup}",
    "library(tidyverse)",
    "library(gt)",
    "",
    "# The report is rendered from results/publication/mi_stability/report.",
    "report_dir <- getwd()",
    "stability_dir <- normalizePath(file.path(report_dir, '..'), mustWork = FALSE)",
    "table_dir <- file.path(stability_dir, 'tables')",
    "figure_dir <- file.path(stability_dir, 'figures')",
    "trajectory_png_report <- file.path(report_dir, 'imputation_stability_trajectories.png')",
    "trajectory_png_source <- file.path(figure_dir, 'imputation_stability_trajectories.png')",
    "stepwise_png_report <- file.path(report_dir, 'imputation_stability_stepwise_change.png')",
    "",
    "# Use the report-local copy for HTML portability. Fall back to the source",
    "# figure path if the report-local copy is missing.",
    "trajectory_png <- if (file.exists(trajectory_png_report)) {",
    "  trajectory_png_report",
    "} else {",
    "  trajectory_png_source",
    "}",
    "```",
    "",
    "# Purpose",
    "",
    "This report checks how posterior summaries change as the number of completed imputations increases. The goal is to assess whether the current number of imputations is sufficient for selected primary parameters.",
    "",
    "This is a practical Monte Carlo stability assessment, not a proof that the missing-data assumptions are correct.",
    "",
    "# Settings",
    "",
    "```{r settings-table, eval=file.exists(file.path(table_dir, 'imputation_stability_settings.csv'))}",
    "settings <- readr::read_csv(",
    "  file.path(table_dir, 'imputation_stability_settings.csv'),",
    "  show_col_types = FALSE",
    ")",
    "gt(settings)",
    "```",
    "",
    "# Final comparison",
    "",
    "The table below compares the final evaluated number of imputations with the previous evaluated batch.",
    "",
    "```{r final-comparison-table, eval=file.exists(file.path(table_dir, 'imputation_stability_final_comparison_display.csv'))}",
    "display <- readr::read_csv(",
    "  file.path(table_dir, 'imputation_stability_final_comparison_display.csv'),",
    "  show_col_types = FALSE",
    ")",
    "gt(display)",
    "```",
    "",
    "# Stepwise quantitative stability",
    "",
    "The table below quantifies how much the selected posterior summaries changed at each increase in the number of imputations.",
    "",
    "```{r stepwise-summary-table, eval=file.exists(file.path(table_dir, 'imputation_stability_stepwise_summary_display.csv'))}",
    "stepwise <- readr::read_csv(",
    "  file.path(table_dir, 'imputation_stability_stepwise_summary_display.csv'),",
    "  show_col_types = FALSE",
    ")",
    "gt(stepwise)",
    "```",
    "",
    "```{r stepwise-change-plot, fig.width=10, fig.height=6, out.width='100%', eval=file.exists(stepwise_png_report)}",
    "knitr::include_graphics(stepwise_png_report)",
    "```",
    "",
    "# Stability plot",
    "",
    "```{r imputation-stability-plot, fig.width=10, fig.height=8, out.width='100%', eval=file.exists(trajectory_png)}",
    "knitr::include_graphics(trajectory_png)",
    "```",
    "",
    "# Interpretation guidance",
    "",
    "A result labelled stable by thresholds means that the change from the previous evaluated batch to the final evaluated batch was within the configured tolerances. The default thresholds are pragmatic and should be replaced by thresholds that match the scientific question whenever possible.",
    "",
    "Suggested reporting language:",
    "",
    "> We used an adaptive multiple-imputation strategy because each imputed-data Bayesian model was computationally expensive. We first fitted an initial set of imputed datasets and assessed the stability of prespecified primary posterior summaries. We increased the number of imputations until the pooled posterior medians, credible intervals, posterior direction probabilities, and substantive conclusions changed negligibly with additional imputations. The final analysis used m = XX imputations.",
    "",
    "# Full numeric results",
    "",
    "Full numeric outputs are available in the tables folder:",
    "",
    "- `imputation_stability_all_batches.csv`",
    "- `imputation_stability_final_comparison_full.csv`",
    "- `imputation_stability_final_comparison_display.csv`",
    "- `imputation_stability_settings.csv`",
    "- `imputation_stability_stepwise_summary.csv`",
    "- `imputation_stability_stepwise_comparison_full.csv`"
  )
  writeLines(report_text, report_qmd)

  if (isTRUE(render_quarto)) {
    quarto_bin <- Sys.which("quarto")
    if (nzchar(quarto_bin)) {
      old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE); setwd(report_dir)
      html_status <- system2(quarto_bin, args = c("render", basename(report_qmd), "--to", "html"))
      docx_status <- system2(quarto_bin, args = c("render", basename(report_qmd), "--to", "docx"))
      log_msg("Quarto HTML render status:", html_status)
      log_msg("Quarto DOCX render status:", docx_status)
    } else {
      log_msg("Quarto not found; wrote .qmd report template without rendering.")
    }
  }

  n_unstable <- sum(!stability_comparison$stable_by_thresholds, na.rm = TRUE)
  log_msg("Imputation stability check completed. Parameters not stable by thresholds:", n_unstable, "of", nrow(stability_comparison))
  log_msg("Outputs written to:", stability_dir)
  gc()

}, analysis_spec)
