# ============================================================
# 08_publication_results.R
# Publication-ready tables, figures, and report template
#
# Robust version:
# - keeps numeric estimates separate from formatted estimates
# - exponentiates only when appropriate
# - supports logistic / log-link models and Gaussian models
# ============================================================

source("00_config.R")
source("00_common_functions.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(gt)
  library(flextable)
  library(officer)
  library(glue)
  library(stringr)
  library(forcats)
  library(brms)
})

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 8: Publication-ready outputs", {

  resolve_rope_range <- function(summary_spec) {
    # Support both old and new config styles:
    #
    # Old/list style:
    #   summary$rope <- list(method = "fixed", fixed_range = c(-0.1, 0.1))
    #
    # Current/simple style:
    #   summary$test <- c("p_direction", "rope")
    #   summary$rope_range <- "auto_logit_5pct"

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
        } else {
          warning(
            "Unknown summary$rope_range value: ",
            rope_range_cfg[[1]],
            ". ROPE values will be shown as not available."
          )
        }
      }
    }

    if (!identical(rope_method, "fixed")) {
      rope_range <- NULL
    }

    rope_range
  }

  format_rope_range <- function(summary_spec) {
    rr <- resolve_rope_range(summary_spec)

    if (is.null(rr) || length(rr) != 2 || any(is.na(rr))) {
      return("not specified")
    }

    paste0(
      round(rr[[1]], 4),
      " to ",
      round(rr[[2]], 4)
    )
  }


  
  # ------------------------------------------------------------
  # Output folders
  # ------------------------------------------------------------
  
  pub_dir <- file.path(paths$results, "publication")
  table_dir <- file.path(pub_dir, "tables")
  figure_dir <- file.path(pub_dir, "figures")
  report_dir <- file.path(pub_dir, "report")
  
  dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------
  
  fmt_num <- function(x, digits = 2) {
    ifelse(
      is.na(x),
      NA_character_,
      formatC(x, format = "f", digits = digits)
    )
  }
  
  fmt_ci <- function(low, high, digits = 2) {
    paste0(
      "[",
      fmt_num(low, digits),
      ", ",
      fmt_num(high, digits),
      "]"
    )
  }
  
  clean_parameter_name <- function(x) {
    x %>%
      str_replace("^b_", "") %>%
      str_replace("^sd_", "SD: ") %>%
      str_replace_all(":", " x ") %>%
      str_replace_all("_z$", " per SD") %>%
      str_replace_all("_", " ")
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
  
  safe_filename <- function(x) {
    x %>%
      stringr::str_replace_all("[^A-Za-z0-9_\\.-]+", "_") %>%
      stringr::str_replace_all("_+", "_") %>%
      stringr::str_replace("^_", "") %>%
      stringr::str_replace("_$", "")
  }
  
  get_report_formula_core <- function(formula) {
    if (inherits(formula, "brmsformula") && !is.null(formula$formula)) {
      return(formula$formula)
    }

    formula
  }

  get_report_formula_text <- function(formula) {
    if (is.null(formula)) {
      return("not available")
    }

    paste(deparse(get_report_formula_core(formula), width.cutoff = 500), collapse = " ")
  }

  compact_report_text <- function(x) {
    x <- paste(x, collapse = " ")
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  get_report_prior_text <- function(prior) {
    if (is.null(prior)) {
      return("not available")
    }

    prior_df <- tryCatch(
      as.data.frame(prior),
      error = function(e) NULL
    )

    if (is.null(prior_df) || nrow(prior_df) == 0) {
      return("not available")
    }

    # Keep the reader-facing prior display compact.  These are the fields
    # most users need in the methods/settings table.
    keep_cols <- intersect(
      c("prior", "class", "coef", "group", "resp", "dpar", "nlpar"),
      names(prior_df)
    )

    prior_df <- prior_df[, keep_cols, drop = FALSE]

    # Convert missing values and empty strings to blanks for clean display.
    prior_df[] <- lapply(prior_df, function(z) {
      z <- as.character(z)
      z[is.na(z)] <- ""
      z
    })

    apply(
      prior_df,
      1,
      function(row_i) {
        parts <- paste0(names(row_i), "=", row_i)
        parts <- parts[nzchar(row_i)]
        paste(parts, collapse = ", ")
      }
    ) |>
      paste(collapse = "; ") |>
      compact_report_text()
  }

  extract_special_term_vars <- function(formula, fun) {
    if (is.null(formula)) return(character(0))
    formula_text <- get_report_formula_text(formula)
    pattern <- paste0("\\b", fun, "\\s*\\(\\s*`?([A-Za-z.][A-Za-z0-9._]*)`?")
    matches <- gregexpr(pattern, formula_text, perl = TRUE)
    hits <- regmatches(formula_text, matches)[[1]]
    if (length(hits) == 0 || identical(hits, character(0))) return(character(0))
    unique(sub(pattern, "\\1", hits, perl = TRUE))
  }
  
  find_col <- function(df, candidates) {
    out <- intersect(candidates, names(df))
    if (length(out) == 0) NA_character_ else out[1]
  }
  
  is_intercept <- function(x) {
    x %in% c("b_Intercept", "Intercept", "(Intercept)")
  }
  
  effect_transform_type <- function(analysis_spec) {
    fam <- analysis_spec$outcome$family
    link <- analysis_spec$outcome$link
    
    if (fam %in% c("bernoulli", "ordinal", "categorical") && link == "logit") {
      return("odds_ratio")
    }
    
    if (fam %in% c("poisson", "negbinomial") && link == "log") {
      return("rate_ratio")
    }
    
    if (fam == "gaussian" && link == "log") {
      return("multiplicative_ratio")
    }
    
    "none"
  }
  
  transform_label <- function(transform_type) {
    switch(
      transform_type,
      odds_ratio = "Odds ratio",
      rate_ratio = "Rate ratio",
      multiplicative_ratio = "Ratio",
      "Transformed estimate"
    )
  }
  
  save_gt_table <- function(gt_obj, filename_base) {
    out_file <- file.path(table_dir, paste0(filename_base, ".html"))
    gt::gtsave(gt_obj, out_file)
    out_file
  }
  
  save_flextable_docx <- function(ft, filename_base, title = NULL) {
    out_file <- file.path(table_dir, paste0(filename_base, ".docx"))
    
    doc <- officer::read_docx()
    
    if (!is.null(title)) {
      doc <- officer::body_add_par(doc, title, style = "heading 1")
    }
    
    doc <- flextable::body_add_flextable(doc, ft)
    
    print(doc, target = out_file)
    
    out_file
  }
  
  # ------------------------------------------------------------
  # Load outputs
  # ------------------------------------------------------------
  
  parameter_summary_file <- file.path(paths$results, "parameter_summary.rds")
  diagnostics_file <- file.path(paths$results, "diagnostics.rds")
  missing_y_summary_file <- file.path(paths$results, "missing_y_summary.rds")
  model_spec_file <- file.path(paths$objects, "model_spec.rds")
  fit_manifest_file <- file.path(paths$objects, "fit_manifest.rds")
  fit_status_file <- file.path(paths$objects, "fit_status.rds")
  parameter_manifest_file <- file.path(paths$objects, "parameter_manifest.rds")
  
  if (!file.exists(parameter_summary_file)) {
    stop("Missing parameter_summary.rds.")
  }
  
  parameter_summary <- readRDS(parameter_summary_file)
  
  diagnostics <- if (file.exists(diagnostics_file)) {
    readRDS(diagnostics_file)
  } else {
    NULL
  }
  
  missing_y_summary <- if (file.exists(missing_y_summary_file)) {
    readRDS(missing_y_summary_file)
  } else {
    NULL
  }
  
  model_spec <- if (file.exists(model_spec_file)) {
    readRDS(model_spec_file)
  } else {
    NULL
  }
  
  fit_manifest <- if (file.exists(fit_manifest_file)) {
    readRDS(fit_manifest_file)
  } else {
    NULL
  }
  
  transform_type <- effect_transform_type(analysis_spec)
  transformed_label <- transform_label(transform_type)
  
  log_msg("Publication transform type:", transform_type)
  
  # ------------------------------------------------------------
  # Standardise parameter summary
  # ------------------------------------------------------------
  
  central_col <- find_col(
    parameter_summary,
    c("Median", "median", "Mean", "mean", "Estimate", "estimate", "MAP")
  )
  
  ci_low_col <- find_col(
    parameter_summary,
    c("CI_low", "CI_low_", "CI_low.", "lower", "Lower")
  )
  
  ci_high_col <- find_col(
    parameter_summary,
    c("CI_high", "CI_high_", "CI_high.", "upper", "Upper")
  )
  
  pd_col <- find_col(
    parameter_summary,
    c("pd", "p_direction", "PD")
  )
  
  rope_col <- find_col(
    parameter_summary,
    c("ROPE_Percentage", "ROPE_Percentage_", "Percentage_in_ROPE")
  )
  
  if (is.na(central_col)) {
    stop("Could not find a central estimate column in parameter_summary.")
  }
  
  if (is.na(ci_low_col) || is.na(ci_high_col)) {
    stop("Could not find CI_low / CI_high columns in parameter_summary.")
  }
  
  main_effect_table <- parameter_summary %>%
    mutate(
      Parameter_raw = Parameter,
      Parameter_Class = if ("Parameter_Class" %in% names(parameter_summary)) {
        as.character(Parameter_Class)
      } else {
        classify_parameter(Parameter)
      },
      Parameter_clean = clean_parameter_name(Parameter),
      Estimate_num = as.numeric(.data[[central_col]]),
      CI_low_num = as.numeric(.data[[ci_low_col]]),
      CI_high_num = as.numeric(.data[[ci_high_col]])
    )
  
  if (!is.na(pd_col)) {
    main_effect_table <- main_effect_table %>%
      mutate(pd_num = as.numeric(.data[[pd_col]]))
  } else {
    main_effect_table <- main_effect_table %>%
      mutate(pd_num = NA_real_)
  }
  
  if (!is.na(rope_col)) {
    main_effect_table <- main_effect_table %>%
      mutate(ROPE_Percentage_num = as.numeric(.data[[rope_col]]))
  } else {
    main_effect_table <- main_effect_table %>%
      mutate(ROPE_Percentage_num = NA_real_)
  }
  
  if (transform_type != "none") {
    main_effect_table <- main_effect_table %>%
      mutate(
        Transformed_num = exp(Estimate_num),
        Transformed_low_num = exp(CI_low_num),
        Transformed_high_num = exp(CI_high_num)
      )
  } else {
    main_effect_table <- main_effect_table %>%
      mutate(
        Transformed_num = NA_real_,
        Transformed_low_num = NA_real_,
        Transformed_high_num = NA_real_
      )
  }
  
  # Keep full numeric version
  saveRDS(
    main_effect_table,
    file.path(table_dir, "main_effect_table_full.rds"),
    compress = FALSE
  )
  
  readr::write_csv(
    main_effect_table,
    file.path(table_dir, "main_effect_table_full.csv")
  )
  
  # ------------------------------------------------------------
  # Special brms parameters from s() and mo()
  # ------------------------------------------------------------
  
  special_parameter_table <- main_effect_table %>%
    dplyr::filter(Parameter_Class %in% c("smooth_sd", "smooth_basis", "monotonic_simplex", "special_brms")) %>%
    dplyr::transmute(
      Parameter = Parameter_clean,
      Parameter_raw = Parameter_raw,
      Parameter_Class = Parameter_Class,
      Estimate = fmt_num(Estimate_num, 3),
      `95% CrI` = fmt_ci(CI_low_num, CI_high_num, 3),
      pd = ifelse(is.na(pd_num), NA_character_, fmt_num(pd_num, 3)),
      `ROPE %` = ifelse(is.na(ROPE_Percentage_num), NA_character_, fmt_num(ROPE_Percentage_num, 1))
    )
  
  special_terms_sentence <- "No smooth-term or monotonic-effect auxiliary parameters were detected in the extracted posterior draws."
  
  if (nrow(special_parameter_table) > 0) {
    saveRDS(special_parameter_table, file.path(table_dir, "special_parameter_table.rds"), compress = FALSE)
    readr::write_csv(special_parameter_table, file.path(table_dir, "special_parameter_table.csv"))
    special_parameter_gt <- special_parameter_table %>% gt() %>%
      tab_header(title = "Supplementary Smooth and Monotonic Parameters", subtitle = "Auxiliary brms parameters from s() and mo() terms")
    save_gt_table(special_parameter_gt, "special_parameter_table")
    special_terms_sentence <- "Supplementary brms parameters from smooth or monotonic terms were detected and written to special_parameter_table.csv. These parameters are usually less directly interpretable than ordinary regression coefficients; conditional-effect plots are recommended for reporting s() and mo() terms."
  }
  
  # ------------------------------------------------------------
  # Display table
  # ------------------------------------------------------------
  
  display_base <- main_effect_table %>%
    filter(Parameter_Class == "fixed") %>%
    filter(!is_intercept(Parameter_raw)) %>%
    transmute(
      Parameter = Parameter_clean,
      Estimate = fmt_num(Estimate_num, 2),
      `95% CrI` = fmt_ci(CI_low_num, CI_high_num, 2),
      pd = ifelse(is.na(pd_num), NA_character_, fmt_num(pd_num, 3)),
      `ROPE %` = ifelse(
        is.na(ROPE_Percentage_num),
        NA_character_,
        fmt_num(ROPE_Percentage_num, 1)
      ),
      Transformed = fmt_num(Transformed_num, 2),
      Transformed_CrI = fmt_ci(Transformed_low_num, Transformed_high_num, 2)
    )
  
  if (transform_type != "none") {
    main_effect_table_display <- display_base %>%
      rename(
        !!transformed_label := Transformed,
        !!paste0("95% CrI, ", transformed_label) := Transformed_CrI
      )
  } else {
    main_effect_table_display <- display_base %>%
      select(-Transformed, -Transformed_CrI)
  }
  
  saveRDS(
    main_effect_table_display,
    file.path(table_dir, "main_effect_table_display.rds"),
    compress = FALSE
  )
  
  readr::write_csv(
    main_effect_table_display,
    file.path(table_dir, "main_effect_table_display.csv")
  )
  
  main_effect_table_center_cols <- setdiff(
    names(main_effect_table_display),
    "Parameter"
  )

  main_effect_gt <- main_effect_table_display %>%
    gt() %>%
    tab_header(
      title = "Posterior Summary of Fixed Effects"
    ) %>%
    cols_align(
      align = "center",
      columns = dplyr::all_of(main_effect_table_center_cols)
    ) %>%
    tab_source_note(
      source_note = "CrI = credible interval; pd = probability of direction; ROPE = region of practical equivalence."
    )
  
  save_gt_table(main_effect_gt, "main_effect_table")
  
  main_effect_ft <- main_effect_table_display %>%
    flextable() %>%
    autofit() %>%
    align(align = "center", part = "all") %>%
    align(j = "Parameter", align = "left", part = "all") %>%
    bold(part = "header") %>%
    set_caption("Posterior Summary of Fixed Effects")
  
  save_flextable_docx(
    main_effect_ft,
    "main_effect_table",
    title = "Posterior Summary of Fixed Effects"
  )
  
  # ------------------------------------------------------------
  # Forest plot
  # ------------------------------------------------------------
  
  forest_data <- main_effect_table %>%
    filter(Parameter_Class == "fixed") %>%
    filter(!is_intercept(Parameter_raw)) %>%
    filter(str_detect(Parameter_raw, "^b_"))
  
  if (nrow(forest_data) > 0) {
    if (transform_type != "none") {
      forest_plot_data <- forest_data %>%
        mutate(
          plot_estimate = Transformed_num,
          plot_low = Transformed_low_num,
          plot_high = Transformed_high_num,
          Parameter_clean = forcats::fct_reorder(Parameter_clean, plot_estimate)
        )
      
      x_label <- transformed_label
      null_value <- 1
      x_scale <- scale_x_log10()
    } else {
      forest_plot_data <- forest_data %>%
        mutate(
          plot_estimate = Estimate_num,
          plot_low = CI_low_num,
          plot_high = CI_high_num,
          Parameter_clean = forcats::fct_reorder(Parameter_clean, plot_estimate)
        )
      
      x_label <- "Estimate"
      null_value <- 0
      x_scale <- scale_x_continuous()
    }
    
    forest_plot <- ggplot(
      forest_plot_data,
      aes(
        x = plot_estimate,
        y = Parameter_clean
      )
    ) +
      geom_vline(xintercept = null_value, linetype = "dashed") +
      geom_pointrange(
        aes(
          xmin = plot_low,
          xmax = plot_high
        )
      ) +
      x_scale +
      labs(
        x = x_label,
        y = NULL,
        title = "Posterior Fixed-Effect Estimates",
        subtitle = "Points indicate posterior medians; intervals indicate 95% credible intervals"
      ) +
      theme_bw(base_size = 12)
    
    ggsave(
      filename = file.path(figure_dir, "forest_plot_fixed_effects.png"),
      plot = forest_plot,
      width = 8,
      height = max(5, 0.30 * nrow(forest_plot_data)),
      dpi = 300
    )
    
    ggsave(
      filename = file.path(figure_dir, "forest_plot_fixed_effects.pdf"),
      plot = forest_plot,
      width = 8,
      height = max(5, 0.30 * nrow(forest_plot_data))
    )
  }
  
  # ------------------------------------------------------------
  # Conditional-effect plots for s() and mo() terms
  # ------------------------------------------------------------
  
  conditional_effects_sentence <- "Conditional-effect plots were not requested or no eligible special terms were detected."
  conditional_effect_manifest <- tibble::tibble()
  reporting_ce <- analysis_spec$reporting$conditional_effects %||% list(enabled = FALSE)
  ce_enabled <- reporting_ce$enabled %||% FALSE
  
  if (isTRUE(ce_enabled) && !is.null(model_spec) && !is.null(fit_manifest)) {
    ce_effects <- reporting_ce$effects %||% "auto"
    if (identical(ce_effects, "auto")) {
      ce_effects <- unique(c(
        extract_special_term_vars(model_spec$formula, "s"),
        extract_special_term_vars(model_spec$formula, "mo")
      ))
    }
    
    if (length(ce_effects) > 0) {
      valid_fit_manifest <- fit_manifest %>%
        dplyr::mutate(fit_valid = purrr::map_lgl(fit_file, rds_ok)) %>%
        dplyr::filter(fit_valid)
      
      if (nrow(valid_fit_manifest) > 0) {
        fit_for_ce <- readRDS(valid_fit_manifest$fit_file[1])
        ce_rows <- list()
        for (effect_i in ce_effects) {
          log_msg("Creating conditional-effect plot for:", effect_i)
          safe_effect <- safe_filename(effect_i)
          png_file <- file.path(figure_dir, paste0("conditional_effect_", safe_effect, ".png"))
          pdf_file <- file.path(figure_dir, paste0("conditional_effect_", safe_effect, ".pdf"))
          ce_result <- tryCatch({
            ce <- brms::conditional_effects(
              fit_for_ce,
              effects = effect_i,
              re_formula = reporting_ce$re_formula %||% NA,
              resolution = reporting_ce$resolution %||% 100
            )
            p <- plot(ce, plot = FALSE)[[1]] +
              ggplot2::theme_bw(base_size = 12) +
              ggplot2::labs(
                title = paste("Conditional effect:", effect_i),
                subtitle = paste0("Generated from representative fit imputation ", valid_fit_manifest$imputation[1], "; posterior tables are pooled across valid fits")
              )
            ggplot2::ggsave(png_file, plot = p, width = 7, height = 5, dpi = 300)
            ggplot2::ggsave(pdf_file, plot = p, width = 7, height = 5)
            tibble::tibble(effect = effect_i, representative_imputation = valid_fit_manifest$imputation[1], png_file = basename(png_file), pdf_file = basename(pdf_file), status = "created", error = NA_character_)
          }, error = function(e) {
            tibble::tibble(effect = effect_i, representative_imputation = valid_fit_manifest$imputation[1], png_file = NA_character_, pdf_file = NA_character_, status = "failed", error = conditionMessage(e))
          })
          ce_rows[[length(ce_rows) + 1]] <- ce_result
        }
        conditional_effect_manifest <- dplyr::bind_rows(ce_rows)
        readr::write_csv(conditional_effect_manifest, file.path(table_dir, "conditional_effects_manifest.csv"))
        n_created <- sum(conditional_effect_manifest$status == "created", na.rm = TRUE)
        conditional_effects_sentence <- glue("Conditional-effect plots were requested for {length(ce_effects)} effect(s); {n_created} plot(s) were created from the first valid representative fit. For multiply imputed analyses, these plots should be interpreted as visual summaries of nonlinear term shape; posterior coefficient summaries remain pooled across valid fits.")
        rm(fit_for_ce)
        gc()
      }
    }
  }
  
  # ------------------------------------------------------------
  # Diagnostics summary
  # ------------------------------------------------------------
  
  diagnostics_sentence <- "Model diagnostics were evaluated using posterior sampler diagnostics."
  
  if (!is.null(diagnostics)) {
    diagnostics_summary <- diagnostics %>%
      summarise(
        n_fits = n(),
        n_fits_with_divergences = sum(divergent > 0, na.rm = TRUE),
        total_divergences = sum(divergent, na.rm = TRUE),
        n_fits_with_treedepth_hits = sum(treedepth_hits > 0, na.rm = TRUE),
        total_treedepth_hits = sum(treedepth_hits, na.rm = TRUE),
        max_divergences_one_fit = max(divergent, na.rm = TRUE),
        max_treedepth_hits_one_fit = max(treedepth_hits, na.rm = TRUE)
      )
    
    readr::write_csv(
      diagnostics_summary,
      file.path(table_dir, "diagnostics_summary.csv")
    )
    
    readr::write_csv(
      diagnostics,
      file.path(table_dir, "diagnostics_detail.csv")
    )
    
    diagnostics_gt <- diagnostics_summary %>%
      gt() %>%
      tab_header(
        title = "Model Diagnostics Across Imputed Datasets"
      )
    
    save_gt_table(diagnostics_gt, "diagnostics_summary")
    
    diagnostics_sentence <- glue(
      "Across {diagnostics_summary$n_fits} fitted imputed datasets, ",
      "there were {diagnostics_summary$total_divergences} total divergent transitions ",
      "and {diagnostics_summary$total_treedepth_hits} total transitions reaching the maximum tree depth."
    )
  }
  
  # ------------------------------------------------------------
  # Missing outcome prediction summary, if available
  # ------------------------------------------------------------
  
  missing_y_sentence <- "Posterior predictions for missing outcome rows were summarised when available."
  
  if (!is.null(missing_y_summary)) {
    readr::write_csv(
      missing_y_summary,
      file.path(table_dir, "missing_y_summary_publication.csv")
    )
    
    median_col <- find_col(
      missing_y_summary,
      c("Median", "median", "Mean", "mean")
    )
    
    if (!is.na(median_col)) {
      missing_y_overall <- missing_y_summary %>%
        summarise(
          n_missing_rows = n(),
          predicted_median_mean = mean(.data[[median_col]], na.rm = TRUE),
          predicted_median_sd = sd(.data[[median_col]], na.rm = TRUE),
          predicted_median_min = min(.data[[median_col]], na.rm = TRUE),
          predicted_median_max = max(.data[[median_col]], na.rm = TRUE)
        )
      
      readr::write_csv(
        missing_y_overall,
        file.path(table_dir, "missing_y_overall_summary.csv")
      )
      
      missing_y_gt <- missing_y_overall %>%
        gt() %>%
        tab_header(
          title = "Posterior Prediction Summary for Missing Outcome Rows"
        ) %>%
        fmt_number(
          columns = where(is.numeric),
          decimals = 3
        )
      
      save_gt_table(missing_y_gt, "missing_y_overall_summary")
      
      missing_y_sentence <- glue(
        "Posterior predictions were summarised for {missing_y_overall$n_missing_rows} rows with missing outcomes."
      )
    }
  }
  
  # ------------------------------------------------------------
  # Report template
  # ------------------------------------------------------------
  
  # Use paths relative to the report .qmd location:
  # results/publication/report/
  #   ../tables/
  #   ../figures/
  #   ../mi_stability/tables/
  #   ../mi_stability/figures/
  #   ../mo_effects/tables/
  #   ../mo_effects/figures/
  table_dir_rel <- "../tables"
  figure_dir_rel <- "../figures"
  mi_stability_table_dir_rel <- "../mi_stability/tables"
  mi_stability_figure_dir_rel <- "../mi_stability/figures"
  mo_effects_table_dir_rel <- "../mo_effects/tables"
  mo_effects_figure_dir_rel <- "../mo_effects/figures"
  
  family_text <- paste0(
    analysis_spec$outcome$family,
    "(",
    analysis_spec$outcome$link,
    ")"
  )
  

  # ------------------------------------------------------------
  # Analysis metadata for the report
  # ------------------------------------------------------------

  fmt_scalar <- function(x, default = "not specified") {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(default)
    }

    paste(as.character(x), collapse = ", ")
  }

  fmt_prior_text <- function(prior_obj) {
    # Use the compact reader-facing formatter defined above.
    # Avoid capture.output(print(...)), because brmsprior objects print as
    # padded tables with many spaces, which looks poor in HTML/DOCX reports.
    out <- get_report_prior_text(prior_obj)

    if (length(out) == 0 || !nzchar(out)) {
      "not available"
    } else {
      out
    }
  }

  formula_text <- if (!is.null(model_spec) && !is.null(model_spec$formula_text)) {
    compact_report_text(model_spec$formula_text)
  } else if (!is.null(model_spec) && !is.null(model_spec$formula)) {
    compact_report_text(get_report_formula_text(model_spec$formula))
  } else {
    "not available"
  }

  family_name <- analysis_spec$outcome$family %||% "not specified"
  link_name <- analysis_spec$outcome$link %||% "not specified"

  family_text <- paste0(
    family_name,
    "(",
    link_name,
    ")"
  )

  family_sentence <- dplyr::case_when(
    family_name == "bernoulli" && link_name == "logit" ~ "Bernoulli regression model with a logit link",
    family_name == "gaussian" && link_name == "identity" ~ "Gaussian regression model with an identity link",
    family_name == "poisson" && link_name == "log" ~ "Poisson regression model with a log link",
    family_name == "negbinomial" && link_name == "log" ~ "negative-binomial regression model with a log link",
    TRUE ~ paste0(family_name, " regression model with a ", link_name, " link")
  )

  prior_text <- if (!is.null(model_spec) && !is.null(model_spec$prior)) {
    fmt_prior_text(model_spec$prior)
  } else {
    "not available"
  }

  outcome_var <- analysis_spec$outcome$y_var %||% "not available"
  id_var <- analysis_spec$data$id_var %||% "not available"
  time_var <- analysis_spec$data$time_var %||% "not available"
  data_structure <- analysis_spec$data$data_structure %||% "not available"
  imputation_strategy <- analysis_spec$imputation$strategy %||% "not available"

  n_imputations_target <- analysis_spec$imputation$m %||% NA
  imputation_maxiter <- analysis_spec$imputation$maxiter %||% NA
  mean_match_k <- analysis_spec$imputation$mean_match_k %||% NA

  n_fitted_models <- if (exists("n_parameter_imputations") && length(n_parameter_imputations) == 1 && !is.na(n_parameter_imputations)) {
    n_parameter_imputations
  } else if (!is.null(diagnostics)) {
    nrow(diagnostics)
  } else {
    NA_integer_
  }

  mcmc_chains <- analysis_spec$model$chains %||% NA
  mcmc_iter <- analysis_spec$model$iter %||% NA
  mcmc_warmup <- analysis_spec$model$warmup %||% NA

  mcmc_sampling <- if (
    is.numeric(mcmc_iter) &&
      is.numeric(mcmc_warmup) &&
      is.finite(mcmc_iter) &&
      is.finite(mcmc_warmup)
  ) {
    mcmc_iter - mcmc_warmup
  } else {
    NA
  }

  mcmc_seed <- analysis_spec$model$seed %||% NA
  adapt_delta <- analysis_spec$model$adapt_delta %||% NA
  max_treedepth <- analysis_spec$model$max_treedepth %||% NA
  run_smoke_fit <- analysis_spec$model$run_smoke_fit %||% NA

  num_impute_threads <- analysis_spec$parallel$num_impute_threads %||% NA
  fit_workers <- analysis_spec$parallel$fit_workers %||% NA
  cores_per_fit <- analysis_spec$parallel$cores_per_fit %||% NA
  future_max_gb <- analysis_spec$parallel$future_globals_maxsize_gb %||% NA

  summary_centrality <- analysis_spec$summary$centrality %||% NA
  summary_ci <- analysis_spec$summary$ci %||% NA
  summary_ci_method <- analysis_spec$summary$ci_method %||% NA

  summary_test <- if (!is.null(analysis_spec$summary$test)) {
    paste(analysis_spec$summary$test, collapse = ", ")
  } else {
    "not specified"
  }

  summary_rope <- if (!is.null(resolve_rope_range(analysis_spec$summary))) {
    paste(resolve_rope_range(analysis_spec$summary), collapse = ", ")
  } else {
    "not specified"
  }

  predictive_draws_text <- analysis_spec$posterior_prediction$ndraws %||% "not specified"

  missing_y_rows_text <- if (!is.null(missing_y_summary)) {
    nrow(missing_y_summary)
  } else {
    0
  }

  analysis_metadata <- tibble::tibble(
    Item = c(
      "Analysis ID",
      "Project label",
      "Data structure",
      "Outcome variable",
      "Subject ID variable",
      "Time variable",
      "Imputation strategy",
      "Model formula",
      "Model family/link",
      "Priors",
      "Target number of imputations",
      "Successfully fitted imputed datasets used in posterior summaries",
      "Imputation iterations",
      "Mean matching candidates",
      "MCMC chains",
      "Total iterations per chain",
      "Warm-up iterations per chain",
      "Post-warm-up sampling iterations per chain",
      "Seed",
      "adapt_delta",
      "max_treedepth",
      "Smoke fit before parallel fitting",
      "Imputation threads",
      "Parallel fit workers",
      "Cores per fit",
      "future.globals.maxSize, GB",
      "Posterior summary centrality",
      "Credible interval",
      "Credible interval method",
      "Posterior tests",
      "ROPE range",
      "Predictive draws for missing outcomes",
      "Rows with missing outcome prediction summaries",
      "Conditional-effect plotting",
      "Conditional-effect plot variables"
    ),
    Value = c(
      fmt_scalar(analysis_spec$analysis_id),
      fmt_scalar(analysis_spec$project_label),
      fmt_scalar(data_structure),
      fmt_scalar(outcome_var),
      fmt_scalar(id_var),
      fmt_scalar(time_var),
      fmt_scalar(imputation_strategy),
      fmt_scalar(formula_text),
      fmt_scalar(family_text),
      fmt_scalar(prior_text),
      fmt_scalar(n_imputations_target),
      fmt_scalar(n_fitted_models),
      fmt_scalar(imputation_maxiter),
      fmt_scalar(mean_match_k),
      fmt_scalar(mcmc_chains),
      fmt_scalar(mcmc_iter),
      fmt_scalar(mcmc_warmup),
      fmt_scalar(mcmc_sampling),
      fmt_scalar(mcmc_seed),
      fmt_scalar(adapt_delta),
      fmt_scalar(max_treedepth),
      fmt_scalar(run_smoke_fit),
      fmt_scalar(num_impute_threads),
      fmt_scalar(fit_workers),
      fmt_scalar(cores_per_fit),
      fmt_scalar(future_max_gb),
      fmt_scalar(summary_centrality),
      fmt_scalar(summary_ci),
      fmt_scalar(summary_ci_method),
      fmt_scalar(summary_test),
      fmt_scalar(summary_rope),
      fmt_scalar(predictive_draws_text),
      fmt_scalar(missing_y_rows_text),
      fmt_scalar(analysis_spec$reporting$conditional_effects$enabled %||% FALSE),
      fmt_scalar(analysis_spec$reporting$conditional_effects$effects %||% "not specified")
    )
  )

  readr::write_csv(
    analysis_metadata,
    file.path(table_dir, "analysis_metadata.csv")
  )

  saveRDS(
    analysis_metadata,
    file.path(table_dir, "analysis_metadata.rds"),
    compress = FALSE
  )

  methods_sentence <- glue(
    "The target number of imputations was {fmt_scalar(n_imputations_target)}. ",
    "Posterior summaries were based on {fmt_scalar(n_fitted_models)} successfully fitted imputed datasets. ",
    "For each fitted model, we used {fmt_scalar(mcmc_chains)} chain(s), ",
    "{fmt_scalar(mcmc_iter)} total iterations per chain, ",
    "including {fmt_scalar(mcmc_warmup)} warm-up iterations and ",
    "{fmt_scalar(mcmc_sampling)} post-warm-up sampling iterations. ",
    "The main sampler-control settings were adapt_delta = {fmt_scalar(adapt_delta)} ",
    "and max_treedepth = {fmt_scalar(max_treedepth)}."
  )

  special_parameter_report_lines <- c(
    "# Supplementary smooth and monotonic parameters",
    "",
    special_terms_sentence,
    "",
    "```{r special-parameter-table, eval=file.exists(file.path(table_dir, 'special_parameter_table.csv'))}",
    "special_parameter_table <- readr::read_csv(file.path(table_dir, 'special_parameter_table.csv'), show_col_types = FALSE)",
    "gt(special_parameter_table)",
    "```"
  )

  conditional_effect_report_lines <- c(
    "# Conditional-effect plots",
    "",
    conditional_effects_sentence,
    "",
    "```{r conditional-effect-manifest, eval=file.exists(file.path(table_dir, 'conditional_effects_manifest.csv'))}",
    "conditional_effects_manifest <- readr::read_csv(file.path(table_dir, 'conditional_effects_manifest.csv'), show_col_types = FALSE)",
    "gt(conditional_effects_manifest)",
    "```",
    "",
    "```{r conditional-effect-images, results='asis', eval=file.exists(file.path(table_dir, 'conditional_effects_manifest.csv'))}",
    "conditional_effects_manifest <- readr::read_csv(file.path(table_dir, 'conditional_effects_manifest.csv'), show_col_types = FALSE)",
    "created_effects <- conditional_effects_manifest |> dplyr::filter(status == 'created')",
    "for (ii in seq_len(nrow(created_effects))) {",
    "  img <- file.path(figure_dir, created_effects$png_file[ii])",
    "  if (file.exists(img)) {",
    "    cat('## Conditional effect: `', created_effects$effect[ii], '`\\n\\n', sep = '')",
    "    cat('![](', img, '){fig-alt=\\\"Conditional effect plot\\\"}\\n\\n', sep = '')",
    "  }",
    "}",
    "```"
  )

  # This chapter pulls in 11_check_imputation_stability.R's outputs by
  # relative path, rather than copying/re-rendering them. Step 11 normally
  # runs after this script in run_all.R, so these files do not exist yet
  # when this .qmd is written; the eval=file.exists(...) checks below are
  # evaluated at render time, once run_all.R renders this report after
  # Step 11 has completed, not when this script runs.
  mi_stability_report_lines <- c(
    "# Imputation-count stability",
    "",
    "This chapter reports whether posterior summaries were stable as the number of imputations increased, generated by 11_check_imputation_stability.R. It assesses numerical Monte Carlo stability as m increases; it does not validate the missing-data mechanism.",
    "",
    "```{r mi-stability-settings, eval=file.exists(file.path(mi_stability_table_dir, 'imputation_stability_settings.csv'))}",
    "mi_stability_settings <- readr::read_csv(file.path(mi_stability_table_dir, 'imputation_stability_settings.csv'), show_col_types = FALSE)",
    "gt(mi_stability_settings)",
    "```",
    "",
    "```{r mi-stability-final-comparison, eval=file.exists(file.path(mi_stability_table_dir, 'imputation_stability_final_comparison_display.csv'))}",
    "mi_stability_final <- readr::read_csv(file.path(mi_stability_table_dir, 'imputation_stability_final_comparison_display.csv'), show_col_types = FALSE)",
    "gt(mi_stability_final)",
    "```",
    "",
    "```{r mi-stability-stepwise-summary, eval=file.exists(file.path(mi_stability_table_dir, 'imputation_stability_stepwise_summary_display.csv'))}",
    "mi_stability_stepwise <- readr::read_csv(file.path(mi_stability_table_dir, 'imputation_stability_stepwise_summary_display.csv'), show_col_types = FALSE)",
    "gt(mi_stability_stepwise)",
    "```",
    "",
    "```{r mi-stability-trajectory-plot, fig.width=10, fig.height=8, eval=file.exists(file.path(mi_stability_figure_dir, 'imputation_stability_trajectories.png'))}",
    "knitr::include_graphics(file.path(mi_stability_figure_dir, 'imputation_stability_trajectories.png'))",
    "```",
    "",
    "```{r mi-stability-stepwise-plot, fig.width=10, fig.height=6, eval=file.exists(file.path(mi_stability_figure_dir, 'imputation_stability_stepwise_change.png'))}",
    "knitr::include_graphics(file.path(mi_stability_figure_dir, 'imputation_stability_stepwise_change.png'))",
    "```",
    "",
    "Full numeric outputs, including every evaluated batch, are available in `results/publication/mi_stability/tables/`."
  )

  # This chapter pulls in 09/10's mo() outputs by relative path, the same
  # way the imputation-stability chapter above pulls in Step 11's outputs.
  # run_all.R runs 09/10 before this script only when the model formula
  # contains mo() terms, so these files may not exist for ordinary models;
  # the eval=file.exists(...) checks below handle that gracefully.
  mo_effects_report_lines <- c(
    "# Monotonic-effect (mo()) results",
    "",
    "This chapter reports category-specific odds ratios derived from monotonic-effect (`mo()`) parameters, generated by `09_check_mo_parameter_columns.R` and `10_publication_mo_results.R`. It is included only for models whose formula contains `mo()` terms.",
    "",
    "```{r mo-effects-cumulative-table, eval=file.exists(file.path(mo_effects_table_dir, 'mo_cumulative_or_table.csv'))}",
    "mo_cumulative_or_table <- readr::read_csv(file.path(mo_effects_table_dir, 'mo_cumulative_or_table.csv'), show_col_types = FALSE)",
    "gt(mo_cumulative_or_table)",
    "```",
    "",
    "```{r mo-effects-adjacent-table, eval=file.exists(file.path(mo_effects_table_dir, 'mo_adjacent_or_table.csv'))}",
    "mo_adjacent_or_table <- readr::read_csv(file.path(mo_effects_table_dir, 'mo_adjacent_or_table.csv'), show_col_types = FALSE)",
    "gt(mo_adjacent_or_table)",
    "```",
    "",
    "```{r mo-effects-simplex-table, eval=file.exists(file.path(mo_effects_table_dir, 'mo_simplex_table.csv'))}",
    "mo_simplex_table <- readr::read_csv(file.path(mo_effects_table_dir, 'mo_simplex_table.csv'), show_col_types = FALSE)",
    "gt(mo_simplex_table)",
    "```",
    "",
    "```{r mo-effects-cumulative-plot, fig.width=12, fig.height=7, eval=file.exists(file.path(mo_effects_figure_dir, 'mo_cumulative_or_plot.png'))}",
    "knitr::include_graphics(file.path(mo_effects_figure_dir, 'mo_cumulative_or_plot.png'))",
    "```",
    "",
    "```{r mo-effects-adjacent-plot, fig.width=12, fig.height=7, eval=file.exists(file.path(mo_effects_figure_dir, 'mo_adjacent_or_plot.png'))}",
    "knitr::include_graphics(file.path(mo_effects_figure_dir, 'mo_adjacent_or_plot.png'))",
    "```",
    "",
    "Full numeric outputs, including average odds ratios and the parameter-column mapping, are available in `results/publication/mo_effects/tables/`."
  )

  missing_y_report_lines <- if (file.exists(file.path(table_dir, "missing_y_overall_summary.csv"))) {
    c(
      "# Missing outcome prediction",
      "",
      missing_y_sentence,
      "",
      "```{r missing-y-summary}",
      "missing_y_overall <- readr::read_csv(file.path(table_dir, 'missing_y_overall_summary.csv'), show_col_types = FALSE)",
      "gt(missing_y_overall)",
      "```"
    )
  } else {
    c(
      "# Missing outcome prediction",
      "",
      missing_y_sentence
    )
  }


  report_lines <- c(
    "---",
    'title: "Bayesian Multiple-Imputation Analysis Report"',
    "format:",
    "  html:",
    "    toc: true",
    "    toc-depth: 3",
    "    number-sections: true",
    "    embed-resources: true",
    "  docx:",
    "    toc: true",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "```{r setup}",
    "library(tidyverse)",
    "library(gt)",
    "library(knitr)",
    glue('table_dir <- "{table_dir_rel}"'),
    glue('figure_dir <- "{figure_dir_rel}"'),
    glue('mi_stability_table_dir <- "{mi_stability_table_dir_rel}"'),
    glue('mi_stability_figure_dir <- "{mi_stability_figure_dir_rel}"'),
    glue('mo_effects_table_dir <- "{mo_effects_table_dir_rel}"'),
    glue('mo_effects_figure_dir <- "{mo_effects_figure_dir_rel}"'),
    "```",
    "",
    "# Overview",
    "",
    glue("This report summarises a Bayesian {family_sentence} using multiple imputation."),
    "",
    "# Model and computational settings",
    "",
    methods_sentence,
    "",
    "The model formula used in the analysis was:",
    "",
    "```text",
    formula_text,
    "```",
    "",
    "The brms family/link specification was:",
    "",
    "```text",
    family_text,
    "```",
    "",
    "The priors used in the fitted model were:",
    "",
    "```text",
    prior_text,
    "```",
    "",
    "The table below records the main analysis, imputation, modelling, parallelisation, and posterior-summary settings used to generate this report.",
    "",
    "```{r analysis-metadata-table}",
    "analysis_metadata <- readr::read_csv(",
    "  file.path(table_dir, 'analysis_metadata.csv'),",
    "  show_col_types = FALSE",
    ")",
    "",
    "gt(analysis_metadata)",
    "```",
    "",
    "# Diagnostics",
    "",
    diagnostics_sentence,
    "",
    "```{r diagnostics-table, eval=file.exists(file.path(table_dir, 'diagnostics_summary.csv'))}",
    "diagnostics_summary <- readr::read_csv(file.path(table_dir, 'diagnostics_summary.csv'), show_col_types = FALSE)",
    "gt(diagnostics_summary)",
    "```",
    "",
    "# Posterior results",
    "",
    "```{r main-effect-table}",
    "main_effect_table <- readr::read_csv(file.path(table_dir, 'main_effect_table_display.csv'), show_col_types = FALSE)",
    "gt(main_effect_table)",
    "```",
    "",
    special_parameter_report_lines,
    "",
    "# Forest plot",
    "",
    "```{r forest-plot, fig.width=8, fig.height=8, eval=file.exists(file.path(figure_dir, 'forest_plot_fixed_effects.png'))}",
    "knitr::include_graphics(file.path(figure_dir, 'forest_plot_fixed_effects.png'))",
    "```",
    "",
    conditional_effect_report_lines,
    "",
    missing_y_report_lines,
    "",
    mi_stability_report_lines,
    "",
    mo_effects_report_lines
  )
  
  report_file <- file.path(report_dir, "bayesian_mi_report_template.qmd")
  
  writeLines(
    report_lines,
    report_file
  )
  
  # ------------------------------------------------------------
  # README for publication outputs
  # ------------------------------------------------------------
  
  readme_lines <- c(
    "# Publication Output Files",
    "",
    glue("Generated on: {Sys.time()}"),
    "",
    "## Tables",
    "",
    "- `tables/main_effect_table_display.csv`",
    "- `tables/main_effect_table_full.csv`",
    "- `tables/main_effect_table.html`",
    "- `tables/main_effect_table.docx`",
    "- `tables/special_parameter_table.csv`, if smooth or monotonic terms are present",
    "- `tables/conditional_effects_manifest.csv`, if conditional-effect plots are created",
    "- `tables/diagnostics_summary.csv`, if diagnostics are available",
    "- `tables/analysis_metadata.csv`",
    "- `tables/analysis_metadata.rds`",
    "",
    "## Figures",
    "",
    "- `figures/forest_plot_fixed_effects.png`",
    "- `figures/forest_plot_fixed_effects.pdf`",
    "- `figures/conditional_effect_*.png`, if conditional-effect plots are created",
    "- `figures/conditional_effect_*.pdf`, if conditional-effect plots are created",
    "",
    "## Report",
    "",
    "- `report/bayesian_mi_report_template.qmd`"
  )
  
  writeLines(
    readme_lines,
    file.path(pub_dir, "README_publication_outputs.md")
  )

  # ------------------------------------------------------------
  # Render the report
  # ------------------------------------------------------------
  #
  # The report's "Imputation-count stability" chapter embeds files from
  # Step 11 (results/publication/mi_stability/...). run_all.R runs Step 11
  # before this script for that reason, so those files already exist by
  # the time we render here. If this script is instead run on its own
  # (e.g. while iterating on Step 8 changes), the embedded chapter renders
  # with whatever Step 11 output already happens to be on disk, or is
  # silently omitted if Step 11 has not been run yet.
  render_quarto <- analysis_spec$publication$render_quarto %||% TRUE

  if (isTRUE(render_quarto)) {
    quarto_bin <- Sys.which("quarto")

    if (nzchar(quarto_bin)) {
      log_msg("Rendering main report.")
      html_status <- system2(quarto_bin, args = c("render", report_file, "--to", "html"))
      docx_status <- system2(quarto_bin, args = c("render", report_file, "--to", "docx"))
      log_msg("Main report HTML render status:", html_status)
      log_msg("Main report DOCX render status:", docx_status)
    } else {
      log_msg("Quarto not found on PATH; wrote", report_file, "without rendering.")
    }
  }

  log_msg("Publication outputs created in:", pub_dir)
  
  gc()
  
  guard_memory("after STEP 8 cleanup")
}, analysis_spec)