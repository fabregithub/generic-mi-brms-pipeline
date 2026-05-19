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
})

init_logging("pipeline")
setup_project_dirs(paths)

safe_step("STEP 8: Publication-ready outputs", {
  
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
  # Display table
  # ------------------------------------------------------------
  
  display_base <- main_effect_table %>%
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
  
  main_effect_gt <- main_effect_table_display %>%
    gt() %>%
    tab_header(
      title = "Posterior Summary of Fixed Effects"
    ) %>%
    cols_align(
      align = "center",
      columns = -Parameter
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
  
  missing_y_sentence <- "Posterior predictions for missing outcome rows were summarized when available."
  
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
        "Posterior predictions were summarized for {missing_y_overall$n_missing_rows} rows with missing outcomes."
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
  table_dir_rel <- "../tables"
  figure_dir_rel <- "../figures"
  
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
    if (is.null(prior_obj)) {
      return("not available")
    }

    out <- tryCatch(
      paste(capture.output(print(prior_obj)), collapse = "; "),
      error = function(e) paste(as.character(prior_obj), collapse = "; ")
    )

    if (length(out) == 0 || !nzchar(out)) {
      "not available"
    } else {
      out
    }
  }

  formula_text <- if (exists("model_spec") && !is.null(model_spec$formula)) {
    paste(deparse(model_spec$formula), collapse = " ")
  } else {
    "not available"
  }

  family_text <- if (exists("model_spec") && !is.null(model_spec$family)) {
    paste(capture.output(print(model_spec$family)), collapse = " ")
  } else {
    paste0(
      analysis_spec$outcome$family,
      "(",
      analysis_spec$outcome$link,
      ")"
    )
  }

  prior_text <- if (exists("model_spec") && !is.null(model_spec$prior)) {
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

  n_fitted_models <- if (!is.na(n_parameter_imputations)) {
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

  summary_rope <- if (!is.null(analysis_spec$summary$rope$fixed_range)) {
    paste(analysis_spec$summary$rope$fixed_range, collapse = ", ")
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
      "Rows with missing outcome prediction summaries"
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
      fmt_scalar(missing_y_rows_text)
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
    "```",
    "",
    "# Overview",
    "",
    glue("This report summarizes a Bayesian `{family_text}` regression analysis using multiple imputation."),
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
    "# Forest plot",
    "",
    "```{r forest-plot, fig.width=8, fig.height=8, eval=file.exists(file.path(figure_dir, 'forest_plot_fixed_effects.png'))}",
    "knitr::include_graphics(file.path(figure_dir, 'forest_plot_fixed_effects.png'))",
    "```",
    "",
    "# Missing outcome prediction",
    "",
    missing_y_sentence,
    "",
    "```{r missing-y-summary, eval=file.exists(file.path(table_dir, 'missing_y_overall_summary.csv'))}",
    "missing_y_overall <- readr::read_csv(file.path(table_dir, 'missing_y_overall_summary.csv'), show_col_types = FALSE)",
    "gt(missing_y_overall)",
    "```"
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
    "- `tables/diagnostics_summary.csv`, if diagnostics are available",
    "- `tables/analysis_metadata.csv`",
    "- `tables/analysis_metadata.rds`",
    "",
    "## Figures",
    "",
    "- `figures/forest_plot_fixed_effects.png`",
    "- `figures/forest_plot_fixed_effects.pdf`",
    "",
    "## Report",
    "",
    "- `report/bayesian_mi_report_template.qmd`"
  )
  
  writeLines(
    readme_lines,
    file.path(pub_dir, "README_publication_outputs.md")
  )
  
  log_msg("Publication outputs created in:", pub_dir)
  
  gc()
  
  guard_memory("after STEP 8 cleanup")
}, analysis_spec)