#!/usr/bin/env Rscript
# ============================================================
# 10_publication_mo_results.R
#
# Publication-ready summaries for brms monotonic effects mo().
#
# v4 fix:
#   - Separates main mo() simplex parameters from time-interaction mo()
#     simplex parameters.
#   - Correctly handles parameter names like:
#       bsp_moC6yincomeidEQC6yincome
#       bsp_moC6yincomeidEQC6yincome:time
#       simo_moC6yincomeidEQC6yincome1[1]
#       simo_moC6yincomeidEQC6yincome:time1[1]
#
# For a model with:
#   time * mo(x)
#
# Adjacent interval j log OR at time t is:
#   main_b * D * main_zeta_j +
#   t * interaction_b * D * interaction_zeta_j
#
# It is NOT:
#   (main_b + t * interaction_b) * D * main_zeta_j
#
# Run:
#   Rscript 10_publication_mo_results.R
#
# Optional:
#   quarto render results/publication/mo_effects/report/mo_effects_report.qmd
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("00_config.R")
source("00_common_functions.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ------------------------------------------------------------
# Settings: discovered from the model formula, overridable via
# analysis_spec$mo_effects in 00_config.R (see 09_check_mo_parameter_columns.R
# for a ready-to-paste config skeleton with the right number of levels).
# ------------------------------------------------------------

parameter_draws_file <- file.path(paths$results, "parameter_draws.rds")
model_spec_file <- file.path(paths$objects, "model_spec.rds")

output_root <- file.path(paths$publication, "mo_effects")
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
report_dir <- file.path(output_root, "report")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(model_spec_file)) {
  stop("File not found: ", model_spec_file, ". Run Step 4 first.")
}
if (!file.exists(parameter_draws_file)) {
  stop("File not found: ", parameter_draws_file, ". Run Step 6 first.")
}

model_spec <- readRDS(model_spec_file)
mo_effects_cfg <- analysis_spec$mo_effects %||% list()

detected_vars <- extract_special_term_vars(model_spec$formula, fun = "mo")

if (length(detected_vars) == 0) {
  message("No mo() terms found in the fitted model's formula. Nothing to do.")
  quit(save = "no", status = 0)
}

# mo_vars: named list of var -> list(label, levels). Falls back to generic
# placeholder labels/levels for any detected variable not configured in
# analysis_spec$mo_effects$vars, so this still runs without manual config,
# just with less informative labels.
configured_vars <- mo_effects_cfg$vars %||% list()

mo_vars <- purrr::map(
  rlang::set_names(detected_vars),
  function(var) {
    cfg <- configured_vars[[var]]

    if (!is.null(cfg)) {
      return(list(label = cfg$label %||% var, levels = cfg$levels))
    }

    NULL # levels resolved later, once the simplex dimension is known
  }
)

# If your model contains time * mo(variable), ORs are calculated at these
# time values. Configure via analysis_spec$mo_effects$time_var/time_values.
time_var <- mo_effects_cfg$time_var %||% NULL
time_values <- mo_effects_cfg$time_values %||% NULL

ci_prob <- analysis_spec$summary$ci %||% 0.95
alpha <- (1 - ci_prob) / 2

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

fmt_ci <- function(low, high, digits = 2) {
  paste0(fmt_num(low, digits), " to ", fmt_num(high, digits))
}

escape_regex <- function(x) {
  gsub("([.\\^$|()\\[\\]{}*+?])", "\\\\\\1", x, perl = TRUE)
}

contains_fixed <- function(x, pattern) {
  stringr::str_detect(x, stringr::fixed(pattern))
}

has_interaction_marker <- function(x) {
  stringr::str_detect(x, ":")
}

extract_simplex_index <- function(x) {
  idx <- stringr::str_match(x, "\\[(\\d+)\\]")[, 2]
  as.integer(idx)
}

select_one <- function(x, label) {
  if (length(x) != 1) {
    stop(
      "Expected exactly one ", label, ", but found ", length(x), ":\n",
      paste(x, collapse = "\n")
    )
  }
  x[[1]]
}

find_mo_parameter_cols <- function(draws, var, time_var = NULL) {
  nms <- names(draws)
  coefficient_cols_all <- nms[stringr::str_detect(nms, "^(b_|bsp_)")]
  simo_cols_all <- nms[stringr::str_detect(nms, "^simo_")]

  # Anchor on the exact brms-mangled token "mo<var>" right after the
  # b_/bsp_/simo_ prefix, followed only by an interaction marker, a
  # simplex-index digit, or the end of the name. This avoids the
  # unanchored substring matching that could previously bind to an
  # unrelated column whose name happened to contain `var` as a substring
  # (e.g. "medu" incorrectly matching "medu2").
  #
  # When mo() is called with an explicit id = "..." argument (used to
  # share a simplex across multiple mo() terms), brms inserts an
  # "idEQ<id>" infix right after "mo<var>", e.g. bsp_momeduidEQmedu or
  # simo_momeduidEQmedu1[1]. The optional group below accounts for that,
  # whether or not id was used.
  token_re <- escape_regex(paste0("mo", var))
  id_infix_re <- "(idEQ[A-Za-z0-9._]+)?"

  related_coef <- coefficient_cols_all[
    stringr::str_detect(coefficient_cols_all, paste0("^(b_|bsp_)", token_re, id_infix_re, "($|:)"))
  ]

  related_simo <- simo_cols_all[
    stringr::str_detect(simo_cols_all, paste0("^simo_", token_re, id_infix_re, "\\d*(:|\\[)"))
  ]

  main_b_candidates <- related_coef[!has_interaction_marker(related_coef)]
  time_b_candidates <- if (is.null(time_var)) {
    character(0)
  } else {
    related_coef[
      has_interaction_marker(related_coef) &
        contains_fixed(related_coef, time_var)
    ]
  }

  main_simo_candidates <- related_simo[!has_interaction_marker(related_simo)]
  time_simo_candidates <- if (is.null(time_var)) {
    character(0)
  } else {
    related_simo[
      has_interaction_marker(related_simo) &
        contains_fixed(related_simo, time_var)
    ]
  }

  if (length(main_b_candidates) != 1) {
    stop(
      "Expected exactly one main coefficient column for monotonic variable ",
      var, ", but found ", length(main_b_candidates), ".\n\n",
      "Candidates found:\n",
      paste(main_b_candidates, collapse = "\n"),
      "\n\nAvailable related coefficient columns are:\n",
      paste(related_coef, collapse = "\n")
    )
  }

  if (length(main_simo_candidates) == 0) {
    stop(
      "No main simo_ columns found for monotonic variable ", var, ".\n\n",
      "Available related simo_ columns are:\n",
      paste(related_simo, collapse = "\n")
    )
  }

  main_idx <- extract_simplex_index(main_simo_candidates)

  if (anyNA(main_idx)) {
    stop("Could not parse main simplex indices from:\n", paste(main_simo_candidates, collapse = "\n"))
  }

  main_simo_candidates <- main_simo_candidates[order(main_idx)]

  has_time_interaction <- length(time_b_candidates) > 0 || length(time_simo_candidates) > 0

  if (has_time_interaction) {
    if (length(time_b_candidates) != 1) {
      stop(
        "Expected exactly one time-interaction coefficient column for monotonic variable ",
        var, ", but found ", length(time_b_candidates), ".\n\n",
        "Candidates found:\n",
        paste(time_b_candidates, collapse = "\n"),
        "\n\nAvailable related coefficient columns are:\n",
        paste(related_coef, collapse = "\n")
      )
    }

    if (length(time_simo_candidates) == 0) {
      stop(
        "A time-interaction coefficient was found for ", var,
        ", but no time-interaction simo_ columns were found.\n\n",
        "Available related simo_ columns are:\n",
        paste(related_simo, collapse = "\n")
      )
    }

    time_idx <- extract_simplex_index(time_simo_candidates)

    if (anyNA(time_idx)) {
      stop("Could not parse interaction simplex indices from:\n", paste(time_simo_candidates, collapse = "\n"))
    }

    time_simo_candidates <- time_simo_candidates[order(time_idx)]

    if (length(time_simo_candidates) != length(main_simo_candidates)) {
      stop(
        "Main and time-interaction simplex lengths differ for ", var, ".\n",
        "Main length: ", length(main_simo_candidates), "\n",
        "Time-interaction length: ", length(time_simo_candidates)
      )
    }
  } else {
    time_b_candidates <- NA_character_
    time_simo_candidates <- character(0)
  }

  list(
    var = var,
    b_main_col = main_b_candidates[[1]],
    b_time_col = if (has_time_interaction) time_b_candidates[[1]] else NA_character_,
    simo_main_cols = main_simo_candidates,
    simo_time_cols = time_simo_candidates,
    has_time_interaction = has_time_interaction
  )
}

compute_interval_log_or <- function(b, zeta) {
  D <- ncol(zeta)
  sweep(zeta * D, 1, b, `*`)
}

compute_mo_or_draws <- function(draws, var, levels, label = var,
                                time_var = NULL, time_values = NULL) {
  cols <- find_mo_parameter_cols(draws, var = var, time_var = time_var)

  b_main <- draws[[cols$b_main_col]]
  zeta_main <- as.matrix(draws[, cols$simo_main_cols, drop = FALSE])

  D <- ncol(zeta_main)
  K <- D + 1L

  if (is.null(levels)) {
    levels <- paste0("Level ", seq_len(K))
  }

  if (length(levels) != K) {
    stop(
      "Variable ", var, " has ", length(levels), " supplied levels, ",
      "but the main simplex has ", D, " intervals, implying ", K, " levels."
    )
  }

  if (cols$has_time_interaction) {
    b_time <- draws[[cols$b_time_col]]
    zeta_time <- as.matrix(draws[, cols$simo_time_cols, drop = FALSE])
  } else {
    b_time <- rep(0, length(b_main))
    zeta_time <- matrix(0, nrow = nrow(zeta_main), ncol = ncol(zeta_main))
  }

  if (!cols$has_time_interaction || is.null(time_values)) {
    effect_grid <- tibble(effect_context = "overall", time_value = NA_real_)
  } else {
    effect_grid <- tibble(
      effect_context = paste0(time_var, "=", time_values),
      time_value = as.numeric(time_values)
    )
  }

  adjacent_contrast <- paste0(levels[1:(K - 1L)], "_to_", levels[2:K])
  adjacent_label <- paste0(levels[1:(K - 1L)], " to ", levels[2:K])
  cumulative_contrast <- paste0(levels[2:K], "_vs_", levels[1])
  cumulative_label <- paste0(levels[2:K], " vs ", levels[1])

  main_interval_log_or <- compute_interval_log_or(b_main, zeta_main)
  time_interval_log_or <- compute_interval_log_or(b_time, zeta_time)

  out_or <- vector("list", nrow(effect_grid))

  for (ii in seq_len(nrow(effect_grid))) {
    context_i <- effect_grid$effect_context[[ii]]
    time_i <- effect_grid$time_value[[ii]]

    log_or_adj <- if (is.na(time_i)) {
      main_interval_log_or
    } else {
      main_interval_log_or + time_i * time_interval_log_or
    }

    adjacent_draws <- as_tibble(log_or_adj, .name_repair = "minimal") %>%
      setNames(adjacent_contrast) %>%
      mutate(.draw_id = row_number()) %>%
      pivot_longer(
        cols = -.draw_id,
        names_to = "contrast",
        values_to = "log_or"
      ) %>%
      mutate(
        variable = var,
        variable_label = label,
        effect_context = context_i,
        time_value = time_i,
        effect_type = "adjacent_increment",
        contrast_label = adjacent_label[match(contrast, adjacent_contrast)],
        or = exp(log_or),
        .before = 1
      )

    log_or_cum <- t(apply(log_or_adj, 1, cumsum))

    cumulative_draws <- as_tibble(log_or_cum, .name_repair = "minimal") %>%
      setNames(cumulative_contrast) %>%
      mutate(.draw_id = row_number()) %>%
      pivot_longer(
        cols = -.draw_id,
        names_to = "contrast",
        values_to = "log_or"
      ) %>%
      mutate(
        variable = var,
        variable_label = label,
        effect_context = context_i,
        time_value = time_i,
        effect_type = "cumulative_vs_lowest",
        contrast_label = cumulative_label[match(contrast, cumulative_contrast)],
        or = exp(log_or),
        .before = 1
      )

    # Average adjacent effect is the average of the interval-specific log ORs.
    avg_log_or <- rowMeans(log_or_adj)

    overall_draws <- tibble(
      variable = var,
      variable_label = label,
      effect_context = context_i,
      time_value = time_i,
      effect_type = "average_adjacent",
      contrast = "average_adjacent_category_increase",
      contrast_label = "Average adjacent category increase",
      .draw_id = seq_along(avg_log_or),
      log_or = avg_log_or,
      or = exp(avg_log_or)
    )

    out_or[[ii]] <- bind_rows(overall_draws, adjacent_draws, cumulative_draws)
  }

  # Save main and interaction simplex proportions separately.
  simplex_main <- as_tibble(zeta_main, .name_repair = "minimal") %>%
    setNames(adjacent_contrast) %>%
    mutate(.draw_id = row_number()) %>%
    pivot_longer(
      cols = -.draw_id,
      names_to = "contrast",
      values_to = "simplex_proportion"
    ) %>%
    mutate(
      variable = var,
      variable_label = label,
      effect_context = "main",
      time_value = NA_real_,
      effect_type = "simplex_proportion",
      contrast_label = adjacent_label[match(contrast, adjacent_contrast)],
      .before = 1
    )

  simplex_list <- list(simplex_main)

  if (cols$has_time_interaction) {
    simplex_time <- as_tibble(zeta_time, .name_repair = "minimal") %>%
      setNames(adjacent_contrast) %>%
      mutate(.draw_id = row_number()) %>%
      pivot_longer(
        cols = -.draw_id,
        names_to = "contrast",
        values_to = "simplex_proportion"
      ) %>%
      mutate(
        variable = var,
        variable_label = label,
        effect_context = paste0(time_var, "_interaction"),
        time_value = NA_real_,
        effect_type = "simplex_proportion",
        contrast_label = adjacent_label[match(contrast, adjacent_contrast)],
        .before = 1
      )

    simplex_list <- c(simplex_list, list(simplex_time))
  }

  list(
    column_map = tibble(
      variable = var,
      variable_label = label,
      b_main_col = cols$b_main_col,
      b_time_col = cols$b_time_col,
      simo_main_cols = paste(cols$simo_main_cols, collapse = "; "),
      simo_time_cols = paste(cols$simo_time_cols, collapse = "; ")
    ),
    or_draws = bind_rows(out_or),
    simplex_draws = bind_rows(simplex_list)
  )
}

summarise_or_draws <- function(df) {
  df %>%
    group_by(variable, variable_label, effect_context, time_value,
             effect_type, contrast, contrast_label) %>%
    summarise(
      n_draws = n(),
      log_or_median = median(log_or, na.rm = TRUE),
      log_or_lower = quantile(log_or, alpha, na.rm = TRUE),
      log_or_upper = quantile(log_or, 1 - alpha, na.rm = TRUE),
      OR_median = median(or, na.rm = TRUE),
      OR_lower = quantile(or, alpha, na.rm = TRUE),
      OR_upper = quantile(or, 1 - alpha, na.rm = TRUE),
      Pr_OR_gt_1 = mean(or > 1, na.rm = TRUE),
      Pr_OR_lt_1 = mean(or < 1, na.rm = TRUE),
      p_direction = pmax(Pr_OR_gt_1, Pr_OR_lt_1),
      .groups = "drop"
    )
}

summarise_simplex_draws <- function(df) {
  df %>%
    group_by(variable, variable_label, effect_context, time_value,
             effect_type, contrast, contrast_label) %>%
    summarise(
      n_draws = n(),
      simplex_median = median(simplex_proportion, na.rm = TRUE),
      simplex_lower = quantile(simplex_proportion, alpha, na.rm = TRUE),
      simplex_upper = quantile(simplex_proportion, 1 - alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

make_or_display_table <- function(summary_df, effect_type_filter) {
  summary_df %>%
    filter(effect_type == effect_type_filter) %>%
    transmute(
      Variable = variable_label,
      Context = effect_context,
      Contrast = contrast_label,
      `OR median` = fmt_num(OR_median, 2),
      `95% CrI` = fmt_ci(OR_lower, OR_upper, 2),
      `Pr(OR > 1)` = fmt_num(Pr_OR_gt_1, 3),
      `Pr(OR < 1)` = fmt_num(Pr_OR_lt_1, 3),
      p_direction = fmt_num(p_direction, 3)
    )
}

make_simplex_display_table <- function(summary_df) {
  summary_df %>%
    transmute(
      Variable = variable_label,
      Context = effect_context,
      Interval = contrast_label,
      `Simplex median` = fmt_num(simplex_median, 3),
      `95% CrI` = fmt_ci(simplex_lower, simplex_upper, 3)
    )
}

save_gt_if_available <- function(df, file_stem, title, subtitle = NULL) {
  if (!requireNamespace("gt", quietly = TRUE)) return(invisible(NULL))
  gt_tbl <- gt::gt(df) %>% gt::tab_header(title = title, subtitle = subtitle)
  gt::gtsave(gt_tbl, filename = file.path(table_dir, paste0(file_stem, ".html")))
  invisible(gt_tbl)
}

save_docx_if_available <- function(tables) {
  if (!requireNamespace("flextable", quietly = TRUE) ||
      !requireNamespace("officer", quietly = TRUE)) {
    return(invisible(NULL))
  }
  doc <- officer::read_docx()
  for (nm in names(tables)) {
    doc <- doc %>%
      officer::body_add_par(nm, style = "heading 1") %>%
      flextable::body_add_flextable(
        flextable::autofit(flextable::flextable(tables[[nm]]))
      ) %>%
      officer::body_add_par("")
  }
  print(doc, target = file.path(output_root, "mo_effect_publication_tables.docx"))
  invisible(NULL)
}

plot_or_summary <- function(summary_df, effect_type_filter, file_stem, title) {
  plot_df <- summary_df %>%
    filter(effect_type == effect_type_filter) %>%
    mutate(contrast_label = factor(contrast_label, levels = rev(unique(contrast_label))))

  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(
    plot_df,
    aes(x = OR_median, y = contrast_label, xmin = OR_lower, xmax = OR_upper)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_pointrange() +
    scale_x_log10() +
    facet_grid(variable_label ~ effect_context, scales = "free_y", space = "free_y") +
    labs(title = title, x = "Odds ratio, log scale", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8)
    )

  ggsave(
    filename = file.path(figure_dir, paste0(file_stem, ".png")),
    plot = p,
    width = 12,
    height = 7,
    dpi = 300
  )

  invisible(p)
}

write_quarto_report <- function() {
  qmd_file <- file.path(report_dir, "mo_effects_report.qmd")
  report_lines <- c(
    "---",
    'title: "Monotonic Ordinal Effect Summaries"',
    "format:",
    "  html:",
    "    toc: true",
    "    toc-depth: 3",
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
    'table_dir <- "../tables"',
    'figure_dir <- "../figures"',
    "```",
    "",
    "# Overview",
    "",
    "This report summarises monotonic ordinal effects from a Bayesian logistic mixed model.",
    "",
    "If a time-by-monotonic-predictor interaction was detected, odds ratios are shown at each configured time value.",
    "",
    "For models with `time * mo(x)`, the interval-specific log odds ratio is calculated as the main monotonic contribution plus the time-specific interaction monotonic contribution.",
    "",
    "# Parameter columns used",
    "",
    "```{r column-map}",
    "column_map <- readr::read_csv(file.path(table_dir, 'mo_parameter_column_map.csv'), show_col_types = FALSE)",
    "gt(column_map)",
    "```",
    "",
    "# Recommended main table: cumulative ORs versus the lowest category",
    "",
    "```{r cumulative-table}",
    "cumulative_or_table <- readr::read_csv(file.path(table_dir, 'mo_cumulative_or_table.csv'), show_col_types = FALSE)",
    "gt(cumulative_or_table)",
    "```",
    "",
    "```{r cumulative-plot, fig.width=12, fig.height=7, eval=file.exists(file.path(figure_dir, 'mo_cumulative_or_plot.png'))}",
    "knitr::include_graphics(file.path(figure_dir, 'mo_cumulative_or_plot.png'))",
    "```",
    "",
    "# Supplementary table: adjacent increment ORs",
    "",
    "```{r adjacent-table}",
    "adjacent_or_table <- readr::read_csv(file.path(table_dir, 'mo_adjacent_or_table.csv'), show_col_types = FALSE)",
    "gt(adjacent_or_table)",
    "```",
    "",
    "```{r adjacent-plot, fig.width=12, fig.height=7, eval=file.exists(file.path(figure_dir, 'mo_adjacent_or_plot.png'))}",
    "knitr::include_graphics(file.path(figure_dir, 'mo_adjacent_or_plot.png'))",
    "```",
    "",
    "# Average adjacent-category effect",
    "",
    "```{r average-table}",
    "average_or_table <- readr::read_csv(file.path(table_dir, 'mo_average_or_table.csv'), show_col_types = FALSE)",
    "gt(average_or_table)",
    "```",
    "",
    "# Supplementary simplex proportions",
    "",
    "```{r simplex-table}",
    "simplex_table <- readr::read_csv(file.path(table_dir, 'mo_simplex_table.csv'), show_col_types = FALSE)",
    "gt(simplex_table)",
    "```",
    "",
    "# Interpretation note",
    "",
    "These are conditional odds ratios on the model's linear predictor scale. They should not be interpreted as risk ratios or marginal probability differences."
  )
  writeLines(report_lines, qmd_file)
  invisible(qmd_file)
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

if (!file.exists(parameter_draws_file)) {
  stop("File not found: ", parameter_draws_file)
}

draws <- readRDS(parameter_draws_file)

mo_results <- purrr::imap(
  mo_vars,
  function(spec, var) {
    compute_mo_or_draws(
      draws = draws,
      var = var,
      levels = spec$levels,
      label = spec$label %||% var,
      time_var = time_var,
      time_values = time_values
    )
  }
)

mo_column_map <- purrr::map_dfr(mo_results, "column_map")
mo_or_draws <- purrr::map_dfr(mo_results, "or_draws")
mo_simplex_draws <- purrr::map_dfr(mo_results, "simplex_draws")

mo_or_summary <- summarise_or_draws(mo_or_draws)
mo_simplex_summary <- summarise_simplex_draws(mo_simplex_draws)

saveRDS(mo_or_draws, file.path(paths$results, "mo_or_draws.rds"), compress = FALSE)
saveRDS(mo_or_summary, file.path(paths$results, "mo_or_summary.rds"), compress = FALSE)
readr::write_csv(mo_or_summary, file.path(paths$results, "mo_or_summary.csv"))

saveRDS(mo_simplex_draws, file.path(paths$results, "mo_simplex_draws.rds"), compress = FALSE)
saveRDS(mo_simplex_summary, file.path(paths$results, "mo_simplex_summary.rds"), compress = FALSE)
readr::write_csv(mo_simplex_summary, file.path(paths$results, "mo_simplex_summary.csv"))

average_or_table <- make_or_display_table(mo_or_summary, "average_adjacent")
adjacent_or_table <- make_or_display_table(mo_or_summary, "adjacent_increment")
cumulative_or_table <- make_or_display_table(mo_or_summary, "cumulative_vs_lowest")
simplex_table <- make_simplex_display_table(mo_simplex_summary)

readr::write_csv(mo_column_map, file.path(table_dir, "mo_parameter_column_map.csv"))
readr::write_csv(average_or_table, file.path(table_dir, "mo_average_or_table.csv"))
readr::write_csv(adjacent_or_table, file.path(table_dir, "mo_adjacent_or_table.csv"))
readr::write_csv(cumulative_or_table, file.path(table_dir, "mo_cumulative_or_table.csv"))
readr::write_csv(simplex_table, file.path(table_dir, "mo_simplex_table.csv"))

saveRDS(mo_column_map, file.path(table_dir, "mo_parameter_column_map.rds"), compress = FALSE)
saveRDS(average_or_table, file.path(table_dir, "mo_average_or_table.rds"), compress = FALSE)
saveRDS(adjacent_or_table, file.path(table_dir, "mo_adjacent_or_table.rds"), compress = FALSE)
saveRDS(cumulative_or_table, file.path(table_dir, "mo_cumulative_or_table.rds"), compress = FALSE)
saveRDS(simplex_table, file.path(table_dir, "mo_simplex_table.rds"), compress = FALSE)

save_gt_if_available(cumulative_or_table, "mo_cumulative_or_table", "Cumulative monotonic-effect odds ratios", "Compared with the lowest category")
save_gt_if_available(adjacent_or_table, "mo_adjacent_or_table", "Adjacent monotonic-effect odds ratios", "Category-to-category increments")
save_gt_if_available(average_or_table, "mo_average_or_table", "Average adjacent-category monotonic-effect odds ratios")
save_gt_if_available(simplex_table, "mo_simplex_table", "Supplementary simplex proportions", "Proportion of the monotonic effect assigned to each adjacent interval")

save_docx_if_available(
  list(
    "Cumulative ORs versus the lowest category" = cumulative_or_table,
    "Adjacent increment ORs" = adjacent_or_table,
    "Average adjacent-category ORs" = average_or_table,
    "Supplementary simplex proportions" = simplex_table,
    "Parameter columns used" = mo_column_map
  )
)

plot_or_summary(mo_or_summary, "cumulative_vs_lowest", "mo_cumulative_or_plot", "Cumulative monotonic-effect odds ratios")
plot_or_summary(mo_or_summary, "adjacent_increment", "mo_adjacent_or_plot", "Adjacent monotonic-effect odds ratios")

qmd_file <- write_quarto_report()

message("Saved monotonic-effect publication outputs:")
message("  ", table_dir)
message("  ", figure_dir)
message("  ", report_dir)
message("")
message("Main table:")
message("  ", file.path(table_dir, "mo_cumulative_or_table.csv"))
message("")
message("Render report with:")
message("  quarto render ", qmd_file)
