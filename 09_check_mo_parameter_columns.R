#!/usr/bin/env Rscript
# ============================================================
# 09_check_mo_parameter_columns.R
#
# Discover monotonic-effect (mo()) parameter columns generically, from
# the fitted model's own formula, rather than from hardcoded variable
# names. Works for any analysis_spec$model$custom_formula containing
# mo() terms.
#
# Run after Step 6 has created results/parameter_draws.rds:
#   Rscript 09_check_mo_parameter_columns.R
# ============================================================

source("00_config.R")
source("00_common_functions.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

escape_regex <- function(x) {
  gsub("([.\\^$|()\\[\\]{}*+?])", "\\\\\\1", x, perl = TRUE)
}

draws_file <- file.path(paths$results, "parameter_draws.rds")
model_spec_file <- file.path(paths$objects, "model_spec.rds")

if (!file.exists(draws_file)) {
  stop("File not found: ", draws_file, ". Run Step 6 first.")
}
if (!file.exists(model_spec_file)) {
  stop("File not found: ", model_spec_file, ". Run Step 4 first.")
}

draws <- readRDS(draws_file)
model_spec <- readRDS(model_spec_file)
nms <- names(draws)

cat("Model formula:", get_brms_formula_text(model_spec$formula), "\n\n")

mo_vars <- extract_special_term_vars(model_spec$formula, fun = "mo")

if (length(mo_vars) == 0) {
  cat("No mo() terms found in the fitted model's formula. Nothing to check.\n")
  quit(save = "no", status = 0)
}

cat("mo() variables found in the model formula:\n")
print(mo_vars)
cat("\nmo()-related parameter columns, by variable:\n\n")

mo_config_lines <- c(
  "analysis_spec$mo_effects <- list(",
  "  vars = list("
)

for (var in mo_vars) {
  # When mo() is called with an explicit id = "..." argument, brms inserts
  # an "idEQ<id>" infix right after "mo<var>" in the parameter name (e.g.
  # bsp_momeduidEQmedu, simo_momeduidEQmedu1[1]). The optional group below
  # matches both with and without an explicit id.
  token_re <- escape_regex(paste0("mo", var))
  id_infix_re <- "(idEQ[A-Za-z0-9._]+)?"

  bsp_cols <- nms[stringr::str_detect(nms, paste0("^bsp_", token_re, id_infix_re, "($|:)"))]
  simo_cols <- nms[stringr::str_detect(nms, paste0("^simo_", token_re, id_infix_re, "\\d*(:|\\[)"))]

  has_interaction <- any(stringr::str_detect(c(bsp_cols, simo_cols), ":"))
  main_simo_cols <- simo_cols[!stringr::str_detect(simo_cols, ":")]
  n_levels <- length(main_simo_cols) + 1L

  cat(sprintf("- %s\n", var))
  cat("    bsp_ columns:  ", paste(bsp_cols, collapse = ", "), "\n")
  cat("    simo_ columns: ", paste(simo_cols, collapse = ", "), "\n")
  cat("    implied number of ordered categories (simplex dim + 1):", n_levels, "\n")
  cat("    time-interaction columns detected:", has_interaction, "\n\n")

  if (length(bsp_cols) == 0 || length(main_simo_cols) == 0) {
    cat(
      "    WARNING: no matching bsp_/simo_ columns found for '", var, "'. ",
      "Check analysis_spec$model$parameter_draw_regex includes bsp_/simo_.\n\n",
      sep = ""
    )
    next
  }

  placeholder_levels <- paste0('"Level ', seq_len(n_levels), '"', collapse = ", ")
  mo_config_lines <- c(
    mo_config_lines,
    sprintf(
      '    %s = list(label = "%s", levels = c(%s)),',
      var, var, placeholder_levels
    )
  )
}

# Drop the trailing comma after the last vars entry so the printed
# snippet is valid R when pasted verbatim.
n_lines <- length(mo_config_lines)
mo_config_lines[n_lines] <- sub(",$", "", mo_config_lines[n_lines])

mo_config_lines <- c(
  mo_config_lines,
  "  ),",
  "  time_var = NULL,    # set this if any variable above showed a time interaction",
  "  time_values = NULL  # e.g. 1:6",
  ")"
)

cat("Paste the following into 00_config.R (edit labels/levels/time settings as needed),\n")
cat("then run 10_publication_mo_results.R:\n\n")
cat(paste(mo_config_lines, collapse = "\n"), "\n")
