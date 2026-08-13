# =============================================================================
# Censored-exposure block-FCS imputation  (PLAN §5 / Phase 4)
# -----------------------------------------------------------------------------
# A congenial imputation path for the special case where a focal EXPOSURE
# (predictor of the outcome) is left-/interval-censored below a reporting limit.
# Imputing such an exposure without the outcome Y biases the exposure-response
# coefficient (PLAN §2); this path handles it with a two-engine block-FCS:
#
#   * Z block (MAR covariates)      -> miceRanger, via the EXISTING pipeline
#                                       helper run_row_level_imputation() (reused,
#                                       not reimplemented), conditioning on the
#                                       current X and Y.
#   * X block (censored exposure)   -> leftcens::impute_censored_conditional()
#                                       (>= 0.9.0): skew-aware, bound-respecting,
#                                       Y-aware, proper-MI conditional draw.
#   * alternate for a few outer sweeps, per completed dataset.
#
# It is opt-in (strategy = "censored_exposure_block_fcs") and returns an
# `imputed_list` in the SAME shape as the other strategies, so 03_impute.R writes
# the standard imputed_###.rds + manifest and Steps 4-12 are untouched.
#
# Input convention (leftcens interval columns): each censored exposure `X` carries
# per-row bounds `X_lo` / `X_hi` (suffixes configurable) of the interval holding
# the true value: X_lo == X_hi where exactly observed; X_lo = -Inf (or 0) with
# X_hi = LOD for a left-censored non-detect; both finite for an interval.
# =============================================================================

#' Build (y, lower, upper) for one exposure from its lo/hi columns, on the
#' modelling scale (log if the variable's dictionary `scale` is "log").
.ce_exposure_bounds <- function(data, exposure, lo_col, hi_col, log_scale) {
  lo <- data[[lo_col]]; hi <- data[[hi_col]]
  if (is.null(lo) || is.null(hi)) {
    stop("Censored exposure '", exposure, "' needs columns '", lo_col,
         "' and '", hi_col, "' in the data.", call. = FALSE)
  }
  if (log_scale) {                       # impute on the log scale (where shash lives)
    lo <- ifelse(lo <= 0, -Inf, log(lo)) # 0/negative lower bound -> left-censored
    hi <- ifelse(is.finite(hi), log(hi), Inf)
  }
  exact <- is.finite(lo) & is.finite(hi) & abs(hi - lo) < 1e-9
  list(
    y     = ifelse(exact, lo, NA_real_),               # observed value where exact
    lower = ifelse(exact, NA_real_, lo),
    upper = ifelse(exact, NA_real_, hi),
    log_scale = log_scale
  )
}

#' Initialise a censored exposure column: observed value where exact, else a
#' point inside the interval (bounded away from an infinite end).
.ce_init_exposure <- function(b) {
  y <- b$y
  cens <- is.na(y)
  lo <- b$lower[cens]; hi <- b$upper[cens]
  fill <- ifelse(is.finite(lo) & is.finite(hi), (lo + hi) / 2,
          ifelse(is.finite(hi), hi - 0.5,                 # left-censored: just below LOD
          ifelse(is.finite(lo), lo + 0.5, 0)))            # right-censored
  y[cens] <- fill
  y
}

#' One block-FCS completed dataset (a single imputation).
.ce_one_imputation <- function(data, analysis_spec, var_dict, ce, seed) {
  set.seed(seed)
  y_var <- analysis_spec$outcome$y_var
  exposures <- ce$exposure_vars
  lo_suffix <- ce$lo_suffix %||% "_lo"
  hi_suffix <- ce$hi_suffix %||% "_hi"
  sweeps <- as.integer(ce$outer_sweeps %||% 5L)
  margin <- ce$margin %||% "shash"

  # Whether to impute on the log scale: an explicit config flag wins; otherwise
  # fall back to the dictionary `scale == "log"`. Keeping this separate from the
  # dictionary scale avoids double-transforming when the model formula already
  # applies log() to the exposure.
  log_for <- function(x) {
    if (!is.null(ce$log_scale)) return(isTRUE(ce$log_scale))
    identical(var_dict$scale[match(x, var_dict$var)], "log")
  }

  # Per-exposure bounds on the modelling scale, and the working column filled in.
  bounds <- list(); work <- data; drop_cols <- character(0)
  for (x in exposures) {
    lo_col <- paste0(x, lo_suffix); hi_col <- paste0(x, hi_suffix)
    b <- .ce_exposure_bounds(data, x, lo_col, hi_col, log_scale = log_for(x))
    bounds[[x]] <- b
    work[[x]] <- .ce_init_exposure(b)                 # complete predictor for the Z block
    drop_cols <- c(drop_cols, lo_col, hi_col)
  }
  # The lo/hi columns have served their purpose; remove them so they never enter
  # the Z-block predictors or the completed dataset (downstream uses the exposure).
  work <- work[, setdiff(names(work), drop_cols), drop = FALSE]

  # Predictor set for the X block: reduced set if given, else auto (outcome +
  # in-model covariates + the OTHER exposures), matching the plan's guidance to
  # avoid conditioning on all ~300 Z.
  auto_preds <- unique(c(
    y_var,
    var_dict$var[var_dict$use_in_model %in% TRUE],
    exposures))
  x_predictors <- ce$predictors %||% setdiff(auto_preds, character(0))

  # Impute Y in the Z block too (MID): the X-block conditions on Y, so Y must be
  # complete for that step; the imputed-Y rows are deleted afterwards. Track the
  # ORIGINAL missing cells so every outer sweep re-imputes them given the current
  # X (miceRanger only targets currently-NA cells — without this the block would
  # stop alternating after the first sweep).
  z_spec_as <- analysis_spec
  z_spec_as$imputation$impute_y <- TRUE
  # Reset only the originally-missing Z/Y cells each sweep — NOT the exposures,
  # whose missing (censored) cells are handled by the X block via the fixed bounds.
  orig_na <- lapply(data[intersect(names(data), names(work))], is.na)
  reset_cols <- setdiff(names(orig_na)[vapply(orig_na, any, logical(1))], exposures)

  for (t in seq_len(sweeps)) {
    # Re-open the originally-missing Z/Y cells so the block re-imputes them.
    for (cc in reset_cols) work[[cc]][orig_na[[cc]]] <- NA
    # --- Z block: miceRanger (reused helper), m = 1, conditioning on current X, Y.
    z_spec <- make_row_level_imputation_spec(work, z_spec_as, var_dict)
    if (length(z_spec$vars) > 0) {
      z_spec$m <- 1L; z_spec$seed <- seed + 1000L + t
      work <- run_row_level_imputation(work, z_spec, analysis_spec)[[1]]
    }
    # --- X block: one congenial censored draw per exposure, given current Z, Y.
    for (x in exposures) {
      b <- bounds[[x]]
      preds <- setdiff(intersect(x_predictors, names(work)), x)
      xmat <- work[, preds, drop = FALSE]
      imp <- leftcens::impute_censored_conditional(
        y = b$y, x = xmat, lower = b$lower, upper = b$upper,
        m = 1L, margin = margin)[, 1]
      work[[x]] <- if (isTRUE(b$log_scale)) exp(imp) else imp   # back to data scale
    }
  }
  work
}

#' Run the censored-exposure block-FCS and return an `imputed_list` (m completed
#' datasets), the same shape 03_impute.R expects from the other strategies.
#'
#' @param data,analysis_spec,var_dict As in the other Step-3 strategies.
#' @param m Number of completed datasets to produce (the new batch size).
#' @param seed Base seed for this batch.
run_censored_exposure_block_fcs <- function(data, analysis_spec, var_dict, m, seed) {
  if (!requireNamespace("leftcens", quietly = TRUE) ||
      !"impute_censored_conditional" %in% getNamespaceExports("leftcens")) {
    stop("strategy 'censored_exposure_block_fcs' requires leftcens >= 0.9.0 ",
         "(exporting impute_censored_conditional()).", call. = FALSE)
  }
  ce <- analysis_spec$imputation$censored_exposure
  if (is.null(ce) || length(ce$exposure_vars) == 0) {
    stop("strategy 'censored_exposure_block_fcs' needs ",
         "analysis_spec$imputation$censored_exposure$exposure_vars.", call. = FALSE)
  }
  y_var <- analysis_spec$outcome$y_var
  mid <- isTRUE(ce$mid_delete_imputed_y %||% TRUE)
  y_missing <- if (!is.null(data[[y_var]])) is.na(data[[y_var]]) else logical(nrow(data))

  log_msg("Censored-exposure block-FCS | exposures:",
          paste(ce$exposure_vars, collapse = ", "),
          "| outer_sweeps:", ce$outer_sweeps %||% 5L,
          "| margin:", ce$margin %||% "shash",
          "| m:", m)

  imputed_list <- vector("list", m)
  for (d in seq_len(m)) {
    work <- .ce_one_imputation(data, analysis_spec, var_dict, ce, seed = seed + d)
    # MID: the imputed-Y rows carry no information for beta; drop them before the
    # brms fit so pooling is over information-bearing rows only (von Hippel 2007).
    if (mid && any(y_missing)) work <- work[!y_missing, , drop = FALSE]
    imputed_list[[d]] <- tibble::as_tibble(work)
    log_msg("  completed dataset", d, "of", m,
            if (mid && any(y_missing)) paste0("(MID dropped ", sum(y_missing), " imputed-Y rows)") else "")
  }
  imputed_list
}
