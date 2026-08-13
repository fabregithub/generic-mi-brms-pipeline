# =============================================================================
# Phase 1 harness -- missingness injection (PLAN §7.2)
# -----------------------------------------------------------------------------
# Left-censor the focal exposure(s) below a per-analyte LOD (MNAR-by-censoring),
# and optionally MCAR a fraction of the covariates.
#
# The LOD is set as the `nd_frac` quantile of the *true* log-exposure, so the
# realised non-detect fraction matches `nd_frac` by construction.
# =============================================================================

#' Left-censor exposures below a per-analyte LOD.
#'
#' @param data A complete data.frame from [simulate_complete()].
#' @param nd_frac Target non-detect fraction (per censored analyte).
#' @param censor_which Integer indices of exposures to censor. Default: only the
#'   focal exposure (X1), which isolates the mechanism under test. Set to
#'   `seq_len(p)` to censor every exposure (PLAN §7.2).
#' @return `data` augmented, per censored analyte `j`, with:
#'   * `logXj`      -- true value where observed, NA where censored
#'   * `logXj_lod`  -- the log LOD (also the reported value for censored rows)
#'   * `logXj_cens` -- brms-style indicator: "none" (observed) / "left" (censored)
#'   plus attribute "lods" (named numeric) and "nd_frac_realised".
inject_left_censoring <- function(data, nd_frac = 0.30, censor_which = 1L) {
  xcols <- grep("^logX[0-9]+$", names(data), value = TRUE)
  lods <- numeric(0); nd_real <- numeric(0)

  for (j in censor_which) {
    xc <- xcols[j]
    lod <- as.numeric(stats::quantile(data[[xc]], probs = nd_frac, names = FALSE))
    is_cens <- data[[xc]] < lod

    data[[paste0(xc, "_lod")]]  <- lod
    data[[paste0(xc, "_cens")]] <- ifelse(is_cens, "left", "none")
    # For brms: censored rows carry the LOD as the (left-)censored response value.
    # For imputation bounds we keep the true value separately via the LOD column;
    # the observed `logXj` is set to NA on censored rows for the MI procedures.
    data[[xc]][is_cens] <- NA_real_

    lods[xc] <- lod
    nd_real[xc] <- mean(is_cens)
  }

  attr(data, "lods") <- lods
  attr(data, "nd_frac_realised") <- nd_real
  attr(data, "censored_cols") <- xcols[censor_which]
  data
}

#' MCAR a fraction of the covariates (PLAN §7.2). Default off in the scaffold.
#'
#' @param data A data.frame (typically after [inject_left_censoring()]).
#' @param mcar_frac Fraction of cells set missing in each covariate column.
#' @param cols Covariate columns to thin. Default: Z1, Z2.
inject_mcar_covariates <- function(data, mcar_frac = 0.0, cols = c("Z1", "Z2")) {
  if (mcar_frac <= 0) return(data)
  cols <- intersect(cols, names(data))
  n <- nrow(data)
  for (cc in cols) {
    idx <- sample.int(n, size = floor(mcar_frac * n))
    data[[cc]][idx] <- NA
  }
  attr(data, "mcar_frac") <- mcar_frac
  data
}
