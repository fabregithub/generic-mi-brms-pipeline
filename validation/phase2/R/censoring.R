# =============================================================================
# Phase 2 harness -- missingness injection (birth-weight mixture)     DESIGN §1
# -----------------------------------------------------------------------------
# Multi-analyte, THREE-TIER left/interval censoring (ND / DNQ / quantified) plus
# MAR covariates, per PROVENANCE.md §7. Emits `leftcens` interval columns
# `<analyte>_lo` / `<analyte>_hi` (concentration scale) so the completed-data
# procedures and the production pipeline consume an identical convention:
#
#   quantified :  x >= LOQ        -> [x , x]     (point observed)
#   DNQ        :  MDL <= x < LOQ  -> [MDL, LOQ]  (interval-censored)
#   ND         :  x <  MDL        -> [0 , MDL]   (left-censored)
#
# The point column `<analyte>` is set to NA wherever not quantified (MI sees NA
# + the interval bounds). Per-analyte MDLs are calibrated so realised detection
# matches the PROVENANCE table; `limit_scale > 1` inflates the limits for the
# increased-censoring scenario.
# =============================================================================

# target detection per analyte (PROVENANCE §1); ND fraction = 1 - detection
DEFAULT_DETECTION <- c(pfos = 0.99, pfoa = 0.998, pfna = 0.93, pfda = 0.92,
                       pfunda = 0.93, pfhxs = 0.92, pftrda = 0.785)

#' Inject three-tier censoring on the analyte concentration columns.
#'
#' @param data complete data from [simulate_complete()].
#' @param analytes analyte column names (default: attr(data,"analytes")).
#' @param detection named target *detection* (above-MDL) fraction per analyte;
#'   MDL is set to its (1-detection) quantile, so ND fraction = 1-detection.
#' @param dnq_frac target DNQ slice just above MDL (LOQ at the (1-detection+dnq_frac)
#'   quantile); a NARROW band, so quantified fraction = detection - dnq_frac.
#' @param mdl,loq optional named explicit limits (ng/mL); override the quantile
#'   calibration for those analytes.
#' @param limit_scale multiplies MDL and LOQ (the increased-censoring axis).
#' @return `data` augmented per analyte with `_lo`, `_hi`, `_tier`
#'   ("quant"/"dnq"/"nd"); the point column set to NA where not quantified; and
#'   attributes "mdl","loq","detected" (above MDL),"quant_frac","dnq_frac","nd_frac".
inject_censoring <- function(data, analytes = attr(data, "analytes"),
                             detection = DEFAULT_DETECTION, dnq_frac = 0.05,
                             mdl = NULL, loq = NULL, limit_scale = 1.0) {
  mdls <- loqs <- detected <- q_frac <- d_frac <- n_frac <-
    stats::setNames(numeric(0), character(0))
  for (a in analytes) {
    x <- data[[a]]
    nd_q  <- 1 - detection[[a]]                       # ND fraction (below MDL)
    loq_q <- min(nd_q + dnq_frac, 0.999)              # cumulative up to LOQ
    m  <- (if (!is.null(mdl) && a %in% names(mdl)) mdl[[a]]
           else as.numeric(stats::quantile(x, nd_q,  names = FALSE))) * limit_scale
    lq <- (if (!is.null(loq) && a %in% names(loq)) loq[[a]]
           else as.numeric(stats::quantile(x, loq_q, names = FALSE))) * limit_scale

    tier <- ifelse(x < m, "nd", ifelse(x < lq, "dnq", "quant"))
    data[[paste0(a, "_lo")]]   <- ifelse(tier == "nd", 0, ifelse(tier == "dnq", m,  x))
    data[[paste0(a, "_hi")]]   <- ifelse(tier == "nd", m, ifelse(tier == "dnq", lq, x))
    data[[paste0(a, "_tier")]] <- tier
    data[[a]][tier != "quant"] <- NA_real_            # point NA unless quantified

    mdls[a] <- m; loqs[a] <- lq
    n_frac[a] <- mean(tier == "nd"); d_frac[a] <- mean(tier == "dnq")
    q_frac[a] <- mean(tier == "quant"); detected[a] <- 1 - n_frac[a]
  }
  attr(data, "analytes")    <- analytes
  attr(data, "mdl")         <- mdls;     attr(data, "loq")        <- loqs
  attr(data, "detected")    <- detected; attr(data, "quant_frac") <- q_frac
  attr(data, "dnq_frac")    <- d_frac;   attr(data, "nd_frac")    <- n_frac
  attr(data, "limit_scale") <- limit_scale
  data
}

#' MAR missingness on covariates: P(missing) depends on an observed variable, so
#' complete-case is biased but MI conditioning on the predictor is valid.
#'
#' @param cols covariate columns to thin.
#' @param frac target overall missing fraction per column.
#' @param depends_on an observed column driving the missingness (default "ga").
#' @param strength log-odds slope on the standardised driver.
inject_mar_covariates <- function(data, cols = c("bmi_z","height_z"),
                                  frac = 0.15, depends_on = "ga",
                                  strength = 0.8) {
  cols <- intersect(cols, names(data)); if (!length(cols) || frac <= 0) return(data)
  z <- as.numeric(scale(data[[depends_on]]))
  n <- nrow(data)
  for (cc in cols) {
    a0 <- stats::qlogis(frac)                       # intercept ~ target frac
    p  <- stats::plogis(a0 + strength * z)
    miss <- stats::runif(n) < p
    data[[cc]][miss] <- NA
  }
  attr(data, "mar_frac") <- frac; attr(data, "mar_driver") <- depends_on
  data
}

#' Convenience: report realised censoring for eyeballing against PROVENANCE.
censoring_summary <- function(data) {
  data.frame(
    analyte  = attr(data, "analytes"),
    mdl      = round(attr(data, "mdl"), 4),
    loq      = round(attr(data, "loq"), 4),
    detected = round(attr(data, "detected"), 3),
    nd_frac  = round(attr(data, "nd_frac"), 3),
    dnq_frac = round(attr(data, "dnq_frac"), 3),
    quant    = round(attr(data, "quant_frac"), 3),
    row.names = NULL)
}
