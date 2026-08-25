# =============================================================================
# Phase 1 harness -- Track V2 robustness axes
# -----------------------------------------------------------------------------
# Track V2 of ../PLAN_pipeline_validation.md §7.
#
# V1 validated the shipped engine in ONE corner of the design space: only X1
# censored, covariates complete, outcome complete, rho = 0.4, skew = 0, n = 800.
# Two consequences noted in FINDINGS_v1.md: MID never fired, and the miceRanger
# Z block was a no-op -- so the block-FCS *alternation itself* was never actually
# exercised under Monte Carlo.
#
# This file adds the missing injectors and defines the scenario grid that drives
# each axis, including the two that make the Z block and MID do real work.
# =============================================================================

#' MCAR a fraction of the OUTCOME (new in V2).
#'
#' The only way to exercise multiple-imputation-then-deletion: the block-FCS
#' imputes Y so it is complete for the X-block predictor step, then deletes those
#' rows before the analysis fit (von Hippel 2007). With Y complete -- as in V1 --
#' `mid_delete_imputed_y` is dead code.
#'
#' @param data A data.frame, typically after [inject_left_censoring()].
#' @param y_frac Fraction of outcome values set missing.
#' @param y_col Outcome column name.
inject_mcar_outcome <- function(data, y_frac = 0.0, y_col = "Y") {
  if (y_frac <= 0) return(data)
  keep <- attributes(data)[c("lods", "nd_frac_realised", "censored_cols")]
  n <- nrow(data)
  idx <- sample.int(n, size = floor(y_frac * n))
  data[[y_col]][idx] <- NA_real_
  for (nm in names(keep)) if (!is.null(keep[[nm]])) attr(data, nm) <- keep[[nm]]
  attr(data, "y_missing_frac") <- y_frac
  data
}

# `inject_mcar_covariates()` (censoring.R) drops the bundle attributes the
# procedures rely on; wrap it so they survive.
.v2_mcar_covariates <- function(data, mcar_frac, cols = c("Z1", "Z2")) {
  if (mcar_frac <= 0) return(data)
  keep <- attributes(data)[c("lods", "nd_frac_realised", "censored_cols")]
  out <- inject_mcar_covariates(data, mcar_frac = mcar_frac, cols = cols)
  for (nm in names(keep)) if (!is.null(keep[[nm]])) attr(out, nm) <- keep[[nm]]
  out
}

# ---- scenario grid -----------------------------------------------------------

#' One scenario = one point in the robustness design space.
#'
#' Fields default to the V1 configuration, so `base` reproduces the V1 additive
#' 40%-ND cell and acts as a calibration check that the V2 machinery did not
#' change the answer.
.v2_scenario <- function(name, axis, ..., erf_form = "additive", nd_frac = 0.4,
                         n = 800L, m = 30L, n_rep = NULL, rho = 0.4, skew = 0.0,
                         censor_all = FALSE, mcar_frac = 0.0, y_frac = 0.0) {
  list(name = name, axis = axis, erf_form = erf_form, nd_frac = nd_frac,
       n = as.integer(n), m = as.integer(m), n_rep = n_rep, rho = rho,
       skew = skew, censor_all = censor_all, mcar_frac = mcar_frac,
       y_frac = y_frac)
}

#' The V2 scenario grid (PLAN §7).
#'
#' @param big_n Sample size for the large-n cell.
#' @param big_n_rep Replications for the large-n cell (it is the expensive one).
v2_scenarios <- function(big_n = 20000L, big_n_rep = 50L) {
  list(
    # Calibration: must reproduce V1's additive ND-0.4 row.
    .v2_scenario("base",            "reference"),

    # Axis 1 -- every exposure censored, not just the focal one. Exercises the
    # X-block's loop over exposures under MC for the first time.
    .v2_scenario("censor_all",      "censoring",  censor_all = TRUE),

    # Axis 2 -- MCAR covariates. Makes the miceRanger Z block do real work, so
    # the block-FCS alternation is finally exercised.
    .v2_scenario("mcar_z20",        "covariates", mcar_frac = 0.2),
    .v2_scenario("mcar_z40",        "covariates", mcar_frac = 0.4),

    # Axis 3 -- missing outcome. The only path that fires MID.
    .v2_scenario("missing_y20",     "outcome",    y_frac = 0.2),

    # The realistic combination: censored exposure + MAR covariates + missing Y.
    .v2_scenario("combined",        "combined",   censor_all = TRUE,
                                                  mcar_frac = 0.2, y_frac = 0.2),

    # Axis 4 -- exposure correlation. High rho is where a reduced X-block
    # predictor set could bite.
    .v2_scenario("rho00",           "correlation", rho = 0.0),
    .v2_scenario("rho80",           "correlation", rho = 0.8),

    # Axis 5 -- skew. Swept for the prototype in FINDINGS.md; confirm the
    # shipped engine inherits the shash margin's robustness.
    .v2_scenario("skew075",         "skew",       skew = 0.75),

    # Axis 6 -- large n. Confirms the scale test at MC level.
    .v2_scenario("large_n",         "scale",      n = big_n, n_rep = big_n_rep)
  )
}

#' Which procedures are DEFINED for a scenario.
#'
#' Not every comparator generalises, and running one outside its design would
#' silently produce a meaningless row:
#'
#'   * `cens_mi_y_shash` hardcodes X1 as the sole censored variable and passes
#'     the other exposures and covariates as complete predictors. Under
#'     censor-all or MCAR covariates those predictors contain NA, so it is not
#'     defined there. It is the V1 prototype, never designed for this.
#'   * `leftcens_prestep` handles multiple censored analytes natively
#'     (`gsimp_mi` over all bounds) and simply listwise-drops NA covariates at
#'     the analysis step -- which is exactly what the status-quo architecture
#'     would do, so it stays in as the comparator throughout.
#'   * `pipeline_block_fcs` is defined everywhere. That is the point of V2.
v2_procedures_for <- function(sc) {
  procs <- c("oracle", "complete_case", "leftcens_prestep", "pipeline_block_fcs")
  prototype_ok <- !isTRUE(sc$censor_all) && sc$mcar_frac <= 0 && sc$y_frac <= 0
  if (prototype_ok) procs <- c(procs, "cens_mi_y_shash")
  procs
}

#' Build one bundle for a scenario (simulate -> censor -> thin).
v2_make_bundle <- function(sc, truth) {
  comp <- simulate_complete(sc$n, truth, rho = sc$rho, skew = sc$skew)
  p <- length(truth$b)
  which_cens <- if (isTRUE(sc$censor_all)) seq_len(p) else 1L

  cens <- inject_left_censoring(comp$data, nd_frac = sc$nd_frac,
                                censor_which = which_cens)
  cens <- .v2_mcar_covariates(cens, sc$mcar_frac)
  cens <- inject_mcar_outcome(cens, sc$y_frac)

  list(complete = comp$data, censored = cens, truth = truth)
}
