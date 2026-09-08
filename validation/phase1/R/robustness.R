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

#' MAR a fraction of the covariates (new in V10).
#'
#' WHY THIS EXISTS. Every scenario in V1-V9 makes covariates missing completely at
#' random. MCAR is the easy case and the unrealistic one: real covariates go
#' missing for reasons, and those reasons are usually visible in the data. MAR --
#' missingness depending on OBSERVED values, not on the missing values themselves
#' -- is the assumption multiple imputation is actually built to handle, and the
#' one every claim in the root README currently rests on WITHOUT having tested.
#'
#' THE MECHANISM. `P(Z missing)` is logistic in the outcome and in an observed,
#' never-censored exposure:
#'
#'     logit P(missing) = a + strength * z(Y) + strength * z(logX_ref)
#'
#' Both drivers are fully observed, so this is MAR and not MNAR -- the missingness
#' does not depend on the value of `Z` itself. That matters: MNAR would be a
#' different (and much harder) question, and conflating the two would make a null
#' result look like a stronger guarantee than it is.
#'
#' `a` is solved numerically so the realised missing fraction matches `mar_frac`,
#' which keeps this cell comparable to the MCAR cell at the same fraction --
#' otherwise a difference in results could just be a difference in how much was
#' missing.
#'
#' THE DRIVER IS DELIBERATELY THE OUTCOME. The pipeline's Z block conditions on Y,
#' so a Y-driven mechanism is one the imputation model can in principle handle.
#' If it still degrades, the problem is the imputer rather than the mechanism
#' being out of scope.
#'
#' @param data A data.frame, typically after [inject_left_censoring()].
#' @param mar_frac Target overall missing fraction per column.
#' @param cols Covariate columns to thin.
#' @param strength Logit coefficient on each standardised driver. 0 reduces to
#'   MCAR (useful as an internal control).
#' @param ref_x Observed exposure used as the second driver.
inject_mar_covariates <- function(data, mar_frac = 0.0, cols = c("Z1", "Z2"),
                                  strength = 1.0, ref_x = "logX2") {
  if (mar_frac <= 0) return(data)
  cols <- intersect(cols, names(data))
  n <- nrow(data)

  zs <- function(v) {
    v <- as.numeric(v)
    sdv <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / sdv
  }
  drv <- zs(data[["Y"]])
  if (ref_x %in% names(data) && !anyNA(data[[ref_x]])) drv <- drv + zs(data[[ref_x]])
  drv[!is.finite(drv)] <- 0

  # Solve the intercept so the expected missing fraction is exactly mar_frac.
  f <- function(a) mean(stats::plogis(a + strength * drv)) - mar_frac
  a <- tryCatch(stats::uniroot(f, interval = c(-50, 50))$root,
                error = function(e) stats::qlogis(mar_frac))

  pr <- stats::plogis(a + strength * drv)
  for (cc in cols) data[[cc]][stats::runif(n) < pr] <- NA

  attr(data, "mar_frac") <- mar_frac
  attr(data, "mar_strength") <- strength
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

# Same attribute-preserving wrapper for the MAR injector.
.v2_mar_covariates <- function(data, mar_frac, strength = 1.0,
                               cols = c("Z1", "Z2")) {
  if (mar_frac <= 0) return(data)
  keep <- attributes(data)[c("lods", "nd_frac_realised", "censored_cols")]
  out <- inject_mar_covariates(data, mar_frac = mar_frac, cols = cols,
                               strength = strength)
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
                         censor_all = FALSE, mcar_frac = 0.0, y_frac = 0.0,
                         z_form = "linear", mar_frac = 0.0, mar_strength = 1.0,
                         y_form = "linear", z_role = "precision", seed_as = NULL) {
  list(name = name, axis = axis, erf_form = erf_form, nd_frac = nd_frac,
       n = as.integer(n), m = as.integer(m), n_rep = n_rep, rho = rho,
       skew = skew, censor_all = censor_all, mcar_frac = mcar_frac,
       y_frac = y_frac, z_form = z_form,
       mar_frac = mar_frac, mar_strength = mar_strength, y_form = y_form,
       # V17: the covariate's causal role -- see make_truth() in dgp.R.
       z_role = z_role,
       # `seed_as` makes this cell draw its data with ANOTHER cell's seed, so the
       # two see byte-identical exposures, covariates and outcome noise and differ
       # only in the feature under test. That turns a cell-vs-cell contrast from
       # unpaired into paired -- worth roughly an 8x cut in Monte-Carlo error
       # here, which is far cheaper than tripling the replications.
       seed_as = seed_as)
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
    .v2_scenario("large_n",         "scale",      n = big_n, n_rep = big_n_rep),

    # --- Track 05: non-linear covariates -------------------------------------
    # The imputation model for Z1 is now genuinely hard (a linear model leaves
    # ~0.39 of explainable R-squared on the table; see R/dgp.R). These mirror the
    # three diagnostic cells so the linear and non-linear designs are directly
    # comparable, and they are what lets a parametric imputer be PENALISED for
    # misspecification rather than only rewarded for calibration.
    .v2_scenario("nl_mcar_z40",     "nonlinear",  mcar_frac = 0.4, z_form = "nonlinear"),
    .v2_scenario("nl_missing_y20",  "nonlinear",  y_frac = 0.2,    z_form = "nonlinear"),
    .v2_scenario("nl_combined",     "nonlinear",  censor_all = TRUE,
                                                  mcar_frac = 0.2, y_frac = 0.2,
                                                  z_form = "nonlinear")
,
    # ---- Axis 6 (V10): MAR covariates ------------------------------------
    # Matched to the MCAR cells at the same fraction, so a difference between
    # `mcar_z40` and `mar_z40` is the MECHANISM and not the amount missing.
    # `mar_z40_mcarctl` sets strength = 0, which reduces the MAR injector to
    # MCAR: it must reproduce `mcar_z40`, and is the control that proves the new
    # injector did not change anything except the mechanism.
    .v2_scenario("mar_z20",        "covariates", mar_frac = 0.2),
    .v2_scenario("mar_z40",        "covariates", mar_frac = 0.4),
    .v2_scenario("mar_z40_strong", "covariates", mar_frac = 0.4, mar_strength = 2.0),
    .v2_scenario("mar_z40_mcarctl","covariates", mar_frac = 0.4, mar_strength = 0.0),
    .v2_scenario("nl_mar_z40",     "covariates", mar_frac = 0.4, z_form = "nonlinear"),
    # Both mechanisms plus a missing outcome -- the analogue of `combined`.
    .v2_scenario("mar_combined",   "covariates", mar_frac = 0.4, y_frac = 0.2),
    .v2_scenario("nl_mar_combined","covariates", mar_frac = 0.4, y_frac = 0.2,
                 z_form = "nonlinear"),
    # ---- Axis 7 (V11): NON-LINEAR OUTCOME in the covariate ----------------
    # Up to V10 the outcome is linear in (logX, Z) in every cell, so a
    # parametric imputation model for Z1 is correctly specified against it and
    # can never be penalised for misspecification. FINDINGS_v7.md names this as
    # the reason parametric arms do well in the outcome-dominated cells -- which
    # makes it a confound in the v1.5.0 `z_imputer = "bart"` decision.
    # These cells add `b_zq * (Z1^2 - 1)` to the outcome, which makes
    # p(Z1 | Y, X) non-linear. The ANALYSIS model gains the matching I(Z1^2)
    # term, so the focal estimand stays recoverable and any bias is the
    # imputation model's.
    # Each is matched to an existing linear cell for a like-for-like contrast.
    .v2_scenario("ynl_mcar_z40",    "y_nonlinear", mcar_frac = 0.4, y_form = "nonlinear", seed_as = "mcar_z40"),
    .v2_scenario("ynl_missing_y20", "y_nonlinear", y_frac = 0.2,    y_form = "nonlinear", seed_as = "missing_y20"),
    .v2_scenario("ynl_combined",    "y_nonlinear", censor_all = TRUE, mcar_frac = 0.2,
                 y_frac = 0.2, y_form = "nonlinear", seed_as = "combined"),
    .v2_scenario("ynl_mar_z40",     "y_nonlinear", mar_frac = 0.4,  y_form = "nonlinear", seed_as = "mar_z40"),

    # ---- V17: the covariate's CAUSAL ROLE --------------------------------------
    # APPENDED, never inserted: task seeds key on a scenario's position in this
    # list, so adding cells at the end leaves every earlier cell's data
    # bit-identical and every published result reproducible.
    #
    # Tracks V0-V16 all ran with Z as a precision covariate or a descendant of
    # the exposures -- never a cause of one. So no track has ever adjusted for a
    # covariate that HAD to be adjusted for, and none has ever adjusted for one
    # that must NOT be. These four cells supply the missing structures; the roles
    # themselves, and the arithmetic that verifies them, are in dgp.R and
    # test_v17_zrole.R.
    #
    # Missingness is held at mcar_z40's (40% non-detects on X1, 40% MCAR on the
    # covariates) so the four differ ONLY in the causal role of Z.
    .v2_scenario("zr_fork",     "z_role", mcar_frac = 0.4, z_role = "fork"),
    .v2_scenario("zr_pipe",     "z_role", mcar_frac = 0.4, z_role = "pipe"),
    .v2_scenario("zr_collider", "z_role", mcar_frac = 0.4, z_role = "collider"),
    .v2_scenario("zr_mixed",    "z_role", mcar_frac = 0.4, z_role = "mixed"),

    # The same four roles under MAR instead of MCAR. NOT redundant with the pair
    # above: `inject_mar_covariates()` drives missingness off **Y** (and logX2),
    # so omitting Y from the Z block does not merely cost information there -- it
    # drops the variable the mechanism depends on, which is the condition MAR
    # validity rests on. Under MCAR the no-Y penalty is a congeniality cost;
    # under MAR-on-Y it should also be a mechanism failure, and the two should
    # therefore differ.
    #
    # `seed_as` pairs each with its MCAR twin: the complete data and the exposure
    # censoring are drawn identically and only the covariate mechanism differs,
    # so MCAR-vs-MAR is a paired contrast rather than an unpaired one.
    .v2_scenario("zrmar_fork",     "z_role", mar_frac = 0.4, z_role = "fork",     seed_as = "zr_fork"),
    .v2_scenario("zrmar_pipe",     "z_role", mar_frac = 0.4, z_role = "pipe",     seed_as = "zr_pipe"),
    .v2_scenario("zrmar_collider", "z_role", mar_frac = 0.4, z_role = "collider", seed_as = "zr_collider"),
    .v2_scenario("zrmar_mixed",    "z_role", mar_frac = 0.4, z_role = "mixed",    seed_as = "zr_mixed"),

    # ---- V22: THE DECIDING CELL ------------------------------------------------
    # A pipe whose X1 -> Z1 arrow is NON-LINEAR, so Z1 is on an X-Y path AND its
    # conditional has a non-linear mean. V21 showed a correctly specified
    # ESTIMATED covariate draw is unbiased where that conditional is
    # linear-Gaussian -- which every cell before this one is, by construction.
    # This is where a linear parametric draw is genuinely misspecified and BART
    # can learn the shape, so it decides whether V21's "a fix exists" is a fix or
    # a different bug.
    #
    # THE EXPOSURE IS FULLY OBSERVED IN THE FIRST TWO CELLS (nd_frac = 0), and
    # that is not a convenience -- it is required. Measured while building this:
    # with logX1 40% censored, the ladder's exact-Z anchor sits at **+17%**
    # instead of ~0, because leftcens draws logX1 from a conditional LINEAR in
    # Z1 while the truth has Z1 = g(logX1). The X block is badly misspecified
    # here, so it swamps the covariate-draw comparison and the ladder has no
    # zero point. With logX1 observed the same anchor is +1.39% +/- 0.89 --
    # consistent with zero, and the comparison is purely about the Z draw.
    #
    # `zr_pipe_nc` is the matched LINEAR comparator (same missingness, same
    # absence of censoring), so the two differ only in the arrow's shape.
    # `zr_pipenl_cens` keeps the censored version as a SECONDARY cell: it has no
    # valid anchor and cannot be decomposed, but it prices the total and flags
    # the X-block interaction as a finding in its own right.
    .v2_scenario("zr_pipe_nc",     "z_role", nd_frac = 0.0, mcar_frac = 0.4, z_role = "pipe"),
    .v2_scenario("zr_pipenl",      "z_role", nd_frac = 0.0, mcar_frac = 0.4, z_role = "pipe_nl",
                 seed_as = "zr_pipe_nc"),
    .v2_scenario("zr_pipenl_cens", "z_role", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl",
                 seed_as = "zr_pipe")
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
  comp <- simulate_complete(sc$n, truth, rho = sc$rho, skew = sc$skew,
                            z_form = sc$z_form %||% "linear")
  p <- length(truth$b)
  which_cens <- if (isTRUE(sc$censor_all)) seq_len(p) else 1L

  cens <- inject_left_censoring(comp$data, nd_frac = sc$nd_frac,
                                censor_which = which_cens)
  cens <- .v2_mcar_covariates(cens, sc$mcar_frac)
  # MAR covariates (V10). Applied after MCAR so a scenario can in principle carry
  # both; in practice a cell sets one or the other, which is what makes the
  # matched MCAR/MAR pair a clean comparison.
  cens <- .v2_mar_covariates(cens, sc$mar_frac %||% 0.0,
                             strength = sc$mar_strength %||% 1.0)
  cens <- inject_mcar_outcome(cens, sc$y_frac)

  list(complete = comp$data, censored = cens, truth = truth)
}
