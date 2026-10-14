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
                         y_form = "linear", z_role = "precision",
                         nl_a = NULL, nl_c = NULL, target = "direct",
                         delta_xz = NULL, delta_zx = NULL, sd_x = 1.0, y_tau = NULL,
                         y_family = "gaussian", seed_as = NULL) {
  list(name = name, axis = axis, erf_form = erf_form, nd_frac = nd_frac,
       n = as.integer(n), m = as.integer(m), n_rep = n_rep, rho = rho,
       skew = skew, censor_all = censor_all, mcar_frac = mcar_frac,
       y_frac = y_frac, z_form = z_form,
       mar_frac = mar_frac, mar_strength = mar_strength, y_form = y_form,
       # V17: the covariate's causal role -- see make_truth() in dgp.R.
       z_role = z_role,
       # V23: the non-linear arrow's shape, swept to test bias ~ u^2.
       nl_a = nl_a, nl_c = nl_c,
       # V29: the exposure's own SD. In a fork Var(logX1) = sd_x^2 + delta_zx^2,
       # so this moves s^2 with the arrow -- and therefore the exposure
       # conditional's Z1 coefficient -- held FIXED. V28 found the chain's
       # derived/measured ratio trending with the arrow in the fork sweep but
       # flat in the pipe sweep, which implicates Var(logX1) rather than the
       # coefficient; only sd_x can test that channel on its own. Default 1.0,
       # the value every pre-V29 cell used, so no existing cell moves.
       sd_x = sd_x,
       # V24: WHICH ESTIMAND. Only a pipe has a choice, and the two answers are
       # different numbers -- so this changes `estimand_true` and the analysis
       # model's adjustment set together, via make_truth(). See dgp.R.
       target = target,
       # V25: the X1 -> Z1 arrow's strength. NULL keeps 0.60; 0 removes the arrow
       # while keeping Z1 present and observed -- the arrow-absent control.
       delta_xz = delta_xz, delta_zx = delta_zx,
       # V26: the outcome's family. Gaussian everywhere before V26.
       y_family = y_family,
       # V30: cumulative-logit cutpoints, carried to make_truth(). Literals on
       # purpose -- see the note in dgp.R's make_truth().
       y_tau = y_tau,
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
                 seed_as = "zr_pipe"),

    # ---- V23: the curvature sweep, testing bias ~ u^2 --------------------------
    # DERIVATION-LED, unlike V22 which found the ~20% defect by surprise. The
    # damage is not "how big is g" but how much of g a LINEAR conditional cannot
    # express where the censored mass sits -- `ef_unrep_curvature()`, written u.
    # And it should enter SQUARED, because that residual is orthogonal to the
    # span of the linear predictors, so the first-order term in the bias
    # expansion vanishes. See ef_unrep_curvature()'s docs in dgp.R.
    #
    # `cv262` is the calibration point (V22's own arrow, u = 0.262). Everything
    # else is a prediction. `cv000` is the null: g == 0 means the arrow is gone
    # and there is nothing for a linear draw to miss.
    .v2_scenario("cv000", "curvature", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl", nl_a = 0.0, nl_c = 0.00),
    .v2_scenario("cv040", "curvature", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl", nl_a = 0.4, nl_c = 0.00),
    .v2_scenario("cv119", "curvature", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl", nl_a = 0.8, nl_c = 0.10),
    .v2_scenario("cv153", "curvature", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl", nl_a = 0.0, nl_c = 0.35),
    .v2_scenario("cv182", "curvature", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl", nl_a = 0.8, nl_c = 0.25),
    .v2_scenario("cv262", "curvature", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl", nl_a = 1.2, nl_c = 0.35),
    .v2_scenario("cv364", "curvature", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe_nl", nl_a = 1.6, nl_c = 0.50),

    # ---- V24: the four cases that need the DAG-factorised exposure draw ------
    # THEORY.md 6b concluded the draw is needed in FOUR cases, differing along
    # TWO axes -- the covariate's role, and (for a pipe) which estimand:
    #
    #   dag_fork        fork      direct   Z1 adjusted (backdoor); Z1 a PARENT
    #   dag_pipe_dir    pipe      direct   Z1 adjusted;            Z1 a CHILD
    #   dag_pipe_tot    pipe      total    Z1 NOT adjusted;        Z1 a CHILD
    #   dag_pnl_dir     pipe_nl   direct   as pipe_dir, arrow non-linear
    #   dag_pnl_tot     pipe_nl   total    as pipe_tot, arrow non-linear
    #   dag_collider    collider  direct   Z1 NOT adjusted;   child of X AND Y
    #
    # Testing all six in one sweep is what makes the DIRECTION claim checkable:
    # the four cases have different bias signs, and a track that ran them
    # separately could not rule out the sign being an artefact of its cell.
    #
    # WHY THE LINEAR `pipe` CELLS ARE HERE and not only pipe_nl. The linear pipe
    # is the control that separates the FACTORISATION from the NON-LINEARITY: its
    # child arrow is linear-Gaussian, so a correct factorisation with a linear
    # child model should be exactly right there. If a pipe_nl cell improves and
    # the linear pipe does not, the gain is not the factorisation.
    #
    # `seed_as` pairs each pipe_nl cell with its linear counterpart, and the two
    # estimand cells with each other, so the contrasts are paired.
    # WHY mcar_frac = 0 HERE, when every zr_* cell uses 0.4. V24's question is
    # about the EXPOSURE draw, and the covariate block was settled by V20 (the
    # block attribution) and V21 (what about the covariate draw has to be right).
    # Leaving Z 40% missing would put V21's ladder back in as a confound, and
    # the DAG factorisation conditions on the covariates in all three of its
    # factors -- so a missing Z1 is not a harder version of the same question,
    # it is a different question. This is the mirror of `zr_pipenl`, which sets
    # nd_frac = 0 to isolate the covariate block. The combined case is open
    # (ROADMAP) and dag_impute_datasets() refuses it rather than returning NA.
    .v2_scenario("dag_fork",     "dag_draw", nd_frac = 0.4, mcar_frac = 0.0,
                 z_role = "fork"),
    .v2_scenario("dag_pipe_dir", "dag_draw", nd_frac = 0.4, mcar_frac = 0.0,
                 z_role = "pipe", target = "direct"),
    .v2_scenario("dag_pipe_tot", "dag_draw", nd_frac = 0.4, mcar_frac = 0.0,
                 z_role = "pipe", target = "total", seed_as = "dag_pipe_dir"),
    .v2_scenario("dag_pnl_dir",  "dag_draw", nd_frac = 0.4, mcar_frac = 0.0,
                 z_role = "pipe_nl", target = "direct", seed_as = "dag_pipe_dir"),
    .v2_scenario("dag_pnl_tot",  "dag_draw", nd_frac = 0.4, mcar_frac = 0.0,
                 z_role = "pipe_nl", target = "total", seed_as = "dag_pipe_dir"),
    .v2_scenario("dag_collider", "dag_draw", nd_frac = 0.4, mcar_frac = 0.0,
                 z_role = "collider"),

    # ---- V25 item 1: the ARROW-ABSENT control -------------------------------
    # THE GAP THIS CLOSES. V24 concluded that discarding an informative
    # covariate is what generates the 6-7% asymptotic bias, and that `u`
    # (unrepresentable curvature) only modifies it. The evidence was that `u` = 0
    # in two situations 6 pp apart: V23's `cv000` (arrow absent) and V24's
    # `dag_pipe_*` (arrow present but linear). But cv000 was measured on the
    # PIPELINE arms and dag_pipe_* on V24's harness draw, so the contrast crossed
    # instruments -- recorded as a caveat in THEORY.md 6b.
    #
    # `delta_xz = 0` closes it inside one instrument. Z1 is still there, still
    # observed, still in the analysis model -- it simply carries no information
    # about the exposure, so p(Z1 | x) = p(Z1) and the child factor has nothing
    # to contribute. Every arm should therefore sit at ~0.
    #
    # AND IT IS GENUINELY PAIRED. Setting delta_xz = 0 leaves the RNG consumption
    # order untouched, so this cell sees BYTE-IDENTICAL exposures to
    # `dag_pipe_dir` (verified: cor(logX1, logX1) = 1.000). Note `seed_as` does
    # NOT achieve that across different roles -- fork and pipe consume the stream
    # in different orders and their exposures correlate 0.057 -- which is why the
    # control had to be built inside the pipe structure rather than as a new role.
    .v2_scenario("dag_pipe_null", "dag_draw", nd_frac = 0.4, mcar_frac = 0.0,
                 z_role = "pipe", target = "direct", delta_xz = 0,
                 seed_as = "dag_pipe_dir"),

    # ---- V26: does the record transfer to a NON-GAUSSIAN estimand? ----------
    # Every track V0-V25 estimated a Gaussian linear coefficient by OLS. The
    # pipeline's real targets are logistic, ordinal, spline and mixed models, so
    # the transfer of 25 tracks' conclusions to those classes is an ASSUMPTION.
    #
    # These cells are the binomial twins of three whose Gaussian answers are
    # already known, so the contrast is family and nothing else: the structural
    # coefficients, the covariate roles, the missingness and the censoring are
    # all unchanged, and the same seed gives BYTE-IDENTICAL exposures (only Y
    # differs, since Y is drawn last). The estimand becomes the focal exposure's
    # log-odds-ratio, still 0.40.
    #
    # Registered comparison, per cell: does the arm ORDERING survive, and does
    # the defect's SIZE in bias/SE units survive? A logistic coefficient is on a
    # different scale, so relative bias is not comparable across families --
    # bias/SE is, which is why V18's currency matters here.
    .v2_scenario("bin_base",  "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 y_family = "binomial", seed_as = "mcar_z40"),
    .v2_scenario("bin_fork",  "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "fork", y_family = "binomial", seed_as = "zr_fork"),
    .v2_scenario("bin_pipe",  "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "pipe", y_family = "binomial", seed_as = "zr_pipe"),

    # ---- V28: what is fork-specific about the scale defect's residual? ------
    # The propagation chain (THEORY.md 1a(ii)) predicts the pipe cell to 1.1%
    # and the fork cell to 15%. The exp-linearisation hypothesis was tested and
    # REFUTED (the agreement improves with s^2, not degrades). Three rebuilt
    # instruments all missed the pipeline's own anchor by 31%-4x, so this sweep
    # uses the REAL pipeline arms -- `pipeline_properZ` against
    # `pipeline_properZ_ds`, which differ only by the pre-v1.6.0 scale defect
    # and run on identical data. The anchor is then hit by construction.
    #
    # THE DESIGN SEPARATES TWO THINGS THE FORK CONFOUNDS. Under a fork,
    # delta_zx moves BOTH Var(logX1) = 1 + delta^2 AND the exposure
    # conditional's Z1 coefficient. Under a pipe, delta_xz moves the Z1-X1
    # relationship while Var(logX1) stays 1. So a ratio that tracks the arrow
    # in the PIPE sweep implicates the coefficient; one that tracks it only in
    # the FORK sweep implicates Var(logX1).
    .v2_scenario("fkd02", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "fork", delta_zx = 0.2),
    .v2_scenario("fkd04", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "fork", delta_zx = 0.4),
    .v2_scenario("fkd06", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "fork", delta_zx = 0.6),
    .v2_scenario("fkd08", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "fork", delta_zx = 0.8),
    .v2_scenario("ppd02", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe", delta_xz = 0.2),
    .v2_scenario("ppd04", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe", delta_xz = 0.4),
    .v2_scenario("ppd06", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe", delta_xz = 0.6),
    .v2_scenario("ppd08", "arrow", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe", delta_xz = 0.8),

    # ---- V29: s^2 ON ITS OWN, at a FIXED arrow -----------------------------
    # V28 (n = 12800, delta >= 0.4) found the chain's derived/measured ratio
    # TRENDING with the arrow in the fork sweep -- 0.944, 0.992, 1.112 -- and
    # FLAT in the pipe sweep -- 0.942, 0.953, 0.952. Under a fork the arrow
    # moves both Var(logX1) = sd_x^2 + delta_zx^2 and the exposure
    # conditional's Z1 coefficient; under a pipe it moves only the coefficient,
    # and Var(logX1) = sd_x^2 regardless. So the flat pipe sweep ACQUITS the
    # coefficient and the trending fork sweep implicates Var(logX1).
    #
    # `sd_x` tests that channel directly, because it moves s^2 in BOTH roles
    # while the arrow stays at 0.60. s^2 = 0.85 / 1.36 / 2.05 (fork) and
    # 0.49 / 1.00 / 1.69 (pipe) -- the fork's middle cell reproduces the
    # canonical `zr_fork` geometry as a within-run anchor.
    #
    # REGISTERED PREDICTION (2026-10-13, before the run):
    #   P1  the ratio rises with sd_x in BOTH roles
    #         -> s^2 is the driver; the chain's linearisation of exp() is what
    #            degrades, and V28's fork/pipe contrast is a consequence of the
    #            fork's arrow moving s^2 rather than of the role.
    #   P2  the ratio rises with sd_x in the FORK ONLY
    #         -> something fork-specific survives, and s^2 is not sufficient.
    #   P3  flat in both -> s^2 is NOT the driver and V28's fork trend needs
    #            another explanation entirely.
    # P1 and P2 are distinguishable because the pipe cells move s^2 over a
    # range (0.49 -> 1.69) WIDER than the fork sweep's own (1.16 -> 1.64),
    # so a pipe that stays flat cannot be dismissed as too small a lever.
    #
    # This re-opens a question THEORY.md records as refuted. The refutation
    # swept s^2 = 0.52 -> 2.05 and found the ratio moving TOWARD 1
    # (0.841 -> 0.908) -- but on a two-block simulation that misses the
    # pipeline's own level by 31%, a limitation THEORY states explicitly. These
    # cells use the real pipeline arms, so they hit the anchor by construction.
    .v2_scenario("sxf07", "sdx", nd_frac = 0.4, mcar_frac = 0.4, z_role = "fork", delta_zx = 0.6, sd_x = 0.7),
    .v2_scenario("sxf10", "sdx", nd_frac = 0.4, mcar_frac = 0.4, z_role = "fork", delta_zx = 0.6, sd_x = 1.0),
    .v2_scenario("sxf13", "sdx", nd_frac = 0.4, mcar_frac = 0.4, z_role = "fork", delta_zx = 0.6, sd_x = 1.3),
    .v2_scenario("sxp07", "sdx", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe", delta_xz = 0.6, sd_x = 0.7),
    .v2_scenario("sxp10", "sdx", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe", delta_xz = 0.6, sd_x = 1.0),
    .v2_scenario("sxp13", "sdx", nd_frac = 0.4, mcar_frac = 0.4, z_role = "pipe", delta_xz = 0.6, sd_x = 1.3),

    # ---- V30: does the record transfer to an ORDINAL outcome? ---------------
    # V26 showed the Gaussian record transfers to a BINARY outcome in direction
    # and ordering but at ~0.474x in bias/SE, and THEORY.md then derived that
    # factor as Fisher information rather than anything about imputation:
    # SE_gauss/SE_fam = sigma * sqrt(E[w]), with w = p(1-p) for a binomial.
    #
    # The cumulative-logit generalisation has NO free parameters:
    #   w_ord(eta) = sum_k [f(tau_k-eta) - f(tau_{k-1}-eta)]^2 / pi_k,  f = F(1-F)
    # It reduces to p(1-p) at K=2 to machine precision, and E[w_ord] -> 1/3 as
    # K -> inf (the location information of a standard logistic -- the latent
    # observed exactly). See derive_ordinal_information.R.
    #
    # sqrt(E[w_ord]) IS NOT THE RATIO. Measuring the complete-data SE ratio
    # BEFORE launching caught two separate errors in it, both one-signed:
    #   (a) the information for beta is E[w*xtilde^2], not E[w]*E[xtilde^2]
    #       (w peaks mid-range, eta rises in X1, so the covariance is negative);
    #   (b) polr ESTIMATES the K-1 cutpoints, and unlike the Gaussian twin's
    #       nuisance parameters they are not orthogonal to beta -- so the
    #       variance is the (1,1) element of the INVERSE full (beta,tau)
    #       information, not 1/I_bb. Flat in n (+0.0073/+0.0061/+0.0071 at
    #       n = 2k/8k/32k), so asymptotic, not finite-sample.
    # With both included the prediction is EXACT: validated against the
    # complete-data SE ratio in six cells, misses -0.0000..+0.0005, mixed
    # signs, no free parameters. See derive_ordinal_information.R.
    #
    # REGISTERED PREDICTIONS (2026-10-14, before the run):
    #   P1 LEVEL  each cell's measured |bias/SE| ratio against its Gaussian
    #             twin matches the exact value:
    #                 K=2  fork 0.4459  pipe 0.4488    (bin_* -- the anchor)
    #                 K=3  fork 0.5091  pipe 0.5120
    #                 K=5  fork 0.5449  pipe 0.5477
    #                 K=7  fork 0.5559  pipe 0.5586
    #   P2 TREND  the ratio RISES with K, concavely, toward 0.5774.
    #   P3 ANCHOR bin_fork/bin_pipe reproduce V26's measured 0.4731 / 0.4527.
    #
    # WHAT WOULD FALSIFY P2, and why it is resolvable: the predicted span
    # K=2 -> K=7 is 0.446 -> 0.559, about 2x V26's across-arm IQR (0.445-0.493).
    # It is resolvable but NOT large, so the gate is on the TREND across K with
    # arms pooled, not on any single cell-arm ratio.
    #
    # NOT A GATE: the cell-to-cell contrast. At K=5 fork and pipe differ by
    # 0.0002, far under any achievable MC error, because E[w_ord] -> 1/3
    # regardless of Var(eta). Registering "reproduces the cell ordering" here
    # would be a criterion that cannot fail for the reason it names -- the sixth
    # such trap in this programme. The two roles are REPLICATION of the K trend,
    # nothing more.
    #
    # tau are the K-quantiles of the latent Y* = eta + Logistic(0,1), so the
    # categories are balanced; they are DGP constants, computed once at
    # N = 2e6 and written here as literals.
    .v2_scenario("ord3_fork", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "fork", y_family = "ordinal", seed_as = "zr_fork",
                 y_tau = c(-0.9732, 0.6698)),
    .v2_scenario("ord5_fork", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "fork", y_family = "ordinal", seed_as = "zr_fork",
                 y_tau = c(-1.7753, -0.6340, 0.3302, 1.4724)),
    .v2_scenario("ord7_fork", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "fork", y_family = "ordinal", seed_as = "zr_fork",
                 y_tau = c(-2.2329, -1.2354, -0.4944, 0.1897, 0.9312, 1.9282)),
    .v2_scenario("ord3_pipe", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "pipe", y_family = "ordinal", seed_as = "zr_pipe",
                 y_tau = c(-0.9740, 0.6719)),
    .v2_scenario("ord5_pipe", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "pipe", y_family = "ordinal", seed_as = "zr_pipe",
                 y_tau = c(-1.7773, -0.6354, 0.3303, 1.4776)),
    .v2_scenario("ord7_pipe", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "pipe", y_family = "ordinal", seed_as = "zr_pipe",
                 y_tau = c(-2.2346, -1.2351, -0.4950, 0.1897, 0.9347, 1.9334)),

    # ---- V30b: UNBALANCED categories, at FIXED K ---------------------------
    # The cells above vary only K, so a measured rise is consistent with the
    # derivation but ALSO with "ratios happen to rise with K". Unbalanced
    # cutpoints move the prediction with K HELD FIXED, which the second story
    # cannot do -- so these are what make the formula earn its keep.
    #
    # The separation is larger than the whole balanced K sweep: at K = 3 the
    # prediction spans 0.509 -> 0.431 -> 0.336 (span 0.173) against 0.446 ->
    # 0.556 for K = 2 -> 7 (span 0.110).
    #
    # And it is COUNTERINTUITIVE, which is the point: a 3-category outcome
    # split 85/10/05 is predicted to carry LESS information than a balanced
    # BINARY one (0.336 vs 0.446). More categories is not more information if
    # they are in the wrong places. Nothing but the formula predicts that.
    #
    # REGISTERED (2026-10-14, before the run):
    #   P4  at K = 3, the ratio FALLS as the split gets more unbalanced, and
    #       the 85/10/05 cell falls BELOW its own K = 2 cell:
    #         balanced   fork 0.5091  pipe 0.5121
    #         70/20/10   fork 0.4306  pipe 0.4327
    #         85/10/05   fork 0.3358  pipe 0.3382   <- below bin_* (0.446/0.449)
    # An unbalanced K = 5 cell (fork 0.5055 / pipe 0.5080) was designed and
    # validated but NOT included: it separates from balanced K = 5 by 0.039,
    # under the +/-0.05 gate, so it could not fail informatively -- two cells
    # of cost for no discrimination. The K = 3 fan carries this test.
    .v2_scenario("ord3u_fork", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "fork", y_family = "ordinal", seed_as = "zr_fork",
                 y_tau = c(0.8557, 2.3818)),
    .v2_scenario("ord3v_fork", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "fork", y_family = "ordinal", seed_as = "zr_fork",
                 y_tau = c(1.8754, 3.1807)),
    .v2_scenario("ord3u_pipe", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "pipe", y_family = "ordinal", seed_as = "zr_pipe",
                 y_tau = c(0.8544, 2.3787)),
    .v2_scenario("ord3v_pipe", "y_family", nd_frac = 0.4, mcar_frac = 0.4,
                 z_role = "pipe", y_family = "ordinal", seed_as = "zr_pipe",
                 y_tau = c(1.8644, 3.1854))
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
                            sd_x = sc$sd_x %||% 1.0,
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
