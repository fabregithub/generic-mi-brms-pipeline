# =============================================================================
# Track V3 harness -- BKMR mixture estimands and their ANALYTIC truth
# -----------------------------------------------------------------------------
# Track V3 of ../PLAN_pipeline_validation.md; roadmap item 04; design-plan
# Phase 2b; requirement R11.
#
# WHY THIS FILE EXISTS. Phase 1 §7.7 showed linear congenial imputation is biased
# on a mixture surface (+19.6%, coverage 0.62) -- but it showed it on a *scaffold*
# estimand: the local `logX1` coefficient of a matched linear model. That is not
# an estimand anyone reports from a BKMR analysis. This file replaces it with the
# estimands that are actually reported, and -- the part that matters -- computes
# each one's truth ANALYTICALLY from the generator.
#
# WHY ANALYTIC AND NOT AN ORACLE FIT. If the reference were "what BKMR recovers
# on complete data", every measured bias would be the sum of the imputation's
# error and the reference's own error, and the two are not separable. The
# generator's surface is known in closed form, so the truth is available exactly.
# The oracle arm then becomes what it should be -- a *test* of the harness (does
# BKMR on complete data recover the analytic truth?) rather than the yardstick.
#
# THE CONSTRUCTION. Every estimand here is a linear contrast over a shared set of
# evaluation points on the exposure surface:
#
#       estimand = sum_k  w_k * h(u_k)
#
# so the whole grid is a matrix `Z` of points plus a contrast matrix `C`. Two
# consequences, both deliberate:
#
#   * The truth is `C %*% h_true(Z)` -- exact, by construction, with no formula
#     retyped by hand and no chance of the truth and the estimator disagreeing
#     about what is being estimated. The closed forms are still asserted in
#     `check_bkmr_truth()`, as a check on the construction rather than a
#     substitute for it.
#   * The posterior is one `SamplePred()` call per fit, then `C %*% draws`.
#     BKMR's h is identified only up to an additive constant, and every contrast
#     here has weights summing to zero, so that constant cancels -- which is why
#     contrasts, not h itself, are the estimands.
#
# THEORETICAL, NOT EMPIRICAL, QUANTILES. `bkmr::OverallRiskSummaries()` and
# friends evaluate at the *empirical* quantiles of the supplied exposure matrix.
# Those move with n, and -- worse here -- they move with the ARM, because an
# imputed exposure column has different empirical quantiles from the true one.
# That would fold a data-dependent shift into every comparison. The points here
# are fixed at the generator's theoretical quantiles, identical for every arm and
# every replication, so an arm's bias is attributable to the arm.
# =============================================================================

# ---- the true surface --------------------------------------------------------

#' The generator's exposure-response surface, evaluated exactly.
#'
#' Mirrors the `eta` construction in [simulate_complete()] for the exposure part
#' only: covariates and intercept are excluded because every estimand is a
#' zero-sum contrast in which they cancel.
#'
#' @param U Numeric matrix, one row per evaluation point, one column per
#'   log-exposure (`logX1 ... logXp`).
#' @param truth A list from [make_truth()].
#' @return Numeric vector of `h(u)` values, one per row of `U`.
h_true <- function(U, truth) {
  U <- as.matrix(U)
  p <- length(truth$b)
  stopifnot(ncol(U) == p)
  out <- as.vector(U %*% truth$b)
  if (identical(truth$erf_form, "mixture")) {
    if (p >= 2L) out <- out + truth$b_int * U[, 1] * U[, 2]
    out <- out + truth$b_quad * U[, 1]^2
  }
  out
}

# ---- the estimand grid -------------------------------------------------------

#' Build the V3 estimand grid: evaluation points, contrasts, and exact truth.
#'
#' The exposures are marginally `N(mu_x, sd_x^2)` by construction in
#' [simulate_complete()], so the theoretical quantiles are available in closed
#' form and do not depend on the realised sample.
#'
#' THE SEVEN ESTIMANDS, and why each is here:
#'
#' | name | what it is | what it isolates |
#' |---|---|---|
#' | `overall_q75_q50` | all exposures at q75 vs the median | the headline "overall mixture effect" |
#' | `overall_q25_q50` | all exposures at q25 vs the median | its mirror -- the ASYMMETRY between the two is a pure curvature signature |
#' | `overall_q75_q25` | all at q75 vs all at q25 | the contrast the design plan names. **Degenerate when `mu_x` = 0** -- see below |
#' | `singvar_X1_q50` | X1 q75 vs q25, others at the median | the focal single-exposure effect |
#' | `singvar_X1_q75` | X1 q75 vs q25, others at q75 | the same effect in a high-mixture background; differs from the above only through the interaction |
#' | `int_X1X2` | (X1 effect with X2 high) − (X1 effect with X2 low) | **pure interaction**: `b_int (q75-q25)^2`, every other term cancels algebraically |
#' | `curv_X1` | second difference of h in X1 at the median background | **pure curvature**: `2 b_quad q75^2`, every other term cancels |
#'
#' A WARNING ABOUT `overall_q75_q25`. With `mu_x` = 0 the quantiles are symmetric
#' (`q25 = -q75`), so the quadratic and interaction terms cancel *exactly* and
#' this estimand collapses to `sum(b) * (q75 - q25)` -- a purely linear quantity.
#' It is kept because the design plan names it, and because an arm that is
#' unbiased here but biased on `int_X1X2` is evidence about *where* the failure
#' lives rather than merely that one exists. It should not be read as "the
#' mixture estimand": on this generator it does not test the mixture at all.
#'
#' `int_X1X2` and `curv_X1` are the sharp tests. Their truths contain nothing but
#' `b_int` and `b_quad`, the two parameters a linear imputation model has no way
#' to represent.
#'
#' @param truth A list from [make_truth()].
#' @param mu_x,sd_x The generator's marginal log-exposure mean and SD -- must
#'   match the values passed to [simulate_complete()], or the points are not the
#'   quantiles they claim to be.
#' @param q_lo,q_hi The low and high evaluation quantiles. Default 0.25 / 0.75.
#'
#'   MOVING THESE INWARD IS A DIAGNOSTIC, NOT A TUNING KNOB. The oracle
#'   calibration sweep found a stable +6-7% bias concentrated on contrasts that
#'   involve the UPPER-tail points, which is what GP behaviour in a sparse region
#'   of the joint exposure distribution would look like: with `rho` = 0.4, all
#'   three exposures simultaneously at their marginal q75 is a corner the data
#'   populates thinly. Narrowing to q30/q70 pulls the evaluation points into
#'   denser territory. If the floor is an extrapolation artifact it should fall;
#'   if it is intrinsic to the estimator it should not.
#'
#'   READ ABSOLUTE AND RELATIVE BIAS TOGETHER when comparing settings. Narrowing
#'   the quantiles SHRINKS every truth (the contrasts span less of the surface),
#'   so a constant absolute error would show up as a LARGER relative bias. Only a
#'   fall in *relative* bias means the gate got easier.
#'
#'   The estimand NAMES keep their `q75`/`q25` spellings whatever these are set
#'   to, so that results align by name across settings. They denote the hi/lo
#'   ROLES, not literal quantiles; `grid$quantiles` records what was actually
#'   used, and the runner prints it.
#' @return A list with `Z` (points x exposures), `C` (estimands x points),
#'   `true` (named numeric), `names`, and the quantiles used.
bkmr_estimand_grid <- function(truth, mu_x = 0, sd_x = 1,
                               q_lo = 0.25, q_hi = 0.75) {
  p <- length(truth$b)
  if (p < 3L) stop("the V3 estimand grid assumes p >= 3 exposures", call. = FALSE)
  if (!(q_lo > 0 && q_lo < 0.5 && q_hi > 0.5 && q_hi < 1)) {
    stop("need 0 < q_lo < 0.5 < q_hi < 1; got q_lo = ", q_lo,
         ", q_hi = ", q_hi, call. = FALSE)
  }

  qs <- c(lo = q_lo, mid = 0.50, hi = q_hi)
  qv <- mu_x + sd_x * stats::qnorm(qs)          # theoretical marginal quantiles
  a <- unname(qv["lo"]); b <- unname(qv["mid"]); cc <- unname(qv["hi"])

  # Evaluation points. Named so the contrasts below read as algebra rather than
  # as row indices; the trailing exposures (4..p, if any) sit at the median
  # throughout, so they never contribute to a zero-sum contrast.
  rest <- function(v) if (p > 3L) rep(b, p - 3L) else numeric(0)
  pt <- function(u1, u2, u3) c(u1, u2, u3, rest())

  pts <- rbind(
    all_hi     = pt(cc, cc, cc),
    all_mid    = pt(b,  b,  b ),
    all_lo     = pt(a,  a,  a ),
    x1hi_mid   = pt(cc, b,  b ),
    x1lo_mid   = pt(a,  b,  b ),
    x1lo_hi    = pt(a,  cc, cc),
    x1hi_x2hi  = pt(cc, cc, b ),
    x1lo_x2hi  = pt(a,  cc, b ),
    x1hi_x2lo  = pt(cc, a,  b ),
    x1lo_x2lo  = pt(a,  a,  b )
  )
  colnames(pts) <- paste0("logX", seq_len(p))

  # Contrast matrix: one row per estimand, weights summing to zero.
  cn <- rownames(pts)
  mkc <- function(...) {
    w <- stats::setNames(rep(0, length(cn)), cn)
    for (kv in list(...)) w[kv[[1]]] <- w[kv[[1]]] + as.numeric(kv[[2]])
    w
  }
  C <- rbind(
    overall_q75_q50 = mkc(list("all_hi", 1), list("all_mid", -1)),
    overall_q25_q50 = mkc(list("all_lo", 1), list("all_mid", -1)),
    overall_q75_q25 = mkc(list("all_hi", 1), list("all_lo", -1)),
    singvar_X1_q50  = mkc(list("x1hi_mid", 1), list("x1lo_mid", -1)),
    singvar_X1_q75  = mkc(list("all_hi", 1), list("x1lo_hi", -1)),
    int_X1X2        = mkc(list("x1hi_x2hi", 1), list("x1lo_x2hi", -1),
                          list("x1hi_x2lo", -1), list("x1lo_x2lo", 1)),
    curv_X1         = mkc(list("x1hi_mid", 1), list("all_mid", -2),
                          list("x1lo_mid", 1))
  )
  stopifnot(all(abs(rowSums(C)) < 1e-12))       # every contrast must be zero-sum

  true <- as.vector(C %*% h_true(pts, truth))
  names(true) <- rownames(C)

  list(Z = pts, C = C, true = true, names = rownames(C),
       quantiles = qv, probs = qs, mu_x = mu_x, sd_x = sd_x, p = p)
}

#' Independent check of the grid's truth against hand-derived closed forms.
#'
#' The grid computes truth as `C %*% h_true(Z)`. That is correct by construction,
#' which is exactly why it deserves a check written a different way: these
#' formulas were derived by hand from the generator and share no code with the
#' matrix construction, so agreement is real evidence rather than a tautology.
#'
#' @return A data.frame with the two values per estimand and their difference;
#'   `stop()`s if any disagreement exceeds `tol`.
check_bkmr_truth <- function(truth, mu_x = 0, sd_x = 1, tol = 1e-10,
                             q_lo = 0.25, q_hi = 0.75) {
  g <- bkmr_estimand_grid(truth, mu_x = mu_x, sd_x = sd_x,
                          q_lo = q_lo, q_hi = q_hi)
  a <- unname(g$quantiles["lo"]); b <- unname(g$quantiles["mid"])
  cc <- unname(g$quantiles["hi"])
  bb <- truth$b; bi <- truth$b_int; bq <- truth$b_quad
  if (!identical(truth$erf_form, "mixture")) { bi <- 0; bq <- 0 }
  sb <- sum(bb)

  closed <- c(
    # all exposures move together: every b_j contributes, and u1*u2 and u1^2
    # both become the square of the common level.
    overall_q75_q50 = sb * (cc - b)  + bi * (cc * cc - b * b) + bq * (cc^2 - b^2),
    overall_q25_q50 = sb * (a  - b)  + bi * (a  * a  - b * b) + bq * (a^2  - b^2),
    overall_q75_q25 = sb * (cc - a)  + bi * (cc * cc - a * a) + bq * (cc^2 - a^2),
    # only X1 moves; X2 sits at the stated background level, so the interaction
    # enters scaled by that level.
    singvar_X1_q50  = bb[1] * (cc - a) + bi * b  * (cc - a) + bq * (cc^2 - a^2),
    singvar_X1_q75  = bb[1] * (cc - a) + bi * cc * (cc - a) + bq * (cc^2 - a^2),
    # pure interaction and pure curvature: linear and nuisance terms cancel.
    int_X1X2        = bi * (cc - a)^2,
    curv_X1         = bq * ((cc - b)^2 + (a - b)^2)
  )

  out <- data.frame(estimand = g$names,
                    grid = unname(g$true[g$names]),
                    closed_form = unname(closed[g$names]),
                    stringsAsFactors = FALSE)
  out$diff <- out$grid - out$closed_form
  bad <- abs(out$diff) > tol
  if (any(bad)) {
    stop("analytic truth disagrees with the closed form for: ",
         paste(out$estimand[bad], collapse = ", "),
         " (max |diff| = ", format(max(abs(out$diff))), ")", call. = FALSE)
  }
  out
}

# ---- posterior side ----------------------------------------------------------

#' Posterior draws of every estimand from one fitted BKMR model.
#'
#' `SamplePred()` returns draws of `h(Znew) + Xnew %*% beta`. `Xnew` is held
#' constant across the evaluation points, so the `beta` part is identical in
#' every column and cancels in each zero-sum contrast -- the returned draws are
#' of `h` contrasts alone, with no assumption that `beta` was estimated well.
#'
#' @param fit A `bkmrfit` from `bkmr::kmbayes()`.
#' @param grid A list from [bkmr_estimand_grid()].
#' @param sel Iterations to keep (post-burn-in). Defaults to the second half.
#' @param x_ref Row vector of covariate values to hold fixed. Its value is
#'   irrelevant to the result (it cancels); zeros keep it inside the data's range.
#' @return A draws x estimands matrix, or NULL on failure.
bkmr_estimand_draws <- function(fit, grid, sel = NULL, x_ref = NULL) {
  n_iter <- tryCatch(fit$iter, error = function(e) NULL)
  if (is.null(n_iter)) n_iter <- nrow(fit$beta)
  if (is.null(sel)) sel <- seq.int(floor(n_iter / 2) + 1L, n_iter)

  n_x <- ncol(fit$X)
  Xnew <- matrix(if (is.null(x_ref)) 0 else x_ref,
                 nrow = nrow(grid$Z), ncol = n_x, byrow = TRUE)
  colnames(Xnew) <- colnames(fit$X)

  Znew <- grid$Z
  colnames(Znew) <- colnames(fit$Z)

  sp <- tryCatch(
    bkmr::SamplePred(fit, Znew = Znew, Xnew = Xnew, sel = sel, type = "link"),
    error = function(e) NULL)
  if (is.null(sp)) return(NULL)

  sp <- as.matrix(sp)
  # COLUMN ORDER, NOT NAMES. `SamplePred()` returns columns named `znew1..znewK`
  # regardless of what `Znew`'s rownames were, so the contrast matrix is applied
  # positionally. `grid$C`'s columns are in `grid$Z`'s row order, and SamplePred
  # preserves the row order of `Znew`; the dimension check below is the only
  # thing standing between that assumption and silently contrasting the wrong
  # points, so it is a hard failure rather than a warning.
  if (ncol(sp) != nrow(grid$Z)) return(NULL)
  draws <- sp %*% t(grid$C)
  colnames(draws) <- rownames(grid$C)
  draws
}

#' Collapse posterior draws to the (estimate, variance) pair Rubin's rules need.
#'
#' Returns the posterior MEAN and VARIANCE per estimand -- not quantiles --
#' because [rubin_pool()] combines moments. For a single-fit (non-MI) arm the
#' caller can use the percentile interval instead; both are provided.
#'
#' @param draws A draws x estimands matrix from [bkmr_estimand_draws()].
#' @return A data.frame with one row per estimand.
bkmr_draws_moments <- function(draws) {
  data.frame(
    estimand = colnames(draws),
    est      = apply(draws, 2, mean),
    var      = apply(draws, 2, stats::var),
    cri_lo   = apply(draws, 2, stats::quantile, probs = 0.025, names = FALSE),
    cri_hi   = apply(draws, 2, stats::quantile, probs = 0.975, names = FALSE),
    stringsAsFactors = FALSE, row.names = NULL
  )
}
