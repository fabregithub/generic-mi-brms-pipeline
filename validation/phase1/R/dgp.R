# =============================================================================
# Phase 1 harness -- data-generating process (DGP)
# -----------------------------------------------------------------------------
# Known-ERF generator for the left-censored-exposure validation study (PLAN §7).
# Everything is simulated, so the exposure-response estimand is known exactly and
# bias is measured against ground truth.
#
# Two ERF forms (PLAN §7.1):
#   * "additive"  : Y = a + Σ b_j logX_j + γ Z + ε           (linear, per-analyte)
#   * "mixture"   : Y = a + h(logX) + γ Z + ε                (non-linear + interaction)
#
# Design note for the SCAFFOLD: exposures are generated on the *log* scale
# (multivariate normal, optionally skewed), which is the scale `leftcens` and the
# outcome model both work on. The focal censored exposure is X1; its coefficient
# `b1` (its main effect on the log scale) is the primary Phase-1 estimand and is
# reported identically by every procedure so they are directly comparable. The
# richer mixture estimands (overall mixture effect, interactions, h at profiles;
# PLAN §7.4) are a documented extension -- see README "Scaffold limitations".
# =============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Draw the true regression coefficients that define the ERF.
#'
#' Returned separately from the data so the "truth" is explicit and auditable.
#'
#' @param p Number of exposures.
#' @param q Number of covariates (here fixed at 2: one continuous, one binary).
#' @param erf_form "additive" or "mixture".
#' @param y_form Shape of the outcome in the COVARIATES. "linear" (default)
#'   reproduces every track up to V10: `Y` is linear in `(logX, Z)`, so the
#'   conditional of `Z1` given `(Y, X)` is linear-Gaussian and a *parametric*
#'   imputation model for `Z1` is correctly specified.
#'
#'   WHY "nonlinear" EXISTS (V11). That linearity is a confound in a shipped
#'   decision. `z_imputer = "bart"` was chosen over `mice pmm` and `forest_boot`
#'   in V6/V7 on a design where the outcome is linear in `Z` by construction --
#'   and FINDINGS_v7.md says so itself: "`Y` is linear in `(logX, Z)` by
#'   construction, which is why parametric arms do well in the outcome-dominated
#'   cells." A linear imputer is never penalised for misspecification against the
#'   outcome, so the comparison could only ever reward calibration.
#'
#'   Under "nonlinear" the outcome gains `b_zq * (Z1^2 - 1)`, which makes
#'   `p(Z1 | Y, X)` non-linear and a linear imputation model for `Z1` genuinely
#'   misspecified. `z_form = "nonlinear"` (V6) is a different axis: it makes `Z1`
#'   a non-linear function of the OTHER PREDICTORS. This one makes the OUTCOME
#'   non-linear in `Z1`, which is what conditions the imputation draw.
#'
#'   The term is centred (`Z1^2 - 1` has mean zero for standard-normal `Z1`) so
#'   the outcome's mean is unchanged, and [dgp_formula()] adds the matching
#'   `I(Z1^2)` term so the ANALYSIS model stays correctly specified. That
#'   separation is the whole point: the focal estimand `b[1]` remains exactly
#'   recoverable, so any bias is attributable to the imputation model and never
#'   to analysis misspecification.
make_truth <- function(p = 3L, q = 2L, erf_form = "additive",
                       y_form = "linear", z_role = "precision",
                       nl_a = NULL, nl_c = NULL, target = "direct",
                       mu_x = 0.0, sd_x = 1.0, delta_xz = NULL) {
  b <- rep(0.0, p)
  b[1] <- 0.40                      # focal exposure main effect (the estimand)
  if (p >= 2L) b[2] <- 0.20
  if (p >= 3L) b[-(1:2)] <- 0.10

  gamma <- c(0.50, -0.30)[seq_len(q)]

  # Interaction / curvature used only by the mixture surface.
  b_int  <- if (erf_form == "mixture" && p >= 2L) 0.25 else 0.0   # logX1 * logX2
  b_quad <- if (erf_form == "mixture") 0.15 else 0.0              # logX1^2

  # Curvature of the outcome in the continuous covariate (V11). Sized comparably
  # to gamma[1] = 0.5 so it is a real feature of the surface rather than a nudge.
  b_zq <- if (identical(y_form, "nonlinear")) 0.40 else 0.0

  # ---- V17: the covariate's CAUSAL ROLE ---------------------------------------
  # Tracks V0-V16 all ran with one covariate structure, in which Z is either
  # independent of the exposures ("precision") or generated FROM them
  # ("descendant", i.e. z_form = "nonlinear"). NEITHER is a confounder, so no
  # track has ever adjusted for a covariate that had to be adjusted for. Drawing
  # V16's DAG is what made that visible. The roles below are the missing ones:
  #
  #   precision   Z1 -> Y only                  (the original design)
  #   descendant  X2, X3 -> Z1 -> Y             (z_form = "nonlinear")
  #   fork        Z1 -> X1,  Z1 -> Y            CONFOUNDER: must be adjusted for
  #   pipe        X1 -> Z1 -> Y                 MEDIATOR: adjusting gives the
  #                                             DIRECT effect, which is b[1]
  #   collider    X1 -> Z1 <- Y                 adjusting INDUCES bias; the
  #                                             correct analysis omits Z1
  #   mixed       Z1 fork AND Z2 pipe           no single covariate is both
  #   pipe_nl     X1 -> Z1 -> Y, arrow NON-LINEAR   THE DECIDING CELL (V22)
  #
  # WHY pipe_nl EXISTS. V21 found a correctly specified ESTIMATED covariate draw
  # is unbiased where the covariate conditional is linear-Gaussian -- which both
  # the `fork` and `pipe` cells are BY CONSTRUCTION. That makes V21's result one
  # half of a trade: R8 adopted BART because a parametric Z block is misspecified
  # when the conditional is NON-linear (V6 measured mice pmm at -4.70% there).
  # Trading a 5% bias under linearity for a 5% bias under non-linearity is not a
  # fix, and no cell could test both sides -- `z_form = "nonlinear"` makes Z1 a
  # DESCENDANT of the exposures, which V17 showed is off-path.
  #
  # `pipe_nl` is the missing cell: Z1 is genuinely ON an X-Y path AND its
  # conditional has a non-linear mean, so a linear imputation model is really
  # misspecified while BART can learn the shape. Verified before adoption:
  # the analysis model still recovers b1 (0.3996 at n = 3e5), the partial
  # correlation of Z1 with logX1 given the rest is +0.601, a linear imputation
  # model's residual SD is 21% worse than the correct form, and the exact
  # conditional of Z1 stays GAUSSIAN (skew 0.001, kurtosis 3.001) -- which is
  # what keeps exact_fork.R's anchor closed-form.
  #
  # TWO STRUCTURAL CONSEQUENCES, both handled here rather than left to the
  # caller. Under "collider" Z1 must NOT cause Y (or it would be a fork as well),
  # so gamma[1] is zeroed; and the correctly-specified analysis model must DROP
  # Z1, which is what `z_in_model` tells dgp_formula(). Everywhere else
  # `z_in_model` is the full set and the formula is unchanged.
  if (!z_role %in% c("precision", "descendant", "fork", "pipe", "pipe_nl",
                     "collider", "mixed")) {
    stop("unknown z_role: ", z_role, call. = FALSE)
  }
  z_in_model <- c("Z1", "Z2")[seq_len(q)]
  if (identical(z_role, "collider")) {
    gamma[1] <- 0.0                       # Z1 is an effect of Y, never a cause
    z_in_model <- setdiff(z_in_model, "Z1")
  }

  # ---- V24: WHICH ESTIMAND, for the roles where there is a choice -----------
  # THEORY.md 1b(i): the adjustment set is not a modelling preference, it is part
  # of the estimand's definition. Under a pipe both answers are correct and they
  # are DIFFERENT numbers -- adjusting for the mediator gives the DIRECT effect
  # b[1], omitting it gives the TOTAL effect. Every track before V24 estimated
  # the direct effect, because `z_in_model` was always the full set. `target`
  # makes the other choice available so the two can be tested side by side.
  #
  # Under a fork there is no choice: Z1 opens a backdoor, so it must be adjusted
  # for and "total" is not a separate estimand. Under a collider Z1 is off every
  # X-Y path, so total and direct coincide and Z1 is out either way. Asking for
  # "total" in those roles is therefore a specification error, not a variant.
  if (!target %in% c("direct", "total"))
    stop("unknown target: ", target, call. = FALSE)
  if (identical(target, "total") && !z_role %in% c("pipe", "pipe_nl")) {
    stop("target = \"total\" is only defined for a pipe; under z_role = \"",
         z_role, "\" the total and direct effects are either not both ",
         "identified (fork) or identical (collider, precision, descendant).",
         call. = FALSE)
  }

  # The total effect is what a model WITHOUT the mediator recovers.
  #
  #   pipe:     Z1 = delta_xz * logX1 + noise, noise independent of X, so
  #             Y = (b1 + delta_xz * gamma1) logX1 + ... exactly.
  #
  #   pipe_nl:  Z1 = g(logX1) + noise with g non-linear, so there is no single
  #             structural total effect -- the derivative varies with x. What the
  #             analysis model targets is the LINEAR PROJECTION, and for jointly
  #             Gaussian exposures Stein's identity makes that projection exact:
  #             Cov(g(X1), Xj) = Cov(X1, Xj) E[g'(X1)], so the coefficient vector
  #             of g(X1) on (X1, X2, X3) is E[g'(X1)] * e1 -- ALL of it lands on
  #             logX1 and none on the other exposures. Hence
  #                 total = b1 + gamma1 * E[g'(X1)],  X1 ~ N(mu_x, sd_x^2).
  #             Verified by simulation at n = 4e5: 0.8269 against 0.8298 here
  #             (1.4 SE), and the logX2 coefficient stays at its structural
  #             0.1996 -- which is the "none on the others" half of the claim.
  #             `mu_x`/`sd_x` must match what is passed to the simulator; the
  #             defaults are the ones every scenario uses.
  # Structural arrow strengths, needed here as well as in the returned list.
  # `delta_xz` is overridable so a scenario can set it to ZERO: that removes the
  # X1 -> Z1 arrow while leaving Z1 present, observed, and consuming the RNG in
  # exactly the same order -- which makes such a cell PAIRED with the ordinary
  # pipe cell on byte-identical exposures. That is the control V24 lacked (see
  # THEORY.md 6b): it separates "an informative covariate is being discarded"
  # from "the covariate's relationship cannot be represented", within one
  # instrument. Both have u = 0; only the second has a p(Z1 | x) that depends on
  # x at all.
  delta_zx <- 0.60; delta_xz <- delta_xz %||% 0.60
  delta_yz <- 0.60; delta_xz2 <- 1.20
  nl_a_v <- nl_a %||% 1.20; nl_b_v <- 1.80; nl_c_v <- nl_c %||% 0.35

  total_effect <- if (identical(z_role, "pipe")) {
    b[1] + delta_xz * gamma[1]
  } else if (identical(z_role, "pipe_nl")) {
    xq <- seq(mu_x - 10 * sd_x, mu_x + 10 * sd_x, length.out = 20001L)
    wq <- stats::dnorm(xq, mu_x, sd_x); wq <- wq / sum(wq)
    gprime <- nl_a_v * nl_b_v / cosh(nl_b_v * xq)^2 + 2 * nl_c_v * xq
    b[1] + gamma[1] * sum(wq * gprime)
  } else b[1]

  if (identical(target, "total")) z_in_model <- setdiff(z_in_model, "Z1")

  list(
    intercept = 0.0,
    b = b, gamma = gamma,
    z_role = z_role, z_in_model = z_in_model,
    # Strengths of the new structural arrows. Sized against gamma[1] = 0.50 and
    # b[1] = 0.40 so each role is a real feature of the DGP, not a nudge.
    delta_zx = delta_zx,    # fork/mixed:     Z1 -> X1
    delta_xz = delta_xz,    # pipe/collider:  X1 -> Z1
    delta_yz = delta_yz,    # collider:        Y -> Z1
    delta_xz2 = delta_xz2,  # mixed:          X1 -> Z2 on the logit scale
    # pipe_nl: Z1 = a*tanh(b*logX1) + c*(logX1^2 - 1) + noise. `tanh` gives
    # saturation and the square gives curvature -- neither representable by a
    # linear conditional, both learnable by BART.
    # V23 sweeps these, so a scenario may override them; NULL keeps V22's values.
    nl_a = nl_a_v, nl_b = nl_b_v, nl_c = nl_c_v, nl_sd = 0.80,
    # Under "pipe" the adjusted analysis estimates the DIRECT effect, which is
    # b[1]. Recorded here so the distinction is in the object rather than only in
    # a comment -- it is not the estimand, and must never be swapped in as one.
    total_effect = total_effect, target = target,
    b_int = b_int, b_quad = b_quad,
    erf_form = erf_form,
    y_form = y_form, b_zq = b_zq,
    # The primary Phase-1 estimand: the focal exposure's main-effect coefficient.
    # For the mixture surface this is the *local* main effect at logX2 = 0,
    # logX1 = 0 (the point about which b_int / b_quad are centred), so the
    # matched analysis model (dgp_formula) recovers it as the `logX1` coefficient.
    estimand_name = "b_logX1",
    # V24: the estimand FOLLOWS the target. Before V24 this was always b[1]
    # because every cell adjusted for everything; a "total" cell that kept
    # b[1] here would score a correct answer as a 75% bias.
    estimand_true = if (identical(target, "total")) total_effect else b[1]
  )
}

#' The non-linear `X1 -> Z1` arrow used by `z_role = "pipe_nl"`.
#'
#' Defined here rather than inline so the DGP and `exact_fork.R`'s exact
#' conditional cannot drift apart -- if they did, the "exact" arm would be wrong
#' and would still return believable numbers.
ef_nl_g <- function(x, truth) {
  truth$nl_a * tanh(truth$nl_b * x) + truth$nl_c * (x^2 - 1)
}

#' The part of the non-linear arrow `g()` that a LINEAR conditional cannot
#' represent over the censored region -- V23's derived predictor.
#'
#' WHY THIS QUANTITY. `leftcens` draws a censored exposure from a conditional
#' LINEAR in its predictors. The true conditional for a censored logX1 contains
#' `p(Z1 | logX1) = N(Z1; g(logX1), s^2)`, whose log contributes
#' `-(Z1 - g(logX1))^2 / 2s^2` -- non-linear in logX1 wherever `g` is curved. So
#' the damage is not "how big is g" but **how much of g a straight line cannot
#' express, where the censored mass actually is**: the density-weighted residual
#' SD of `g` after its best linear fit below the LOD.
#'
#' WHY IT SHOULD ENTER SQUARED. That residual is by construction ORTHOGONAL (in
#' the density-weighted L2 sense) to the span of the linear predictors. A
#' first-order expansion of the bias functional pairs the misspecification with
#' the score, and orthogonality kills that term -- so the leading contribution is
#' second order. Hence `bias ~ u^2` rather than `~ u`. This is a heuristic
#' argument, not a proof, and V23 is its test.
#'
#' @param truth A list from `make_truth()` with `nl_a`, `nl_b`, `nl_c`.
#' @param nd_frac The censored fraction, which sets where the LOD falls.
ef_unrep_curvature <- function(truth, nd_frac = 0.40, n_grid = 4000L) {
  L  <- stats::qnorm(nd_frac)
  xs <- seq(L - 5, L, length.out = n_grid)
  wt <- stats::dnorm(xs); wt <- wt / sum(wt)
  g  <- truth$nl_a * tanh(truth$nl_b * xs) + truth$nl_c * (xs^2 - 1)
  if (stats::sd(g) == 0) return(0)
  sqrt(sum(wt * stats::residuals(stats::lm(g ~ xs, weights = wt))^2))
}

#' Simulate one complete (uncensored, fully observed) dataset from the ERF.
#'
#' @param n Sample size.
#' @param truth A list from [make_truth()].
#' @param rho Exchangeable correlation among the log-exposures.
#' @param sd_x Marginal SD of each log-exposure.
#' @param mu_x Marginal mean of each log-exposure (sets the censoring geometry).
#' @param skew sinh-arcsinh skewness applied to the log-exposures (0 = none).
#' @param sigma_y Residual SD of the outcome.
#' @return A list: `data` (data.frame with Y, logX1.., Z1, Z2), `truth`.
#' @param z_form Structure of the covariates. "linear" (default) reproduces the
#'   original design: Z1 ~ N(0,1) and Z2 ~ Bern(0.5), independent of the
#'   exposures. "nonlinear" makes Z1 a non-linear, interacting function of the
#'   NON-focal exposures and Z2.
#'
#'   WHY "nonlinear" EXISTS (Track 05). Every covariate in the original design is
#'   linear and Gaussian, so a *parametric* imputation model for Z1 is correctly
#'   specified and is never penalised for misspecification. That makes it
#'   impossible to compare a flexible imputer (random forest) against a
#'   parametric one on equal terms -- the comparison can only reward calibration,
#'   never punish getting the conditional mean wrong. Under "nonlinear", a linear
#'   imputation model for Z1 IS misspecified while a forest can capture the
#'   structure, so the two families finally face a fair test.
#'
#'   The non-linearity is deliberately built from the non-focal exposures
#'   (logX2, logX3) and Z2 -- never the focal logX1 -- so the focal estimand's
#'   geometry is untouched. The OUTCOME model stays exactly linear in
#'   (logX, Z), so `dgp_formula()` remains correctly specified and `b[1]` is
#'   still recovered without bias by every procedure. Only the imputation model
#'   for Z1 becomes hard.
simulate_complete <- function(n, truth, rho = 0.4, sd_x = 1.0, mu_x = 0.0,
                              skew = 0.0, sigma_y = 1.0, z_form = "linear",
                              z_role = NULL) {
  p <- length(truth$b)
  q <- length(truth$gamma)

  # ---- V17: causal roles for the covariate -----------------------------------
  # These need a DIFFERENT GENERATION ORDER from the original design (a fork
  # draws Z1 before the exposures; a collider draws it after Y), so they live in
  # their own branch. The default path below is byte-identical to the pre-V17
  # code -- same calls, same order, same RNG stream -- which is what keeps every
  # earlier track's results reproducible from their recorded seeds.
  z_role <- z_role %||% truth$z_role %||% "precision"
  if (z_role %in% c("fork", "pipe", "pipe_nl", "collider", "mixed")) {
    return(.simulate_causal_z(n, truth, rho = rho, sd_x = sd_x, mu_x = mu_x,
                              skew = skew, sigma_y = sigma_y, z_role = z_role))
  }

  # --- log-exposures: exchangeable-correlation MVN, optional skew --------------
  Sigma <- matrix(rho, p, p); diag(Sigma) <- 1
  Z0 <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  if (!is.matrix(Z0)) Z0 <- matrix(Z0, ncol = p)
  if (skew != 0) Z0 <- sinh(asinh(Z0) + skew)          # sinh-arcsinh skew
  logX <- mu_x + sd_x * Z0
  colnames(logX) <- paste0("logX", seq_len(p))

  # --- covariates: one continuous, one binary ---------------------------------
  z2 <- stats::rbinom(n, 1, 0.5)

  if (identical(z_form, "nonlinear")) {
    # Z1 as a non-linear, interacting function of the NON-focal exposures and Z2.
    # tanh gives saturation, the square gives curvature, the product gives an
    # interaction -- none of which a linear imputation model can represent.
    x2 <- if (p >= 2L) logX[, 2] else rep(0, n)
    x3 <- if (p >= 3L) logX[, 3] else rep(0, n)
    z1_raw <- 0.9 * tanh(1.8 * x2) +
              0.45 * (x3^2 - 1) +
              -0.7 * x2 * z2 +
              0.4 * z2 +
              stats::rnorm(n, 0, 0.6)
    # Standardise so Z1 keeps a unit-SD scale and gamma retains its meaning,
    # making the linear and non-linear designs comparable.
    z1 <- as.vector(scale(z1_raw))
  } else {
    z1 <- stats::rnorm(n)
  }

  Z <- data.frame(Z1 = z1, Z2 = z2)
  Zmat <- as.matrix(Z[, seq_len(q), drop = FALSE])

  # --- outcome ----------------------------------------------------------------
  eta <- truth$intercept + as.vector(logX %*% truth$b) + as.vector(Zmat %*% truth$gamma)
  if (truth$erf_form == "mixture") {
    if (p >= 2L) eta <- eta + truth$b_int * logX[, 1] * logX[, 2]
    eta <- eta + truth$b_quad * logX[, 1]^2
  }
  # V11: outcome curvature in the continuous covariate. Centred, so the mean of
  # Y is unchanged and the linear and non-linear designs stay comparable.
  if (identical(truth$y_form %||% "linear", "nonlinear") && q >= 1L) {
    eta <- eta + (truth$b_zq %||% 0) * (Zmat[, 1]^2 - 1)
  }

  Y <- eta + stats::rnorm(n, 0, sigma_y)

  data <- data.frame(Y = Y, logX, Z, check.names = FALSE)
  list(data = data, truth = truth)
}

#' Simulate a complete dataset in which the covariate has a CAUSAL ROLE (V17).
#'
#' Separate from `simulate_complete()`'s main body because these designs need a
#' different generation order -- a confounder is drawn before the exposures, a
#' collider after the outcome -- and mixing the orders into one code path would
#' change the RNG stream of every pre-V17 scenario.
#'
#' The exposures keep their exchangeable correlation and the outcome keeps its
#' linear form, so `dgp_formula()` stays correctly specified and `b[1]` stays the
#' estimand in every role. What changes is only which arrows exist:
#'
#'   fork      Z1 -> logX1 (delta_zx),  Z1 -> Y (gamma[1])
#'             Adjusting for Z1 is REQUIRED; omitting it biases b[1].
#'   pipe      logX1 -> Z1 (delta_xz),  Z1 -> Y (gamma[1])
#'             Adjusting for Z1 blocks the indirect path, so the `logX1`
#'             coefficient is the DIRECT effect -- which is exactly b[1]. The
#'             total effect, b[1] + delta_xz * gamma[1], is reported in `truth`
#'             as `total_effect` but is NOT the estimand.
#'   collider  logX1 -> Z1 <- Y (delta_xz, delta_yz), gamma[1] forced to 0.
#'             Adjusting for Z1 induces bias; make_truth() therefore drops Z1
#'             from `z_in_model` and the analysis model omits it.
#'   mixed     Z1 is a fork and Z2 is a pipe (binary, on the logit scale), so a
#'             single adjustment set cannot be right for both roles at once.
.simulate_causal_z <- function(n, truth, rho = 0.4, sd_x = 1.0, mu_x = 0.0,
                               skew = 0.0, sigma_y = 1.0, z_role = "fork") {
  p <- length(truth$b); q <- length(truth$gamma)
  Sigma <- matrix(rho, p, p); diag(Sigma) <- 1

  draw_x <- function(shift1) {
    Z0 <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    if (!is.matrix(Z0)) Z0 <- matrix(Z0, ncol = p)
    if (skew != 0) Z0 <- sinh(asinh(Z0) + skew)
    lx <- mu_x + sd_x * Z0
    lx[, 1] <- lx[, 1] + shift1
    colnames(lx) <- paste0("logX", seq_len(p))
    lx
  }
  eta_x <- function(lx) truth$intercept + as.vector(lx %*% truth$b)

  if (z_role == "fork") {
    z1   <- stats::rnorm(n)                                   # Z1 first
    z2   <- stats::rbinom(n, 1, 0.5)
    logX <- draw_x(truth$delta_zx * z1)                       # Z1 -> logX1
    Zm   <- cbind(Z1 = z1, Z2 = z2)[, seq_len(q), drop = FALSE]
    Y    <- eta_x(logX) + as.vector(Zm %*% truth$gamma) + stats::rnorm(n, 0, sigma_y)

  } else if (z_role == "pipe") {
    z2   <- stats::rbinom(n, 1, 0.5)
    logX <- draw_x(0)
    z1   <- truth$delta_xz * logX[, 1] + stats::rnorm(n, 0, 0.8)   # logX1 -> Z1
    Zm   <- cbind(Z1 = z1, Z2 = z2)[, seq_len(q), drop = FALSE]
    Y    <- eta_x(logX) + as.vector(Zm %*% truth$gamma) + stats::rnorm(n, 0, sigma_y)

  } else if (z_role == "collider") {
    z2   <- stats::rbinom(n, 1, 0.5)
    logX <- draw_x(0)
    # gamma[1] is 0 here (make_truth enforces it), so Z1 is absent from Y.
    Zm0  <- cbind(Z1 = 0, Z2 = z2)[, seq_len(q), drop = FALSE]
    Y    <- eta_x(logX) + as.vector(Zm0 %*% truth$gamma) + stats::rnorm(n, 0, sigma_y)
    z1   <- truth$delta_xz * logX[, 1] + truth$delta_yz * Y +      # X1 -> Z1 <- Y
            stats::rnorm(n, 0, 0.8)

  } else if (z_role == "pipe_nl") {
    z2   <- stats::rbinom(n, 1, 0.5)
    logX <- draw_x(0)
    z1   <- ef_nl_g(logX[, 1], truth) + stats::rnorm(n, 0, truth$nl_sd)
    Zm   <- cbind(Z1 = z1, Z2 = z2)[, seq_len(q), drop = FALSE]
    Y    <- eta_x(logX) + as.vector(Zm %*% truth$gamma) + stats::rnorm(n, 0, sigma_y)

  } else {                                                    # mixed: fork + pipe
    z1   <- stats::rnorm(n)                                   # fork, drawn first
    logX <- draw_x(truth$delta_zx * z1)
    z2   <- stats::rbinom(n, 1, stats::plogis(truth$delta_xz2 * logX[, 1]))  # pipe
    Zm   <- cbind(Z1 = z1, Z2 = z2)[, seq_len(q), drop = FALSE]
    Y    <- eta_x(logX) + as.vector(Zm %*% truth$gamma) + stats::rnorm(n, 0, sigma_y)
  }

  data <- data.frame(Y = Y, logX, Z1 = z1, Z2 = z2, check.names = FALSE)
  list(data = data, truth = truth)
}

#' The matched analysis-model formula for a given ERF form.
#'
#' Every procedure fits *this* model, so all differences are attributable to how
#' missingness was handled, not to model form.
dgp_formula <- function(truth) {
  p <- length(truth$b); q <- length(truth$gamma)
  xterms <- paste0("logX", seq_len(p))
  # V17: a collider must be left OUT of the analysis model -- adjusting for it
  # induces bias that no imputation method can undo. `z_in_model` defaults to the
  # full set, so every pre-V17 truth object produces the same formula as before.
  zterms <- truth$z_in_model %||% c("Z1", "Z2")[seq_len(q)]
  rhs <- c(xterms, zterms)
  if (truth$erf_form == "mixture") {
    if (p >= 2L) rhs <- c(rhs, "logX1:logX2")
    rhs <- c(rhs, "I(logX1^2)")
  }
  # V11: keep the ANALYSIS model correctly specified when the outcome is
  # non-linear in Z1, so the focal estimand stays recoverable and any bias is
  # attributable to the imputation model alone.
  if (identical(truth$y_form %||% "linear", "nonlinear") && q >= 1L &&
      "Z1" %in% zterms) {
    rhs <- c(rhs, "I(Z1^2)")
  }
  stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")))
}
