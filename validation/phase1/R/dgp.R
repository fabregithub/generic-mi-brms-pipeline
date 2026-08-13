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

#' Draw the true regression coefficients that define the ERF.
#'
#' Returned separately from the data so the "truth" is explicit and auditable.
#'
#' @param p Number of exposures.
#' @param q Number of covariates (here fixed at 2: one continuous, one binary).
#' @param erf_form "additive" or "mixture".
make_truth <- function(p = 3L, q = 2L, erf_form = "additive") {
  b <- rep(0.0, p)
  b[1] <- 0.40                      # focal exposure main effect (the estimand)
  if (p >= 2L) b[2] <- 0.20
  if (p >= 3L) b[-(1:2)] <- 0.10

  gamma <- c(0.50, -0.30)[seq_len(q)]

  # Interaction / curvature used only by the mixture surface.
  b_int  <- if (erf_form == "mixture" && p >= 2L) 0.25 else 0.0   # logX1 * logX2
  b_quad <- if (erf_form == "mixture") 0.15 else 0.0              # logX1^2

  list(
    intercept = 0.0,
    b = b, gamma = gamma,
    b_int = b_int, b_quad = b_quad,
    erf_form = erf_form,
    # The primary Phase-1 estimand: the focal exposure's main-effect coefficient.
    # For the mixture surface this is the *local* main effect at logX2 = 0,
    # logX1 = 0 (the point about which b_int / b_quad are centred), so the
    # matched analysis model (dgp_formula) recovers it as the `logX1` coefficient.
    estimand_name = "b_logX1",
    estimand_true = b[1]
  )
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
simulate_complete <- function(n, truth, rho = 0.4, sd_x = 1.0, mu_x = 0.0,
                              skew = 0.0, sigma_y = 1.0) {
  p <- length(truth$b)
  q <- length(truth$gamma)

  # --- log-exposures: exchangeable-correlation MVN, optional skew --------------
  Sigma <- matrix(rho, p, p); diag(Sigma) <- 1
  Z0 <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  if (!is.matrix(Z0)) Z0 <- matrix(Z0, ncol = p)
  if (skew != 0) Z0 <- sinh(asinh(Z0) + skew)          # sinh-arcsinh skew
  logX <- mu_x + sd_x * Z0
  colnames(logX) <- paste0("logX", seq_len(p))

  # --- covariates: one continuous, one binary ---------------------------------
  Z <- data.frame(Z1 = stats::rnorm(n), Z2 = stats::rbinom(n, 1, 0.5))
  Zmat <- as.matrix(Z[, seq_len(q), drop = FALSE])

  # --- outcome ----------------------------------------------------------------
  eta <- truth$intercept + as.vector(logX %*% truth$b) + as.vector(Zmat %*% truth$gamma)
  if (truth$erf_form == "mixture") {
    if (p >= 2L) eta <- eta + truth$b_int * logX[, 1] * logX[, 2]
    eta <- eta + truth$b_quad * logX[, 1]^2
  }
  Y <- eta + stats::rnorm(n, 0, sigma_y)

  data <- data.frame(Y = Y, logX, Z, check.names = FALSE)
  list(data = data, truth = truth)
}

#' The matched analysis-model formula for a given ERF form.
#'
#' Every procedure fits *this* model, so all differences are attributable to how
#' missingness was handled, not to model form.
dgp_formula <- function(truth) {
  p <- length(truth$b); q <- length(truth$gamma)
  xterms <- paste0("logX", seq_len(p))
  zterms <- c("Z1", "Z2")[seq_len(q)]
  rhs <- c(xterms, zterms)
  if (truth$erf_form == "mixture") {
    if (p >= 2L) rhs <- c(rhs, "logX1:logX2")
    rhs <- c(rhs, "I(logX1^2)")
  }
  stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")))
}
