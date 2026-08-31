# =============================================================================
# Track V9 harness -- substantive-model-compatible (SMC) censored imputation
# -----------------------------------------------------------------------------
# WHAT THIS IS FOR. V3 concluded that the shipped X block destroys curvature
# because it draws the censored exposure from a conditional that is LINEAR in
# the predictors, and so imposes linearity on exactly the rows where curvature
# would appear. That conclusion is an INFERENCE from a pattern (curvature dies,
# interaction survives, LOD/sqrt(2) beats the congenial engine), not a
# demonstration. Nobody has shown that replacing the linear conditional fixes it.
#
# This file is that experiment. It draws the censored exposure from the CORRECT
# conditional under the generator's own surface. If curvature is restored, the
# V3 diagnosis is confirmed causally and roadmap item 07 has a defined target and
# a measured upper bound. If it is not restored, the diagnosis is wrong and 07
# would have been aimed at the wrong thing.
#
# THIS IS A HARNESS INSTRUMENT, NOT PIPELINE CODE. It does not call, and is not
# called by, anything in `00_censored_exposure.R`. V8 recorded what happens when
# an instrument is mistaken for the shipped path; this one is labelled so it
# cannot be. What it measures is what SMC imputation COULD achieve here, not what
# the pipeline does.
#
# THE CONDITIONAL. For a censored row, the exposure to be drawn is `x`, the
# others are held at their current values, and the target is
#
#     p(x | Y, other exposures, Z, x < LOD)
#       propto  p(Y | x, rest) * p(x | other exposures) * 1{x < LOD}
#
# The outcome term is Gaussian in Y with a mean that is QUADRATIC in x -- the
# generator's surface contains `b_quad * x1^2` and `b_int * x1 * x2`, and each is
# at most quadratic in any single coordinate. Writing that mean as
#
#     mu(x) = A + B x + C x^2
#
# makes the log target a QUARTIC polynomial in x plus the Gaussian prior term:
#
#     log p  =  -(Y - A - Bx - Cx^2)^2 / (2 sigma_y^2)  -  (x - m)^2 / (2 s^2)
#
# which has no closed-form sampler but is one-dimensional and smooth. A grid
# inverse-CDF over the truncated support is exact to grid resolution and costs
# nothing. `test_v9_smc.R` checks it against numerical integration.
#
# A, B and C are not hand-coded per exposure. Because the surface is quadratic in
# each coordinate, they are recovered by evaluating the mean function at x = 0,
# 1, -1 and solving:
#
#     A = mu(0);  B = (mu(1) - mu(-1))/2;  C = (mu(1) + mu(-1) - 2 mu(0))/2
#
# That is general, works for whichever exposure is being drawn, and works
# unchanged for the plug-in arm where the coefficients are estimated rather than
# known.
# =============================================================================

#' Recover (A, B, C) in mu(x) = A + B x + C x^2 for one exposure, one row set.
#'
#' @param eta_fun A function of a numeric vector `x` returning the outcome mean
#'   for the rows in question, with every other predictor held fixed.
#' @return A list of numeric vectors `A`, `B`, `C`, one entry per row.
.smc_quadratic_coefs <- function(eta_fun, n) {
  e0 <- eta_fun(rep(0, n))
  ep <- eta_fun(rep(1, n))
  em <- eta_fun(rep(-1, n))
  list(A = e0, B = (ep - em) / 2, C = (ep + em - 2 * e0) / 2)
}

#' Conditional mean and SD of one log-exposure given the others, under the
#' generator's exchangeable-correlation MVN.
#'
#' @param j Index of the exposure being drawn.
#' @param other A matrix of the OTHER exposures' current values (n x (p-1)).
#' @param rho,sd_x,mu_x Generator parameters.
.smc_exposure_prior <- function(j, other, p, rho, sd_x, mu_x) {
  if (p == 1L) {
    return(list(m = rep(mu_x, nrow(other)), s = sd_x))
  }
  S <- matrix(rho, p, p); diag(S) <- 1          # correlation of the log-exposures
  s11 <- S[j, j]; s12 <- S[j, -j, drop = FALSE]; S22 <- S[-j, -j, drop = FALSE]
  Ainv <- solve(S22)
  beta <- as.vector(s12 %*% Ainv)               # regression of x_j on the others
  cond_var <- as.numeric(s11 - s12 %*% Ainv %*% t(s12))
  ctr <- sweep(other, 2, mu_x, "-") / sd_x      # standardise the conditioning set
  m <- mu_x + sd_x * as.vector(ctr %*% beta)
  list(m = m, s = sd_x * sqrt(max(cond_var, 1e-12)))
}

#' Draw from p(x | Y, rest) truncated to x < upper, by grid inverse-CDF.
#'
#' Vectorised over rows: each row gets its own grid spanning its own prior, so a
#' row whose conditional sits far from the others is not sampled on a grid built
#' for someone else.
#'
#' @param Y Outcome values (n).
#' @param A,B,C Quadratic coefficients of the outcome mean in x (each length n).
#' @param sigma_y Residual SD of the outcome.
#' @param m,s Prior conditional mean (n) and SD (scalar) of x.
#' @param upper Upper truncation point per row (n) -- the log LOD.
#' @param n_grid Grid points per row.
#' @return A numeric vector of draws (n).
.smc_draw_truncated <- function(Y, A, B, C, sigma_y, m, s, upper, n_grid = 512L) {
  n <- length(Y)
  if (n == 0L) return(numeric(0))

  # Grid: from well below the prior mean up to the truncation point. The prior is
  # the only thing guaranteed to bound the target from the left, so span it
  # generously; the right edge is the LOD by construction.
  lo <- pmin(m - 8 * s, upper - 8 * s)
  hi <- upper
  # Degenerate rows (upper at or below lo) would give a zero-width grid.
  bad <- !is.finite(lo) | !is.finite(hi) | (hi <= lo)
  lo[bad] <- hi[bad] - 1e-6

  u <- seq(0, 1, length.out = n_grid)
  X <- outer(lo, rep(1, n_grid)) + outer(hi - lo, u)      # n x n_grid

  mu <- A + B * X + C * X^2                                # broadcast over columns
  logp <- -(Y - mu)^2 / (2 * sigma_y^2) - (X - m)^2 / (2 * s^2)
  logp <- logp - apply(logp, 1, max)                       # stabilise
  dens <- exp(logp)

  # Trapezoidal CDF along each row.
  w <- (X[, 2] - X[, 1])                                   # equal spacing per row
  cdf <- t(apply(dens, 1, cumsum)) - dens / 2
  cdf <- cdf * w
  tot <- cdf[, n_grid]
  ok <- is.finite(tot) & tot > 0
  cdf <- cdf / ifelse(tot > 0, tot, 1)

  # Inverse-CDF by linear interpolation, row by row.
  out <- numeric(n)
  target <- stats::runif(n)
  for (i in seq_len(n)) {
    if (!ok[i]) { out[i] <- m[i]; next }                   # degenerate: fall back to the prior mean
    out[i] <- stats::approx(x = cdf[i, ], y = X[i, ], xout = target[i],
                            rule = 2, ties = "ordered")$y
  }
  out
}

#' Produce `m` completed datasets by SMC imputation of the censored exposures.
#'
#' @param d The censored data.frame (from `inject_left_censoring`).
#' @param truth The generator's truth list; used directly by the "oracle" mode.
#' @param m Number of completed datasets.
#' @param mode "oracle" uses the generator's true coefficients -- the upper bound
#'   on what a correctly specified SMC imputation can achieve. "plugin" estimates
#'   them from the current completed data each sweep, drawing them from their
#'   posterior so the imputation stays PROPER (R8's standard, applied here).
#' @param sweeps FCS sweeps over the censored exposures.
#' @param rho,sd_x,mu_x,sigma_y Generator parameters for the exposure prior and
#'   the outcome residual. In "plugin" mode `sigma_y` is re-estimated.
#' @return A list of `m` completed data.frames.
smc_impute_datasets <- function(d, truth, m = 10L, mode = c("oracle", "plugin"),
                                sweeps = 3L, rho = 0.4, sd_x = 1, mu_x = 0,
                                sigma_y = 1, n_grid = 512L) {
  mode <- match.arg(mode)
  p <- length(truth$b)
  xcols <- paste0("logX", seq_len(p))
  zcols <- c("Z1", "Z2")
  cens_cols <- attr(d, "censored_cols")
  if (is.null(cens_cols)) cens_cols <- "logX1"

  is_cens <- lapply(stats::setNames(cens_cols, cens_cols), function(cc) is.na(d[[cc]]))
  upper   <- lapply(stats::setNames(cens_cols, cens_cols), function(cc) d[[paste0(cc, "_lod")]])

  # The matched analysis-model formula: the SAME surface the generator used. This
  # is what makes the imputation substantive-model-COMPATIBLE.
  fm <- dgp_formula(truth)

  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d
    # Initialise censored cells just below the LOD, so the first sweep starts
    # somewhere legal rather than at NA.
    for (cc in cens_cols) w[[cc]][is_cens[[cc]]] <- upper[[cc]][is_cens[[cc]]] - 0.5

    for (t in seq_len(sweeps)) {
      # --- parameters for this sweep -------------------------------------------
      if (mode == "oracle") {
        coefs <- NULL; sig <- sigma_y
      } else {
        fit <- stats::lm(fm, data = w)
        # Proper draw: sigma^2 from its scaled-inverse-chi-square posterior, then
        # beta | sigma^2 from its Gaussian posterior. Without this the arm would
        # be improper and its intervals would be too narrow for reasons that have
        # nothing to do with the functional form under test.
        df <- fit$df.residual
        sig <- sqrt(sum(stats::residuals(fit)^2) / stats::rchisq(1, df))
        V <- stats::summary.lm(fit)$cov.unscaled
        beta <- as.vector(stats::coef(fit) +
                          t(chol(V)) %*% stats::rnorm(length(stats::coef(fit))) * sig)
        names(beta) <- names(stats::coef(fit))
        coefs <- beta
      }

      # --- draw each censored exposure in turn ---------------------------------
      for (cc in cens_cols) {
        idx <- which(is_cens[[cc]])
        if (!length(idx)) next
        j <- match(cc, xcols)

        sub <- w[idx, , drop = FALSE]
        # eta as a function of THIS exposure, others held at current values.
        eta_fun <- function(xv) {
          s2 <- sub; s2[[cc]] <- xv
          if (mode == "oracle") {
            U <- as.matrix(s2[, xcols, drop = FALSE])
            Zm <- as.matrix(s2[, zcols, drop = FALSE])
            as.vector(U %*% truth$b) + as.vector(Zm %*% truth$gamma) +
              truth$intercept +
              (if (identical(truth$erf_form, "mixture"))
                 truth$b_int * U[, 1] * U[, 2] + truth$b_quad * U[, 1]^2 else 0)
          } else {
            MM <- stats::model.matrix(stats::delete.response(stats::terms(fm)), data = s2)
            as.vector(MM %*% coefs[colnames(MM)])
          }
        }

        qc <- .smc_quadratic_coefs(eta_fun, length(idx))
        other <- as.matrix(sub[, setdiff(xcols, cc), drop = FALSE])
        pr <- .smc_exposure_prior(j, other, p, rho, sd_x, mu_x)

        w[[cc]][idx] <- .smc_draw_truncated(
          Y = sub$Y, A = qc$A, B = qc$B, C = qc$C, sigma_y = sig,
          m = pr$m, s = pr$s, upper = upper[[cc]][idx], n_grid = n_grid)
      }
    }
    out[[i]] <- w
  }
  out
}
