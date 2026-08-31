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


#' Outcome mean for a set of rows, under either the true coefficients ("oracle")
#' or the per-sweep posterior draw ("plugin"). Shared by the X and Z blocks so
#' the two cannot drift apart in what surface they assume.
.smc_eta <- function(dat, truth, mode, coefs, fm, xcols, zcols) {
  if (mode == "oracle") {
    U  <- as.matrix(dat[, xcols, drop = FALSE])
    Zm <- as.matrix(dat[, zcols, drop = FALSE])
    out <- as.vector(U %*% truth$b) + as.vector(Zm %*% truth$gamma) + truth$intercept
    if (identical(truth$erf_form, "mixture")) {
      if (length(truth$b) >= 2L) out <- out + truth$b_int * U[, 1] * U[, 2]
      out <- out + truth$b_quad * U[, 1]^2
    }
    if (identical(truth$y_form %||% "linear", "nonlinear"))
      out <- out + (truth$b_zq %||% 0) * (Zm[, 1]^2 - 1)
    out
  } else {
    MM <- stats::model.matrix(stats::delete.response(stats::terms(fm)), data = dat)
    as.vector(MM %*% coefs[colnames(MM)])
  }
}

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
#' @param z_mode NULL (default) leaves covariates untouched -- correct when they
#'   are complete, as in V9's cells. "exact" or "gaussian" additionally imputes
#'   missing covariates, which is what V12 contrasts: both use the same correct
#'   conditional mean, and differ only in whether the DRAW respects the true
#'   spread and skew.
smc_impute_datasets <- function(d, truth, m = 10L, mode = c("oracle", "plugin"),
                                sweeps = 3L, rho = 0.4, sd_x = 1, mu_x = 0,
                                sigma_y = 1, n_grid = 512L, z_mode = NULL) {
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

    # Which covariates need imputing (V12 only; NULL z_mode leaves them alone).
    zmiss <- if (is.null(z_mode)) character(0) else
      zcols[vapply(zcols, function(z) anyNA(d[[z]]), logical(1))]
    z_na  <- lapply(stats::setNames(zmiss, zmiss), function(z) is.na(d[[z]]))
    for (z in zmiss) {                       # start from the observed mean/mode
      w[[z]][z_na[[z]]] <- if (length(unique(stats::na.omit(d[[z]]))) <= 2)
        as.numeric(stats::median(d[[z]], na.rm = TRUE))
      else mean(d[[z]], na.rm = TRUE)
    }

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

      # --- Z block (V12): draw each missing covariate from its conditional -----
      for (z in zmiss) {
        idx <- which(z_na[[z]]); if (!length(idx)) next
        sub <- w[idx, , drop = FALSE]
        eta_z <- function(zv) { s2 <- sub; s2[[z]] <- zv; .smc_eta(s2, truth, mode, coefs, fm, xcols, zcols) }
        if (length(unique(stats::na.omit(d[[z]]))) <= 2) {
          w[[z]][idx] <- .smc_draw_z_binary(sub$Y, eta_z(rep(0, length(idx))),
                                            eta_z(rep(1, length(idx))), sig)
        } else {
          qc <- .smc_quadratic_coefs(eta_z, length(idx))
          w[[z]][idx] <- .smc_draw_z(sub$Y, qc$A, qc$B, qc$C, sig,
                                     m = rep(0, length(idx)), s = 1,
                                     mode = z_mode, n_grid = n_grid)
        }
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
          .smc_eta(s2, truth, mode, coefs, fm, xcols, zcols)
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


# =============================================================================
# V12: the Z-block draw -- shape versus mean
# -----------------------------------------------------------------------------
# WHAT V11 LEFT UNEXPLAINED. A non-linear outcome cost ~3 pp of bias for EVERY
# imputer -- BART, mice pmm, forest and bootstrapped forest alike, spread only
# 0.44 pp. That is odd if the problem were the conditional MEAN, because those
# arms differ enormously in how flexibly they model it.
#
# THE DIAGNOSIS. Every one of them draws Z as
#
#     (fitted conditional mean)  +  homoscedastic Gaussian noise
#
# and under a linear outcome that is EXACTLY right: p(Z1 | Y, X) is Gaussian
# with constant SD 0.894 and zero skew at every Y. Under a non-linear outcome
# the same conditional becomes heteroscedastic and skewed -- SD ranging 0.63 to
# 1.22 across Y, skew reaching -1.34. The imputers share the part that is wrong
# (the noise) and differ only in the part that is not (the mean), which is
# exactly the pattern V11 measured.
#
# THE EXPERIMENT. Two draws that use the SAME, CORRECT conditional mean and
# differ ONLY in shape:
#
#   z_mode = "gaussian"  mean + homoscedastic Gaussian noise -- what every
#                        shipped imputer effectively does, but with the mean
#                        handed to it exactly, so no mean-modelling error remains
#   z_mode = "exact"     a draw from the true conditional, spread and skew
#                        included
#
# If "exact" removes the ~3 pp and "gaussian" does not, the defect is the SHAPE
# of the draw, and the fix is a shape-aware Z-block sampler rather than a better
# mean model. If neither helps, the cause is elsewhere again.
#
# NOTE this is a DIFFERENT defect from the mixture failure of V3/V9. That one is
# the X block drawing the exposure from a conditional with the wrong functional
# FORM. This one is the Z block drawing a covariate with the right mean and the
# wrong SHAPE. They are not the same bug and one fix will not close both.
# =============================================================================

#' Grid-based draw from p(Z | Y, rest), or its mean-matched Gaussian counterfeit.
#'
#' The outcome mean is quadratic in Z (`A + B Z + C Z^2`), so the log target is a
#' quartic plus the Gaussian prior -- no closed form, but one-dimensional.
#'
#' @param mode "exact" samples the true conditional. "gaussian" computes that
#'   conditional's mean exactly and then draws `N(mean, sd_pooled)` with a single
#'   pooled SD, mimicking a mean-model imputer that adds homoscedastic noise.
#' @return Numeric vector of draws.
.smc_draw_z <- function(Y, A, B, C, sigma_y, m, s, mode = c("exact","gaussian"),
                        n_grid = 512L) {
  mode <- match.arg(mode)
  n <- length(Y)
  if (n == 0L) return(numeric(0))
  lo <- m - 8 * s; hi <- m + 8 * s
  u <- seq(0, 1, length.out = n_grid)
  X <- outer(lo, rep(1, n_grid)) + outer(hi - lo, u)
  mu <- A + B * X + C * X^2
  logp <- -(Y - mu)^2 / (2 * sigma_y^2) - (X - m)^2 / (2 * s^2)
  logp <- logp - apply(logp, 1, max)
  dens <- exp(logp)
  w <- X[, 2] - X[, 1]
  tot <- rowSums(dens) * w
  ok <- is.finite(tot) & tot > 0

  # Conditional mean and SD on the grid -- both modes need the mean; only
  # "gaussian" discards everything else about the shape.
  pw <- dens / ifelse(rowSums(dens) > 0, rowSums(dens), 1)
  cmean <- rowSums(pw * X)
  csd   <- sqrt(pmax(rowSums(pw * (X - cmean)^2), 1e-12))

  if (mode == "gaussian") {
    # One pooled SD for all rows: that is precisely the homoscedasticity
    # assumption every shipped imputer makes.
    sd_pool <- sqrt(mean(csd^2))
    out <- stats::rnorm(n, cmean, sd_pool)
    out[!ok] <- m[!ok]
    return(out)
  }

  cdf <- t(apply(dens, 1, cumsum)) - dens / 2
  cdf <- cdf * w
  cdf <- cdf / ifelse(tot > 0, tot, 1)
  out <- numeric(n); target <- stats::runif(n)
  for (i in seq_len(n)) {
    if (!ok[i]) { out[i] <- m[i]; next }
    out[i] <- stats::approx(cdf[i, ], X[i, ], xout = target[i],
                            rule = 2, ties = "ordered")$y
  }
  out
}

#' Exact Bernoulli conditional for the binary covariate.
#'
#' Not the quantity under test -- Z2 is drawn exactly in both modes so that any
#' difference between them is attributable to the CONTINUOUS covariate's shape.
.smc_draw_z_binary <- function(Y, eta0, eta1, sigma_y, p_prior = 0.5) {
  l1 <- -(Y - eta1)^2 / (2 * sigma_y^2) + log(p_prior)
  l0 <- -(Y - eta0)^2 / (2 * sigma_y^2) + log(1 - p_prior)
  p <- 1 / (1 + exp(l0 - l1))
  stats::rbinom(length(Y), 1L, pmin(pmax(p, 0), 1))
}


#' V12 arm: SMC imputation, then the matched scalar model.
#'
#' Same shape as `proc_pipeline_block_fcs` but with both blocks drawn from the
#' generator's own conditionals, and with the Z-block draw switchable between
#' the true shape and a mean-matched Gaussian counterfeit. Everything downstream
#' -- the matched `lm`, Rubin pooling -- is shared with every other arm, so a
#' difference between `smc_zexact` and `smc_zgauss` is the shape of the Z draw
#' and nothing else.
#'
#' A HARNESS INSTRUMENT, not pipeline code. It measures what a shape-aware
#' Z-block sampler could achieve; it is not wired into `00_common_functions.R`.
proc_smc_scalar <- function(bundle, m = 30L, seed = NULL, sweeps = 3L,
                            z_mode = c("exact", "gaussian"), mode = "oracle",
                            rho = 0.4, sd_x = 1, mu_x = 0, sigma_y = 1,
                            label = NULL) {
  z_mode <- match.arg(z_mode)
  label  <- label %||% paste0("smc_z", if (z_mode == "exact") "exact" else "gauss")
  if (!is.null(seed)) set.seed(seed)
  truth <- bundle$truth

  imputed <- tryCatch(
    smc_impute_datasets(bundle$censored, truth, m = as.integer(m), mode = mode,
                        sweeps = sweeps, rho = rho, sd_x = sd_x, mu_x = mu_x,
                        sigma_y = sigma_y, z_mode = z_mode),
    error = function(e) structure(list(), err = conditionMessage(e)))
  if (!length(imputed))
    return(one_row(label, NA, NA, NA, NA,
                   substr(paste("smc failed:", attr(imputed, "err") %||% "none"), 1, 120)))

  ests <- vars <- rep(NA_real_, length(imputed))
  for (i in seq_along(imputed)) {
    e <- fit_lm_estimand(as.data.frame(imputed[[i]]), truth)
    ests[i] <- e["est"]; vars[i] <- e["se"]^2
  }
  pl <- rubin_pool(ests, vars)
  n_bad <- sum(!is.finite(ests))
  note <- sprintf("smc-%s; z=%s; m=%d; sweeps=%d", mode, z_mode, m, sweeps)
  if (n_bad) note <- paste0(note, sprintf("; %d/%d fits failed", n_bad, length(ests)))
  one_row(label, pl["est"], pl["se"], pl["ci_lo"], pl["ci_hi"], note,
          ubar = pl["ubar"], b = pl["b"], fmi = pl["fmi"])
}
