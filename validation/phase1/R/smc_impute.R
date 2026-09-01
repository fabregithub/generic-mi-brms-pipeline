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
                                sigma_y = 1, n_grid = 512L, z_mode = NULL,
                                x_mode = c("exact", "shipped", "grid"), margin = "shash") {
  mode   <- match.arg(mode)
  x_mode <- match.arg(x_mode)
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

        if (x_mode == "grid") {
          # ITEM 07's candidate: the general grid draw. Unlike "exact" it does
          # NOT know the generator's coefficients -- it is handed the analysis
          # FORMULA and estimates the rest, which is what a shipped version
          # would have. V9 showed that suffices (plug-in matched oracle to 0.2 pp).
          ll <- function(dat) stats::dnorm(
            dat$Y, .smc_eta(dat, truth, mode, coefs, fm, xcols, zcols), sig, log = TRUE)
          gr <- .ce_smc_x_grid(w, cc, cen = is_cens[[cc]], upper = upper[[cc]],
                               preds = setdiff(c(xcols, zcols), cc),
                               loglik = ll, n_grid = n_grid, margin = margin)
          if (isTRUE(gr$ok) && length(gr$x)) w[[cc]][idx] <- gr$x
        } else if (x_mode == "shipped") {
          # The SHIPPED exposure draw: leftcens's linear conditional, exactly as
          # `.ce_one_imputation()` calls it. Holding this fixed while the Z draw
          # varies is what separates the X block's contribution from the Z
          # block's -- V12 could not, because both differed at once.
          preds <- setdiff(c("Y", xcols, zcols), cc)
          preds <- intersect(preds, names(w))
          # leftcens's convention, per `.ce_exposure_bounds()` in
          # 00_censored_exposure.R: an OBSERVED row supplies `y` and has NA
          # bounds; a CENSORED row supplies bounds and has y = NA. Passing
          # filled values *and* point bounds (the intuitive reading) silently
          # produces a different, much worse draw -- caught by a smoke test in
          # which this arm was 5 pp off in a cell where it should have matched.
          cen <- is_cens[[cc]]
          yv <- w[[cc]]; yv[cen] <- NA_real_
          lo <- rep(NA_real_, nrow(w)); hi <- rep(NA_real_, nrow(w))
          lo[cen] <- -Inf
          hi[cen] <- upper[[cc]][cen]
          dr <- tryCatch(
            leftcens::impute_censored_conditional(
              y = yv, x = w[, preds, drop = FALSE], lower = lo, upper = hi,
              m = 1L, margin = margin)[, 1],
            error = function(e) NULL)
          if (!is.null(dr)) w[[cc]][idx] <- dr[idx]
        } else {
          qc <- .smc_quadratic_coefs(eta_fun, length(idx))
          other <- as.matrix(sub[, setdiff(xcols, cc), drop = FALSE])
          pr <- .smc_exposure_prior(j, other, p, rho, sd_x, mu_x)
          w[[cc]][idx] <- .smc_draw_truncated(
            Y = sub$Y, A = qc$A, B = qc$B, C = qc$C, sigma_y = sig,
            m = pr$m, s = pr$s, upper = upper[[cc]][idx], n_grid = n_grid)
        }
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
                            x_mode = c("exact", "shipped", "grid"), margin = "shash",
                            rho = 0.4, sd_x = 1, mu_x = 0, sigma_y = 1,
                            label = NULL) {
  z_mode <- match.arg(z_mode); x_mode <- match.arg(x_mode)
  label  <- label %||% if (x_mode == "grid") "smc_xgrid" else
    paste0("smc_z", if (z_mode == "exact") "exact" else "gauss",
           if (x_mode == "shipped") "_xship" else "")
  if (!is.null(seed)) set.seed(seed)
  truth <- bundle$truth

  imputed <- tryCatch(
    smc_impute_datasets(bundle$censored, truth, m = as.integer(m), mode = mode,
                        sweeps = sweeps, rho = rho, sd_x = sd_x, mu_x = mu_x,
                        sigma_y = sigma_y, z_mode = z_mode, x_mode = x_mode,
                        margin = margin),
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
  note <- sprintf("smc-%s; z=%s; x=%s; m=%d; sweeps=%d", mode, z_mode, x_mode, m, sweeps)
  if (n_bad) note <- paste0(note, sprintf("; %d/%d fits failed", n_bad, length(ests)))
  one_row(label, pl["est"], pl["se"], pl["ci_lo"], pl["ci_hi"], note,
          ubar = pl["ubar"], b = pl["b"], fmi = pl["fmi"])
}


# =============================================================================
# Item 07: a GENERAL substantive-model-compatible exposure draw
# -----------------------------------------------------------------------------
# WHAT V9 AND V13 ESTABLISHED. Replacing the X block's linear conditional with
# the correct one restores the mixture estimands (V9: 89-100% of the gap) and
# accounts for 3.1-3.4 pp of the ~3.9 pp non-linear-outcome penalty (V13). Both
# used a sampler that KNEW the generator's surface. A shippable version cannot.
#
# WHY IMPORTANCE RESAMPLING AND NOT A GRID. The V9/V12 samplers evaluate the
# target on a grid, which needs the exposure prior's DENSITY. `leftcens` exposes
# draws (`impute_censored_conditional`) but no density, and its shash margin is
# not optional -- V1 found a Gaussian conditional gives coverage 0.47 on skewed
# exposures. Importance resampling needs only draws, so it can keep that margin:
#
#   1. PROPOSE  K candidates per censored cell from leftcens, conditioning on the
#               other exposures and covariates but NOT on Y -- so the proposal is
#               the (skew-aware, bound-respecting) exposure model alone.
#   2. WEIGHT   each candidate by the substantive model's likelihood
#               p(Y_i | x = cand, rest_i; theta).
#   3. RESAMPLE one candidate per cell with probability proportional to weight.
#
# The result is a draw from p(x | Y, rest) that never needs that density in
# closed form, and it is GENERAL: any analysis formula works, because the
# likelihood is only ever EVALUATED at candidate values. Splines, mo() terms and
# interactions all pass through untouched.
#
# THE FAILURE MODE IS WEIGHT DEGENERACY. If Y is very informative about x, the
# prior is a poor proposal, one candidate takes nearly all the weight, and the
# draw collapses to a point -- which understates between-imputation variance
# exactly as the improper imputation of V4 did. `ess` is therefore returned, not
# optional: it is the diagnostic that says whether K was large enough.
# =============================================================================

#' Substantive-model-compatible draw for one censored exposure, by importance
#' resampling against a leftcens proposal.
#'
#' @param w Current completed data (one imputation, mid-sweep).
#' @param cc Name of the exposure column being drawn.
#' @param cen Logical vector: which rows are censored.
#' @param upper Upper bound (log LOD) per row.
#' @param preds Predictor names for the PROPOSAL. `Y` must be excluded -- the
#'   outcome enters through the weights, not the proposal, and including it both
#'   places would double-count it.
#' @param loglik Function(dat) returning the per-row log-likelihood of the
#'   outcome under the substantive model. Receives a data frame with `cc` set to
#'   candidate values.
#' @param K Number of proposals per censored cell.
#' @param margin leftcens margin (`"shash"` keeps skew handling).
#' @param propose Optional `function(n_rows, K)` returning an `n_rows x K` matrix
#'   of candidates, used INSTEAD of leftcens. Exists so the weighting and
#'   resampling can be tested against a proposal whose target is known in closed
#'   form -- validating the machinery separately from the proposal model, which
#'   is the only way to tell which of the two is at fault when a result looks
#'   wrong.
#' @return A list: `x` (one draw per censored row) and `ess` (mean effective
#'   sample size, out of K -- the degeneracy diagnostic).
.smc_x_importance <- function(w, cc, cen, upper, preds, loglik, K = 50L,
                              margin = "shash", propose = NULL) {
  idx <- which(cen)
  if (!length(idx)) return(list(x = numeric(0), ess = NA_real_))
  if (is.null(propose)) stopifnot(!("Y" %in% preds))

  # 1. PROPOSE: K candidates per row, from the exposure model without Y.
  if (!is.null(propose)) {
    cand <- propose(length(idx), as.integer(K))
  } else {
    yv <- w[[cc]]; yv[cen] <- NA_real_
    lo <- rep(NA_real_, nrow(w)); hi <- rep(NA_real_, nrow(w))
    lo[cen] <- -Inf; hi[cen] <- upper[cen]
    prop <- tryCatch(
      leftcens::impute_censored_conditional(
        y = yv, x = w[, intersect(preds, names(w)), drop = FALSE],
        lower = lo, upper = hi, m = as.integer(K), margin = margin),
      error = function(e) NULL)
    if (is.null(prop)) return(list(x = rep(NA_real_, length(idx)), ess = NA_real_))
    cand <- as.matrix(prop)[idx, , drop = FALSE]        # length(idx) x K
  }

  # 2. WEIGHT: the substantive model's likelihood at each candidate.
  sub <- w[idx, , drop = FALSE]
  lw <- matrix(NA_real_, nrow = length(idx), ncol = K)
  for (k in seq_len(K)) {
    s2 <- sub; s2[[cc]] <- cand[, k]
    lw[, k] <- loglik(s2)
  }
  lw[!is.finite(lw)] <- -Inf
  lw <- lw - apply(lw, 1, max)
  wt <- exp(lw)
  rs <- rowSums(wt)
  ok <- is.finite(rs) & rs > 0
  wt[ok, ] <- wt[ok, , drop = FALSE] / rs[ok]

  # 3. RESAMPLE one candidate per row.
  out <- numeric(length(idx))
  for (i in seq_along(idx)) {
    if (!ok[i]) { out[i] <- cand[i, 1L]; next }         # degenerate: keep the proposal
    out[i] <- cand[i, sample.int(K, 1L, prob = wt[i, ])]
  }
  # Effective sample size: 1/sum(w^2). K means the outcome added nothing; 1 means
  # a single candidate took all the weight and the draw is effectively a point.
  ess <- mean(1 / rowSums(wt[ok, , drop = FALSE]^2), na.rm = TRUE)
  list(x = out, ess = ess)
}


# =============================================================================
# Item 07: the GRID exposure draw -- general, and robust where IS is not
# -----------------------------------------------------------------------------
# WHY NOT IMPORTANCE SAMPLING. The obvious approach -- propose from the exposure
# model, weight by the outcome likelihood -- fails, and fails quietly. With a
# non-linear exposure-response the equation mu(x) = Y can have a second root far
# from the linear solution: for mu(x) = 0.1 + 0.4x + 0.15x^2 and Y = 1.5 the
# roots are +2.0 and -4.667, and left-censoring admits only -4.667. As the
# outcome tightens the target migrates there (exact mean -0.76 -> -3.86 -> -4.62
# as sigma_y goes 1.0 -> 0.3 -> 0.1) while a Y-conditional proposal piles up near
# the bound. Measured: the draw was off by 4.6 while ESS/K reported 0.999.
# **ESS detects weight concentration, not proposal misplacement**, so it cannot
# be relied on to catch this.
#
# A GRID HAS NO PROPOSAL TO MISPLACE. It covers the admissible range by
# construction, so a distant or multimodal target is found automatically. The
# cost is one grid per censored cell per sweep -- comparable to the K = 200
# proposals importance sampling needed, and robust instead of fragile.
#
# THE PRIOR IS THE PIPELINE'S OWN EXPOSURE MODEL, not an approximation of it.
# `leftcens::impute_censored_conditional()` fits a shash margin, transforms to a
# standard-normal `z` scale via `x_to_z()`, fits an interval-censored Gaussian
# AFT (`survreg`) of `z` on the predictors, and draws the censored `z` from a
# truncated normal. This function reuses those same steps, so the only thing that
# changes is *what the draw is conditioned on*: leftcens conditions on Y
# linearly, and this conditions on Y through the substantive model's actual
# likelihood.
#
# WORKING IN z IS WHY THIS IS SIMPLE. On the z scale the prior is a plain
# truncated normal. Sampling z from `phi(z; mu_lin, s) * L(x(z))` and mapping
# back through `z_to_x()` yields exactly x ~ p(x | preds) * L(x): the Jacobians
# cancel between the two parameterisations, so no derivative of the shash
# transform is needed anywhere.
# =============================================================================

#' Substantive-model-compatible draw for one censored exposure, by grid
#' inverse-CDF on the shash-transformed scale.
#'
#' @param w Current completed data for one imputation, mid-sweep.
#' @param cc Exposure column being drawn.
#' @param cen Logical: which rows are censored.
#' @param upper Upper bound (log LOD, data scale) per censored row.
#' @param preds Predictor names for the EXPOSURE model. `Y` must be excluded --
#'   the outcome enters through `loglik`, and including it in both places would
#'   condition on it twice.
#' @param loglik `function(dat)` returning the per-row log-likelihood of the
#'   outcome under the substantive model, with `cc` set to candidate values.
#' @param n_grid Grid points per censored row.
#' @param margin Passed to the shash margin fit.
#' @param proper Draw the exposure model's parameters from their posterior, as
#'   the shipped path does. Leave TRUE: without it the imputation is improper and
#'   between-imputation variance is understated (the V4 defect).
#' @return list(`x` = one draw per censored row, `ok` = did the fit succeed).
.ce_smc_x_grid <- function(w, cc, cen, upper, preds, loglik, n_grid = 512L,
                           margin = "shash", proper = TRUE) {
  idx <- which(cen)
  if (!length(idx)) return(list(x = numeric(0), ok = TRUE))
  stopifnot(!("Y" %in% preds))
  if (!requireNamespace("survival", quietly = TRUE)) return(list(x = NULL, ok = FALSE))

  n  <- nrow(w)
  xv <- w[[cc]]; xv[cen] <- NA_real_          # pipeline convention: NA where censored
  lo <- rep(NA_real_, n); hi <- rep(NA_real_, n)
  lo[cen] <- -Inf; hi[cen] <- upper[cen]
  obs <- !cen

  # --- 1. the exposure model, exactly as the shipped path fits it -------------
  mfit <- tryCatch(leftcens:::fit_shash_margin(xv[obs], lo[cen], hi[cen]),
                   error = function(e) NULL)
  if (is.null(mfit)) return(list(x = NULL, ok = FALSE))
  mp <- if (isTRUE(proper)) leftcens::draw_margin(mfit)
        else list(mu = mfit$mu, sigma = mfit$sigma, eps = mfit$eps)

  z_obs <- leftcens::x_to_z(xv, mp$mu, mp$sigma, mp$eps)
  z_lo  <- leftcens::x_to_z(lo,  mp$mu, mp$sigma, mp$eps)
  z_hi  <- leftcens::x_to_z(hi,  mp$mu, mp$sigma, mp$eps)

  pn <- intersect(preds, names(w))
  fm <- stats::as.formula(paste(
    "survival::Surv(.t1, .t2, type = 'interval2') ~",
    if (length(pn)) paste(sprintf("`%s`", pn), collapse = " + ") else "1"))
  t1 <- z_obs; t2 <- z_obs
  t1[cen] <- ifelse(is.finite(z_lo[cen]), z_lo[cen], NA_real_)
  t2[cen] <- ifelse(is.finite(z_hi[cen]), z_hi[cen], NA_real_)
  dat <- cbind(data.frame(.t1 = t1, .t2 = t2), w[, pn, drop = FALSE])
  sr <- tryCatch(survival::survreg(fm, data = dat, dist = "gaussian"),
                 error = function(e) NULL)
  if (is.null(sr)) return(list(x = NULL, ok = FALSE))

  beta <- stats::coef(sr); k <- length(beta)
  if (isTRUE(proper)) {
    V <- stats::vcov(sr); par <- c(beta, log(sr$scale))
    drawn <- tryCatch(as.numeric(par + t(chol(V)) %*% stats::rnorm(length(par))),
                      error = function(e) par)
    beta_i <- drawn[seq_len(k)]; s_i <- exp(drawn[k + 1L])
  } else { beta_i <- beta; s_i <- sr$scale }
  mu_lin <- as.vector(stats::model.matrix(sr) %*% beta_i)

  # --- 2. grid over the ADMISSIBLE z range, per censored row ------------------
  ml <- mu_lin[idx]; zh <- z_hi[idx]
  zh[!is.finite(zh)] <- ml[!is.finite(zh)] + 8 * s_i
  lo_z <- pmin(ml - 8 * s_i, zh - 8 * s_i)
  u  <- seq(0, 1, length.out = n_grid)
  Zg <- outer(lo_z, rep(1, n_grid)) + outer(zh - lo_z, u)

  # --- 3. target = truncated-normal prior x substantive-model likelihood ------
  sub <- w[idx, , drop = FALSE]
  logL <- matrix(NA_real_, nrow = length(idx), ncol = n_grid)
  for (g in seq_len(n_grid)) {
    s2 <- sub
    s2[[cc]] <- leftcens::z_to_x(Zg[, g], mp$mu, mp$sigma, mp$eps)
    logL[, g] <- loglik(s2)
  }
  logp <- stats::dnorm(Zg, ml, s_i, log = TRUE) + logL
  logp[!is.finite(logp)] <- -Inf
  logp <- logp - apply(logp, 1, max)
  dens <- exp(logp)

  wgt <- Zg[, 2] - Zg[, 1]
  cdf <- t(apply(dens, 1, cumsum)) - dens / 2
  cdf <- cdf * wgt
  tot <- cdf[, n_grid]
  good <- is.finite(tot) & tot > 0
  cdf <- cdf / ifelse(tot > 0, tot, 1)

  zdraw <- numeric(length(idx)); target <- stats::runif(length(idx))
  for (i in seq_along(idx)) {
    if (!good[i]) { zdraw[i] <- ml[i]; next }
    zdraw[i] <- stats::approx(cdf[i, ], Zg[i, ], xout = target[i],
                              rule = 2, ties = "ordered")$y
  }
  list(x = leftcens::z_to_x(zdraw, mp$mu, mp$sigma, mp$eps), ok = TRUE)
}
