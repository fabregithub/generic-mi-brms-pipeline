# =============================================================================
# V20 -- EXACT conditional draws for the causal-role DGPs (fork / pipe / precision)
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS RATHER THAN REUSING smc_impute.R.
# V12/V13's exact Z draw (`.smc_draw_z`) computes
#
#     p(Z1 | Y)  ~  p(Y | Z1, rest) * N(Z1 | 0, 1)
#
# -- an N(0,1) prior for Z1 and the outcome model as the only likelihood term.
# That is correct in the PRECISION structure every earlier track used, where Z1
# is exogenous and independent of the exposures. **Under a fork it is wrong**:
# Z1 -> logX1 with delta = 0.60, so logX1 is a CHILD of Z1 and carries real
# information about it. The correct conditional needs a p(logX1 | Z1) factor
# that `.smc_draw_z` does not have. `.smc_exposure_prior()` has the mirror-image
# gap: it derives p(logX_j | other X) from the exchangeable Sigma alone, ignoring
# the Z1 -> logX1 shift.
#
# Using those draws under a fork would have produced a plausible decomposition
# from an incorrect "exact" arm -- the same failure mode as the four leftcens
# incidents in FINDINGS_v13.md, where a wrong conditional returned believable
# numbers instead of an error.
#
# THE APPROACH. Given Z2, the fork, pipe and precision DGPs are all JOINTLY
# GAUSSIAN in (Z1, logX, Y). So the exact conditionals are available in closed
# form by conditioning a multivariate normal -- no grid, no importance weights,
# and exact rather than approximate. Each role differs only in the loading matrix
# that maps the independent primitives to (Z1, logX, Y), which is what
# `.ef_joint()` builds.
#
# Everything here is verified against known answers in test_v20_exact.R: because
# the joint is Gaussian, the true conditional mean IS the population regression
# and the true conditional SD IS its residual SD, so a large simulated sample
# checks the algebra directly.
# =============================================================================

#' Loading matrix and mean for (Z1, logX1..p, Y) given Z2, by causal role.
#'
#' Writes V = M %*% L + const with L a vector of INDEPENDENT primitives whose
#' covariance is `Lcov` (block diagonal). Cov(V) is then M %*% Lcov %*% t(M),
#' which is easier to get right than hand-deriving each covariance.
#'
#' @param truth A list from `make_truth()` (supplies b, gamma, deltas, z_role).
#' @param rho,sd_x,mu_x,sigma_y DGP knobs, matching `simulate_complete()`.
#' @return list(M, Lcov, const_fixed, gamma2) -- `const_fixed` excludes the
#'   Z2 term, which is row-specific and applied by the callers.
.ef_joint <- function(truth, rho = 0.4, sd_x = 1, mu_x = 0, sigma_y = 1) {
  p <- length(truth$b); b <- truth$b
  role <- truth$z_role %||% "precision"
  g1 <- truth$gamma[1]; g2 <- if (length(truth$gamma) >= 2) truth$gamma[2] else 0
  Sig <- matrix(rho, p, p); diag(Sig) <- 1
  e1 <- c(1, rep(0, p - 1))
  I  <- diag(p)

  if (role == "fork") {
    d <- truth$delta_zx
    # L = (Z1, U[1..p], eps)
    M <- rbind(c(1,            rep(0, p),  0),
               cbind(d * e1,   sd_x * I,   0),
               c(b[1] * d + g1, sd_x * b,  1))
    Lcov <- rbind(cbind(1, matrix(0, 1, p), 0),
                  cbind(matrix(0, p, 1), Sig, matrix(0, p, 1)),
                  cbind(0, matrix(0, 1, p), sigma_y^2))
    const <- c(0, rep(mu_x, p), truth$intercept + sum(b) * mu_x)

  } else if (role == "pipe") {
    d <- truth$delta_xz; s_nu <- 0.8
    # L = (U[1..p], nu, eps);  Z1 = d*logX1 + nu
    M <- rbind(c(d * sd_x * e1,               1,  0),
               cbind(sd_x * I,                0,  0),
               c(sd_x * b + g1 * d * sd_x * e1, g1, 1))
    Lcov <- rbind(cbind(Sig, matrix(0, p, 1), matrix(0, p, 1)),
                  cbind(matrix(0, 1, p), s_nu^2, 0),
                  cbind(matrix(0, 1, p), 0, sigma_y^2))
    const <- c(d * mu_x, rep(mu_x, p),
               truth$intercept + sum(b) * mu_x + g1 * d * mu_x)

  } else if (role == "precision") {
    M <- rbind(c(1,   rep(0, p), 0),
               cbind(rep(0, p), sd_x * I, 0),
               c(g1, sd_x * b,  1))
    Lcov <- rbind(cbind(1, matrix(0, 1, p), 0),
                  cbind(matrix(0, p, 1), Sig, matrix(0, p, 1)),
                  cbind(0, matrix(0, 1, p), sigma_y^2))
    const <- c(0, rep(mu_x, p), truth$intercept + sum(b) * mu_x)

  } else if (role == "pipe_nl") {
    # NO GAUSSIAN JOINT EXISTS for pipe_nl: Z1's mean is non-linear in logX1, so
    # (Z1, logX, Y) is not jointly normal. The exact Z1 conditional is still
    # closed form -- see ef_draw_z1_pipenl() -- but it is a product of two
    # Gaussians in Z1, not a marginalisation of a joint. Anything that needs the
    # joint (notably the exact logX1 draw) is therefore NOT available here, and
    # refusing is the point: a fabricated joint would give a wrong "exact" arm
    # that returned believable numbers.
    stop("z_role 'pipe_nl' has no Gaussian joint: use ef_draw_z1_pipenl() for ",
         "the Z1 conditional, and hold the X block at the shipped draw. The ",
         "exact logX1 draw is not available (its conditional is non-Gaussian).",
         call. = FALSE)
  } else {
    stop("exact_fork.R supports z_role fork / pipe / pipe_nl / precision; got ",
         role, call. = FALSE)
  }

  V <- M %*% Lcov %*% t(M)
  dimnames(V) <- list(c("Z1", paste0("logX", seq_len(p)), "Y"),
                      c("Z1", paste0("logX", seq_len(p)), "Y"))
  names(const) <- rownames(V)
  list(V = V, const = const, gamma2 = g2, p = p)
}

#' Conditional mean and SD of one component of a Gaussian, given the others.
#'
#' @param J From `.ef_joint()`.
#' @param target Name of the variable to draw.
#' @param given Named matrix/data.frame of the conditioning variables (one row
#'   per observation), columns named as in the joint.
#' @param z2 Vector of Z2 values (its effect enters Y's mean, row by row).
.ef_cond <- function(J, target, given, z2) {
  gn <- colnames(given)
  mu <- J$const
  mu_row <- matrix(mu, nrow = nrow(given), ncol = length(mu), byrow = TRUE,
                   dimnames = list(NULL, names(mu)))
  mu_row[, "Y"] <- mu_row[, "Y"] + J$gamma2 * z2      # Z2 shifts Y's mean only

  S11 <- J$V[target, target]
  S12 <- J$V[target, gn, drop = FALSE]
  S22 <- J$V[gn, gn, drop = FALSE]
  A   <- S12 %*% solve(S22)                            # 1 x k
  dev <- as.matrix(given) - mu_row[, gn, drop = FALSE]
  cm  <- as.vector(mu_row[, target] + dev %*% t(A))
  cv  <- as.numeric(S11 - A %*% t(S12))
  list(mean = cm, sd = sqrt(max(cv, 1e-12)))
}

#' Exact draw of the missing Z1 values, given the exposures, Z2 and Y.
ef_draw_z1 <- function(w, truth, idx, J) {
  if (!length(idx)) return(w$Z1)
  xc <- paste0("logX", seq_len(J$p))
  gv <- as.matrix(w[idx, c(xc, "Y"), drop = FALSE])
  cc <- .ef_cond(J, "Z1", gv, z2 = w$Z2[idx])
  out <- w$Z1
  out[idx] <- stats::rnorm(length(idx), cc$mean, cc$sd)
  out
}

#' Exact draw of the left-censored logX1 values: the same Gaussian conditional,
#' truncated above at each row's LOD.
ef_draw_x1 <- function(w, truth, idx, lod, J) {
  if (identical(truth$z_role, "pipe_nl"))
    stop("the exact logX1 draw is not available under z_role 'pipe_nl': its ",
         "conditional is non-Gaussian because g(logX1) enters Z1's prior. Hold ",
         "the X block at the shipped draw.", call. = FALSE)
  if (!length(idx)) return(w$logX1)
  xc <- paste0("logX", seq_len(J$p))
  gn <- c("Z1", setdiff(xc, "logX1"), "Y")
  gv <- as.matrix(w[idx, gn, drop = FALSE])
  cc <- .ef_cond(J, "logX1", gv, z2 = w$Z2[idx])
  out <- w$logX1
  out[idx] <- leftcens::rnorm_trunc(length(idx), cc$mean, cc$sd,
                                    lower = -Inf, upper = lod[idx])
  out
}

#' Exact draw of the missing Z1 under `z_role = "pipe_nl"`.
#'
#' The joint is NOT Gaussian -- Z1's mean is non-linear in logX1 -- but the
#' CONDITIONAL of Z1 still is, because the two factors it multiplies are both
#' Gaussian in Z1:
#'
#'   prior       Z1 | logX  ~  N( g(logX1), nl_sd^2 )
#'   likelihood  Y | Z1, .  ~  N( eta + gamma1 * Z1, sigma_y^2 )
#'
#' so the product is Gaussian with
#'
#'   precision = 1/nl_sd^2 + gamma1^2/sigma_y^2
#'   mean      = [ g(logX1)/nl_sd^2 + gamma1 * (Y - eta)/sigma_y^2 ] / precision
#'
#' where eta is the outcome mean with Z1 removed. `g()` is `ef_nl_g()` from
#' dgp.R -- the SAME function the generator uses, deliberately shared so the
#' exact draw cannot drift away from the truth it is supposed to be exact for.
ef_draw_z1_pipenl <- function(w, truth, idx, sigma_y = 1) {
  if (!length(idx)) return(w$Z1)
  p  <- length(truth$b)
  xc <- paste0("logX", seq_len(p))
  g1 <- truth$gamma[1]
  g2 <- if (length(truth$gamma) >= 2) truth$gamma[2] else 0
  sub <- w[idx, , drop = FALSE]
  eta <- truth$intercept +
    drop(as.matrix(sub[, xc, drop = FALSE]) %*% truth$b) + g2 * sub$Z2
  prior_mu <- ef_nl_g(sub[[xc[1]]], truth)
  prec <- 1 / truth$nl_sd^2 + g1^2 / sigma_y^2
  mu   <- (prior_mu / truth$nl_sd^2 + g1 * (sub$Y - eta) / sigma_y^2) / prec
  out <- w$Z1
  out[idx] <- stats::rnorm(length(idx), mu, sqrt(1 / prec))
  out
}

#' Exact draw of the missing binary Z2, given everything else.
#'
#' Z2 is independent of Z1 and the exposures in every role here, so its
#' conditional comes from the outcome model alone: a Bernoulli whose logit is the
#' log-likelihood ratio between Z2 = 1 and Z2 = 0, plus the prior logit.
ef_draw_z2 <- function(w, truth, idx, J, sigma_y = 1) {
  if (!length(idx)) return(w$Z2)
  xc <- paste0("logX", seq_len(J$p))
  mu0 <- truth$intercept +
    drop(as.matrix(w[idx, xc, drop = FALSE]) %*% truth$b) +
    truth$gamma[1] * w$Z1[idx]
  ll <- function(z2) stats::dnorm(w$Y[idx], mu0 + J$gamma2 * z2, sigma_y, log = TRUE)
  p1 <- stats::plogis(ll(1) - ll(0) + stats::qlogis(0.5))
  out <- w$Z2
  out[idx] <- stats::rbinom(length(idx), 1, p1)
  out
}

#' V20's instrument: block-FCS with EACH BLOCK SWITCHABLE between the shipped
#' draw and the exact conditional.
#'
#' Mirrors 00_censored_exposure.R's alternation (Z block, then X block, repeated)
#' so the only thing that varies is which draw each block uses:
#'
#'   z_src = "bart"    the shipped v1.5.0 Z block's draw (.v7_bart_draw)
#'   z_src = "exact"   the true conditional (ef_draw_z1 / ef_draw_z2)
#'   x_src = "shipped" leftcens::impute_censored_conditional, as the pipeline
#'   x_src = "exact"   the true truncated conditional (ef_draw_x1)
#'
#' The (bart, shipped) combination is the CONTROL: it must reproduce
#' `pipeline_bartMI`. If it does not, this instrument is not the pipeline and no
#' attribution from the other three cells means anything.
ef_impute_datasets <- function(d, truth, m = 30L, sweeps = 3L,
                               z_src = c("bart", "exact"),
                               x_src = c("shipped", "exact"),
                               margin = "shash", rho = 0.4, sd_x = 1,
                               mu_x = 0, sigma_y = 1) {
  z_src <- match.arg(z_src); x_src <- match.arg(x_src)
  J <- .ef_joint(truth, rho = rho, sd_x = sd_x, mu_x = mu_x, sigma_y = sigma_y)
  p <- J$p; xc <- paste0("logX", seq_len(p))
  cens <- attr(d, "censored_cols") %||% "logX1"
  if (!identical(as.character(cens), "logX1"))
    stop("V20's instrument covers focal-only censoring; got: ",
         paste(cens, collapse = ","), call. = FALSE)
  is_c <- is.na(d$logX1)
  lod  <- d$logX1_lod
  z1_na <- if (is.null(d$Z1)) logical(nrow(d)) else is.na(d$Z1)
  z2_na <- if (is.null(d$Z2)) logical(nrow(d)) else is.na(d$Z2)

  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d
    w$logX1[is_c] <- lod[is_c] - 0.5
    if (any(z1_na)) w$Z1[z1_na] <- mean(d$Z1, na.rm = TRUE)
    if (any(z2_na)) w$Z2[z2_na] <- round(mean(d$Z2, na.rm = TRUE))

    for (t in seq_len(sweeps)) {
      # --- Z block ------------------------------------------------------------
      if (z_src == "exact") {
        w$Z1 <- ef_draw_z1(w, truth, which(z1_na), J)
        w$Z2 <- ef_draw_z2(w, truth, which(z2_na), J, sigma_y = sigma_y)
      } else {
        for (zz in c("Z1", "Z2")) {
          na <- if (zz == "Z1") z1_na else z2_na
          if (!any(na)) next
          pr <- c(xc, setdiff(c("Z1", "Z2"), zz), "Y")
          # .v7_bart_draw returns a LIST of `ndraw` draws, not a vector.
          dr <- .v7_bart_draw(
            y_obs = w[[zz]][!na],
            x_obs = as.data.frame(w[!na, pr, drop = FALSE]),
            x_mis = as.data.frame(w[na,  pr, drop = FALSE]),
            binary = (zz == "Z2"))
          if (is.null(dr) || !length(dr))
            stop("BART Z draw failed for ", zz, call. = FALSE)
          w[[zz]][na] <- dr[[1L]]
        }
      }
      # --- X block ------------------------------------------------------------
      if (!any(is_c)) next
      if (x_src == "exact") {
        w$logX1 <- ef_draw_x1(w, truth, which(is_c), lod, J)
      } else {
        pr <- c("Y", setdiff(xc, "logX1"), "Z1", "Z2")
        yv <- w$logX1; yv[is_c] <- NA_real_
        lo <- rep(NA_real_, nrow(w)); hi <- rep(NA_real_, nrow(w))
        lo[is_c] <- -Inf; hi[is_c] <- lod[is_c]
        w$logX1 <- leftcens::impute_censored_conditional(
          y = yv, x = w[, pr, drop = FALSE], lower = lo, upper = hi,
          m = 1L, margin = margin)[, 1]
      }
    }
    out[[i]] <- w
  }
  out
}

#' One V20 arm: run the instrument and pool.
proc_ef_arm <- function(bundle, m = 30L, seed = NULL, sweeps = 3L,
                        z_src = "bart", x_src = "shipped", margin = "shash",
                        label = NULL, ...) {
  label <- label %||% paste0("ef_", z_src, "_", x_src)
  if (!is.null(seed)) set.seed(seed)
  truth <- bundle$truth
  imp <- tryCatch(ef_impute_datasets(bundle$censored, truth, m = as.integer(m),
                                     sweeps = sweeps, z_src = z_src,
                                     x_src = x_src, margin = margin),
                  error = function(e) structure(list(), err = conditionMessage(e)))
  if (!length(imp))
    return(one_row(label, NA, NA, NA, NA,
                   substr(paste("ef failed:", attr(imp, "err") %||% "none"), 1, 120)))
  es <- vs <- rep(NA_real_, length(imp))
  for (i in seq_along(imp)) {
    e <- fit_lm_estimand(as.data.frame(imp[[i]]), truth)
    es[i] <- e["est"]; vs[i] <- e["se"]^2
  }
  pl <- rubin_pool(es, vs)
  one_row(label, pl["est"], pl["se"], pl["ci_lo"], pl["ci_hi"],
          sprintf("z=%s; x=%s; m=%d; sweeps=%d", z_src, x_src, m, sweeps),
          ubar = pl["ubar"], b = pl["b"], fmi = pl["fmi"])
}

proc_ef_bart_ship  <- function(bundle, ...) proc_ef_arm(bundle, z_src = "bart",  x_src = "shipped", label = "ef_bart_ship",  ...)
proc_ef_exact_ship <- function(bundle, ...) proc_ef_arm(bundle, z_src = "exact", x_src = "shipped", label = "ef_exact_ship", ...)
proc_ef_bart_exact <- function(bundle, ...) proc_ef_arm(bundle, z_src = "bart",  x_src = "exact",   label = "ef_bart_exact", ...)
proc_ef_exact_exact<- function(bundle, ...) proc_ef_arm(bundle, z_src = "exact", x_src = "exact",   label = "ef_exact_exact", ...)

# =============================================================================
# V21 -- WHAT about a covariate draw has to be right?
# -----------------------------------------------------------------------------
# V20 located the bias in the covariate draw: replacing it with the true
# conditional removes 85-104%. V19 showed that is NOT simply "use a parametric
# imputer" -- `micePmm` and `properZ` are correctly specified IN FORM for this
# Z1 and still carry +4 to +6% bias that does not shrink with n.
#
# So something between "the exact conditional" (unbiased, +-0.5%) and "a
# correctly specified parametric draw" (+4 to +6%) is doing the damage. In THIS
# DGP the separable candidates are few, because the true conditional of Z1 is
# Gaussian with CONSTANT variance -- so V12's heteroscedastic-shape effect
# cannot arise here, and what remains is:
#
#   estimation error   fitting the conditional from data rather than knowing it
#   properness         propagating parameter uncertainty, or not
#   donor matching     pmm returns an observed neighbour, not a density draw
#   smoothing          BART's nonparametric fit of a linear truth
#   the inner-FCS loop present in shipped code, absent from V20's instrument
#
# EACH HAS A DIFFERENT n-SIGNATURE, which is what makes the ladder decisive
# rather than merely descriptive:
#
#   estimation error  -> decays as n^-1/2   (parametric rate, correct form)
#   smoothing         -> decays as n^-1/3   (V18/V19/V20 all measured this)
#   donor matching    -> does NOT decay     (a fixed discretisation)
#   properness        -> affects the INTERVAL, not the point estimate (V4)
#
# The X block is held at the shipped `leftcens` draw throughout, because V20
# measured its contribution at 0.1-11.7% -- removing that dimension makes every
# arm cheap and changes nothing being asked.
# =============================================================================

#' Fit the CORRECT linear conditional for Z1 and draw from it.
#'
#' @param proper TRUE draws coefficients and sigma from their posterior (proper
#'   MI); FALSE plugs in the MLE, which is what isolates properness.
#' @param pmm TRUE returns an observed donor whose predicted value is nearest
#'   (k = 5), which is what isolates donor matching. Overrides `proper`, since
#'   pmm's randomness comes from the donor pool rather than the parameters.
.ef_fit_draw_z1 <- function(w, idx, obs, proper = TRUE, pmm = FALSE, k = 5L) {
  xc <- grep("^logX[0-9]+$", names(w), value = TRUE)
  fo <- stats::as.formula(paste("Z1 ~", paste(c(xc, "Z2", "Y"), collapse = "+")))
  fit <- stats::lm(fo, data = w[obs, , drop = FALSE])
  mm  <- stats::model.matrix(fo, data = w)
  bh  <- stats::coef(fit); sg <- stats::summary.lm(fit)$sigma

  if (isTRUE(pmm)) {
    # Predictive mean matching on the FITTED means -- no parameter draw, no
    # density draw: the imputed value is always some other subject's observed Z1.
    yh_obs <- as.vector(mm[obs, , drop = FALSE] %*% bh)
    yh_mis <- as.vector(mm[idx, , drop = FALSE] %*% bh)
    zo <- w$Z1[obs]
    out <- w$Z1
    out[idx] <- vapply(yh_mis, function(v) {
      d <- abs(yh_obs - v)
      zo[sample(order(d)[seq_len(min(k, length(d)))], 1L)]
    }, 0)
    return(out)
  }

  if (isTRUE(proper)) {
    df <- fit$df.residual
    s2 <- sum(stats::residuals(fit)^2) / stats::rchisq(1, df)
    V  <- stats::summary.lm(fit)$cov.unscaled
    bh <- as.vector(bh + t(chol(V)) %*% stats::rnorm(length(bh)) * sqrt(s2))
    sg <- sqrt(s2)
  }
  mu <- as.vector(mm[idx, , drop = FALSE] %*% bh)
  out <- w$Z1
  out[idx] <- stats::rnorm(length(idx), mu, sg)
  out
}

#' V21's instrument: the Z-draw ladder, X block held at the shipped draw.
#' @param z2_same When TRUE, Z2 is drawn by the SAME method as Z1 instead of by
#'   its exact conditional. This is V22's probe of the V19 discrepancy:
#'   `micePmm` and `properZ` carried +2.4 to +6.3% ASYMPTOTIC bias while V21's
#'   parametric arms sat at +-0.6% -- and the registered explanation was that
#'   V21 draws Z1 alone while mice imputes Z1, Z2 and Y together. Holding Z2
#'   exact versus drawing it the same way isolates the first half of that.
ef21_impute_datasets <- function(d, truth, m = 30L, sweeps = 3L,
                                 z_src = c("exact", "fit_proper", "fit_improper",
                                           "pmm", "bart"),
                                 inner_iter = 1L, z2_same = FALSE,
                                 margin = "shash",
                                 rho = 0.4, sd_x = 1, mu_x = 0, sigma_y = 1) {
  z_src <- match.arg(z_src)
  nl <- identical(truth$z_role, "pipe_nl")
  # pipe_nl has no Gaussian joint, and needs none: the ladder holds the X block
  # at the shipped draw, so the joint was only ever used for the Z conditionals.
  # For pipe_nl, `J` is only a small parameter bundle -- ef_draw_z2() needs `p`
  # and `gamma2`, and nothing here needs the joint. Populating both rather than
  # leaving a stub: a NULL gamma2 silently collapses Z2's conditional mean to
  # numeric(0), which surfaces as "subscript out of bounds" three calls later.
  J <- if (nl) list(p = length(truth$b),
                    gamma2 = if (length(truth$gamma) >= 2) truth$gamma[2] else 0) else
    .ef_joint(truth, rho = rho, sd_x = sd_x, mu_x = mu_x, sigma_y = sigma_y)
  xc <- paste0("logX", seq_len(J$p))
  is_c <- is.na(d$logX1); lod <- d$logX1_lod
  z1_na <- is.na(d$Z1); z2_na <- is.na(d$Z2)
  obs1 <- which(!z1_na); idx1 <- which(z1_na)

  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d
    w$logX1[is_c] <- lod[is_c] - 0.5
    if (any(z1_na)) w$Z1[z1_na] <- mean(d$Z1, na.rm = TRUE)
    if (any(z2_na)) w$Z2[z2_na] <- round(mean(d$Z2, na.rm = TRUE))

    for (t in seq_len(sweeps)) {
      # The Z block, repeated `inner_iter` times -- the shipped engine's inner
      # FCS loop, which V20's instrument omitted (and ran 0.26-0.80 pp better).
      for (q in seq_len(inner_iter)) {
        if (length(idx1)) {
          w$Z1 <- switch(z_src,
            exact        = if (nl) ef_draw_z1_pipenl(w, truth, idx1, sigma_y = sigma_y)
                           else    ef_draw_z1(w, truth, idx1, J),
            fit_proper   = .ef_fit_draw_z1(w, idx1, obs1, proper = TRUE),
            fit_improper = .ef_fit_draw_z1(w, idx1, obs1, proper = FALSE),
            pmm          = .ef_fit_draw_z1(w, idx1, obs1, pmm = TRUE),
            bart         = {
              pr <- c(xc, "Z2", "Y")
              dr <- .v7_bart_draw(y_obs = w$Z1[obs1],
                                  x_obs = as.data.frame(w[obs1, pr, drop = FALSE]),
                                  x_mis = as.data.frame(w[idx1, pr, drop = FALSE]))
              if (is.null(dr) || !length(dr)) stop("BART Z1 draw failed", call. = FALSE)
              z <- w$Z1; z[idx1] <- dr[[1L]]; z
            })
        }
        # Z2: exact by default, so the ladder is about Z1 alone. With
        # z2_same = TRUE it is drawn by the same method as Z1 -- see the arg doc.
        if (any(z2_na)) {
          if (!isTRUE(z2_same) || z_src == "exact") {
            w$Z2 <- ef_draw_z2(w, truth, which(z2_na), J, sigma_y = sigma_y)
          } else if (z_src == "bart") {
            pr <- c(xc, "Z1", "Y")
            o2 <- which(!z2_na)
            dr <- .v7_bart_draw(y_obs = w$Z2[o2],
                                x_obs = as.data.frame(w[o2, pr, drop = FALSE]),
                                x_mis = as.data.frame(w[which(z2_na), pr, drop = FALSE]),
                                binary = TRUE)
            if (is.null(dr) || !length(dr)) stop("BART Z2 draw failed", call. = FALSE)
            w$Z2[z2_na] <- dr[[1L]]
          } else {
            # Parametric counterpart: a logistic fit on the correct predictor
            # set, drawn with parameter uncertainty unless the arm is improper.
            o2 <- which(!z2_na)
            fo <- stats::as.formula(paste("Z2 ~", paste(c(xc, "Z1", "Y"), collapse = "+")))
            fit <- suppressWarnings(stats::glm(fo, data = w[o2, , drop = FALSE],
                                               family = stats::binomial()))
            bh <- stats::coef(fit)
            if (z_src == "fit_proper" || z_src == "pmm") {
              V <- stats::vcov(fit)
              bh <- tryCatch(as.vector(bh + t(chol(V)) %*% stats::rnorm(length(bh))),
                             error = function(e) bh)
            }
            mm <- stats::model.matrix(fo, data = w)
            pr1 <- stats::plogis(as.vector(mm[which(z2_na), , drop = FALSE] %*% bh))
            w$Z2[z2_na] <- stats::rbinom(sum(z2_na), 1, pr1)
          }
        }
      }
      # --- X block: always the shipped leftcens draw --------------------------
      if (!any(is_c)) next
      pr <- c("Y", setdiff(xc, "logX1"), "Z1", "Z2")
      yv <- w$logX1; yv[is_c] <- NA_real_
      lo <- rep(NA_real_, nrow(w)); hi <- rep(NA_real_, nrow(w))
      lo[is_c] <- -Inf; hi[is_c] <- lod[is_c]
      w$logX1 <- leftcens::impute_censored_conditional(
        y = yv, x = w[, pr, drop = FALSE], lower = lo, upper = hi,
        m = 1L, margin = margin)[, 1]
    }
    out[[i]] <- w
  }
  out
}

proc_ef21_arm <- function(bundle, m = 30L, seed = NULL, sweeps = 3L,
                          z_src = "exact", inner_iter = 1L, z2_same = FALSE,
                          margin = "shash", label = NULL, ...) {
  label <- label %||% paste0("z21_", z_src, if (inner_iter > 1L) paste0("_inner", inner_iter),
                             if (isTRUE(z2_same)) "_z2same")
  if (!is.null(seed)) set.seed(seed)
  truth <- bundle$truth
  imp <- tryCatch(ef21_impute_datasets(bundle$censored, truth, m = as.integer(m),
                                       sweeps = sweeps, z_src = z_src,
                                       inner_iter = as.integer(inner_iter),
                                       z2_same = z2_same, margin = margin),
                  error = function(e) structure(list(), err = conditionMessage(e)))
  if (!length(imp))
    return(one_row(label, NA, NA, NA, NA,
                   substr(paste("z21 failed:", attr(imp, "err") %||% "none"), 1, 120)))
  es <- vs <- rep(NA_real_, length(imp))
  for (i in seq_along(imp)) {
    e <- fit_lm_estimand(as.data.frame(imp[[i]]), truth)
    es[i] <- e["est"]; vs[i] <- e["se"]^2
  }
  pl <- rubin_pool(es, vs)
  one_row(label, pl["est"], pl["se"], pl["ci_lo"], pl["ci_hi"],
          sprintf("z=%s; inner=%d; m=%d; sweeps=%d", z_src, inner_iter, m, sweeps),
          ubar = pl["ubar"], b = pl["b"], fmi = pl["fmi"])
}

proc_z21_exact        <- function(b, ...) proc_ef21_arm(b, z_src = "exact",        label = "z21_exact", ...)
proc_z21_fit_proper   <- function(b, ...) proc_ef21_arm(b, z_src = "fit_proper",   label = "z21_fit_proper", ...)
proc_z21_fit_improper <- function(b, ...) proc_ef21_arm(b, z_src = "fit_improper", label = "z21_fit_improper", ...)
proc_z21_pmm          <- function(b, ...) proc_ef21_arm(b, z_src = "pmm",          label = "z21_pmm", ...)
proc_z21_bart         <- function(b, ...) proc_ef21_arm(b, z_src = "bart",         label = "z21_bart", ...)
proc_z21_bart_inner3  <- function(b, ...) proc_ef21_arm(b, z_src = "bart", inner_iter = 3L, label = "z21_bart_inner3", ...)
# V22: the two Z2-same-method arms, for the V19 discrepancy.
proc_z21_fit_proper_z2same <- function(b, ...) proc_ef21_arm(b, z_src = "fit_proper", z2_same = TRUE, label = "z21_fit_proper_z2same", ...)
proc_z21_pmm_z2same        <- function(b, ...) proc_ef21_arm(b, z_src = "pmm",        z2_same = TRUE, label = "z21_pmm_z2same", ...)
