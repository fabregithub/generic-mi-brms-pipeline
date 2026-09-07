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

  } else {
    stop("exact_fork.R supports z_role fork / pipe / precision; got ", role,
         call. = FALSE)
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
