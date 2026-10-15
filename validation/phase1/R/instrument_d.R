# =============================================================================
# Recording the imputation error `d`, and the first-order bias it predicts.
# -----------------------------------------------------------------------------
# V30 showed a defect's ABSOLUTE bias is not family-invariant (0.26x-2.2x), and
# that no reduction of the transfer to a function of the information ratio
# survives (FINDINGS_v30 §4a, both rejected). The failing step is factoring
# E[w] out of E[w x d'] -- which cannot be done from the design, because it
# depends on the imputation error distribution `d` itself. So `d` has to be
# recorded.
#
# THE QUANTITY. The analysis solves an M-estimating equation; the imputation
# delivers x_hat = x + d. Writing g = dlogL/deta for the score in the linear
# predictor and theta = (beta, tau) for all estimated parameters (tau empty
# except for a cumulative-logit fit, where the cutpoints are estimated jointly
# and are NOT orthogonal to beta -- see derive_ordinal_information.R):
#
#     [I_bb I_bt; I_tb I_tt] [b_b; b_t]
#         = [E[d g]; 0] - E[J (d'beta0)] - [E[w d d'] beta0; 0]
#     J = sum_k (1/pi_k) g_k (dpi_k/deta)
#
# THE THIRD TERM IS SECOND ORDER IN d AND CANNOT BE DROPPED. A first-order
# expansion predicts the BORROWING channel well (self-test ratio 1.01-1.05 at
# small |d|, all three families) and the ATTENUATION channel not at all
# (0.003-0.15, sign-flipped under binomial) -- because for an error with
# E[d|x] = 0 the attenuation is exactly the term that first order kills. It is
# the classical measurement-error result: contaminating x_j with independent
# noise of variance s_d^2 biases beta_j by -beta_j s_d^2 / Var(x_j | rest),
# which is O(d^2). Keeping E[d g_eta d'beta0] = -E[w d d'] beta0 restores it.
#
# The beta block of J is x * w, recovering E[w x d'] beta0; the tau block is the
# cross term that the known-cutpoint version drops.
#
# WHY IT IS BUILT HERE AND NOT AS A STANDALONE SIMULATION. Three rebuilt
# instruments in this programme missed the pipeline's own level by 31%-4x
# (THEORY.md's caution on the s^2 sweep). The error has to come off the real
# arms.
# =============================================================================

#' Per-observation score in eta and the information blocks, for one family.
#'
#' @return list(g, J, I) -- g: n-vector score in eta; J: n x q cross matrix;
#'   I: q x q expected information (already averaged over rows). q = ncol(X)
#'   plus the number of cutpoints.
.d_family_parts <- function(X, eta, y, truth) {
  fam <- truth$y_family %||% "gaussian"
  n <- nrow(X); p <- ncol(X)
  if (identical(fam, "gaussian")) {
    g <- as.numeric(y - eta)                    # sigma^2 = 1 throughout
    w <- rep(1, n)
    return(list(g = g, J = X * w, I = crossprod(X * sqrt(w)) / n))
  }
  if (identical(fam, "binomial")) {
    mu <- stats::plogis(eta); w <- mu * (1 - mu)
    return(list(g = as.numeric(y - mu), J = X * w,
                I = crossprod(X * sqrt(w)) / n))
  }
  if (identical(fam, "ordinal")) {
    tau <- truth$y_tau; K <- length(tau) + 1L
    Fm <- vapply(tau, function(t) stats::plogis(t - eta), numeric(n))
    fm <- vapply(tau, function(t) stats::dlogis(t - eta), numeric(n))
    if (is.null(dim(Fm))) { Fm <- matrix(Fm, n); fm <- matrix(fm, n) }
    G  <- cbind(0, Fm, 1); fz <- cbind(0, fm, 0)
    pi <- G[, -1, drop = FALSE] - G[, -(K + 1), drop = FALSE]
    de <- -(fz[, -1, drop = FALSE] - fz[, -(K + 1), drop = FALSE])   # dpi/deta
    q  <- p + K - 1L
    I  <- matrix(0, q, q); J <- matrix(0, n, q)
    for (k in seq_len(K)) {
      gt <- matrix(0, n, K - 1L)
      if (k <= K - 1L) gt[, k]     <- gt[, k]     + fm[, k]
      if (k >= 2L)     gt[, k - 1] <- gt[, k - 1] - fm[, k - 1]
      gk <- cbind(de[, k] * X, gt)                 # dpi_k/dtheta
      iw <- 1 / pmax(pi[, k], 1e-300)
      I  <- I + crossprod(gk * sqrt(iw))
      J  <- J + gk * (iw * de[, k])
    }
    idx <- cbind(seq_len(n), as.integer(y))        # the realised category
    g <- de[idx] / pmax(pi[idx], 1e-300)           # dlogpi_Y/deta
    return(list(g = g, J = J, I = I / n))
  }
  stop("instrument_d: unsupported y_family ", fam, call. = FALSE)
}

#' First-order predicted bias in the focal coefficient, from a recorded `d`.
#'
#' @param comp True complete data (the DGP's own, before censoring/missingness).
#' @param di   One completed (imputed) dataset, on the ANALYSIS scale.
#' @param truth make_truth() output; `truth$y_tau` for an ordinal outcome.
#' @param focal Name of the focal covariate (default "logX1").
#' @return named numeric: pred_b (predicted bias in the focal coefficient, to
#'   second order in `d`),
#'   plus the two components -- borrow (the E[d g] term) and atten (the
#'   E[J d'beta0] term) -- so an arm's MIX of them is recorded, which is what
#'   V30 showed the transfer depends on. NA if the rows do not align.
d_predicted_bias <- function(comp, di, truth, focal = "logX1") {
  na <- c(pred_b = NA_real_, borrow = NA_real_, atten = NA_real_, d_rms = NA_real_)
  if (is.null(comp) || is.null(di) || nrow(comp) != nrow(di)) return(na)
  fam <- truth$y_family %||% "gaussian"
  rhs <- all.vars(dgp_formula(truth))[-1]
  if (!all(rhs %in% names(comp)) || !all(rhs %in% names(di))) return(na)

  Xt <- as.matrix(comp[, rhs, drop = FALSE])
  Xh <- as.matrix(di[,   rhs, drop = FALSE])
  if (anyNA(Xt) || anyNA(Xh)) return(na)
  # An ordinal fit has no intercept -- the cutpoints absorb it, and the columns
  # must be centred to match how derive_ordinal_information.R builds them.
  if (identical(fam, "ordinal")) {
    ctr <- colMeans(Xt); Xt <- sweep(Xt, 2, ctr); Xh <- sweep(Xh, 2, ctr)
  } else {
    Xt <- cbind(`(Intercept)` = 1, Xt); Xh <- cbind(`(Intercept)` = 1, Xh)
  }
  d <- Xh - Xt
  if (all(abs(d) < 1e-12)) return(c(pred_b = 0, borrow = 0, atten = 0, d_rms = 0))

  # theta0 on the design's own scale: the DGP's linear predictor.
  b0 <- rep(0, ncol(Xt)); names(b0) <- colnames(Xt)
  xn <- paste0("logX", seq_along(truth$b)); zn <- c("Z1", "Z2")[seq_along(truth$gamma)]
  b0[intersect(xn, names(b0))] <- truth$b[match(intersect(xn, names(b0)), xn)]
  b0[intersect(zn, names(b0))] <- truth$gamma[match(intersect(zn, names(b0)), zn)]
  if ("(Intercept)" %in% names(b0)) b0[["(Intercept)"]] <- truth$intercept %||% 0
  eta <- as.numeric(Xt %*% b0)
  if (identical(fam, "ordinal")) eta <- eta + (truth$intercept %||% 0)

  fp <- .d_family_parts(Xt, eta, comp[["Y"]], truth)
  q  <- ncol(fp$I); n <- nrow(Xt)
  dbeta <- as.numeric(d %*% b0)                      # d'beta0

  rhs_borrow <- rep(0, q); rhs_borrow[seq_len(ncol(Xt))] <- colMeans(d * fp$g)
  rhs_atten  <- -colMeans(fp$J * dbeta)
  # second-order attenuation: -E[w d d'] beta0, in the beta block only
  # J's beta block is X * w in every family, so this recovers w without
  # special-casing the intercept (absent under ordinal, where tau absorbs it).
  wd <- rowSums(fp$J[, seq_len(ncol(Xt)), drop = FALSE] * Xt) /
        pmax(rowSums(Xt * Xt), 1e-300)
  rhs_atten[seq_len(ncol(Xt))] <- rhs_atten[seq_len(ncol(Xt))] -
    as.numeric(crossprod(d * wd, d) %*% b0) / n
  Iinv <- tryCatch(solve(fp$I), error = function(e) NULL)
  if (is.null(Iinv)) return(na)
  j <- match(focal, colnames(Xt))
  bor <- as.numeric(Iinv[j, ] %*% rhs_borrow)
  att <- as.numeric(Iinv[j, ] %*% rhs_atten)
  c(pred_b = bor + att, borrow = bor, atten = att,
    d_rms = sqrt(mean(d[, j]^2)))
}
