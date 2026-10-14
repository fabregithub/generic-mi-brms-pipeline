# =============================================================================
# V30 -- the attenuation factor for an ORDINAL outcome, derived.
# -----------------------------------------------------------------------------
# V26 measured defects reading ~0.474x as large in `bias/SE` under a binary
# outcome, and THEORY.md then derived that factor as FISHER INFORMATION rather
# than anything about imputation: for a GLM the information in the eta
# direction is w = (dmu/deta)^2 / Var(Y), so
#
#     SE_gauss / SE_fam = sigma * sqrt(E[w]),   w = p(1-p) for a binomial.
#
# This script generalises that to a cumulative-logit (proportional-odds)
# outcome and -- importantly -- corrects it.
#
# STEP 1. The cumulative-logit information in the eta direction.
#   P(Y<=k|eta) = F(tau_k - eta),  pi_k = F(tau_k-eta) - F(tau_{k-1}-eta)
#   dpi_k/deta  = -[f(tau_k-eta) - f(tau_{k-1}-eta)],   f = F(1-F)
#   w_ord(eta)  = sum_k (dpi_k/deta)^2 / pi_k
#
# Two checks it has to pass, and does:
#   * K = 2 reduces to p(1-p) to MACHINE PRECISION (diff 0.00e+00).
#   * E[w_ord] -> 1/3 as K -> inf, monotonically (0.33328 at K=129) -- the
#     location information of a standard logistic, i.e. the latent observed
#     exactly. So an ordinal outcome is strictly more informative than a binary
#     one and strictly less than the latent, which is what it must be.
#
# STEP 2. sqrt(E[w]) IS NOT THE RATIO, for two separate reasons, and both were
# caught by measuring the complete-data SE ratio before running anything.
#
#   (a) The information for beta is E[w * xtilde^2], not E[w] * E[xtilde^2]:
#       w peaks mid-range and eta increases in X1, so Cov(w, xtilde^2) < 0.
#       Using E[w] over-predicts. Measured miss: -0.0068 .. -0.0112, one-signed.
#
#   (b) THE CUTPOINTS ARE ESTIMATED. The derivation above treats tau as known;
#       polr estimates K-1 of them jointly with beta, and -- unlike the
#       Gaussian twin's nuisance parameters -- they are NOT orthogonal to beta.
#       So Var(beta-hat) is the (1,1) element of the INVERSE of the full
#       (beta, tau) information, not 1/I_bb. This is an asymptotic effect, not
#       a finite-sample one: measured FLAT at +0.0073 / +0.0061 / +0.0071 for
#       n = 2000 / 8000 / 32000, i.e. 14 sigma at the largest and no decay.
#
# With both corrected the prediction is exact -- see the validation at the
# bottom: six cells, misses -0.0000 .. +0.0005 with mixed signs, no free
# parameters.
#
# NOTE FOR V26. Its E[p(1-p)] form is the known-threshold approximation, and a
# logistic regression estimates its intercept too. V26's misses (+0.013,
# +0.016, -0.004) are consistent with the same (a)+(b) effect; the exact form
# gives bin_fork 0.4459 and bin_pipe 0.4488 against V26's measured 0.4731 and
# 0.4527. Worth re-reading THEORY.md's logistic table against this.
#
#   Rscript derive_ordinal_information.R
# =============================================================================
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f))})

#' Per-observation Fisher information in the eta direction (cumulative logit).
w_ord <- function(eta, tau) {
  n <- length(eta); K <- length(tau) + 1L
  cum <- matrix(vapply(tau, function(t) plogis(t - eta), numeric(n)), nrow = n)
  den <- matrix(vapply(tau, function(t) dlogis(t - eta), numeric(n)), nrow = n)
  G <- cbind(0, cum, 1); f <- cbind(0, den, 0)
  pi <- G[, -1, drop = FALSE] - G[, -(K + 1), drop = FALSE]
  d  <- -(f[, -1, drop = FALSE] - f[, -(K + 1), drop = FALSE])
  rowSums(d^2 / pmax(pi, 1e-300))
}

#' EXACT asymptotic SE ratio, Gaussian twin vs proportional-odds fit.
#' Builds the full (beta, tau) expected information and inverts it, so both
#' corrections (a) and (b) above are included by construction.
ord_se_ratio <- function(X, eta, tau, sigma = 1) {
  n <- nrow(X); p <- ncol(X); K <- length(tau) + 1L
  F <- sapply(tau, function(t) plogis(t - eta))
  f <- sapply(tau, function(t) dlogis(t - eta))
  if (is.null(dim(F))) { F <- matrix(F, n); f <- matrix(f, n) }
  G  <- cbind(0, F, 1); fz <- cbind(0, f, 0)
  pi <- G[, -1, drop = FALSE] - G[, -(K + 1), drop = FALSE]
  de <- -(fz[, -1, drop = FALSE] - fz[, -(K + 1), drop = FALSE])
  I  <- matrix(0, p + K - 1L, p + K - 1L)
  for (k in seq_len(K)) {
    w  <- 1 / pmax(pi[, k], 1e-300)
    gt <- matrix(0, n, K - 1L)
    if (k <= K - 1L) gt[, k]     <- gt[, k]     + f[, k]
    if (k >= 2L)     gt[, k - 1] <- gt[, k - 1] - f[, k - 1]
    I <- I + crossprod(cbind(de[, k] * X, gt) * sqrt(w))
  }
  I <- I / n
  c(ratio  = (sigma * sqrt(solve(crossprod(X) / n)[1, 1] / n)) /
             sqrt(solve(I)[1, 1] / n))
}

.cells <- list(
  bin_fork  = list("fork", 0),
  bin_pipe  = list("pipe", 0),
  ord3_fork = list("fork", c(-0.9732, 0.6698)),
  ord5_fork = list("fork", c(-1.7753, -0.6340, 0.3302, 1.4724)),
  ord7_fork = list("fork", c(-2.2329, -1.2354, -0.4944, 0.1897, 0.9312, 1.9282)),
  ord3_pipe = list("pipe", c(-0.9740, 0.6719)),
  ord5_pipe = list("pipe", c(-1.7773, -0.6354, 0.3303, 1.4776)),
  ord7_pipe = list("pipe", c(-2.2346, -1.2351, -0.4950, 0.1897, 0.9347, 1.9334)))

#' The registered predictions, recomputed from scratch.
v30_predictions <- function(N = 300000L) {
  do.call(rbind, lapply(names(.cells), function(nm) {
    role <- .cells[[nm]][[1]]; tau <- .cells[[nm]][[2]]
    set.seed(101)
    tr <- make_truth(p = 3L, erf_form = "additive", z_role = role, target = "direct")
    d  <- simulate_complete(N, tr, rho = 0.4, skew = 0.0, sigma_y = 0)$data
    rhs <- setdiff(all.vars(dgp_formula(tr))[-1], "logX1")
    X <- scale(cbind(logX1 = d$logX1, as.matrix(d[, rhs, drop = FALSE])), TRUE, FALSE)
    data.frame(cell = nm, K = length(tau) + 1L,
               predicted = unname(ord_se_ratio(X, d$Y, tau)))
  }))
}

if (sys.nframe() == 0L) {
  set.seed(1); e <- rnorm(400000)
  cat("CHECK 1  K=2 reduces to p(1-p):\n")
  for (t1 in c(-0.8, 0, 0.3)) {
    p <- plogis(t1 - e)
    cat(sprintf("   tau=%+.1f  diff = %.2e\n", t1, mean(w_ord(e, t1)) - mean(p * (1 - p))))
  }
  cat("\nCHECK 2  K -> inf approaches the logistic location information 1/3:\n")
  set.seed(2); lat <- e + rlogis(length(e))
  for (K in c(2, 3, 5, 9, 33, 129)) {
    tau <- as.numeric(quantile(lat, seq_len(K - 1) / K))
    cat(sprintf("   K=%4d  E[w_ord] = %.5f\n", K, mean(w_ord(e, tau))))
  }
  cat(sprintf("   limit          = %.5f\n", 1 / 3))
  cat("\nREGISTERED PREDICTIONS\n")
  print(format(v30_predictions(), digits = 4), row.names = FALSE)
}
