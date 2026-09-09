# EXACT finite-sample decomposition of the imputation-induced bias.
#
# Truth:  Y = b0*X + g0'Z + eps.   Impute X~ for censored rows; X~ = X elsewhere.
# Let D = X~ - X (zero on observed rows). Then X = X~ - D, so with W = (X~, Z):
#
#   b_hat = (W'W)^-1 W'Y
#         = (b0, g0)' - b0 * (W'W)^-1 W'D + (W'W)^-1 W'eps
#
# so the bias in the X~ coefficient is EXACTLY the sum of two terms:
#
#   T1 = -b0 * [(W'W)^-1 W'D]_1     imputation error   -> classical attenuation
#   T2 = +     [(W'W)^-1 W'eps]_1   outcome borrowing  -> nonzero only if X~ | eps
#
# T2 vanishes when the imputation ignores Y. It does NOT vanish when the
# imputation conditions on Y and gets the conditional wrong -- and its sign is
# what decides whether the net bias attenuates or inflates.
`%||%` <- function(a,b) if (is.null(a)) b else a
for (f in c("dgp.R","censoring.R","metrics.R","procedures.R","robustness.R",
            "bart_impute.R","exact_fork.R")) source(file.path("R",f))
suppressMessages({library(MASS); library(survival)})

one <- function(role, nl_a = NULL, nl_c = NULL, with_y = TRUE, NR = 120L, n = 2000L) {
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role,
                   nl_a = nl_a, nl_c = nl_c)
  b0 <- tr$b[1]; out <- matrix(NA_real_, NR, 3)
  for (r in seq_len(NR)) {
    set.seed(5000 + r)
    d  <- simulate_complete(n, tr)$data
    cn <- inject_left_censoring(d, nd_frac = 0.40, censor_which = 1L)
    isc <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
    eps <- d$Y - (as.matrix(d[, paste0("logX",1:3)]) %*% tr$b +
                  as.matrix(d[, c("Z1","Z2")]) %*% tr$gamma + tr$intercept)
    # the SHIPPED exposure draw: leftcens, linear in its predictors
    pr <- c(if (with_y) "Y", "logX2","logX3","Z1","Z2")
    yv <- d$logX1; yv[isc] <- NA_real_
    lo <- rep(NA_real_, n); hi <- rep(NA_real_, n); lo[isc] <- -Inf; hi[isc] <- lod
    xt <- leftcens::impute_censored_conditional(y = yv, x = d[, pr, drop = FALSE],
             lower = lo, upper = hi, m = 1L, margin = "shash")[, 1]
    D  <- xt - d$logX1
    W  <- cbind(xt, as.matrix(d[, c("logX2","logX3","Z1","Z2")]), 1)
    XtXi <- solve(crossprod(W))
    T1 <- -b0 * (XtXi %*% crossprod(W, D))[1]
    T2 <-        (XtXi %*% crossprod(W, eps))[1]
    bh <- (XtXi %*% crossprod(W, d$Y))[1]
    out[r, ] <- c(bh - b0, T1, T2)
  }
  m <- colMeans(out)
  cat(sprintf("  %-26s bias %+7.4f = T1 %+7.4f (attenuate) + T2 %+7.4f (borrow)  | check %+.2e\n",
              sprintf("%s%s", role, if (!with_y) " [no Y]" else ""),
              m[1], m[2], m[3], m[1] - m[2] - m[3]))
  invisible(m)
}
cat("\n  b0 = 0.40, n = 2000, 40% left-censored, 120 reps, shipped leftcens draw\n\n")
one("precision")
one("precision", with_y = FALSE)
one("fork")
one("pipe")
one("pipe_nl", nl_a = 1.2, nl_c = 0.35)
one("pipe_nl", nl_a = 1.2, nl_c = 0.35, with_y = FALSE)
