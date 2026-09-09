# TESTABLE CONSEQUENCE of the decomposition: T2 is OUTCOME BORROWING. So if the
# imputation's loading on Y is scaled by a factor k, T2 should scale roughly with
# k while T1 should not. k = 1 is the correctly specified draw.
`%||%` <- function(a,b) if (is.null(a)) b else a
for (f in c("dgp.R","censoring.R")) source(file.path("R",f))
suppressMessages({library(MASS); library(survival)})
tr <- make_truth(p = 3L, erf_form = "additive", z_role = "precision")
b0 <- tr$b[1]; n <- 3000L; NR <- 100L
cat("\n  precision cell, correctly specified draw, Y-loading scaled by k:\n\n")
cat(sprintf("  %5s %10s %10s %10s\n", "k", "bias", "T1", "T2"))
for (k in c(0, 0.5, 1, 1.5, 2)) {
  acc <- matrix(NA_real_, NR, 3)
  for (r in seq_len(NR)) {
    set.seed(7000 + r)
    d <- simulate_complete(n, tr)$data
    cn <- inject_left_censoring(d, nd_frac = 0.40, censor_which = 1L)
    isc <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
    eps <- d$Y - (as.matrix(d[, paste0("logX",1:3)]) %*% tr$b +
                  as.matrix(d[, c("Z1","Z2")]) %*% tr$gamma + tr$intercept)
    # fit the correct interval-censored conditional, then SCALE its Y coefficient
    dd <- d; dd$event <- as.integer(!isc); dd$resp <- ifelse(isc, lod, d$logX1)
    sr <- survreg(Surv(resp, event, type = "left") ~ Y + logX2 + logX3 + Z1 + Z2,
                  data = dd, dist = "gaussian")
    cf <- coef(sr); cf[["Y"]] <- k * cf[["Y"]]
    mm <- model.matrix(~ Y + logX2 + logX3 + Z1 + Z2, data = d)
    mu <- as.vector(mm %*% cf)
    xt <- d$logX1
    xt[isc] <- leftcens::rnorm_trunc(sum(isc), mu[isc], sr$scale, -Inf, lod)
    D <- xt - d$logX1
    W <- cbind(xt, as.matrix(d[, c("logX2","logX3","Z1","Z2")]), 1)
    Q <- solve(crossprod(W))
    acc[r, ] <- c((Q %*% crossprod(W, d$Y))[1] - b0,
                  -b0 * (Q %*% crossprod(W, D))[1],
                  (Q %*% crossprod(W, eps))[1])
  }
  m <- colMeans(acc)
  cat(sprintf("  %5.1f %+10.4f %+10.4f %+10.4f\n", k, m[1], m[2], m[3]))
}
