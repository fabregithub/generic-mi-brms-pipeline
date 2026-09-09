# DO WE LOSE ANYTHING BY COMPARING POINT ESTIMATES + COVERAGE INSTEAD OF WHOLE
# POSTERIORS? For a Gaussian linear model with a flat prior the posterior for a
# coefficient is closed form -- t_{n-p}(beta_hat, s^2 (X'X)^-1) -- so the whole
# pooled posterior can be built exactly, with no MCMC, and its SHAPE inspected.
# Shape is the thing bias / emp_se / coverage cannot see.
`%||%` <- function(a,b) if (is.null(a)) b else a
for (f in c("dgp.R","censoring.R","metrics.R","procedures.R")) source(file.path("R",f))
suppressMessages({library(MASS); library(survival); library(parallel)})
n <- 2000L; NR <- 60L; M <- 20L; NDRAW <- 400L

post_draws <- function(d, fo) {                    # exact conjugate posterior draws
  f <- lm(fo, data = d); df <- f$df.residual
  s2 <- sum(residuals(f)^2) / rchisq(NDRAW, df)
  se <- sqrt(vcov(f)["logX1","logX1"] / summary(f)$sigma^2)
  coef(f)[["logX1"]] + rnorm(NDRAW, 0, se * sqrt(s2))
}
run <- function(role) {
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role)
  fo <- Y ~ logX1 + logX2 + logX3 + Z1 + Z2
  out <- do.call(rbind, mclapply(seq_len(NR), function(r) {
    set.seed(9500 + r)
    d  <- simulate_complete(n, tr)$data
    cn <- inject_left_censoring(d, nd_frac = 0.40, censor_which = 1L)
    isc <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
    yv <- d$logX1; yv[isc] <- NA_real_
    lo <- rep(NA_real_, n); hi <- rep(NA_real_, n); lo[isc] <- -Inf; hi[isc] <- lod
    xm <- leftcens::impute_censored_conditional(y = yv,
            x = d[, c("Y","logX2","logX3","Z1","Z2")], lower = lo, upper = hi,
            m = M, margin = "shash")
    pooled <- unlist(lapply(seq_len(M), function(i) {
      dd <- d; dd$logX1 <- xm[, i]; post_draws(dd, fo) }))
    full <- post_draws(d, fo)                       # complete-data posterior
    zc <- function(v) (v - mean(v))/sd(v)
    c(mean(pooled) - mean(full), sd(pooled)/sd(full),
      mean(zc(pooled)^3), mean(zc(pooled)^4), mean(zc(full)^3), mean(zc(full)^4))
  }, mc.cores = 20))
  m <- colMeans(out)
  cat(sprintf("  %-9s location shift %+7.4f | scale ratio %5.2f | pooled skew %+6.3f kurt %5.2f | complete-data skew %+6.3f kurt %5.2f\n",
              role, m[1], m[2], m[3], m[4], m[5], m[6]))
}
cat(sprintf("\n  Exact conjugate posteriors (no MCMC), m = %d, %d reps, %d draws each\n\n", M, NR, NDRAW))
run("fork"); run("pipe"); run("pipe_nl")
