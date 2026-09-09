# Is the TOTAL effect still immune when the exposure draw is BADLY misspecified?
# pipe_nl is the cell where the ~20% asymptotic defect lives for the direct
# effect. If the total effect survives there too, the estimand choice is a
# genuine mitigation, not an artefact of a well-behaved cell.
`%||%` <- function(a,b) if (is.null(a)) b else a
for (f in c("dgp.R","censoring.R","metrics.R","procedures.R")) source(file.path("R",f))
suppressMessages({library(MASS); library(survival); library(parallel)})
n <- 2000L; NR <- 200L; M <- 20L

run <- function(role, target, nd = 0.40) {
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role)
  b0 <- tr$b[1]
  # total effect: for pipe_nl, Z1 = g(X1) + nu, so the induced X1 -> Y path is
  # gamma1 * dg/dx averaged over the population -- get it by simulation rather
  # than assuming linearity.
  set.seed(1); big <- simulate_complete(400000L, tr)$data
  tau <- coef(lm(Y ~ logX1 + logX2 + logX3 + Z2, data = big))[["logX1"]]
  tv <- if (target == "total") tau else b0
  fo <- if (target == "total") Y ~ logX1 + logX2 + logX3 + Z2
        else                   Y ~ logX1 + logX2 + logX3 + Z1 + Z2
  res <- do.call(rbind, mclapply(seq_len(NR), function(r) {
    set.seed(8800 + r)
    d  <- simulate_complete(n, tr)$data
    cn <- inject_left_censoring(d, nd_frac = nd, censor_which = 1L)
    isc <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
    yv <- d$logX1; yv[isc] <- NA_real_
    lo <- rep(NA_real_, n); hi <- rep(NA_real_, n); lo[isc] <- -Inf; hi[isc] <- lod
    xm <- leftcens::impute_censored_conditional(y = yv,
             x = d[, c("Y","logX2","logX3","Z1","Z2")], lower = lo, upper = hi,
             m = M, margin = "shash")
    es <- vs <- numeric(M)
    for (i in seq_len(M)) {
      dd <- d; dd$logX1 <- xm[, i]; f <- lm(fo, data = dd)
      es[i] <- coef(f)[["logX1"]]; vs[i] <- vcov(f)["logX1","logX1"]
    }
    p <- rubin_pool(es, vs)
    c(p[["est"]], p[["se"]], as.numeric(p[["ci_lo"]] <= tv && tv <= p[["ci_hi"]]))
  }, mc.cores = 20))
  cat(sprintf("  %-9s %-6s  truth %.3f  est %.4f  rel bias %+7.2f%%  bias/SE %5.2f  cov %.3f\n",
              role, target, tv, mean(res[,1]), (mean(res[,1])-tv)/tv*100,
              abs(mean(res[,1])-tv)/sd(res[,1]), mean(res[,3])))
}
cat(sprintf("\n  40%% censored, m = %d, %d reps, Z1 always IN the imputation\n\n", M, NR))
run("pipe",    "direct"); run("pipe",    "total")
run("pipe_nl", "direct"); run("pipe_nl", "total")
