# Is conditioning the exposure draw on Z LINEARLY worse than not conditioning on
# it at all, when the true X->Z arrow is non-linear? The grid comparison hinted
# so (+6.50% ignoring Z vs leftcens +22.75% using it linearly), but those differ
# in more than Z. Test it inside leftcens, changing only its predictor set.
`%||%` <- function(a,b) if (is.null(a)) b else a
for (f in c("dgp.R","censoring.R","metrics.R","procedures.R")) source(file.path("R",f))
suppressMessages({library(MASS); library(survival); library(parallel)})
n <- 2000L; NR <- 200L; M <- 20L

run <- function(role, z_in, target) {
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role)
  set.seed(1); big <- simulate_complete(400000L, tr)$data
  tv <- if (target=="total") coef(lm(Y~logX1+logX2+logX3+Z2, data=big))[["logX1"]] else tr$b[1]
  fo <- if (target=="total") Y~logX1+logX2+logX3+Z2 else Y~logX1+logX2+logX3+Z1+Z2
  res <- do.call(rbind, mclapply(seq_len(NR), function(r) {
    set.seed(9300 + r)
    d <- simulate_complete(n, tr)$data
    cn <- inject_left_censoring(d, nd_frac = 0.40, censor_which = 1L)
    isc <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
    pr <- c("Y","logX2","logX3","Z2", if (z_in) "Z1")
    yv <- d$logX1; yv[isc] <- NA_real_
    lo <- rep(NA_real_, n); hi <- rep(NA_real_, n); lo[isc] <- -Inf; hi[isc] <- lod
    xm <- leftcens::impute_censored_conditional(y = yv, x = d[, pr, drop=FALSE],
             lower = lo, upper = hi, m = M, margin = "shash")
    es <- vs <- numeric(M)
    for (i in seq_len(M)) {
      dd <- d; dd$logX1 <- xm[, i]; f <- lm(fo, data = dd)
      es[i] <- coef(f)[["logX1"]]; vs[i] <- vcov(f)["logX1","logX1"]
    }
    q <- rubin_pool(es, vs)
    c(q[["est"]], as.numeric(q[["ci_lo"]] <= tv && tv <= q[["ci_hi"]]))
  }, mc.cores = 20))
  cat(sprintf("  %-8s Z1 %-3s in draw  %-6s  rel %+7.2f%%  bias/SE %5.2f  cov %.3f\n",
              role, if (z_in) "IN" else "OUT", target,
              (mean(res[,1])-tv)/tv*100, abs(mean(res[,1])-tv)/sd(res[,1]), mean(res[,2])))
}
cat("\n  leftcens, 40% censored, m = 20, 200 reps -- only the predictor set changes\n\n")
run("pipe_nl", TRUE,  "direct"); run("pipe_nl", FALSE, "direct")
run("pipe_nl", TRUE,  "total");  run("pipe_nl", FALSE, "total")
run("fork",    TRUE,  "direct"); run("fork",    FALSE, "direct")
