# THE DAG-FACTORISED TARGET, tested as a DESIGN.
#
# Bayes for a censored X, factorised by the DAG rather than by regression habit:
#
#   p(X | rest, X<=L)  ∝  p(X | Pa(X))                      parents
#                       * PROD_{C in Ch(X)} p(C | X, Pa(C)\X)  CHILDREN  <-- missing today
#                       * 1{X <= L}
#
#   fork      Pa(X)={Z}, Ch(X)={Y}      ∝ p(X|Z) p(Y|X,Z)
#   pipe      Pa(X)={},  Ch(X)={Z,Y}    ∝ p(X)   p(Z|X) p(Y|X,Z)   <-- p(Z|X) omitted
#   collider  Pa(X)={},  Ch(X)={Z,Y}    ∝ p(X)   p(Z|X,Y) p(Y|X)
#
# `leftcens` regresses X on (Y, all covariates): that IS the fork form, and it is
# exact when everything is jointly Gaussian-linear. It has NO child factor, so it
# is wrong for a pipe whose X->Z arrow is non-linear -- precisely V22/V23.
#
# Tested here with g KNOWN, to separate "is the design right" from "can g be
# estimated" (the latter is the non-identifiability of THEORY 0b).
`%||%` <- function(a,b) if (is.null(a)) b else a
for (f in c("dgp.R","censoring.R","metrics.R","procedures.R")) source(file.path("R",f))
suppressMessages({library(MASS); library(survival); library(parallel)})

draw_dag <- function(d, isc, lod, tr, child = TRUE, ngrid = 400L) {
  n <- nrow(d); p <- 3L
  xc <- paste0("logX", seq_len(p))
  idx <- which(isc)
  # p(X | other exposures): exchangeable MVN, exact
  Sig <- matrix(0.4, p, p); diag(Sig) <- 1
  s11 <- Sig[1,1]; s12 <- Sig[1,-1,drop=FALSE]; S22 <- Sig[-1,-1]
  A <- s12 %*% solve(S22)
  mu_x <- as.vector(as.matrix(d[idx, xc[-1], drop=FALSE]) %*% t(A))
  sd_x <- sqrt(as.numeric(s11 - A %*% t(s12)))
  # grid on (-inf, L]
  G <- outer(rep(1, length(idx)), seq(-6, 0, length.out = ngrid)) * sd_x + lod
  eta0 <- tr$intercept +
    as.matrix(d[idx, xc[-1], drop=FALSE]) %*% tr$b[-1] +
    tr$gamma[1] * d$Z1[idx] + tr$gamma[2] * d$Z2[idx]
  lw <- dnorm(G, mu_x, sd_x, log = TRUE) +                      # p(X | Pa)
        dnorm(d$Y[idx], as.vector(eta0) + tr$b[1] * G, 1, log = TRUE)  # p(Y | X, Z)
  if (child)                                                     # p(Z | X)  <-- the fix
    lw <- lw + dnorm(d$Z1[idx], ef_nl_g(G, tr), tr$nl_sd, log = TRUE)
  lw <- lw - apply(lw, 1, max); W <- exp(lw)
  cdf <- t(apply(W, 1, cumsum)); cdf <- cdf / cdf[, ngrid]
  u <- runif(length(idx))
  out <- d$logX1
  out[idx] <- vapply(seq_along(idx), function(i)
    approx(cdf[i, ], G[i, ], xout = u[i], rule = 2, ties = "ordered")$y, 0)
  out
}

run <- function(child, target, NR = 200L, M = 20L, n = 2000L) {
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = "pipe_nl")
  set.seed(1); big <- simulate_complete(400000L, tr)$data
  tv <- if (target == "total") coef(lm(Y ~ logX1+logX2+logX3+Z2, data=big))[["logX1"]] else tr$b[1]
  fo <- if (target == "total") Y ~ logX1+logX2+logX3+Z2 else Y ~ logX1+logX2+logX3+Z1+Z2
  res <- do.call(rbind, mclapply(seq_len(NR), function(r) {
    set.seed(9100 + r)
    d <- simulate_complete(n, tr)$data
    cn <- inject_left_censoring(d, nd_frac = 0.40, censor_which = 1L)
    isc <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
    es <- vs <- numeric(M)
    for (i in seq_len(M)) {
      dd <- d; dd$logX1 <- draw_dag(d, isc, lod, tr, child = child)
      f <- lm(fo, data = dd); es[i] <- coef(f)[["logX1"]]; vs[i] <- vcov(f)["logX1","logX1"]
    }
    q <- rubin_pool(es, vs)
    c(q[["est"]], as.numeric(q[["ci_lo"]] <= tv && tv <= q[["ci_hi"]]))
  }, mc.cores = 20))
  cat(sprintf("  child factor %-3s  %-6s  truth %.3f  est %.4f  rel %+7.2f%%  bias/SE %5.2f  cov %.3f\n",
              if (child) "ON" else "OFF", target, tv, mean(res[,1]),
              (mean(res[,1])-tv)/tv*100, abs(mean(res[,1])-tv)/sd(res[,1]), mean(res[,2])))
}
cat("\n  pipe_nl, 40% censored, m = 20, 200 reps. g KNOWN (design test, not a diagnostic)\n")
cat("  For reference: shipped leftcens gave direct +22.75% (cov 0.600), total +17.85% (cov 0.100)\n\n")
run(FALSE, "direct"); run(TRUE, "direct")
run(FALSE, "total");  run(TRUE, "total")
