#!/usr/bin/env Rscript
# =============================================================================
# V9 -- self-test for the SMC censored-exposure sampler (no BKMR, ~2 s)
# -----------------------------------------------------------------------------
# The whole V9 experiment rests on `.smc_draw_truncated()` actually sampling from
#
#     p(x) propto exp( -(Y - A - Bx - Cx^2)^2 / (2 sigma^2) - (x-m)^2 / (2 s^2) ),  x < upper
#
# If it does not, a "curvature restored" result would be an artifact of a broken
# sampler and a "not restored" result would be meaningless. These checks compare
# the sampler against that density evaluated by fine numerical integration --
# computed independently of the sampler's own grid code.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.here, "R", "dgp.R"))
source(file.path(.here, "R", "smc_impute.R"))

fails <- 0L
ok <- function(lbl, pass, detail="") {
  cat(sprintf("  [%s] %s%s\n", if (pass) "PASS" else "FAIL", lbl,
              if (nzchar(detail)) paste0("  --  ", detail) else ""))
  if (!pass) fails <<- fails + 1L
}

# Reference density by numerical integration on a very fine independent grid.
ref_moments <- function(Y, A, B, C, sig, m, s, upper, N = 200001L) {
  lo <- min(m - 12*s, upper - 12*s)
  x  <- seq(lo, upper, length.out = N)
  lp <- -(Y - A - B*x - C*x^2)^2/(2*sig^2) - (x-m)^2/(2*s^2)
  w  <- exp(lp - max(lp))
  w  <- w/sum(w)
  mu <- sum(w*x); v <- sum(w*(x-mu)^2)
  cdf <- cumsum(w)
  q <- sapply(c(.1,.25,.5,.75,.9), function(p) x[which.min(abs(cdf-p))])
  list(mean = mu, sd = sqrt(v), q = q)
}

cat("\n=== 1. Sampler vs numerical integration, across regimes ===\n")
cat("    (100k draws per case; compares mean, SD and five quantiles)\n\n")
set.seed(7)
cases <- list(
  list(lbl="strong curvature (C>0)",  Y= 0.8, A=0.1, B=0.40, C= 0.15, sig=1.0, m= 0.0, s=1.0, up= 0.2),
  list(lbl="negative curvature",      Y=-0.5, A=0.0, B=0.30, C=-0.25, sig=1.0, m=-0.2, s=0.9, up= 0.0),
  list(lbl="no curvature (C=0)",      Y= 0.3, A=0.0, B=0.50, C= 0.00, sig=1.0, m= 0.0, s=1.0, up=-0.2),
  list(lbl="tight outcome (sig=0.3)", Y= 0.6, A=0.1, B=0.40, C= 0.15, sig=0.3, m= 0.0, s=1.0, up= 0.5),
  list(lbl="hard truncation",         Y= 0.2, A=0.0, B=0.40, C= 0.15, sig=1.0, m= 0.0, s=1.0, up=-1.5)
)
for (cs in cases) {
  n <- 100000L
  dr <- .smc_draw_truncated(Y=rep(cs$Y,n), A=rep(cs$A,n), B=rep(cs$B,n), C=rep(cs$C,n),
                            sigma_y=cs$sig, m=rep(cs$m,n), s=cs$s, upper=rep(cs$up,n))
  rf <- ref_moments(cs$Y, cs$A, cs$B, cs$C, cs$sig, cs$m, cs$s, cs$up)
  dm <- abs(mean(dr)-rf$mean); ds <- abs(sd(dr)-rf$sd)
  dq <- max(abs(quantile(dr, c(.1,.25,.5,.75,.9), names=FALSE) - rf$q))
  # MC error of the mean at 100k draws is ~sd/316; require agreement well inside 0.01.
  ok(sprintf("%-24s", cs$lbl), dm < 0.01 && ds < 0.01 && dq < 0.02,
     sprintf("d_mean=%.4f d_sd=%.4f max_dq=%.4f", dm, ds, dq))
  ok(sprintf("%-24s respects truncation", cs$lbl), max(dr) <= cs$up + 1e-9,
     sprintf("max draw = %.4f, upper = %.4f", max(dr), cs$up))
}

cat("\n=== 2. Quadratic-coefficient recovery from the generator's own surface ===\n")
cat("    A, B, C are recovered by evaluating eta at x = 0, 1, -1. If the surface\n")
cat("    were not quadratic in that coordinate this would silently be wrong.\n\n")
truth <- make_truth(3L, 2L, "mixture")
set.seed(3)
sim <- simulate_complete(n = 5, truth = truth, rho = .4, sd_x = 1, mu_x = 0, sigma_y = 0)
d <- sim$data
for (j in 1:3) {
  cc <- paste0("logX", j)
  eta_fun <- function(xv) {
    s2 <- d; s2[[cc]] <- xv
    U <- as.matrix(s2[, paste0("logX",1:3)]); Zm <- as.matrix(s2[, c("Z1","Z2")])
    as.vector(U %*% truth$b) + as.vector(Zm %*% truth$gamma) + truth$intercept +
      truth$b_int*U[,1]*U[,2] + truth$b_quad*U[,1]^2
  }
  qc <- .smc_quadratic_coefs(eta_fun, nrow(d))
  # Check the recovered quadratic reproduces eta at fresh points it never saw.
  xt <- c(-2.3, 0.7, 1.9)
  err <- max(sapply(xt, function(v) max(abs(eta_fun(rep(v,nrow(d))) - (qc$A + qc$B*v + qc$C*v^2)))))
  expC <- if (j==1) truth$b_quad else 0
  ok(sprintf("logX%d: quadratic reproduces eta off-grid", j), err < 1e-10,
     sprintf("max err = %.3g", err))
  ok(sprintf("logX%d: C = %.2f as the surface implies", j, expC),
     max(abs(qc$C - expC)) < 1e-12)
}

cat("\n=== 3. Exposure prior matches the generator's MVN conditional ===\n\n")
set.seed(11)
big <- simulate_complete(n = 200000, truth = truth, rho = .4, sd_x = 1, mu_x = 0, sigma_y = 1)$data
U <- as.matrix(big[, paste0("logX",1:3)])
pr <- .smc_exposure_prior(1L, U[, -1, drop=FALSE], p = 3L, rho = .4, sd_x = 1, mu_x = 0)
resid <- U[,1] - pr$m
ok("conditional mean is unbiased", abs(mean(resid)) < 0.01, sprintf("mean resid = %.4f", mean(resid)))
ok("conditional SD matches empirical", abs(sd(resid) - pr$s) < 0.01,
   sprintf("empirical %.4f vs claimed %.4f", sd(resid), pr$s))

cat(sprintf("\n%s  (%d check%s failed)\n\n",
            if (fails==0L) "ALL CHECKS PASSED" else "FAILURES PRESENT",
            fails, if (fails==1L) "" else "s"))
quit(status = if (fails==0L) 0L else 1L)
