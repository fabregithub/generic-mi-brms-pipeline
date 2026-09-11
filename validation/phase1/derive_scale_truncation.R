# LOD TRUNCATION AND lambda.
# Derived: the X block draws L | rest ~ N(m, sx^2) truncated above at c = LOD.
# With alpha = (c-m)/sx and h = phi/Phi,  E[L|L<=c] = m - sx*h(alpha),  and
# h' = -alpha*h - h^2, so
#     dE[L|L<=c]/dm = 1 - alpha*h - h^2 = Var(L|L<=c)/sx^2
# i.e. the response to a shift in the mean is the truncated variance RATIO,
# always < 1. So lambda_eff = lambda * mean(variance ratio over censored rows).
#
# Rather than argue the sign again, trace the predicted bias across lambda and
# see (a) whether it is monotone, (b) which direction truncation moves it, and
# (c) whether the lambda that would reproduce the measurement is reachable.
suppressMessages(for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f)))
`%||%` <- function(a,b) if (is.null(a)) b else a
oth <- c("logX2","logX3","Z2","Y"); N <- 300000L
meas <- c(zr_fork = 4.40, zr_pipe = 5.36)
for (cell in c("zr_fork","zr_pipe")) {
  sc <- Filter(function(s) s$name == cell, v2_scenarios())[[1]]
  tr <- make_truth(p = 3L, erf_form = sc$erf_form, z_role = sc$z_role)
  set.seed(23); d <- simulate_complete(N, tr)$data
  L <- d$logX1; s2 <- var(L)
  fL <- lm(reformulate(oth, response="L"), data=transform(d, L=L))
  sigt2 <- sum(residuals(fL)^2)/N; tau2 <- s2 - sigt2
  f_ret <- sigt2/(exp(s2) - 1 - tau2)
  fa <- lm(reformulate(c("L", oth), response="Z1"), data=transform(d, L=L))
  beta <- coef(fa)[["L"]]; a0 <- -beta*(1 - f_ret)
  fx <- lm(reformulate(c("Y","logX2","logX3","Z1","Z2"), response="L"),
           data = transform(d, L=L))
  lam0 <- coef(fx)[["Z1"]]; sx <- sigma(fx)
  lod <- as.numeric(quantile(L, sc$nd_frac)); x_cen <- L < lod
  # the derived attenuation factor, averaged over censored rows
  m  <- fitted(fx)[x_cen]
  al <- (lod - m)/sx
  h  <- dnorm(al)/pnorm(al)
  vr <- 1 - al*h - h^2                       # = Var(L|L<=c)/sx^2
  att <- mean(vr)
  lam_eff <- lam0 * att
  set.seed(99); z_mis <- runif(N) < sc$mcar_frac
  Lt <- residuals(fL)
  fit_t <- lm(Y ~ L + Z1 + logX2 + logX3 + Z2, data=transform(d, L=L))
  b1_t <- coef(fit_t)[["L"]]
  pred <- function(lam) {
    D  <- ifelse(z_mis, a0 * Lt * ifelse(x_cen, 1/(1 - a0*lam), 1), 0)
    Dx <- ifelse(x_cen & z_mis, lam * D, 0)
    dd <- transform(d, Lh = L + Dx, Zh = d$Z1 + D)
    100*(coef(lm(Y ~ Lh + Zh + logX2 + logX3 + Z2, data=dd))[["Lh"]] - b1_t)/tr$estimand_true
  }
  cat(sprintf("\n=== %s ===  lambda_untrunc %.4f | attenuation %.4f | lambda_eff %.4f\n",
              cell, lam0, att, lam_eff))
  cat(sprintf("%10s %12s\n", "lambda", "pred bias%"))
  for (lam in c(0, lam_eff, lam0, 2*lam0, 4*lam0))
    cat(sprintf("%10.4f %+12.2f%s\n", lam, pred(lam),
                if (isTRUE(all.equal(lam, lam_eff))) "   <- truncation-corrected"
                else if (isTRUE(all.equal(lam, lam0))) "   <- untruncated" else ""))
  cat(sprintf("   measured %+.2f  -> needs lambda %s\n", meas[[cell]],
              if (pred(0) > meas[[cell]]) "LARGER than any tried (monotone decreasing?)" else "in range"))
}
