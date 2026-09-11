# THE CORRECT STOCHASTIC MODEL OF THE IMPUTATION.
# An imputation REPLACES the value: on missing rows Zhat1 = m_wrong + e_wrong,
# with residual variance sigma_wrong^2. It is NOT "the truth plus noise" -- that
# is the classical measurement-error structure, and writing Zhat1 = Z1 + D + noise
# keeps the true residual AND adds another, roughly doubling the imputed
# covariate's residual variance. That is why the earlier noise term read +28%.
#
# PREDICTION: building Zhat1 as a genuine replacement draw must land close to the
# mean-shift-only figure, because the misspecified model's residual SD is only
# 1.048x the correct one -- so there is very little excess spread to attenuate
# anything.
suppressMessages(for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f)))
`%||%` <- function(a,b) if (is.null(a)) b else a
oth <- c("logX2","logX3","Z2","Y"); N <- 300000L
meas <- c(zr_fork = 4.40, zr_pipe = 5.36)
for (cell in c("zr_fork","zr_pipe")) {
  sc <- Filter(function(s) s$name == cell, v2_scenarios())[[1]]
  tr <- make_truth(p=3L, erf_form=sc$erf_form, z_role=sc$z_role)
  set.seed(23); d <- simulate_complete(N, tr)$data
  L <- d$logX1; s2 <- var(L)
  fL <- lm(reformulate(oth, response="L"), data=transform(d, L=L))
  sigt2 <- sum(residuals(fL)^2)/N; tau2 <- s2 - sigt2
  f_ret <- sigt2/(exp(s2)-1-tau2)
  fa <- lm(reformulate(c("L", oth), response="Z1"), data=transform(d, L=L))
  fb <- lm(reformulate(c("R", oth), response="Z1"), data=transform(d, R=exp(L)))
  s_r <- sigma(fa); s_w <- sigma(fb)
  fx <- lm(reformulate(c("Y","logX2","logX3","Z1","Z2"), response="L"), data=transform(d, L=L))
  lam0 <- coef(fx)[["Z1"]]; sx <- sigma(fx)
  lod <- as.numeric(quantile(L, sc$nd_frac)); x_cen <- L < lod
  mm <- fitted(fx)[x_cen]; al <- (lod-mm)/sx; h <- dnorm(al)/pnorm(al)
  lam <- lam0 * mean(1 - al*h - h^2)
  set.seed(99); z_mis <- runif(N) < sc$mcar_frac
  b1_t <- coef(lm(Y ~ L + Z1 + logX2 + logX3 + Z2, data=transform(d, L=L)))[["L"]]
  # genuine replacement draw, with the feedback applied to the exposure
  # averaged over draws: a single realisation carries its own Monte-Carlo error
  set.seed(5)
  ps <- replicate(20L, {
    Zh <- d$Z1
    Zh[z_mis] <- fitted(fb)[z_mis] + rnorm(sum(z_mis), 0, s_w)
    D  <- Zh - d$Z1
    Dx <- ifelse(x_cen & z_mis, lam * D, 0)
    dd <- transform(d, Lh = L + Dx, Zh = Zh)
    100*(coef(lm(Y ~ Lh + Zh + logX2 + logX3 + Z2, data=dd))[["Lh"]] - b1_t)/tr$estimand_true
  })
  p <- mean(ps); pse <- sd(ps)/sqrt(length(ps))
  cat(sprintf("%-9s sigma_w/sigma_r %.3f | derived %+6.2f%% +/- %.2f | measured %+6.2f%% | ratio %.3f\n",
              cell, s_w/s_r, p, pse, meas[[cell]], p/meas[[cell]]))
}
