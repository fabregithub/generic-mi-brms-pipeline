# THE FORK RESIDUAL. Derived +5.05% vs measured +4.40% (1.147); the pipe closed
# at 1.011. The fork's distinguishing feature is Var(logX) = 1.36 vs 1.00, so
# e^{s^2} = 3.90 vs 2.72 -- the exp-nonlinearity my chain LINEARISES is much
# stronger there. Link 1's miss was already +0.004 in the same direction.
#
# REGISTERED PREDICTION: the ratio derived/measured is a LINEARISATION error, so
# it must approach 1 as s^2 -> 0 and grow with s^2. If it is flat in s^2 the
# hypothesis is wrong and something else is fork-specific.
#
# Measured here by a FAITHFUL two-block simulation -- misspecified Z model on
# exp(L), leftcens-style truncated exposure draw, sweeps, m imputations, pooled
# -- so the comparison isolates the closed form from the pipeline's other
# machinery.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f))})
`%||%` <- function(a,b) if (is.null(a)) b else a
oth <- c("logX2","logX3","Z2","Y")

closed_form <- function(d, tr, mcar, ndf) {
  L <- d$logX1; N <- nrow(d); s2 <- var(L)
  fL <- lm(reformulate(oth, response="L"), data=transform(d, L=L))
  sigt2 <- sum(residuals(fL)^2)/N; tau2 <- s2 - sigt2
  f <- sigt2/(exp(s2) - 1 - tau2)
  fa <- lm(reformulate(c("L", oth), response="Z1"), data=transform(d, L=L))
  fb <- lm(reformulate(c("R", oth), response="Z1"), data=transform(d, R=exp(L)))
  a0 <- -coef(fa)[["L"]]*(1-f); s_w <- sigma(fb)
  fx <- lm(reformulate(c("Y","logX2","logX3","Z1","Z2"), response="L"), data=transform(d, L=L))
  lam0 <- coef(fx)[["Z1"]]; sx <- sigma(fx)
  lod <- as.numeric(quantile(L, ndf)); xc <- L < lod
  mm <- fitted(fx)[xc]; al <- (lod-mm)/sx; h <- dnorm(al)/pnorm(al)
  lam <- lam0*mean(1 - al*h - h^2)
  set.seed(99); zm <- runif(N) < mcar
  b1t <- coef(lm(Y ~ L + Z1 + logX2 + logX3 + Z2, data=transform(d, L=L)))[["L"]]
  set.seed(5)
  ps <- replicate(10L, {
    Zh <- d$Z1; Zh[zm] <- fitted(fb)[zm] + rnorm(sum(zm), 0, s_w)
    Dx <- ifelse(xc & zm, lam*(Zh - d$Z1), 0)
    coef(lm(Y ~ Lh + Zh + logX2 + logX3 + Z2,
            data=transform(d, Lh = L + Dx, Zh = Zh)))[["Lh"]] })
  list(pred = 100*(mean(ps)-b1t)/tr$estimand_true, s2 = s2, f = f, lam = lam)
}

# faithful two-block FCS with the WRONG scale in the Z block
faithful <- function(d, tr, mcar, ndf, m = 20L, sweeps = 3L) {
  N <- nrow(d); L <- d$logX1
  lod <- as.numeric(quantile(L, ndf)); xc <- L < lod
  set.seed(99); zm <- runif(N) < mcar
  est <- numeric(m)
  for (i in seq_len(m)) {
    w <- d; w$logX1[xc] <- lod - 0.5; w$Z1[zm] <- mean(d$Z1[!zm])
    for (t in seq_len(sweeps)) {
      # Z block: regress on exp(logX1) -- the defect -- fitted on observed rows
      ww <- transform(w, R = exp(w$logX1))
      fz <- lm(reformulate(c("R", oth), response="Z1"), data=ww[!zm, , drop=FALSE])
      w$Z1[zm] <- predict(fz, newdata=ww[zm, , drop=FALSE]) + rnorm(sum(zm), 0, sigma(fz))
      # X block: INTERVAL-CENSORED fit on ALL rows, then a truncated draw -- what
      # leftcens does. An lm() on the above-LOD rows is a response-truncated
      # regression (sigma ~24% low) and injects a SECOND defect on top of the
      # scale defect under test; the first version of this script did exactly
      # that and read +17% where the pipeline reads +4.4%.
      lo <- w$logX1; hi <- w$logX1; lo[xc] <- -Inf; hi[xc] <- lod
      fxx <- survreg(Surv(lo, hi, type="interval2") ~ Y + logX2 + logX3 + Z1 + Z2,
                     data = cbind(w, lo = lo, hi = hi), dist = "gaussian")
      mu <- predict(fxx, newdata=w[xc, , drop=FALSE], type="response"); sg <- fxx$scale
      u <- runif(sum(xc))*pnorm(lod, mu, sg)
      w$logX1[xc] <- qnorm(pmax(u, 1e-300), mu, sg)
    }
    est[i] <- coef(lm(Y ~ logX1 + Z1 + logX2 + logX3 + Z2, data=w))[["logX1"]]
  }
  b1t <- coef(lm(Y ~ logX1 + Z1 + logX2 + logX3 + Z2, data=d))[["logX1"]]
  100*(mean(est)-b1t)/tr$estimand_true
}

cat(sprintf("%-6s %7s %7s %7s | %9s %9s %8s\n",
            "sd_x","s^2","f","lambda","derived","faithful","ratio"))
for (sdx in c(0.4, 0.7, 1.0, 1.3)) {
  tr <- make_truth(p=3L, erf_form="additive", z_role="fork")
  set.seed(23); d <- simulate_complete(120000L, tr, sd_x = sdx)$data
  cf <- closed_form(d, tr, 0.4, 0.4)
  fa <- faithful(d, tr, 0.4, 0.4)
  cat(sprintf("%-6.1f %7.3f %7.4f %7.4f | %+9.2f %+9.2f %8.3f\n",
              sdx, cf$s2, cf$f, cf$lam, cf$pred, fa, cf$pred/fa))
}
