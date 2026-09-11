# DERIVATION, then test. Retained fraction of Z1's explainable variance when the
# exposure enters the covariate block as R = e^L instead of L = logX1, GIVEN the
# other predictors W = (logX2, logX3, Z2, Y).
#
#   L = mu_W + Ltilde,  tau^2 = Var(mu_W),  sigt^2 = Var(Ltilde),  s^2 = tau^2 + sigt^2
#   R = e^L = A*B,  A = e^{mu_W},  B = e^{Ltilde},  A independent of B
#
#   Cov(Ltilde, Rtilde) = E[A]E[B] sigt^2            (Stein on Cov(X, e^X))
#   Var(Rtilde)         = E[A]^2 E[B]^2 (e^{s^2} - 1 - tau^2)
#   => retained = Corr(L, R | W)^2 = sigt^2 / (e^{s^2} - 1 - tau^2)
#
# Nests the marginal form: tau^2 = 0 gives s^2/(e^{s^2} - 1).
#
# TWO REGISTERED PREDICTIONS, stated before looking:
#   P1  retained = sigt^2 / (e^{s^2} - 1 - tau^2), against the measured partial-R2 ratio
#   P2  the marginal form s^2/(e^{s^2}-1) must OVER-predict, since tau^2 > 0 shrinks
#       the numerator and the -tau^2 in the denominator does not compensate
suppressMessages(for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f)))
`%||%` <- function(a,b) if (is.null(a)) b else a
oth <- c("logX2","logX3","Z2","Y")
cat(sprintf("%-9s %7s %7s %7s | %9s %9s | %9s %8s\n",
            "cell","s^2","tau^2","sigt^2","P1 deriv","measured","miss","marginal"))
for (cell in c("zr_fork","zr_pipe","zr_collider")) {
  sc <- Filter(function(s) s$name == cell, v2_scenarios())
  if (!length(sc)) next
  sc <- sc[[1]]
  tr <- make_truth(p = 3L, erf_form = sc$erf_form, z_role = sc$z_role)
  set.seed(21); d <- simulate_complete(400000L, tr)$data
  L <- d$logX1; s2 <- var(L)
  # tau^2 / sigt^2 from the projection of L on W
  fL <- lm(reformulate(oth, response = "L"), data = transform(d, L = L))
  sigt2 <- sum(residuals(fL)^2)/nrow(d); tau2 <- s2 - sigt2
  pred  <- sigt2 / (exp(s2) - 1 - tau2)
  marg  <- s2 / (exp(s2) - 1)
  # measured: partial R^2 of the exposure term, R vs L, other predictors held in
  f0 <- lm(reformulate(oth, response = "Z1"), data = d)
  fa <- lm(reformulate(c("L", oth), response = "Z1"), data = transform(d, L = L))
  fb <- lm(reformulate(c("R", oth), response = "Z1"), data = transform(d, R = exp(L)))
  pr <- function(f) 1 - sum(residuals(f)^2)/sum(residuals(f0)^2)
  meas <- pr(fb)/pr(fa)
  cat(sprintf("%-9s %7.3f %7.3f %7.3f | %9.4f %9.4f | %+9.4f %8.4f\n",
              cell, s2, tau2, sigt2, pred, meas, meas - pred, marg))
}
