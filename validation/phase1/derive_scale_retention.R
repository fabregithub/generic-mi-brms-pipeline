# DERIVE the scale defect, then check the derivation -- not the other way round.
#
# For a fork, logX1 is Gaussian and Z1 is jointly Gaussian with it, so the true
# conditional E[Z1 | logX1, ...] is LINEAR in L = logX1. The pipeline regressed
# on R = exp(L). For L ~ N(m, s^2):
#     Cov(L, R) = s^2 e^{m + s^2/2},   Var(R) = e^{2m + s^2}(e^{s^2} - 1)
#  => Corr(L, R)^2 = s^2 / (e^{s^2} - 1)
# which is the fraction of the explainable variance a linear predictor in R can
# retain. Everything else is lost into the residual, so the imputed Z1 is
# attenuated toward its marginal by exactly that factor.
#
# PREDICTIONS, before looking:
#   (1) retained R^2 ratio = s^2/(e^{s^2}-1)
#   (2) the Z1 residual SD inflates by 1/sqrt(1 - rho^2(1 - that ratio))
suppressMessages(for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f)))
`%||%` <- function(a,b) if (is.null(a)) b else a
for (cell in c("zr_fork","zr_pipe")) {
  sc <- Filter(function(s) s$name == cell, v2_scenarios())[[1]]
  tr <- make_truth(p = 3L, erf_form = sc$erf_form, z_role = sc$z_role)
  set.seed(21); d <- simulate_complete(400000L, tr)$data
  L <- d$logX1; s2 <- var(L)
  pred_ratio <- s2 / (exp(s2) - 1)
  # measured: variance of Z1 explained by R, relative to that explained by L,
  # holding the other predictors fixed in both models
  oth <- c("logX2", "logX3", "Z2", "Y")
  f_L <- lm(reformulate(c("L", oth), response = "Z1"), data = transform(d, L = L))
  f_R <- lm(reformulate(c("R", oth), response = "Z1"), data = transform(d, R = exp(L)))
  # partial R^2 of the exposure term in each
  f_0 <- lm(reformulate(oth, response = "Z1"), data = d)
  pr <- function(f) 1 - sum(residuals(f)^2) / sum(residuals(f_0)^2)
  meas_ratio <- pr(f_R) / pr(f_L)
  cat(sprintf("%-8s Var(logX1) = %.3f\n", cell, s2))
  cat(sprintf("   predicted retained fraction  s^2/(e^{s^2}-1) = %.4f\n", pred_ratio))
  cat(sprintf("   MEASURED  partial-R2(R) / partial-R2(L)      = %.4f   (miss %+.4f)\n",
              meas_ratio, meas_ratio - pred_ratio))
  cat(sprintf("   predicted residual-SD inflation = %.4f, measured = %.4f\n\n",
              sqrt((1 - pr(f_L) * pred_ratio) / (1 - pr(f_L))),
              sigma(f_R) / sigma(f_L)))
}
