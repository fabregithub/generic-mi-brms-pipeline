# LINK 1: L-loading of the Z-block error D = -beta*(1-f), where beta is the true
#         partial coefficient of L in E[Z1 | L, W] and f the retained fraction.
# LINK 2: the X block re-draws L given the corrupted Z1, so an error D in Z1
#         produces Delta_x = lambda*D on censored rows, with lambda the Z1
#         coefficient in the X-block conditional E[L | Y, X_other, Z1, Z2].
# Each link is checked on its own -- a final number that happens to match is not
# evidence that the chain is right.
suppressMessages(for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f)))
`%||%` <- function(a,b) if (is.null(a)) b else a
oth <- c("logX2","logX3","Z2","Y")
cat(sprintf("%-9s %8s %8s | %10s %10s %8s | %8s\n",
            "cell","beta","f","L-load pred","L-load meas","miss","lambda"))
for (cell in c("zr_fork","zr_pipe","zr_collider")) {
  sc <- Filter(function(s) s$name == cell, v2_scenarios())
  if (!length(sc)) next
  sc <- sc[[1]]
  tr <- make_truth(p = 3L, erf_form = sc$erf_form, z_role = sc$z_role)
  set.seed(23); d <- simulate_complete(300000L, tr)$data
  L <- d$logX1; s2 <- var(L)
  fL <- lm(reformulate(oth, response="L"), data=transform(d, L=L))
  sigt2 <- sum(residuals(fL)^2)/nrow(d); tau2 <- s2 - sigt2
  f_ret <- sigt2/(exp(s2) - 1 - tau2)
  # beta: true partial coefficient of L in E[Z1 | L, W]
  fa <- lm(reformulate(c("L", oth), response="Z1"), data=transform(d, L=L))
  beta <- coef(fa)[["L"]]
  # measured L-loading of D, partialled on W (the derivation is a partial one)
  fb <- lm(reformulate(c("R", oth), response="Z1"), data=transform(d, R=exp(L)))
  D  <- fitted(fb) - fitted(fa)
  Lt <- residuals(fL)                      # L residualised on W
  meas <- coef(lm(D ~ Lt))[["Lt"]]
  pred <- -beta * (1 - f_ret)
  # LINK 2: lambda from the X-block conditional
  fx <- lm(reformulate(c("Y","logX2","logX3","Z1","Z2"), response="L"),
           data = transform(d, L = L))
  lambda <- coef(fx)[["Z1"]]
  cat(sprintf("%-9s %8.4f %8.4f | %+10.5f %+10.5f %+8.5f | %+8.4f\n",
              cell, beta, f_ret, pred, meas, meas-pred, lambda))
}
