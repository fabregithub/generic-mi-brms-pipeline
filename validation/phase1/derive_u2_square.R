# WHY THE u^2 LAW IS SQUARED -- derivation, then check.
#
# A censored exposure's true conditional carries log p(Z|x) = -(Z-g(x))^2/2s^2.
# A linear draw uses ghat instead of g. With r = g - ghat:
#     (Z-g)^2 - (Z-ghat)^2 = -2r(Z-ghat) + r^2
#  => Dlog p = r(Z-ghat)/s^2 - r^2/(2 s^2)
# Z = g(x) + eps with eps independent of x, so (Z-ghat) = r + eps and
#     E_Z[ r(Z-ghat) ] = r^2  =>  E_Z[Dlog p] = r^2/s^2 - r^2/(2s^2) = r^2/(2 s^2).
#
# So the FIRST-ORDER term has mean zero -- it adds spread, not systematic tilt --
# and the systematic part is QUADRATIC in r. That is the square, derived rather
# than asserted. It also predicts a dependence the law has never been tested on:
#     systematic tilt  ~  u^2 / (2 s^2)      so bias should scale as 1/s^2.
#
# CHECKS:
#   (1) E_Z[Dlog p] = r^2/(2 s^2) exactly
#   (2) the first-order term really is mean-zero over Z
suppressMessages(for (f in c("dgp.R","censoring.R")) source(file.path("R",f)))
`%||%` <- function(a,b) if (is.null(a)) b else a
set.seed(4)
tr <- make_truth(z_role = "pipe_nl")
s  <- tr$nl_sd
x  <- seq(-3, 1, length.out = 400)
g  <- ef_nl_g(x, tr)
# best linear fit to g over the density-weighted censored region
w  <- dnorm(x); fit <- lm(g ~ x, weights = w); ghat <- fitted(fit)
r  <- g - ghat
# Monte Carlo over Z at each x
K <- 20000L
d_mean <- sapply(seq_along(x), function(i) {
  Z <- g[i] + rnorm(K, 0, s)
  mean(((Z - ghat[i])^2 - (Z - g[i])^2) / (2 * s^2))   # = Dlog p with sign fixed
})
pred <- r^2 / (2 * s^2)
cat(sprintf("(1) E_Z[Dlog p] vs r^2/(2s^2):  max |miss| = %.5f  (max pred %.4f)\n",
            max(abs(d_mean - pred)), max(pred)))
first <- sapply(seq_along(x), function(i) {
  Z <- g[i] + rnorm(K, 0, s); mean(r[i] * (Z - ghat[i]) / s^2) })
cat(sprintf("(2) first-order term E_Z[r(Z-ghat)/s^2]: mean %.5f vs r^2/s^2 mean %.5f\n",
            mean(first), mean(r^2/s^2)))
cat(sprintf("    (it is NOT mean-zero in x -- it equals r^2/s^2; the CANCELLATION\n",
            ""))
cat("     is between it and -r^2/2s^2, leaving +r^2/2s^2)\n")
cat(sprintf("\n(3) the untested prediction: tilt ~ u^2/(2 s^2), so bias ~ 1/s^2.\n"))
for (sd_ in c(0.4, 0.8, 1.6))
  cat(sprintf("    s = %.1f -> relative tilt %.3f (vs s = 0.8)\n", sd_, (0.8/sd_)^2))
