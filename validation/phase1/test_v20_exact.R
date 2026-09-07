#!/usr/bin/env Rscript
# =============================================================================
# V20 -- are the "exact" conditionals actually exact?
# -----------------------------------------------------------------------------
# THE WHOLE DECOMPOSITION RESTS ON THIS. V20 attributes bias by swapping the
# shipped draw for the true conditional one block at a time. If the "exact" draw
# is not exact, the attribution is meaningless -- and it would still return
# believable numbers, which is exactly how the four leftcens incidents in
# FINDINGS_v13.md happened.
#
# THE KNOWN ANSWER. Given Z2, these DGPs are jointly Gaussian, so the true
# conditional mean IS the population regression of the target on the
# conditioning set, and the true conditional SD IS that regression's residual
# SD. A large simulated sample therefore checks the algebra directly, with no
# appeal to the algebra being checked.
#
# It also checks the thing that motivated the file: that V12/V13's
# `.smc_draw_z` -- correct in the precision structure -- is WRONG under a fork.
# If that comparison ever stops failing, either .smc_draw_z gained the missing
# p(logX1 | Z1) factor or this test stopped testing anything.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
for (f in c("dgp.R", "censoring.R", "exact_fork.R")) source(file.path(.here, "R", f))
`%||%` <- function(a, b) if (is.null(a)) b else a

fails <- 0L
ok <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "FAIL", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- fails + 1L
}
N <- 300000L

for (role in c("precision", "fork", "pipe")) {
  cat(sprintf("\n=== %s ===\n", role))
  set.seed(4242)
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role)
  d  <- simulate_complete(N, tr)$data
  J  <- .ef_joint(tr)

  # --- the analytic joint covariance must match the empirical one -------------
  emp <- cov(d[, c("Z1", "logX1", "logX2", "logX3", "Y")])
  # Z2 is conditioned out of the analytic joint, so add back its contribution
  # to Y's variance and to Cov(., Y) before comparing.
  g2 <- J$gamma2; vz2 <- 0.25
  ana <- J$V
  ana["Y", "Y"] <- ana["Y", "Y"] + g2^2 * vz2
  # Compare on an ABSOLUTE scale. A relative bar is the wrong instrument here:
  # several true covariances are exactly 0 (Z1 vs logX2/logX3 off a fork), so
  # dividing sampling noise of order 1/sqrt(N) by a near-zero denominator
  # manufactures a large "relative" error out of a correct entry. Monte-Carlo
  # noise on a covariance at N = 3e5 is ~0.002, so 0.02 is a wide, honest bar.
  ok(sprintf("%-9s analytic Cov matches empirical (max |diff| < 0.02)", role),
     max(abs(ana - emp)) < 0.02,
     sprintf("max |diff| %.4f", max(abs(ana - emp))))

  # --- conditional mean == population regression; SD == residual SD ----------
  # Z1 | logX, Y, Z2
  fit <- lm(Z1 ~ logX1 + logX2 + logX3 + Y + Z2, data = d)
  cc  <- .ef_cond(J, "Z1", as.matrix(d[1:5000, c("logX1","logX2","logX3","Y")]),
                  z2 = d$Z2[1:5000])
  pred <- predict(fit, newdata = d[1:5000, ])
  ok(sprintf("%-9s exact Z1 conditional MEAN matches the regression", role),
     max(abs(cc$mean - pred)) < 0.02,
     sprintf("max |diff| = %.5f", max(abs(cc$mean - pred))))
  ok(sprintf("%-9s exact Z1 conditional SD matches the residual SD", role),
     abs(cc$sd - summary(fit)$sigma) < 0.01,
     sprintf("analytic %.4f vs residual %.4f", cc$sd, summary(fit)$sigma))

  # logX1 | Z1, other X, Y, Z2
  fx <- lm(logX1 ~ Z1 + logX2 + logX3 + Y + Z2, data = d)
  cx <- .ef_cond(J, "logX1", as.matrix(d[1:5000, c("Z1","logX2","logX3","Y")]),
                 z2 = d$Z2[1:5000])
  px <- predict(fx, newdata = d[1:5000, ])
  ok(sprintf("%-9s exact logX1 conditional MEAN matches the regression", role),
     max(abs(cx$mean - px)) < 0.02, sprintf("max |diff| = %.5f", max(abs(cx$mean - px))))
  ok(sprintf("%-9s exact logX1 conditional SD matches the residual SD", role),
     abs(cx$sd - summary(fx)$sigma) < 0.01,
     sprintf("analytic %.4f vs residual %.4f", cx$sd, summary(fx)$sigma))

  # --- the draws themselves recover the conditional -------------------------
  set.seed(9)
  w <- d[1:60000, ]
  z <- ef_draw_z1(w, tr, seq_len(nrow(w)), J)
  ok(sprintf("%-9s Z1 draws have the right marginal SD", role),
     abs(sd(z) - sd(d$Z1)) < 0.03, sprintf("draw %.3f vs true %.3f", sd(z), sd(d$Z1)))
}

cat("\n=== the gap that motivated this file: .smc_draw_z under a fork ===\n")
if (file.exists(file.path(.here, "R", "smc_impute.R"))) {
  source(file.path(.here, "R", "smc_impute.R"))
  set.seed(11)
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = "fork")
  d  <- simulate_complete(40000L, tr)$data
  J  <- .ef_joint(tr)
  fit <- lm(Z1 ~ logX1 + logX2 + logX3 + Y + Z2, data = d)
  truth_sd <- summary(fit)$sigma
  ef <- .ef_cond(J, "Z1", as.matrix(d[, c("logX1","logX2","logX3","Y")]), z2 = d$Z2)
  # What .smc_draw_z assumes: prior N(0,1) for Z1, likelihood p(Y|Z1,rest) only.
  eta_z <- function(zv) { s2 <- d; s2$Z1 <- zv; .smc_eta(s2, tr, "oracle", NULL,
                                                         dgp_formula(tr),
                                                         paste0("logX",1:3), c("Z1","Z2")) }
  qc <- .smc_quadratic_coefs(eta_z, nrow(d))
  old <- .smc_draw_z(d$Y, qc$A, qc$B, qc$C, 1, m = rep(0, nrow(d)), s = 1,
                     mode = "exact")
  cat(sprintf("  true conditional SD of Z1 | X, Y, Z2  = %.4f\n", truth_sd))
  cat(sprintf("  exact_fork.R analytic SD              = %.4f\n", ef$sd))
  cat(sprintf("  .smc_draw_z implied SD (draws)        = %.4f\n", sd(old - ef$mean)))
  ok("exact_fork.R matches the true conditional SD", abs(ef$sd - truth_sd) < 0.01)
  # The MARGINAL SD of the old draws is ~1 and looks fine -- which is the point,
  # and why this had to be checked deliberately. The diagnostic is the spread
  # AROUND THE TRUE CONDITIONAL MEAN: if the draw were exact that spread would
  # equal the true conditional SD.
  ok(".smc_draw_z is WRONG under a fork (it ignores p(logX1 | Z1))",
     abs(sd(old - ef$mean) - truth_sd) > 0.10,
     sprintf("spread around the true conditional mean %.4f vs true %.4f (%.0f%% too wide); its marginal SD %.3f looks fine",
             sd(old - ef$mean), truth_sd,
             100 * (sd(old - ef$mean) / truth_sd - 1), sd(old)))
} else ok("smc_impute.R present", FALSE)

cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL PASS" else "FAILURES",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
