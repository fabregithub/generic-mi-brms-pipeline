#!/usr/bin/env Rscript
# =============================================================================
# Track V3 -- fast self-test of the estimand machinery (no MCMC, runs in ~1 s)
# -----------------------------------------------------------------------------
# Everything V3 reports is measured against `h_true()`. If that function and the
# generator's own `eta` construction ever drift apart, every bias in the study is
# wrong by the same silent amount and nothing in the output would show it. These
# checks are the link between the two, plus the algebraic properties the contrast
# construction relies on.
#
#   Rscript test_v3_harness.R
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.here, "R", "dgp.R"))
source(file.path(.here, "R", "estimands_bkmr.R"))

fails <- 0L
ok <- function(label, pass, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (pass) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  --  ", detail) else ""))
  if (!pass) fails <<- fails + 1L
  invisible(pass)
}

cat("\n=== 1. h_true() must reproduce the GENERATOR's exposure term exactly ===\n")
cat("    Rebuilds eta from simulate_complete()'s own output and strips the\n")
cat("    covariate and noise parts, leaving only what h_true() claims to be.\n\n")
for (form in c("mixture", "additive")) {
  set.seed(99)
  truth <- make_truth(p = 3L, q = 2L, erf_form = form)
  s <- simulate_complete(n = 500, truth = truth, rho = 0.4, sd_x = 1, mu_x = 0,
                         sigma_y = 0, z_form = "linear")   # sigma_y = 0: no noise
  d <- s$data
  U <- as.matrix(d[, paste0("logX", 1:3)])
  Zm <- as.matrix(d[, c("Z1", "Z2")])
  # With sigma_y = 0, Y IS eta exactly. Removing the covariate contribution and
  # the intercept must leave h(logX) and nothing else.
  h_from_generator <- d$Y - as.vector(Zm %*% truth$gamma) - truth$intercept
  mx <- max(abs(h_from_generator - h_true(U, truth)))
  ok(sprintf("%-9s h_true() == generator eta - gamma'Z - intercept", form),
     mx < 1e-10, sprintf("max |diff| = %.3g over 500 rows", mx))
}

cat("\n=== 2. Contrast algebra ===\n")
truth <- make_truth(p = 3L, q = 2L, erf_form = "mixture")
g <- bkmr_estimand_grid(truth)
ok("every contrast is zero-sum (so BKMR's additive non-identifiability cancels)",
   all(abs(rowSums(g$C)) < 1e-12),
   sprintf("max |rowSum| = %.3g", max(abs(rowSums(g$C)))))
ok("contrast columns align with the grid's evaluation points",
   ncol(g$C) == nrow(g$Z), sprintf("%d x %d vs %d points", nrow(g$C), ncol(g$C), nrow(g$Z)))
ok("every estimand has a name and a finite truth",
   length(g$names) == length(g$true) && all(is.finite(g$true)))

cat("\n=== 3. Analytic truth vs independently derived closed forms ===\n")
for (cfg in list(list(0, 1), list(0.5, 1.3), list(-0.8, 0.7))) {
  r <- tryCatch(check_bkmr_truth(truth, mu_x = cfg[[1]], sd_x = cfg[[2]]),
                error = function(e) e)
  ok(sprintf("mu_x = %+.1f, sd_x = %.1f", cfg[[1]], cfg[[2]]),
     is.data.frame(r),
     if (is.data.frame(r)) sprintf("max |diff| = %.3g", max(abs(r$diff)))
     else conditionMessage(r))
}

cat("\n=== 4. The additive generator must zero the non-linear estimands ===\n")
ga <- bkmr_estimand_grid(make_truth(p = 3L, q = 2L, erf_form = "additive"))
ok("int_X1X2 == 0 with b_int = 0", abs(ga$true["int_X1X2"]) < 1e-12)
ok("curv_X1  == 0 with b_quad = 0", abs(ga$true["curv_X1"]) < 1e-12)
ok("overall_q75_q25 is IDENTICAL on both generators (the documented degeneracy)",
   abs(ga$true["overall_q75_q25"] - g$true["overall_q75_q25"]) < 1e-12,
   "confirms it does not test the mixture at all")

cat("\n=== 5. The sharp estimands must move with -- and only with -- their parameter ===\n")
t2 <- truth; t2$b_int <- t2$b_int * 2
g2 <- bkmr_estimand_grid(t2)
ok("doubling b_int doubles int_X1X2",
   abs(g2$true["int_X1X2"] - 2 * g$true["int_X1X2"]) < 1e-12)
ok("doubling b_int leaves curv_X1 unchanged",
   abs(g2$true["curv_X1"] - g$true["curv_X1"]) < 1e-12)
t3 <- truth; t3$b_quad <- t3$b_quad * 2
g3 <- bkmr_estimand_grid(t3)
ok("doubling b_quad doubles curv_X1",
   abs(g3$true["curv_X1"] - 2 * g$true["curv_X1"]) < 1e-12)
ok("doubling b_quad leaves int_X1X2 unchanged",
   abs(g3$true["int_X1X2"] - g$true["int_X1X2"]) < 1e-12)
t4 <- truth; t4$b <- t4$b * 3
g4 <- bkmr_estimand_grid(t4)
ok("tripling the LINEAR coefficients leaves int_X1X2 unchanged",
   abs(g4$true["int_X1X2"] - g$true["int_X1X2"]) < 1e-12)
ok("tripling the LINEAR coefficients leaves curv_X1 unchanged",
   abs(g4$true["curv_X1"] - g$true["curv_X1"]) < 1e-12)

cat(sprintf("\n%s  (%d check%s failed)\n\n",
            if (fails == 0L) "ALL CHECKS PASSED" else "FAILURES PRESENT",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
