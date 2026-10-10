#!/usr/bin/env Rscript
# =============================================================================
# Step 6's MI pooling machinery -- the first test it has ever had
# -----------------------------------------------------------------------------
# WHY THIS EXISTS. Twenty-six validation tracks have exercised the imputation
# and never once exercised the code that COMBINES the results. Step 6 does not
# concatenate draws: it weights them 1/(m*K_i), applies a finite-m variance
# correction on a support-respecting transform, and gates that correction on a
# bimodality diagnostic. None of those three paths has ever run in validation,
# because no validation cell produces a non-Gaussian pooled posterior
# (THEORY.md 6b). So the machinery is unverified shipped code.
#
# WHAT IS CHECKED AGAINST WHAT. Nothing is checked against Step 6's own
# arithmetic. The targets are Rubin's rules, computed independently here:
#
#     Ubar = mean within-imputation variance
#     B    = between-imputation variance of the per-imputation means
#     T    = Ubar + (1 + 1/m) B          <- the total the correction must reach
#
# and the inputs are synthetic posteriors with known moments, so every target
# is known before the code runs.
#
# THE HELPERS ARE EXTRACTED FROM THE REAL FILE BY MARKER, not copied. If
# 06_posterior_summary.R is refactored this test fails loudly rather than
# quietly testing a stale duplicate.
# =============================================================================
suppressMessages({library(dplyr); library(tibble)})
root <- dirname(dirname(normalizePath(
  sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))))
# The pooling logic now lives in 00_common_functions.R as mi_pool_draws(), the
# single source of truth shared by Step 6 (model parameters) and Step 10
# (derived monotonic odds ratios). Step 6's pool_one_parameter() is still
# extracted from its own file by marker so BOTH layers are covered: the shared
# arithmetic, and Step 6's use of it.
paths <- list(cache = tempdir(), logs = tempdir(), objects = tempdir())
suppressMessages(suppressWarnings(source(file.path(root, "00_common_functions.R"))))
src <- readLines(file.path(root, "06_posterior_summary.R"))


i0 <- grep("^  weighted_mean_ +<- mi_weighted_mean", src)
i1 <- grep("^  pool_one_parameter <- function", src)
if (length(i0) != 1L || length(i1) != 1L)
  stop("cannot locate Step 6's pooling helpers -- 06_posterior_summary.R has been ",
       "refactored and this test must be updated rather than skipped. ",
       "(This firing is the test working: it refuses to silently check a stale copy.)")
# pool_one_parameter ends at the first line that is exactly "  }" after the
# tibble it returns; find it by brace balance from i1.
depth <- 0L; i2 <- NA_integer_
for (k in i1:length(src)) {
  depth <- depth + lengths(regmatches(src[k], gregexpr("\\{", src[k]))) -
                   lengths(regmatches(src[k], gregexpr("\\}", src[k])))
  if (k > i1 && depth == 0L) { i2 <- k; break }
}
if (is.na(i2)) stop("could not find the end of pool_one_parameter()")
classify_parameter <- function(p) "fixed"
log_msg <- function(...) invisible(NULL)
eval(parse(text = paste(src[i0:i2], collapse = "\n")), envir = globalenv())

fails <- 0L
ok <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "FAIL", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- fails + 1L
}
cat(sprintf("extracted Step 6 helpers from lines %d-%d\n", i0, i2))

#' Synthetic draws: m imputations, imputation i centred at mu[i] with SD `sw`,
#' and DELIBERATELY UNEQUAL draw counts so the 1/(m*K_i) weighting is exercised.
mk <- function(mu, sw = 1, K = NULL, seed = 1) {
  set.seed(seed); m <- length(mu)
  K <- K %||% rep(400L, m)
  do.call(rbind, lapply(seq_len(m), function(i) {
    z <- stats::rnorm(K[i]); z <- (z - mean(z)) / stats::sd(z)   # exact moments
    tibble(imputation = i, b_x = mu[i] + sw * z)
  }))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
rubin <- function(d) {
  p <- d %>% group_by(imputation) %>%
    summarise(Q = mean(b_x), U = var(b_x), .groups = "drop")
  m <- nrow(p); Ubar <- mean(p$U); B <- var(p$Q)
  c(Ubar = Ubar, B = B, T = Ubar + (1 + 1/m) * B, m = m)
}

cat("\n=== the finite-m correction reaches Rubin's total ===\n")
# Unimodal, spans zero -> identity transform, so the corrected draws' variance
# is directly comparable to T with no transform in the way.
d <- mk(mu = c(-0.30, -0.10, 0.10, 0.30), sw = 1, K = c(200L, 400L, 600L, 800L))
r <- rubin(d)
res <- pool_one_parameter(d, "b_x", 0.95, 0.025, NULL, as.integer(r[["m"]]))
ok("the correction is applied on a unimodal, zero-spanning parameter",
   isTRUE(res$variance_corrected) && res$transform_used == "identity",
   sprintf("transform=%s, BC=%.3f", res$transform_used, res$bimodality_coef))
ok("between_var matches an independent Rubin calculation",
   abs(res$between_var - r[["B"]]) < 1e-9, sprintf("%.6f vs %.6f", res$between_var, r[["B"]]))
ok("within_var matches an independent Rubin calculation",
   abs(res$within_var - r[["Ubar"]]) < 1e-9, sprintf("%.6f vs %.6f", res$within_var, r[["Ubar"]]))
# THE CENTRAL CLAIM: the corrected draws carry Rubin's total variance.
ok("corrected SD^2 equals Ubar + (1 + 1/m)B",
   abs(res$SD^2 / r[["T"]] - 1) < 0.02,
   sprintf("achieved %.5f vs target %.5f (%+.2f%%)", res$SD^2, r[["T"]],
           100*(res$SD^2/r[["T"]] - 1)))
# and it must be strictly LARGER than the uncorrected mixture
ok("and is larger than the uncorrected mixture variance",
   res$SD^2 > r[["Ubar"]] + r[["B"]],
   sprintf("%.5f > %.5f", res$SD^2, r[["Ubar"]] + r[["B"]]))

cat("\n=== every imputation contributes 1/m, whatever its draw count ===\n")
# K_i varies 4x above. If the weighting were per-draw, the 800-draw imputation
# would dominate and the pooled mean would be pulled toward its centre.
ok("the pooled mean is the unweighted mean of the per-imputation means",
   abs(res$Mean - mean(c(-0.30, -0.10, 0.10, 0.30))) < 0.01,
   sprintf("%.5f vs 0.00000", res$Mean))
d_skew <- mk(mu = c(0, 0, 0, 2), K = c(100L, 100L, 100L, 3000L), seed = 3)
res_skew <- pool_one_parameter(d_skew, "b_x", 0.95, 0.025, NULL, 4L)
ok("a 30x draw-count imbalance does not tilt the pooled mean",
   abs(res_skew$Mean - 0.5) < 0.05,
   sprintf("%.4f (per-draw weighting would give ~1.8)", res_skew$Mean))

cat("\n=== the bimodality gate ===\n")
d_bi <- mk(mu = c(-4, -4, 4, 4), sw = 0.4, seed = 7)
res_bi <- pool_one_parameter(d_bi, "b_x", 0.95, 0.025, NULL, 4L)
ok("a clearly bimodal pooled posterior is NOT corrected",
   !isTRUE(res_bi$variance_corrected) && res_bi$bimodality_coef > 0.555,
   sprintf("BC=%.3f (threshold 0.555)", res_bi$bimodality_coef))
ok("a unimodal one IS corrected", isTRUE(res$variance_corrected),
   sprintf("BC=%.3f", res$bimodality_coef))
ok("m = 1 is never corrected (B is undefined)",
   {r1 <- pool_one_parameter(mk(mu = 0, K = 500L), "b_x", 0.95, 0.025, NULL, 1L)
    !isTRUE(r1$variance_corrected)})

cat("\n=== the support-respecting transform ===\n")
d_pos <- mk(mu = c(1.6, 1.8, 2.0, 2.2), sw = 0.25, seed = 11)
res_pos <- pool_one_parameter(d_pos, "b_x", 0.95, 0.025, NULL, 4L)
ok("a strictly positive parameter uses the log transform",
   res_pos$transform_used == "log", res_pos$transform_used)
ok("and stays strictly positive after correction", res_pos$CI_low > 0,
   sprintf("CI_low=%.4f", res_pos$CI_low))

# A PROBABILITY-SCALE parameter: every draw in (0,1). The documented intent is
# a logit transform, so the correction cannot push the interval through 1.
d_prob <- mk(mu = c(0.88, 0.90, 0.92, 0.94), sw = 0.03, seed = 13)
d_prob$b_x <- pmin(pmax(d_prob$b_x, 1e-6), 1 - 1e-6)
res_prob <- pool_one_parameter(d_prob, "b_x", 0.95, 0.025, NULL, 4L)
ok("a (0,1) parameter uses the LOGIT transform",
   res_prob$transform_used == "logit",
   sprintf("got '%s'", res_prob$transform_used))
ok("and the corrected interval stays below 1",
   res_prob$CI_high < 1, sprintf("CI_high=%.5f", res_prob$CI_high))

cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL PASS" else "FAILURES",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
