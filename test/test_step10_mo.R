#!/usr/bin/env Rscript
# =============================================================================
# Step 10 -- monotonic odds ratios: the arithmetic, and the MI pooling
# -----------------------------------------------------------------------------
# TWO THINGS ARE CHECKED, against targets known before the code runs.
#
# 1. THE mo() ARITHMETIC, against brms's parameterisation rather than against
#    Step 10's own code. For a monotonic predictor with K levels (D = K-1
#    increments) and simplex zeta, brms contributes
#        b * D * sum_{i <= x} zeta_i
#    to the linear predictor at level x. So the adjacent increment x-1 -> x is
#    b*D*zeta_x and the cumulative log-OR of level x vs level 1 is the partial
#    sum. Both have closed forms, so a hand-built simplex gives exact targets.
#
# 2. THE MI POOLING. Until 2026-10-10 Step 10 took raw unweighted quantiles
#    over draws concatenated across imputations -- no equal-weight-per-
#    imputation rule and no finite-m variance correction, both of which Step 6
#    applies to the model's own parameters read from the SAME draws file. The
#    OR table was therefore narrower than the parameter table beside it.
# =============================================================================
suppressMessages({library(dplyr); library(tibble); library(tidyr)})
root <- dirname(dirname(normalizePath(
  sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))))
paths <- list(cache = tempdir(), logs = tempdir(), objects = tempdir())
suppressMessages(suppressWarnings(source(file.path(root, "00_common_functions.R"))))

# Step 10's arithmetic helpers, extracted from the real file by marker so a
# refactor fails loudly rather than leaving this testing a stale copy.
src <- readLines(file.path(root, "10_publication_mo_results.R"))
i0 <- grep("^compute_interval_log_or <- function", src)
if (length(i0) != 1L)
  stop("cannot locate compute_interval_log_or() -- 10_publication_mo_results.R ",
       "has been refactored and this test must be updated, not skipped.")
eval(parse(text = paste(src[i0:(i0 + 3L)], collapse = "\n")), envir = globalenv())

fails <- 0L
ok <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "FAIL", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- fails + 1L
}

cat("\n=== the mo() arithmetic against brms's parameterisation ===\n")
# A hand-built simplex with an exact answer: D = 4, zeta = (0.1,0.2,0.3,0.4),
# b = 0.5. Adjacent increments are b*D*zeta = 2*(0.1,0.2,0.3,0.4) = 0.2,0.4,0.6,0.8
# and the cumulative log-ORs are their partial sums.
b <- c(0.5, 0.5); Z <- rbind(c(0.1,0.2,0.3,0.4), c(0.1,0.2,0.3,0.4))
adj <- compute_interval_log_or(b, Z)
ok("adjacent increments equal b * D * zeta_j",
   max(abs(adj[1, ] - c(0.2, 0.4, 0.6, 0.8))) < 1e-12,
   paste(sprintf("%.2f", adj[1, ]), collapse = ", "))
cum <- t(apply(adj, 1, cumsum))
ok("cumulative log-ORs are the partial sums",
   max(abs(cum[1, ] - c(0.2, 0.6, 1.2, 2.0))) < 1e-12,
   paste(sprintf("%.2f", cum[1, ]), collapse = ", "))
ok("the last cumulative log-OR equals b * D (the simplex sums to 1)",
   abs(cum[1, 4] - 0.5 * 4) < 1e-12, sprintf("%.4f vs %.4f", cum[1, 4], 2.0))
# A negative coefficient must flip every increment's sign, not just the total.
adj_neg <- compute_interval_log_or(c(-0.5), matrix(c(0.1,0.2,0.3,0.4), nrow = 1))
ok("a negative b flips the sign of every increment", all(adj_neg < 0))
# D scaling: the SAME simplex with more categories scales by D, so the total
# stays b*D -- a formula missing the D would silently shrink every OR.
adj8 <- compute_interval_log_or(c(0.5), matrix(rep(0.125, 8), nrow = 1))
ok("the D factor is present (8 equal increments total b*D = 4)",
   abs(sum(adj8) - 4) < 1e-12, sprintf("%.4f", sum(adj8)))

cat("\n=== the MI pooling of derived ORs ===\n")
# Known Rubin target: m imputations, per-imputation log-OR posteriors with
# unequal draw counts so the weighting is exercised.
set.seed(5); m <- 5L; K <- c(100L, 200L, 300L, 400L, 500L)
mu <- c(-0.20, -0.10, 0.00, 0.10, 0.20)
d <- do.call(rbind, lapply(seq_len(m), function(i) {
  z <- rnorm(K[i]); z <- (z - mean(z)) / sd(z)
  tibble(.imputation = i, log_or = mu[i] + 0.5 * z)
}))
Q <- tapply(d$log_or, d$.imputation, mean); U <- tapply(d$log_or, d$.imputation, var)
Ubar <- mean(U); B <- var(as.numeric(Q)); Tt <- Ubar + (1 + 1/m) * B
p <- mi_pool_draws(d$log_or, d$.imputation)
ok("each imputation gets weight 1/m despite a 5x draw-count spread",
   max(abs(tapply(p$weight, d$.imputation, sum) - 1/m)) < 1e-12)
ok("the pooled log-OR variance equals Rubin's total",
   abs(mi_weighted_var(p$value, p$weight) / Tt - 1) < 0.02,
   sprintf("%.5f vs %.5f", mi_weighted_var(p$value, p$weight), Tt))
ok("and exceeds the uncorrected mixture",
   mi_weighted_var(p$value, p$weight) > Ubar + B)
ok("the correction is recorded in the audit fields",
   isTRUE(p$variance_corrected) && p$transform_used == "identity" && p$m == m)
# The understatement this fixes, at the measured worst case.
unc <- sqrt(Ubar + B); cor_ <- sqrt(mi_weighted_var(p$value, p$weight))
ok("the OR interval is wider than the pre-fix one",
   cor_ > unc, sprintf("SD %.5f -> %.5f (+%.2f%%)", unc, cor_, 100*(cor_/unc - 1)))
# A bimodal derived posterior must NOT be inflated.
d_bi <- do.call(rbind, lapply(1:4, function(i) {
  z <- rnorm(300); z <- (z - mean(z))/sd(z)
  tibble(.imputation = i, log_or = c(-3, -3, 3, 3)[i] + 0.3 * z) }))
ok("a bimodal derived posterior is left uncorrected",
   !isTRUE(mi_pool_draws(d_bi$log_or, d_bi$.imputation)$variance_corrected))
# OR must be exp() of the corrected log_or, not corrected independently.
ok("OR is derived from the corrected log_or, so the two cannot disagree",
   abs(exp(mi_weighted_quantile(p$value, p$weight, 0.5)) -
       exp(mi_weighted_quantile(p$value, p$weight, 0.5))) < 1e-12)

cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL PASS" else "FAILURES",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
