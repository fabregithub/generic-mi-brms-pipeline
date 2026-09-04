#!/usr/bin/env Rscript
# =============================================================================
# V18 -- DERIVE the n-scaling prediction. Closed form, no simulation.
# -----------------------------------------------------------------------------
# THE CLAIM UNDER TEST. V17 found the shipped default biased +5% to +10.4% where
# the covariate is a confounder or a mediator, with coverage still 0.93-0.97 and
# honestly-sized intervals. The proposed explanation was that the concealment is
# a SAMPLE-SIZE ACCIDENT: relative bias is ~constant in n while the SE falls as
# 1/sqrt(n), so bias/SE grows and coverage must eventually collapse.
#
# That explanation is an assumption, not a measurement, and it has a competitor:
# the bias could be FINITE-SAMPLE -- an artefact of imputing from n = 800 rows --
# in which case it shrinks with n and coverage holds.
#
# The two make opposite predictions, so the run discriminates them. This script
# derives the first one exactly.
#
# THE ARITHMETIC. With the estimate ~ N(truth + b, s^2) and the pooled interval
# half-width ~ 1.96 * w * s (w = the measured width/SE ratio),
#
#     coverage = Phi(1.96w - b/s) - Phi(-1.96w - b/s)
#
# and under the constant-relative-bias hypothesis b/s scales as sqrt(n/800).
# Calibrating at n = 800 against V17's MEASURED coverage is the check that this
# formula describes the harness at all -- if it cannot reproduce the row we have,
# it cannot be trusted for the rows we do not.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
setwd(.here)

s <- read.csv("results/v17a_summary.csv", stringsAsFactors = FALSE)
s <- subset(s, procedure == "pipeline_bartMI" &
              scenario %in% c("zr_fork", "zr_pipe", "zrmar_fork", "zrmar_pipe"))
TRUE_B <- 0.40; N0 <- 800L
NS <- c(800L, 3200L, 12800L)

cov_at <- function(r, w) stats::pnorm(1.96 * w - r) - stats::pnorm(-1.96 * w - r)

cat("\n  CALIBRATION at n = 800 -- does the formula reproduce V17's coverage?\n")
cat(sprintf("  %-14s %8s %8s %9s %9s %8s\n",
            "cell", "bias/se", "width/se", "predicted", "measured", "diff"))
ok <- TRUE
for (i in seq_len(nrow(s))) {
  r <- abs(s$rel_bias[i] * TRUE_B) / s$emp_se[i]; w <- s$width_se_ratio[i]
  p <- cov_at(r, w); m <- s$coverage[i]
  if (abs(p - m) > 0.02) ok <- FALSE
  cat(sprintf("  %-14s %8.3f %8.3f %9.3f %9.3f %+8.3f\n",
              s$scenario[i], r, w, p, m, p - m))
}
cat(sprintf("  -> %s\n", if (ok)
  "every cell within 0.02; the formula describes the harness." else
  "SOME CELL OFF BY >0.02 -- do not trust the projections below."))

cat("\n  PREDICTION under CONSTANT RELATIVE BIAS (bias/se grows as sqrt(n/800))\n")
cat(sprintf("  %-14s %10s %10s %10s\n", "cell", "n=800", "n=3200", "n=12800"))
pred <- list()
for (i in seq_len(nrow(s))) {
  r0 <- abs(s$rel_bias[i] * TRUE_B) / s$emp_se[i]; w <- s$width_se_ratio[i]
  cv <- vapply(NS, function(n) cov_at(r0 * sqrt(n / N0), w), 0)
  pred[[s$scenario[i]]] <- cv
  cat(sprintf("  %-14s %10.3f %10.3f %10.3f\n", s$scenario[i], cv[1], cv[2], cv[3]))
}
cat("\n  THE COMPETING HYPOTHESIS -- bias is finite-sample and shrinks as 1/sqrt(n):\n")
cat("  then bias/se is CONSTANT, and coverage stays at its n = 800 value at every n.\n")
cat(sprintf("  The two differ by %.2f to %.2f coverage points at n = 12800 --\n",
            min(vapply(pred, function(c) c[1] - c[3], 0)),
            max(vapply(pred, function(c) c[1] - c[3], 0))))
cat("  far larger than the ~0.03 MC error of coverage at 200 reps.\n")

cat("\n  ALSO PREDICTED, and the sharper test of the two hypotheses:\n")
cat("  relative bias itself. Constant-bias says it holds; finite-sample says it\n")
cat("  falls as 1/sqrt(n). At n = 12800 (16x) that is a 4x drop:\n")
cat(sprintf("  %-14s %10s %10s %10s\n", "cell", "n=800", "constant", "finite-sample"))
for (i in seq_len(nrow(s)))
  cat(sprintf("  %-14s %9.2f%% %9.2f%% %13.2f%%\n", s$scenario[i],
              100 * s$rel_bias[i], 100 * s$rel_bias[i], 100 * s$rel_bias[i] / 4))
cat("\n  Reading relative bias needs no coverage model and no normality\n")
cat("  assumption, so it is the primary outcome; coverage is the consequence.\n")

cat("\n  ASSUMPTIONS, declared:\n")
cat("   1. width/SE is taken as constant in n. It was 0.98-1.13 at n = 800; if it\n")
cat("      drifts, the coverage projections move but the bias test does not.\n")
cat("   2. Normality of the pooled estimate -- safe at these n.\n")
cat("   3. The oracle stays unbiased at every n. It is the control; if it drifts,\n")
cat("      the DGP or the estimand is n-dependent and nothing else is readable.\n")
saveRDS(list(pred = pred, ns = NS), "results/predictions_v18.rds")
