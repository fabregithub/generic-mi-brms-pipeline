#!/usr/bin/env Rscript
# V35, second view: bias against ANALYTIC TRUTH rather than against each arm's own oracle.
#
# analyze_v35.R scores the three registered predictions, all of which are stated as
# excess over the arm's OWN oracle -- the right currency for "what does the imputation
# cost this analysis model". It is the wrong currency for "what does the shipped engine
# report to a user", because it divides out whatever the analysis model does on complete
# data. V35 measured that term for the first time (+13 to +17%, flat in f), so the
# vs-truth view is now computable and is not the same curve.
#
# Nothing here is a registered gate. It is a re-read of a completed run.
suppressPackageStartupMessages(library(stats))
tag <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "v35"
s <- read.csv(sprintf("results/%s_summary.csv", tag), stringsAsFactors = FALSE)
fmap <- c(mixf10=0.10, nd20=0.20, mixf30=0.30, nd40=0.40, mixf50=0.50, mixf60=0.60)
s$f <- fmap[s$scenario]
s <- s[!is.na(s$f), ]
s$rel_pct  <- 100 * s$rel_bias
s$mcse_pct <- 100 * s$mcse_bias / abs(s$true)
s$bias_se  <- s$bias / s$emp_se          # the PLAN section 11b currency; gate is 0.3

arms  <- c("pipeline_bkmr","oracle_bkmr","pipeline_lm","oracle_lm")
ests  <- unique(s$estimand)
fs    <- sort(unique(s$f))

cat("=============== V35 vs ANALYTIC TRUTH ===============\n")
cat("rel% (mcse) [bias/SE] cov   -- 100 reps/cell\n")
for (e in ests) {
  cat(sprintf("\n-- %s --\n", e))
  cat(sprintf("%6s", "f"))
  for (a in arms) cat(sprintf(" %26s", a))
  cat("\n")
  for (ff in fs) {
    cat(sprintf("%6.2f", ff))
    for (a in arms) {
      r <- s[s$f==ff & s$estimand==e & s$procedure==a, ]
      if (!nrow(r)) { cat(sprintf(" %26s", "-")); next }
      cat(sprintf(" %+8.1f (%4.1f) [%+5.2f] %.2f",
                  r$rel_pct, r$mcse_pct, r$bias_se, r$coverage))
    }
    cat("\n")
  }
}

# --- sign changes in the SHIPPED arm, with MC error --------------------------
cat("\n=============== SIGN CHANGES (pipeline_bkmr vs truth) ===============\n")
cat("A sign that flips with a design parameter is PLAN section 11b's worst category.\n")
for (e in ests) {
  d <- s[s$procedure=="pipeline_bkmr" & s$estimand==e, ]
  d <- d[order(d$f), ]
  if (nrow(d) < 2) next
  # Count a flip only between cells that are each individually significant, and in
  # EITHER direction -- an earlier version tested only positive-then-negative and so
  # silently passed estimands that run negative-then-positive.
  z   <- d$rel_pct / d$mcse_pct
  sig <- abs(z) > 2
  if (sum(sig) < 2 || length(unique(sign(d$rel_pct[sig]))) < 2) {
    cat(sprintf("  %-16s no significant flip (range %+.1f to %+.1f)\n",
                e, min(d$rel_pct), max(d$rel_pct)))
    next
  }
  si <- which(sig)
  # Declare the flip from the significant cells; LOCATE the zero from the adjacent
  # bracketing pair, so an intermediate non-significant cell still informs the crossing.
  adj   <- which(sign(d$rel_pct[-1]) != sign(d$rel_pct[-nrow(d)]))
  parts <- character(0)
  for (k in adj) {
    f0 <- d$f[k]; f1 <- d$f[k+1]; y0 <- d$rel_pct[k]; y1 <- d$rel_pct[k+1]
    fc <- f0 + (f1-f0) * y0/(y0-y1)
    lo <- max(si[si <= k]); hi <- min(si[si >= k+1])
    parts <- c(parts, sprintf("%+.1f%% (%.1f sd) at f=%.2f -> %+.1f%% (%.1f sd) at f=%.2f [zero near f=%.2f]",
                              d$rel_pct[lo], abs(z[lo]), d$f[lo],
                              d$rel_pct[hi], abs(z[hi]), d$f[hi], fc))
  }
  cat(sprintf("  %-16s FLIPS: %s\n", e, paste(parts, collapse = "; ")))
}

# --- is the complete-data term flat in f? ------------------------------------
cat("\n=============== BKMR's COMPLETE-DATA TERM (oracle_bkmr vs truth) ===============\n")
cat("If this is flat in f, the vs-truth curve is 'flat offset + growing attenuation',\n")
cat("and the near-zero cell is a CANCELLATION, not accuracy.\n")
for (e in ests) {
  d <- s[s$procedure=="oracle_bkmr" & s$estimand==e, ]; d <- d[order(d$f), ]
  if (nrow(d) < 3) next
  fit <- lm(rel_pct ~ f, data = d, weights = 1/d$mcse_pct^2)
  ci  <- confint(fit)[2, ]
  cat(sprintf("  %-16s mean %+6.1f%%  | trend %+7.1f %%/unit f [%+.1f, %+.1f] %s\n",
              e, mean(d$rel_pct), coef(fit)[2], ci[1], ci[2],
              if (ci[1] <= 0 && ci[2] >= 0) "FLAT" else "trending"))
}

# --- where does the shipped arm cross the project's own 0.3 bias/SE bar? -----
cat("\n=============== bias/SE vs the 0.3 bar (pipeline_bkmr) ===============\n")
for (e in ests) {
  d <- s[s$procedure=="pipeline_bkmr" & s$estimand==e, ]; d <- d[order(d$f), ]
  ok <- d$f[abs(d$bias_se) <= 0.3]
  cat(sprintf("  %-16s within 0.3 at f = %s\n", e,
              if (length(ok)) paste(sprintf("%.2f", ok), collapse=", ") else "NONE"))
}
