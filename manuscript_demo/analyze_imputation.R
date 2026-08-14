#!/usr/bin/env Rscript
# =============================================================================
# Single-pollutant demonstration: the four censored-exposure methods vs oracle,
# on the paper-grounded JECS PFAS -> Kawasaki disease data. The key contrast is
# impute_noY vs impute_congenial (the Y-inclusion / congeniality test).
#
# Start SMALL (N ~ 2000, few reps) to verify; scale N and reps for the study.
# Env: N (2000), M (10), REPS (5), LIMIT_SCALE (1), CORRELATION (real),
#      METHODS (oracle,substitution,impute_noY,impute_congenial).
# Writes results/single_pollutant_<tag>.csv.
# =============================================================================
suppressWarnings(suppressMessages(library(leftcens)))
.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "." })
source(file.path(.here, "R", "make_demo_data.R"))
source(file.path(.here, "R", "analysis_methods.R"))

geti <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }
gets <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
N     <- as.integer(geti("N", 2000)); M <- as.integer(geti("M", 10)); REPS <- as.integer(geti("REPS", 5))
scale <- geti("LIMIT_SCALE", 1); corr <- gets("CORRELATION", "real")
outcome <- gets("OUTCOME", "kd")                                  # "kd" (rare) or "asthma" (common)
methods <- strsplit(gets("METHODS", "oracle,substitution,impute_noY,impute_congenial"), ",")[[1]]
tag <- sprintf("%s_scale-%.2g_cor-%s_n%d", outcome, scale, corr, N)

cat(sprintf("Single-pollutant study | %s | outcome=%s m=%d reps=%d methods=%s\n",
            tag, outcome, M, REPS, paste(methods, collapse=",")))

acc <- list()
for (r in seq_len(REPS)) {
  dat <- make_demo_data(n = N, seed = 1000 + r, limit_scale = scale, correlation = corr)
  # oracle estimate for this replicate = the scoring reference.
  oracle_fit <- fit_single_pollutant(impute_pfas(dat$data, dat$truth, "oracle", outcome = outcome), dat$truth, outcome = outcome)
  oref <- setNames(oracle_fit$beta, oracle_fit$pfas)
  for (meth in methods) {
    comp <- impute_pfas(dat$data, dat$truth, meth, m = M, seed = 7 * r, outcome = outcome)
    fit  <- fit_single_pollutant(comp, dat$truth, outcome = outcome)
    fit$rep <- r; fit$method <- meth
    fit$beta_oracle <- oref[fit$pfas]
    fit$err_vs_oracle <- fit$beta - fit$beta_oracle
    fit$covers_oracle <- fit$lo <= fit$beta_oracle & fit$beta_oracle <= fit$hi
    acc[[length(acc) + 1]] <- fit
  }
  cat(".")
}
cat("\n")
raw <- do.call(rbind, acc)

# --- aggregate: bias vs oracle + oracle-coverage, per method x pfas ----------
agg <- aggregate(cbind(err_vs_oracle, covers_oracle, or) ~ method + pfas, raw,
                 function(x) mean(x, na.rm = TRUE))
agg <- agg[order(agg$pfas, match(agg$method, methods)), ]

cat("\n=== mean error vs oracle (beta per log2), and oracle-coverage ===\n")
for (j in unique(agg$pfas)) {
  cat(sprintf("-- %s --\n", j))
  s <- agg[agg$pfas == j, ]
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-18s err=%+.3f  cover=%.2f  OR=%.3f\n",
                s$method[i], s$err_vs_oracle[i], s$covers_oracle[i], s$or[i]))
}
# headline: mean |error| across PFAS, by method (the Y-inclusion contrast)
head <- aggregate(abs(err_vs_oracle) ~ method, raw, mean)
cat("\n=== mean |error vs oracle| across all PFAS (lower = better) ===\n")
for (i in order(match(head$method, methods)))
  cat(sprintf("  %-18s %.3f\n", head$method[i], head$`abs(err_vs_oracle)`[i]))

dir.create(file.path(.here, "results"), showWarnings = FALSE)
utils::write.csv(raw, file.path(.here, "results", sprintf("single_pollutant_%s.csv", tag)), row.names = FALSE)
cat(sprintf("\nWrote results/single_pollutant_%s.csv\n", tag))
