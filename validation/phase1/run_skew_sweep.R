#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 — skew sweep (PLAN §7.5 design factor)
# -----------------------------------------------------------------------------
# Question: the Phase-1 gold-standard reference `cens_mi_y` uses a Gaussian
# (tobit-flavoured) conditional model. The main run used skew = 0 (Gaussian
# log-exposures), where that is correctly specified. This sweep adds right-skew
# to the log-exposures to size how much the tobit reference degrades — and
# whether the copula-based `leftcens_prestep` (skew-aware margin) holds up.
#
# Additive DGP only (the case where cens_mi_y was the valid oracle at skew=0, so
# any degradation is attributable to the margin, not the §7.7 mixture confound).
#
# Env: N_REP (default 150), M (20), N (800), NCORES (14), SKEWS ("0,0.75").
# Writes results/skew_sweep.{rds,csv}.
# =============================================================================
suppressWarnings(suppressMessages(library(leftcens)))
.here <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "."
})
for (f in c("dgp.R", "censoring.R", "procedures.R", "metrics.R")) source(file.path(.here, "R", f))

getenv_i <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.integer(v) else d }
n_rep  <- getenv_i("N_REP", 150L); m <- getenv_i("M", 20L); n <- getenv_i("N", 800L)
ncores <- getenv_i("NCORES", 14L)
skews  <- { v <- Sys.getenv("SKEWS"); if (nzchar(v)) as.numeric(strsplit(v, ",")[[1]]) else c(0, 0.75) }
nd_fracs <- c(0.2, 0.4)
procs  <- c("oracle", "complete_case", "leftcens_prestep", "cens_mi_y", "cens_mi_y_shash")
base_seed <- 20260813L

raw <- list(); k <- 0L
for (sk in skews) for (nd in nd_fracs) {
  truth <- make_truth(p = 3L, erf_form = "additive")
  cat(sprintf("[skew=%.2f nd=%.1f] n_rep=%d ", sk, nd, n_rep))
  for (r in seq_len(n_rep)) {
    set.seed(base_seed + round(sk * 1000) * 1000000L + as.integer(nd * 10) * 100000L + r)
    comp <- simulate_complete(n, truth, rho = 0.4, skew = sk)
    cens <- inject_left_censoring(comp$data, nd_frac = nd, censor_which = 1L)
    bundle <- list(complete = comp$data, censored = cens, truth = truth)
    res <- run_procedures(bundle, which = procs, m = m, seed = base_seed + r, n_cores = ncores)
    res$rep <- r; res$erf_form <- "additive"; res$nd_frac <- nd; res$skew <- sk
    res$estimand_true <- truth$estimand_true
    k <- k + 1L; raw[[k]] <- res
    if (r %% max(1L, n_rep %/% 10L) == 0L) cat(".")
  }
  cat("\n")
}
raw <- do.call(rbind, raw)

# Aggregate per skew × nd × procedure.
grp <- interaction(raw$skew, raw$nd_frac, raw$procedure, drop = TRUE)
summ <- do.call(rbind, lapply(split(raw, grp), function(df) {
  t <- df$estimand_true[1]; e <- df$estimate
  data.frame(skew = df$skew[1], nd_frac = df$nd_frac[1], procedure = df$procedure[1],
             n_ok = sum(is.finite(e)), mean_est = mean(e, na.rm = TRUE),
             rel_bias = mean(e - t, na.rm = TRUE) / t,
             rmse = sqrt(mean((e - t)^2, na.rm = TRUE)),
             coverage = mean(df$ci_lo <= t & t <= df$ci_hi, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
ord <- c("oracle", "complete_case", "leftcens_prestep", "cens_mi_y", "cens_mi_y_shash")
summ$procedure <- factor(summ$procedure, levels = ord)
summ <- summ[order(summ$skew, summ$nd_frac, summ$procedure), ]

cat("\n============== Phase 1 skew sweep (additive, b_logX1, true=0.40) ==============\n")
cat(paste(sprintf("skew=%.2f nd=%.1f  %-16s | n_ok=%3d est=%.3f bias=%+5.1f%% rmse=%.3f cover=%.2f",
    summ$skew, summ$nd_frac, summ$procedure, summ$n_ok, summ$mean_est,
    100 * summ$rel_bias, summ$rmse, summ$coverage), collapse = "\n"), "\n")

results_dir <- file.path(.here, "results"); dir.create(results_dir, showWarnings = FALSE)
saveRDS(list(raw = raw, summary = summ,
             meta = list(n_rep = n_rep, m = m, n = n, skews = skews, ts = Sys.time())),
        file.path(results_dir, "skew_sweep.rds"), compress = FALSE)
utils::write.csv(summ, file.path(results_dir, "skew_sweep_summary.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s/skew_sweep.{rds,csv}\n", results_dir))
