#!/usr/bin/env Rscript
# =============================================================================
# Phase 2 -- realistic-profile CENSORING SWEEP of the mixture-benchmark trap
# -----------------------------------------------------------------------------
# Maps how the recovered mixture-HBEB benchmark bias grows as censoring on the
# real JECS detection profile is scaled up (`limit_scale` = 1.0 -> 3.0). Answers
# DESIGN §7: *when* is surface-congenial imputation essential? Expected story:
# at good detection (ls=1.0) all imputation recovers the benchmark; as detection
# degrades, substitution and the outcome-blind / linear-congenial arms diverge
# from the oracle -- the regime where the fix (procedure 5) will earn its keep.
#
# Reuses the validated machinery in run_trap_demo.R. The ORACLE is
# censoring-independent, so it is computed ONCE and reused across the sweep.
# Saves incrementally (results/trap_sweep_*.rds) so a long run is never lost.
#
# Config via env:
#   LIMITS  comma list of limit_scale values   (default 1.0,1.6,2.2,3.0)
#   N ITER M NKNOTS SWEEPS REPS ERF SEED REL    (see run_trap_demo.R defaults)
#
# Time: ~1 bkmr fit/min. Per limit_scale = (1 + M + M) fits x REPS; plus REPS for
#   the one-off oracle. Default (4 limits, M=5, REPS=3) ~= 135 fits ~= 2-2.5 h.
#   To go faster: fewer LIMITS, M=3, or N=600.
# =============================================================================

suppressWarnings(suppressMessages({library(MASS); library(leftcens); library(bkmr)}))

.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "R"
}
source(file.path(.here, "run_trap_demo.R"))     # brings in benchmark/dgp/censoring + procedures

getenv <- function(k, d) { v <- Sys.getenv(k, NA); if (is.na(v) || !nzchar(v)) d else v }

#' Mean benchmark (theta + all analyte concentrations) over a list of completed
#' datasets -- one arm's estimate. NA row if every fit failed.
arm_bench <- function(dfs, truth, iter, nknots, rel) {
  bs <- lapply(dfs, function(df)
    tryCatch(bench_of(df, truth, iter, nknots, rel), error = function(e) NULL))
  bs <- bs[!vapply(bs, is.null, logical(1))]
  if (!length(bs))
    return(stats::setNames(rep(NA_real_, length(truth$analytes) + 1L), c("theta", truth$analytes)))
  rowMeans(do.call(cbind, bs))
}

row_of <- function(limit_scale, rep, method, est)
  cbind(data.frame(limit_scale = limit_scale, rep = rep, method = method),
        as.data.frame(as.list(est)))

run_sweep <- function(limits, n, iter, m, nknots, sweeps, erf_form, reps, seed, rel) {
  truth <- make_truth(erf_form, "real"); an <- truth$analytes
  outdir <- file.path(.here, "..", "results"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  tag <- sprintf("trap_sweep_%s_n%d", erf_form, n)

  # ground-truth benchmark (causal surface; independent of censoring)
  tbp <- true_benchmark(truth, rel_targets = rel, n_pop = 1e5, seed = 99L)
  tbp <- tbp[tbp$direction == "prop", ]
  truth_vec <- c(theta = tbp$theta, unlist(tbp[, an]))
  cat(sprintf("\n=== CENSORING SWEEP  erf=%s n=%d iter=%d m=%d reps=%d (%d%%-reduction) ===\n",
              erf_form, n, iter, m, reps, round(rel * 100)))
  cat(sprintf("TRUTH: theta*=%.4f  PFOS=%.3f  PFTrDA=%.4f ng/mL\n", truth_vec[["theta"]],
              truth_vec[["pfos"]], truth_vec[["pftrda"]]))

  # --- ORACLE once (censoring-independent) -----------------------------------
  cat("computing oracle (once, reused across the sweep)...\n")
  rows <- list()
  for (rp in seq_len(reps)) {
    dcomplete <- simulate_complete(n, truth, seed = seed + rp)$data
    rows[[length(rows) + 1]] <- row_of(NA_real_, rp, "oracle",
                                       arm_bench(list(dcomplete), truth, iter, nknots, rel))
  }
  res <- do.call(rbind, rows)
  oracle_mean <- colMeans(res[res$method == "oracle", c("theta", an), drop = FALSE], na.rm = TRUE)

  # --- sweep over limit_scale -------------------------------------------------
  detection <- list()
  for (ls in limits) {
    dprobe <- inject_censoring(simulate_complete(20000L, truth, seed = 1L)$data, limit_scale = ls)
    det <- attr(dprobe, "detected"); detection[[sprintf("%.2f", ls)]] <- det
    cat(sprintf("\n--- limit_scale=%.2f  mean detection %.2f  (PFTrDA %.2f, PFHxS %.2f) ---\n",
                ls, mean(det), det[["pftrda"]], det[["pfhxs"]]))
    for (rp in seq_len(reps)) {
      t0 <- Sys.time()
      dcomplete <- simulate_complete(n, truth, seed = seed + rp)$data
      dcens <- inject_censoring(dcomplete, limit_scale = ls)
      arms <- list(
        substitution  = proc_substitution(dcens, truth),
        impute_blind  = block_fcs(dcens, truth, use_y = FALSE, m = m, sweeps = sweeps, seed = seed + rp),
        impute_linear = block_fcs(dcens, truth, use_y = TRUE,  m = m, sweeps = sweeps, seed = seed + rp))
      for (nm in names(arms))
        res <- rbind(res, row_of(ls, rp, nm, arm_bench(arms[[nm]], truth, iter, nknots, rel)))
      cat(sprintf("  rep %d/%d (%.0fs)\n", rp, reps, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    saveRDS(list(res = res, truth = truth_vec, oracle_mean = oracle_mean,
                 detection = detection, config = list(n = n, iter = iter, m = m,
                 reps = reps, erf_form = erf_form, rel = rel, limits = limits)),
            file.path(outdir, paste0(tag, ".rds")))                 # incremental
  }

  # --- summary ---------------------------------------------------------------
  summarise_sweep(res, truth_vec, oracle_mean, detection, an)
  utils::write.csv(res, file.path(outdir, paste0(tag, ".csv")), row.names = FALSE)
  cat(sprintf("\nsaved: results/%s.rds + .csv\n", tag))
  invisible(res)
}

#' Print the "does the trap grow with censoring" tables (bias vs oracle).
summarise_sweep <- function(res, truth_vec, oracle_mean, detection, an) {
  cens <- res[!is.na(res$limit_scale), ]
  ls_levels <- sort(unique(cens$limit_scale))
  methods <- c("substitution", "impute_blind", "impute_linear")
  meanby <- function(metric) {
    m <- matrix(NA_real_, length(ls_levels), length(methods),
                dimnames = list(sprintf("ls=%.1f", ls_levels), methods))
    for (i in seq_along(ls_levels)) for (j in seq_along(methods)) {
      v <- cens[[metric]][cens$limit_scale == ls_levels[i] & cens$method == methods[j]]
      m[i, j] <- mean(v, na.rm = TRUE)
    }
    m
  }
  det_mean <- vapply(ls_levels, function(l) mean(detection[[sprintf("%.2f", l)]]), numeric(1))
  pftrda_det <- vapply(ls_levels, function(l) detection[[sprintf("%.2f", l)]][["pftrda"]], numeric(1))

  for (metric in c("pfos", "pftrda", "theta")) {
    mm <- meanby(metric)
    ref <- oracle_mean[[metric]]
    bias <- 100 * (mm - ref) / ref
    cat(sprintf("\n### %s benchmark bias vs ORACLE (%%)   [oracle %s = %.4f] ###\n",
                toupper(metric), metric, ref))
    tab <- cbind(`det%` = round(100 * det_mean, 1), `PFTrDA det%` = round(100 * pftrda_det, 1),
                 round(bias, 1))
    print(tab)
  }
  cat("\n(Read down each column: does |bias| grow as detection falls? ",
      "Compare impute_linear vs impute_blind vs substitution.)\n")
}

if (sys.nframe() == 0) {
  lim <- as.numeric(strsplit(getenv("LIMITS", "1.0,1.6,2.2,3.0"), ",")[[1]])
  run_sweep(limits = lim,
            n = as.integer(getenv("N", 800)), iter = as.integer(getenv("ITER", 2500)),
            m = as.integer(getenv("M", 5)), nknots = as.integer(getenv("NKNOTS", 50)),
            sweeps = as.integer(getenv("SWEEPS", 3)), erf_form = getenv("ERF", "mixture"),
            reps = as.integer(getenv("REPS", 3)), seed = as.integer(getenv("SEED", 1)),
            rel = as.numeric(getenv("REL", 0.25)))
}
