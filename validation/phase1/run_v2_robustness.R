#!/usr/bin/env Rscript
# =============================================================================
# Track V2 runner -- robustness sweep for the shipped censored-exposure engine
# -----------------------------------------------------------------------------
# Track V2 of ../PLAN_pipeline_validation.md §7.
#
# V1 validated one corner of the design space. This sweeps the axes that corner
# left untouched -- crucially the two that were *dead code* in V1: MCAR
# covariates (which makes the miceRanger Z block real, so the block-FCS
# alternation is finally exercised) and a missing outcome (which fires MID).
#
# PARALLELISM. run_phase1.R parallelises INSIDE procedures and runs replications
# serially, which leaves most of a 24-core box idle whenever the inner work is
# small. This runner inverts that: replications are the unit of parallelism
# (mclapply over scenario x rep), and every procedure is told n_cores = 1.
# Rationale:
#   * There are thousands of independent tasks, so scaling is near-linear.
#   * Task cost varies ~100x across scenarios (large_n vs base), so
#     mc.preschedule = FALSE load-balances instead of chunking naively.
#   * The pipeline is sourced ONCE in the parent; forks inherit it copy-on-write
#     rather than each re-sourcing ~57 KB of R.
#   * Forking inside a fork is what must be avoided -- hence n_cores = 1 in the
#     procedures. Nested forks are the failure mode, not forks as such.
#
# RNG. Each task seeds explicitly from (scenario, rep), so results are
# reproducible regardless of how the scheduler interleaves tasks -- which
# relying on parallel RNG streams alone would not guarantee. The scenario part of
# the seed is its position in the CANONICAL (unfiltered) grid, not in whatever
# subset SCENARIOS selected: otherwise `SCENARIOS=combined` would silently
# generate different data than the same scenario in a full run, and the two could
# not be compared. That matters directly -- e.g. re-running one scenario at a
# larger M to test whether under-coverage is a small-M artefact requires the data
# to be identical across the two arms.
#
# Sourceable (defines run_v2()) and runnable as a script.
#
# Config via environment variables:
#   NCORES     fork workers                     (default: cores - 2)
#   N_REP      replications per scenario         (default 300)
#   M          imputations                       (default 30)
#   SEED       base seed                         (default 20260825)
#   SCENARIOS  comma list to subset              (default: all)
#   BIG_N      n for the large_n scenario        (default 20000)
#   BIG_N_REP  reps for the large_n scenario     (default 50)
#   SWEEPS     block-FCS outer sweeps            (default 3)
#   MARGIN     shash | gaussian                  (default shash)
#   PIPELINE_ROOT  where 00_censored_exposure.R lives
# =============================================================================

suppressWarnings(suppressMessages({
  library(leftcens)
  library(parallel)
}))

# --- locate and source the component files -----------------------------------
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "validation/phase1"
}
for (f in c("dgp.R", "censoring.R", "procedures.R", "procedures_pipeline.R",
            "metrics.R", "robustness.R")) {
  source(file.path(.here, "R", f))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

getenv <- function(key, default) {
  v <- Sys.getenv(key, unset = NA); if (is.na(v) || !nzchar(v)) default else v
}
split_csv <- function(x) trimws(strsplit(x, ",")[[1]])

# ---- summarising -------------------------------------------------------------

#' Summarise raw V2 results by scenario x procedure.
#'
#' `summarise_phase1()` groups by erf_form x nd_frac, which is the wrong key
#' here: V2 varies six axes at once, so the scenario name is the grouping unit.
summarise_v2 <- function(raw) {
  grp <- interaction(raw$scenario, raw$procedure, drop = TRUE)
  parts <- split(raw, grp)

  rows <- lapply(parts, function(df) {
    truth <- df$estimand_true[1]
    est <- df$estimate
    covered <- df$ci_lo <= truth & truth <= df$ci_hi
    data.frame(
      scenario   = df$scenario[1],
      axis       = df$axis[1],
      procedure  = df$procedure[1],
      n_rep      = nrow(df),
      n_ok       = sum(is.finite(est)),
      true       = truth,
      mean_est   = mean(est, na.rm = TRUE),
      bias       = mean(est - truth, na.rm = TRUE),
      rel_bias   = mean(est - truth, na.rm = TRUE) / truth,
      rmse       = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      emp_se     = stats::sd(est, na.rm = TRUE),
      mean_ci_w  = mean(df$ci_hi - df$ci_lo, na.rm = TRUE),
      coverage   = mean(covered, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  proc_order <- c("oracle", "complete_case", "leftcens_prestep",
                  "cens_mi_y_shash", "pipeline_block_fcs")
  out$procedure <- factor(out$procedure, levels = proc_order)
  out <- out[order(out$scenario, out$procedure), ]
  rownames(out) <- NULL
  out
}

#' Console print, rounded, one block per scenario.
print_v2 <- function(s) {
  for (sc in unique(s$scenario)) {
    x <- s[s$scenario == sc, ]
    cat(sprintf("\n--- %s  (axis: %s, true = %.2f) ---\n", sc, x$axis[1], x$true[1]))
    cat(sprintf("  %-20s %5s %9s %8s %8s %8s %7s\n",
                "procedure", "n_ok", "est", "rel%", "rmse", "ci_w", "cov"))
    for (i in seq_len(nrow(x))) {
      cat(sprintf("  %-20s %5d %9.4f %8.2f %8.4f %8.4f %7.3f\n",
                  as.character(x$procedure[i]), x$n_ok[i], x$mean_est[i],
                  100 * x$rel_bias[i], x$rmse[i], x$mean_ci_w[i], x$coverage[i]))
    }
  }
  invisible(s)
}

# ---- one task ----------------------------------------------------------------

#' Run every applicable procedure on one (scenario, rep) task.
#'
#' Returns a data.frame of procedure rows, or NULL if the task died outright.
#' A task failure must never kill the sweep -- it is recorded and skipped.
.v2_run_task <- function(task, m, base_seed, sweeps, margin, project_root) {
  sc <- task$sc
  r  <- task$rep

  # Deterministic per-task seed, keyed on the CANONICAL scenario index so that
  # subsetting via SCENARIOS does not change the data a scenario sees.
  set.seed(base_seed + task$canon_index * 100000L + r)

  t0 <- Sys.time()
  res <- tryCatch({
    truth  <- make_truth(p = 3L, erf_form = sc$erf_form)
    bundle <- v2_make_bundle(sc, truth)

    out <- run_procedures(
      bundle,
      which        = v2_procedures_for(sc),
      m            = m,
      seed         = base_seed + r,
      n_cores      = 1L,                       # rep-level forking owns the cores
      ce_control   = list(outer_sweeps = sweeps, margin = margin,
                          project_root = project_root)
    )
    out$estimand_true <- truth$estimand_true
    out
  }, error = function(e) {
    data.frame(procedure = "TASK_ERROR", estimate = NA_real_, se = NA_real_,
               ci_lo = NA_real_, ci_hi = NA_real_,
               note = substr(conditionMessage(e), 1, 120),
               estimand_true = NA_real_, stringsAsFactors = FALSE)
  })

  res$secs     <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res$rep      <- r
  res$scenario <- sc$name
  res$axis     <- sc$axis
  res$n        <- sc$n
  res$nd_frac  <- sc$nd_frac
  res$rho      <- sc$rho
  res$skew     <- sc$skew
  res$mcar     <- sc$mcar_frac
  res$y_frac   <- sc$y_frac
  res$cens_all <- sc$censor_all
  res
}

# ---- runner ------------------------------------------------------------------

#' Run the V2 robustness sweep.
#'
#' @return list(raw, summary, meta)
run_v2 <- function(n_rep = 300L, m = 30L, base_seed = 20260825L,
                   n_cores = NULL, scenarios = NULL, which_scenarios = NULL,
                   sweeps = 3L, margin = "shash", project_root = NULL,
                   big_n = 20000L, big_n_rep = 50L, verbose = TRUE) {

  n_cores <- n_cores %||% max(1L, parallel::detectCores() - 2L)
  if (.Platform$OS.type == "windows") n_cores <- 1L   # mclapply forks: unix only

  scs <- scenarios %||% v2_scenarios(big_n = big_n, big_n_rep = big_n_rep)

  # Canonical ordering, captured BEFORE any filtering -- this is what seeds key on.
  canon_names <- vapply(scs, `[[`, "", "name")

  if (!is.null(which_scenarios)) {
    keep <- vapply(scs, function(s) s$name %in% which_scenarios, logical(1))
    if (!any(keep)) stop("No scenarios matched: ",
                         paste(which_scenarios, collapse = ", "), call. = FALSE)
    scs <- scs[keep]
  }

  # Source the pipeline ONCE here, in the parent. Forked workers inherit it
  # copy-on-write instead of re-sourcing it thousands of times.
  if (verbose) cat("Loading pipeline once in parent ...\n")
  invisible(v1_load_pipeline(project_root, quiet = TRUE))

  # Flatten to a single task list so one mclapply call load-balances across ALL
  # scenarios -- otherwise a cheap scenario finishes and leaves cores idle while
  # an expensive one is still chunked across too few workers.
  tasks <- list(); k <- 0L
  for (si in seq_along(scs)) {
    sc <- scs[[si]]
    reps <- sc$n_rep %||% n_rep
    for (r in seq_len(reps)) {
      k <- k + 1L
      tasks[[k]] <- list(sc = sc, sc_index = si, rep = r,
                         canon_index = match(sc$name, canon_names))
    }
  }

  if (verbose) {
    cat(sprintf("Scenarios: %s\n", paste(vapply(scs, `[[`, "", "name"), collapse = ", ")))
    cat(sprintf("Tasks: %d  |  workers: %d  |  m = %d  |  sweeps = %d\n",
                length(tasks), n_cores, m, sweeps))
  }

  t0 <- Sys.time()
  out <- parallel::mclapply(
    tasks, .v2_run_task,
    m = m, base_seed = base_seed, sweeps = sweeps, margin = margin,
    project_root = project_root,
    mc.cores = n_cores,
    mc.preschedule = FALSE      # load-balance: task costs vary ~100x
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # A worker that dies returns a try-error rather than a data.frame; drop those
  # loudly instead of letting rbind fail or silently shrink the grid.
  bad <- vapply(out, function(z) !is.data.frame(z), logical(1))
  if (any(bad)) {
    warning(sprintf("%d of %d tasks failed at the worker level and were dropped.",
                    sum(bad), length(bad)), call. = FALSE, immediate. = TRUE)
    out <- out[!bad]
  }

  raw <- do.call(rbind, out)

  n_task_err <- sum(raw$procedure == "TASK_ERROR")
  if (n_task_err > 0) {
    warning(sprintf("%d task(s) errored inside the harness; see note column.",
                    n_task_err), call. = FALSE, immediate. = TRUE)
  }
  raw <- raw[raw$procedure != "TASK_ERROR", ]

  summary <- summarise_v2(raw)

  # Per-scenario wall-clock, so the cost of a bigger run is predictable rather
  # than guessed. One row per task carries the same `secs`, hence the [1].
  per_task <- raw[!duplicated(raw[, c("scenario", "rep")]), c("scenario", "secs")]
  cost <- aggregate(secs ~ scenario, per_task, function(z)
    c(n = length(z), mean = mean(z), total = sum(z)))
  cost <- do.call(data.frame, cost)
  names(cost) <- c("scenario", "n_tasks", "mean_secs", "total_secs")
  cost <- cost[order(-cost$mean_secs), ]
  rownames(cost) <- NULL
  meta <- list(n_rep = n_rep, m = m, base_seed = base_seed, n_cores = n_cores,
               sweeps = sweeps, margin = margin, big_n = big_n,
               big_n_rep = big_n_rep, n_tasks = length(tasks),
               n_worker_failures = sum(bad), n_task_errors = n_task_err,
               elapsed_sec = elapsed,
               scenarios = vapply(scs, `[[`, "", "name"),
               cost = cost,
               timestamp = Sys.time())

  if (verbose) {
    cat(sprintf("\nDone in %.1f min (%.2f s/task effective across %d workers)\n",
                elapsed / 60, elapsed / max(1L, length(tasks)), n_cores))
    cat("\nPer-scenario cost (single-worker seconds per replication):\n")
    print(data.frame(scenario = cost$scenario, n = cost$n_tasks,
                     mean_secs = round(cost$mean_secs, 1)), row.names = FALSE)
  }

  list(raw = raw, summary = summary, meta = meta)
}

# --- script entry point ------------------------------------------------------
# sys.nframe() == 0 only when run as a script; sourcing does NOT trigger a run.
if (sys.nframe() == 0L) {
  n_rep     <- as.integer(getenv("N_REP", 300L))
  m         <- as.integer(getenv("M", 30L))
  seed      <- as.integer(getenv("SEED", 20260825L))
  ncores    <- as.integer(getenv("NCORES", max(1L, parallel::detectCores() - 2L)))
  sweeps    <- as.integer(getenv("SWEEPS", 3L))
  margin    <- getenv("MARGIN", "shash")
  big_n     <- as.integer(getenv("BIG_N", 20000L))
  big_n_rep <- as.integer(getenv("BIG_N_REP", 50L))
  which_sc  <- if (nzchar(Sys.getenv("SCENARIOS"))) split_csv(Sys.getenv("SCENARIOS")) else NULL
  proot     <- if (nzchar(Sys.getenv("PIPELINE_ROOT"))) Sys.getenv("PIPELINE_ROOT") else NULL

  results_dir <- file.path(.here, "results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  res <- run_v2(n_rep = n_rep, m = m, base_seed = seed, n_cores = ncores,
                which_scenarios = which_sc, sweeps = sweeps, margin = margin,
                project_root = proot, big_n = big_n, big_n_rep = big_n_rep)

  cat("\n================ V2 robustness summary (estimand b_logX1) ================\n")
  print_v2(res$summary)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  saveRDS(res, file.path(results_dir, "v2_latest.rds"), compress = FALSE)
  saveRDS(res, file.path(results_dir, sprintf("v2_%s.rds", stamp)), compress = FALSE)
  utils::write.csv(res$summary, file.path(results_dir, "v2_summary.csv"), row.names = FALSE)
  utils::write.csv(res$raw, file.path(results_dir, "v2_raw.csv"), row.names = FALSE)
  cat(sprintf("\nWrote results to %s (v2_latest.rds, v2_summary.csv, v2_raw.csv)\n",
              results_dir))
}
