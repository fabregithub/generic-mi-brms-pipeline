#!/usr/bin/env Rscript
# =============================================================================
# Phase 1 runner -- confirm the ERF bias and test the joint-model fix
# -----------------------------------------------------------------------------
# PLAN §10, Phase 1. Runs the four "no new code" procedures across a small grid
# of ERF forms × non-detect fractions, for N_REP replications each, and writes a
# bias/coverage table. This is the acceptance test for everything downstream:
#
#   H1 -- does the independent leftcens pre-step (no Y) bias the focal estimand,
#         and does the congenial brms joint model remove that bias?
#
# Sourceable (defines run_phase1()) and runnable as a script.
#
# Config via environment variables:
#   CONFIG   quick | full          (default quick)
#   N_REP    replications per cell  (default 50 quick / 300 full)
#   M        imputations for leftcens_prestep (default 10 / 30)
#   N        sample size            (default 400 / 800)
#   PROCS    comma list of procedures (default excludes brms; add brms_joint to include)
#   ERF      comma list: additive,mixture  (default both)
#   ND       comma list of nd fractions     (default 0.2,0.4)
#   SEED     base seed               (default 20260813)
#
# brms is EXCLUDED by default because each fit compiles/samples Stan and dominates
# runtime. Add it explicitly, e.g.:
#   PROCS=oracle,complete_case,leftcens_prestep,brms_joint N_REP=20 Rscript run_phase1.R
# =============================================================================

suppressWarnings(suppressMessages({
  library(leftcens)
}))

# --- locate and source the component files -----------------------------------
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "validation/phase1"
}
for (f in c("dgp.R", "censoring.R", "procedures.R", "metrics.R")) {
  source(file.path(.here, "R", f))
}

getenv <- function(key, default) {
  v <- Sys.getenv(key, unset = NA); if (is.na(v) || !nzchar(v)) default else v
}
split_csv <- function(x) trimws(strsplit(x, ",")[[1]])

#' Run the Phase-1 study. Returns a list(raw, summary, meta).
run_phase1 <- function(config = "quick", n_rep = NULL, m = NULL, n = NULL,
                       procs = NULL, erf_forms = NULL, nd_fracs = NULL,
                       base_seed = 20260813L, p = 3L, n_cores = 1L,
                       rho = 0.4, skew = 0.0, sigma_y = 1.0,
                       brms_control = list(), verbose = TRUE) {
  full <- identical(config, "full")
  n_rep     <- n_rep     %||% (if (full) 300L else 50L)
  m         <- m         %||% (if (full) 30L else 10L)
  n         <- n         %||% (if (full) 800L else 400L)
  procs     <- procs     %||% c("oracle", "complete_case", "leftcens_prestep", "cens_mi_y_shash")
  erf_forms <- erf_forms %||% c("additive", "mixture")
  nd_fracs  <- nd_fracs  %||% c(0.2, 0.4)

  # brms fits are SERIAL and one-per-rep; warn before a long unattended run.
  if ("brms_joint" %in% procs && "additive" %in% erf_forms && n_rep > 60L) {
    n_fits <- n_rep * sum(erf_forms == "additive")
    warning(sprintf(
      paste0("brms_joint with n_rep=%d means ~%d serial Stan fits (mixture cells ",
             "skip brms). This can take hours. Consider a separate run with ",
             "N_REP<=50 for the joint-model reference."),
      n_rep, n_fits), call. = FALSE, immediate. = TRUE)
  }

  grid <- expand.grid(erf_form = erf_forms, nd_frac = nd_fracs,
                      stringsAsFactors = FALSE)
  raw <- list(); k <- 0L
  ckpt_dir <- getOption("phase1.results_dir", NULL)

  # One compiled Stan model reused across the whole grid (the additive joint
  # model is structurally identical across cells; only the data differs, which
  # update(recompile = FALSE) handles). Mixture cells skip brms entirely.
  brms_cache <- new.env()

  for (gi in seq_len(nrow(grid))) {
    erf_form <- grid$erf_form[gi]; nd_frac <- grid$nd_frac[gi]
    truth <- make_truth(p = p, erf_form = erf_form)

    if (verbose) cat(sprintf("[cell %d/%d] erf=%s nd=%.2f  n_rep=%d\n",
                             gi, nrow(grid), erf_form, nd_frac, n_rep))

    for (r in seq_len(n_rep)) {
      set.seed(base_seed + gi * 100000L + r)
      comp <- simulate_complete(n, truth, rho = rho, skew = skew, sigma_y = sigma_y)
      cens <- inject_left_censoring(comp$data, nd_frac = nd_frac, censor_which = 1L)
      bundle <- list(complete = comp$data, censored = cens, truth = truth)

      res <- run_procedures(bundle, which = procs, m = m,
                            seed = base_seed + r, n_cores = n_cores,
                            brms_control = brms_control, brms_cache = brms_cache)
      res$rep <- r; res$erf_form <- erf_form; res$nd_frac <- nd_frac
      res$estimand_true <- truth$estimand_true
      k <- k + 1L; raw[[k]] <- res
      if (verbose && r %% max(1L, n_rep %/% 10L) == 0L) cat(".")
    }
    if (verbose) cat("\n")

    # Checkpoint after each cell so a later failure never loses finished work.
    if (!is.null(ckpt_dir)) {
      tryCatch(
        saveRDS(do.call(rbind, raw), file.path(ckpt_dir, "checkpoint_raw.rds"),
                compress = FALSE),
        error = function(e) if (verbose) cat("  (checkpoint save failed)\n")
      )
    }
  }

  raw <- do.call(rbind, raw)
  summary <- summarise_phase1(raw)
  meta <- list(config = config, n_rep = n_rep, m = m, n = n, p = p,
               n_cores = n_cores, procs = procs, erf_forms = erf_forms,
               nd_fracs = nd_fracs, rho = rho, skew = skew, sigma_y = sigma_y,
               base_seed = base_seed, timestamp = Sys.time())
  list(raw = raw, summary = summary, meta = meta)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- script entry point ------------------------------------------------------
# sys.nframe() == 0 only when run as a script (Rscript); when this file is
# source()'d to reuse the functions it is > 0, so sourcing does NOT trigger a run.
if (sys.nframe() == 0L) {
  config <- getenv("CONFIG", "quick")
  procs  <- if (!is.na(Sys.getenv("PROCS", unset = NA)) && nzchar(Sys.getenv("PROCS")))
    split_csv(Sys.getenv("PROCS")) else NULL
  erf    <- if (nzchar(Sys.getenv("ERF"))) split_csv(Sys.getenv("ERF")) else NULL
  nd     <- if (nzchar(Sys.getenv("ND")))  as.numeric(split_csv(Sys.getenv("ND"))) else NULL
  n_rep  <- if (nzchar(Sys.getenv("N_REP"))) as.integer(Sys.getenv("N_REP")) else NULL
  m      <- if (nzchar(Sys.getenv("M")))     as.integer(Sys.getenv("M")) else NULL
  n      <- if (nzchar(Sys.getenv("N")))     as.integer(Sys.getenv("N")) else NULL
  ncores <- if (nzchar(Sys.getenv("NCORES"))) as.integer(Sys.getenv("NCORES")) else 1L
  seed   <- if (nzchar(Sys.getenv("SEED")))  as.integer(Sys.getenv("SEED")) else 20260813L

  results_dir <- file.path(.here, "results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  options(phase1.results_dir = results_dir)        # enables per-cell checkpointing

  res <- run_phase1(config = config, n_rep = n_rep, m = m, n = n, n_cores = ncores,
                    procs = procs, erf_forms = erf, nd_fracs = nd, base_seed = seed)

  cat("\n================ Phase 1 summary (focal estimand b_logX1) ================\n")
  print_phase1(res$summary)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  saveRDS(res, file.path(results_dir, "latest.rds"), compress = FALSE)
  saveRDS(res, file.path(results_dir, sprintf("phase1_%s.rds", stamp)), compress = FALSE)
  utils::write.csv(res$summary, file.path(results_dir, "phase1_summary.csv"), row.names = FALSE)
  utils::write.csv(res$raw, file.path(results_dir, "phase1_raw.csv"), row.names = FALSE)
  cat(sprintf("\nWrote results to %s (latest.rds, phase1_summary.csv, phase1_raw.csv)\n",
              results_dir))
}
