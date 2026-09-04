#!/usr/bin/env Rscript
# =============================================================================
# Track V4 runner -- attributing the V2 interval under-coverage
# -----------------------------------------------------------------------------
# Track V4 of ../PLAN_pipeline_validation.md §9.
#
# V2 established: the shipped engine's point estimates are sound everywhere
# (|bias| <= 2.71%) but its 95% intervals are ~17% too narrow when covariate AND
# outcome missingness are both substantial (coverage 0.893). A paired m=30 vs
# m=100 test on identical data showed raising `m` changes nothing, so this is a
# systematic shortfall in the variance estimator, not Monte-Carlo noise.
#
# Two mechanisms remain in play, and the V2 scenario pattern rules out either
# acting alone:
#
#   (a) The miceRanger Z block is IMPROPER MI -- random forest + predictive mean
#       matching does not draw imputation-model parameters from a posterior,
#       which understates between-imputation variance. Supports: mcar_z40 has a
#       complete outcome (no MID) yet width/SE = 0.90.
#   (b) Naive Rubin pooling on MULTIPLE-IMPUTATION-THEN-DELETION data understates
#       variance (von Hippel 2007 derives a corrected estimator; this pipeline
#       uses the naive one). Supports: missing_y20 has complete covariates (Z
#       block idle) yet width/SE = 0.96.
#
# THE DESIGN. Four arms per scenario, each differing from the shipped engine in
# exactly ONE respect, so a shift in width/SE attributes the shortfall:
#
#   pipeline_block_fcs        as shipped                      (reference)
#   pipeline_properZ          Z block -> proper Bayesian draw (isolates (a))
#   pipeline_noMID            MID off, imputed-Y rows kept    (isolates (b))
#   pipeline_properZ_noMID    both                            (tests additivity)
#
# Plus the variance DECOMPOSITION, now returned by rubin_pool(): `ubar`
# (within-imputation), `b` (between-imputation) and `fmi`. A shortfall in `b`
# points at the imputation being improper; a shortfall in `ubar` points somewhere
# else entirely. The pooled SE alone -- all V2 stored -- cannot distinguish them.
#
# Scenarios: the three diagnostic cells plus `base` as a negative control (the Z
# block is idle there, so all four arms should agree and sit at width/SE ~1.0).
#
# Config via environment variables:
#   NCORES     fork workers                    (default: cores - 2)
#   N_REP      replications per scenario        (default 150)
#   M          imputations                      (default 30)
#   SEED       base seed                        (default 20260825, matching V2)
#   SCENARIOS  comma list                       (default base,mcar_z40,missing_y20,combined)
#   ARMS       comma list of arms               (default: all four)
#   SWEEPS / MARGIN / PIPELINE_ROOT             as V2
#   BART_NTREE number of trees for the BART arm  (default 50)
#   BART_NSKIP BART burn-in iterations           (default 100)
#   OUT_TAG    prefix for result/checkpoint files, so several runs can be chained
#              in one session without overwriting each other (default "v4")
#
# SEED note: the default matches V2 and the canonical scenario indexing is shared,
# so arm `pipeline_block_fcs` here sees byte-identical data to V2's corresponding
# scenario. That makes this run directly comparable to FINDINGS_v2.md rather than
# an independent sample.
# =============================================================================

suppressWarnings(suppressMessages({
  library(leftcens)
  library(parallel)
}))

.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "validation/phase1"
}
for (f in c("dgp.R", "censoring.R", "procedures.R", "procedures_pipeline.R",
            "metrics.R", "robustness.R", "proper_impute.R", "mice_impute.R",
            "bart_impute.R", "smc_impute.R")) {
  source(file.path(.here, "R", f))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
getenv <- function(key, default) {
  v <- Sys.getenv(key, unset = NA); if (is.na(v) || !nzchar(v)) default else v
}
split_csv <- function(x) trimws(strsplit(x, ",")[[1]])

V4_ARMS <- c("pipeline_block_fcs", "pipeline_properBoot", "pipeline_micePmm",
             "pipeline_bartMI", "pipeline_bartMI_iter3", "pipeline_bartHarness",
             "pipeline_properZ",
             "pipeline_noMID", "pipeline_properZ_noMID",
             # V12: Z-block draw shape (harness instruments, not pipeline code)
             "smc_zexact", "smc_zgauss", "smc_zexact_xship", "smc_zgauss_xship",
             "smc_xgrid",
             # V16: the root claim -- Y dropped from one imputation block at a
             # time, against pipeline_bartMI as the both-blocks reference.
             "pipeline_noYx", "pipeline_noYz", "pipeline_noYboth",
             # V17: auxiliary covariate, shipped vs documented behaviour.
             "pipeline_auxZ_shipped", "pipeline_auxZ_asdoc")

V4_SCENARIOS <- c("base", "mcar_z40", "missing_y20", "combined")

# ---- summarising -------------------------------------------------------------

#' Summarise by scenario x arm, with the variance decomposition.
#'
#' `width_se_ratio` is the headline diagnostic: (mean CI width / 3.92) divided by
#' the empirical SE of the estimates. 1.0 means the intervals match the true
#' sampling spread; below 1.0 they are too narrow.
#'
#' `b_share` is the fraction of total pooled variance coming from the
#' between-imputation term. If mechanism (a) holds, the improper arms should show
#' a *smaller* `b_share` than the proper ones on the same data.
summarise_v4 <- function(raw) {
  grp <- interaction(raw$scenario, raw$procedure, drop = TRUE)
  parts <- split(raw, grp)

  rows <- lapply(parts, function(df) {
    truth <- df$estimand_true[1]
    est <- df$estimate
    emp_se <- stats::sd(est, na.rm = TRUE)
    ci_w <- mean(df$ci_hi - df$ci_lo, na.rm = TRUE)
    ubar <- mean(df$ubar, na.rm = TRUE)
    b <- mean(df$b, na.rm = TRUE)
    tot <- ubar + b
    data.frame(
      scenario       = df$scenario[1],
      procedure      = df$procedure[1],
      n_rep          = nrow(df),
      n_ok           = sum(is.finite(est)),
      mean_est       = mean(est, na.rm = TRUE),
      rel_bias       = (mean(est, na.rm = TRUE) - truth) / truth,
      emp_se         = emp_se,
      mean_ci_w      = ci_w,
      width_se_ratio = (ci_w / 3.92) / emp_se,
      ubar           = ubar,
      b              = b,
      b_share        = if (is.finite(tot) && tot > 0) b / tot else NA_real_,
      fmi            = mean(df$fmi, na.rm = TRUE),
      # What the pooled SE claims, against what the sampling distribution shows.
      claimed_se     = mean(df$se, na.rm = TRUE),
      se_deficit     = 1 - mean(df$se, na.rm = TRUE) / emp_se,
      coverage       = mean(df$ci_lo <= truth & truth <= df$ci_hi, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$procedure <- factor(out$procedure, levels = c("oracle", V4_ARMS))
  out <- out[order(out$scenario, out$procedure), ]
  rownames(out) <- NULL
  out
}

print_v4 <- function(s) {
  for (sc in unique(s$scenario)) {
    x <- s[s$scenario == sc, ]
    cat(sprintf("\n--- %s ---\n", sc))
    cat(sprintf("  %-24s %5s %7s %8s %8s %8s %8s %7s\n",
                "arm", "n_ok", "rel%", "width/se", "ubar", "b", "b_share", "cov"))
    for (i in seq_len(nrow(x))) {
      cat(sprintf("  %-24s %5d %7.2f %8.3f %8.5f %8.5f %8.3f %7.3f\n",
                  as.character(x$procedure[i]), x$n_ok[i], 100 * x$rel_bias[i],
                  x$width_se_ratio[i], x$ubar[i], x$b[i], x$b_share[i],
                  x$coverage[i]))
    }
  }
  invisible(s)
}

# ---- one task ----------------------------------------------------------------

.v4_run_task <- function(task, m, base_seed, sweeps, margin, project_root, arms) {
  sc <- task$sc
  r  <- task$rep

  # Same seeding rule as V2, keyed on the canonical scenario index, so the
  # reference arm reproduces V2's data exactly.
  set.seed(base_seed + task$canon_index * 100000L + r)

  t0 <- Sys.time()
  res <- tryCatch({
    # y_form defaults to "linear", so every pre-V11 scenario is unchanged.
    truth  <- make_truth(p = 3L, erf_form = sc$erf_form,
                         y_form = sc$y_form %||% "linear",
                         z_role = sc$z_role %||% "precision")
    bundle <- v2_make_bundle(sc, truth)

    out <- run_procedures(
      bundle,
      which      = c("oracle", arms),
      m          = m,
      seed       = base_seed + r,
      n_cores    = 1L,
      ce_control = list(outer_sweeps = sweeps, margin = margin,
                        project_root = project_root)
    )
    out$estimand_true <- truth$estimand_true
    out
  }, error = function(e) {
    data.frame(procedure = "TASK_ERROR", estimate = NA_real_, se = NA_real_,
               ci_lo = NA_real_, ci_hi = NA_real_,
               note = substr(conditionMessage(e), 1, 120),
               ubar = NA_real_, b = NA_real_, fmi = NA_real_,
               estimand_true = NA_real_, stringsAsFactors = FALSE)
  })

  res$secs     <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res$rep      <- r
  res$scenario <- sc$name
  res
}

# ---- runner ------------------------------------------------------------------

run_v4 <- function(n_rep = 150L, m = 30L, base_seed = 20260825L, n_cores = NULL,
                   which_scenarios = V4_SCENARIOS, arms = V4_ARMS,
                   sweeps = 3L, margin = "shash", project_root = NULL,
                   resume = FALSE, verbose = TRUE, n_obs = NULL) {

  n_cores <- n_cores %||% max(1L, parallel::detectCores() - 2L)
  if (.Platform$OS.type == "windows") n_cores <- 1L

  bad_arms <- setdiff(arms, V4_ARMS)
  if (length(bad_arms)) {
    stop("Unknown arm(s): ", paste(bad_arms, collapse = ", "),
         ". Valid: ", paste(V4_ARMS, collapse = ", "), call. = FALSE)
  }

  scs_all <- v2_scenarios()
  canon_names <- vapply(scs_all, `[[`, "", "name")
  keep <- canon_names %in% which_scenarios
  if (!any(keep)) stop("No scenarios matched.", call. = FALSE)
  scs <- scs_all[keep]

  # V18: override every selected scenario's n. The scenarios' own `n` is 800
  # everywhere, so an n-scaling track would otherwise need a duplicate cell per
  # (cell, n) pair -- and duplicates drift. Overriding keeps ONE definition of
  # each cell and varies only n.
  #
  # Off unless asked for, so no earlier run is affected. When on it is announced
  # loudly, because results at a different n must never be mistaken for the
  # n = 800 ones they will sit next to in results/.
  if (!is.null(n_obs)) {
    for (i in seq_along(scs)) scs[[i]]$n <- as.integer(n_obs)
    if (verbose) cat(sprintf("*** N_OBS OVERRIDE: every cell runs at n = %d (default 800) ***\n",
                             as.integer(n_obs)))
  }

  # Load BOTH pipeline environments (shipped + proper-Z) once in the parent, so
  # forks inherit them copy-on-write rather than re-sourcing per task.
  if (verbose) cat("Loading pipeline environments once in parent ...\n")
  invisible(v1_load_pipeline(project_root, quiet = TRUE, proper_z = FALSE))
  if (any(grepl("properZ", arms))) {
    invisible(v1_load_pipeline(project_root, quiet = TRUE, proper_z = TRUE))
  }

  tasks <- list(); k <- 0L
  for (sc in scs) {
    for (r in seq_len(n_rep)) {
      k <- k + 1L
      # Seed on `seed_as` when a cell is deliberately paired with another (V11):
      # the pair then sees byte-identical data and differs only in the feature
      # under test, making their contrast paired rather than unpaired.
      tasks[[k]] <- list(sc = sc, rep = r,
                         canon_index = match(sc$seed_as %||% sc$name, canon_names))
    }
  }

  if (verbose) {
    cat(sprintf("Scenarios: %s\n", paste(vapply(scs, `[[`, "", "name"), collapse = ", ")))
    cat(sprintf("Arms: %s\n", paste(arms, collapse = ", ")))
    cat(sprintf("Tasks: %d  |  workers: %d  |  m = %d  |  sweeps = %d\n",
                length(tasks), n_cores, m, sweeps))
  }

  # CHECKPOINTED EXECUTION. A single mclapply over all tasks writes nothing until
  # it returns, so any interruption -- and the first V4 attempt died unexplained
  # 18 min in -- loses the whole run. Tasks are therefore processed in chunks and
  # the accumulated raw results are flushed to disk after each one. Worst case a
  # kill costs one chunk, not the run.
  ckpt_file <- getOption("v4.checkpoint", NULL)

  # RESUME IS OFF BY DEFAULT AND CURRENTLY NOT TRUSTED. Empirically, resuming
  # from an accumulated checkpoint killed the parent process right after the
  # first chunk -- reproducibly, in both detached and foreground runs -- while
  # the identical configuration with a FRESH checkpoint ran fine. The cause is
  # not yet understood (the checkpoint's structure looks correct: right columns,
  # classes and procedures), so resume must be opted into explicitly.
  #
  # This does NOT make a crash costly: the checkpoint is still written after every
  # chunk, and partial results can be read and summarised directly --
  #     summarise_v4(readRDS("results/v4_checkpoint.rds"))
  # -- which recovers everything completed without going through resume at all.
  if (isTRUE(resume)) {
    warning("resume = TRUE is not currently trusted; see the note in ",
            "run_v4_variance.R. Prefer summarising the checkpoint directly.",
            call. = FALSE, immediate. = TRUE)
  }

  done <- list(); k0 <- 0L
  if (!is.null(ckpt_file) && isTRUE(resume) && file.exists(ckpt_file)) {
    prev <- tryCatch(readRDS(ckpt_file), error = function(e) NULL)
    if (is.data.frame(prev) && nrow(prev)) {
      seen <- unique(paste(prev$scenario, prev$rep))
      keep <- !(vapply(tasks, function(t) paste(t$sc$name, t$rep), "") %in% seen)
      k0 <- sum(!keep)
      # MUST be wrapped in a list: `done` is later flattened with
      # unlist(recursive = FALSE), and unlisting a bare data.frame explodes it
      # into its columns -- which silently discarded every resumed row and then
      # reported them as "worker failures".
      done[[1]] <- list(prev)
      tasks <- tasks[keep]
      if (verbose) cat(sprintf("Resuming: %d task(s) already done, %d to go\n",
                               k0, length(tasks)))
    }
  }

  chunk_size <- max(n_cores, 2L * n_cores)
  chunks <- split(seq_along(tasks), ceiling(seq_along(tasks) / chunk_size))

  t0 <- Sys.time()
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    res_chunk <- parallel::mclapply(
      tasks[idx], .v4_run_task,
      m = m, base_seed = base_seed, sweeps = sweeps, margin = margin,
      project_root = project_root, arms = arms,
      mc.cores = n_cores, mc.preschedule = FALSE)
    done[[length(done) + 1L]] <- res_chunk

    if (!is.null(ckpt_file)) {
      flat <- unlist(done, recursive = FALSE)
      dfs <- Filter(is.data.frame, flat)
      if (length(dfs)) {
        tryCatch(saveRDS(do.call(rbind, dfs), ckpt_file, compress = FALSE),
                 error = function(e) if (verbose) cat("  (checkpoint failed)\n"))
      }
    }
    if (verbose) {
      cat(sprintf("  chunk %d/%d done (%d/%d tasks, %.1f min elapsed)\n",
                  ci, length(chunks), k0 + max(idx), k0 + length(tasks),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      utils::flush.console()
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  out <- unlist(done, recursive = FALSE)
  bad <- vapply(out, function(z) !is.data.frame(z), logical(1))
  if (any(bad)) {
    warning(sprintf("%d of %d tasks failed at the worker level and were dropped.",
                    sum(bad), length(bad)), call. = FALSE, immediate. = TRUE)
    out <- out[!bad]
  }
  raw <- do.call(rbind, out)

  n_task_err <- sum(raw$procedure == "TASK_ERROR")
  if (n_task_err > 0) {
    warning(sprintf("%d task(s) errored; see note column.", n_task_err),
            call. = FALSE, immediate. = TRUE)
  }
  raw <- raw[raw$procedure != "TASK_ERROR", ]

  summary <- summarise_v4(raw)
  meta <- list(n_rep = n_rep, m = m, base_seed = base_seed, n_cores = n_cores,
               sweeps = sweeps, margin = margin, arms = arms,
               bart_ntree = getOption("v7.bart.ntree", 50L),
               bart_nskip = getOption("v7.bart.nskip", 100L),
               scenarios = which_scenarios, n_tasks = length(tasks),
               n_worker_failures = sum(bad), n_task_errors = n_task_err,
               elapsed_sec = elapsed, timestamp = Sys.time())

  if (verbose) {
    cat(sprintf("\nDone in %.1f min\n", elapsed / 60))
  }

  list(raw = raw, summary = summary, meta = meta)
}

# --- script entry point ------------------------------------------------------
if (sys.nframe() == 0L) {
  n_rep  <- as.integer(getenv("N_REP", 150L))
  m      <- as.integer(getenv("M", 30L))
  seed   <- as.integer(getenv("SEED", 20260825L))
  ncores <- as.integer(getenv("NCORES", max(1L, parallel::detectCores() - 2L)))
  sweeps <- as.integer(getenv("SWEEPS", 3L))
  margin <- getenv("MARGIN", "shash")
  n_obs  <- if (nzchar(Sys.getenv("N_OBS"))) as.integer(Sys.getenv("N_OBS")) else NULL
  scs    <- if (nzchar(Sys.getenv("SCENARIOS"))) split_csv(Sys.getenv("SCENARIOS")) else V4_SCENARIOS
  arms   <- if (nzchar(Sys.getenv("ARMS")))      split_csv(Sys.getenv("ARMS"))      else V4_ARMS
  proot  <- if (nzchar(Sys.getenv("PIPELINE_ROOT"))) Sys.getenv("PIPELINE_ROOT") else NULL

  # BART tuning, exposed so a tuning-sensitivity run needs no file edit.
  if (nzchar(Sys.getenv("BART_NTREE")))
    options(v7.bart.ntree = as.integer(Sys.getenv("BART_NTREE")))
  if (nzchar(Sys.getenv("BART_NSKIP")))
    options(v7.bart.nskip = as.integer(Sys.getenv("BART_NSKIP")))

  tag <- getenv("OUT_TAG", "v4")

  results_dir <- file.path(.here, "results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  options(v4.checkpoint = file.path(results_dir, paste0(tag, "_checkpoint.rds")))

  res <- run_v4(n_rep = n_rep, m = m, base_seed = seed, n_cores = ncores,
                which_scenarios = scs, arms = arms, sweeps = sweeps,
                margin = margin, project_root = proot, n_obs = n_obs)

  cat("\n================ V4 variance attribution (estimand b_logX1) ================\n")
  print_v4(res$summary)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  saveRDS(res, file.path(results_dir, paste0(tag, "_latest.rds")), compress = FALSE)
  saveRDS(res, file.path(results_dir, sprintf("%s_%s.rds", tag, stamp)), compress = FALSE)
  utils::write.csv(res$summary, file.path(results_dir, paste0(tag, "_summary.csv")), row.names = FALSE)
  utils::write.csv(res$raw, file.path(results_dir, paste0(tag, "_raw.csv")), row.names = FALSE)
  cat(sprintf("\nWrote results to %s (%s_latest.rds, %s_summary.csv, %s_raw.csv)\n",
              results_dir, tag, tag, tag))
}
