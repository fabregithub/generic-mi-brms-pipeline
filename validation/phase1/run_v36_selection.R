#!/usr/bin/env Rscript
# =============================================================================
# Track V36 -- does censored-exposure imputation AMPLIFY outcome-dependent
#              selection bias, or are the two merely additive?
# -----------------------------------------------------------------------------
# Design, registered predictions and the target-population argument:
#   DESIGN_v36_selection.md   (written 2026-11-09, BEFORE this driver ran)
#
# Its own driver rather than a scenario added to run_v4_variance.R, because that
# file seeds on a CANONICAL SCENARIO INDEX and inserting cells there would move
# the RNG stream of every earlier V4-family track.
#
# THE QUANTITY. Not selection bias (textbook) and not imputation bias (V17/V18),
# but their INTERACTION:
#     interaction = (arm_sel - oracle_sel) - (arm_none - oracle_none)
# positive = imputation makes selection bias worse.
#
# ENV:
#   N_REP   reps per cell            (default 100)
#   N_OBS   rows per rep             (default 800)
#   M_IMP   imputations              (default 30)
#   ND_FRAC non-detect fraction      (default 0.30)
#   CELLS   comma list               (default all seven)
#   ARMS    comma list               (default all five)
#   NCORES  fork workers             (default cores - 2)
#   OUT_TAG result prefix            (default "v36")
#   RESUME  1 = skip checkpointed    (default off)
# =============================================================================
suppressWarnings(suppressPackageStartupMessages({
  source("R/dgp.R"); source("R/censoring.R"); source("R/procedures.R")
  source("R/procedures_pipeline.R"); source("R/robustness.R")
  source("R/smc_impute.R"); source("R/bart_impute.R"); source("R/mice_impute.R")
  source("R/proper_impute.R"); source("R/metrics.R"); source("R/selection.R")
}))
`%||%` <- function(a, b) if (is.null(a)) b else a
getenv <- function(k, d) { v <- Sys.getenv(k, unset = NA)
  if (is.na(v) || !nzchar(v)) d else v }
split_csv <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])

V36_CELLS <- c("none", "descnull", "desc", "desc2",
               "outcome", "outcome2", "collider")
V36_ARMS  <- c("complete_case", "sub_lod2", "pipeline_block_fcs", "smc_xgrid")
# `oracle` is added automatically and is the per-cell reference.

.v36_run_task <- function(task, n_obs, nd_frac, m, base_seed, sweeps, margin,
                          project_root, arms) {
  sc <- task$sc; r <- task$rep
  set.seed(base_seed + task$cell_index * 100000L + r)
  t0 <- Sys.time()
  res <- tryCatch({
    truth  <- make_truth(p = 3L, erf_form = "additive")
    bundle <- v36_make_bundle(
      list(n = n_obs, nd_frac = nd_frac, sel = sc$sel, y_frac = 0.0), truth)
    std <- setdiff(arms, "sub_lod2")
    out <- run_procedures(
      bundle, which = c("oracle", std), m = m, seed = base_seed + r,
      n_cores = 1L,
      ce_control = list(outer_sweeps = sweeps, margin = margin,
                        project_root = project_root))
    if ("sub_lod2" %in% arms) out <- rbind(out, proc_sub_lod2_lm(bundle))
    out$estimand_true <- truth$estimand_true
    out$n_kept        <- bundle$n_kept
    out$nd_realised   <- bundle$nd_realised
    out
  }, error = function(e) {
    data.frame(procedure = "TASK_ERROR", estimate = NA_real_, se = NA_real_,
               ci_lo = NA_real_, ci_hi = NA_real_,
               note = substr(conditionMessage(e), 1, 160),
               ubar = NA_real_, b = NA_real_, fmi = NA_real_,
               pred_b = NA_real_, borrow = NA_real_, atten = NA_real_,
               d_rms = NA_real_, estimand_true = NA_real_,
               n_kept = NA_real_, nd_realised = NA_real_,
               stringsAsFactors = FALSE)
  })
  res$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res$rep  <- r
  res$cell <- sc$name
  res$scenario <- sc$name     # alias, so the shared RESUME logic applies
  res
}

v36_summarise <- function(raw) {
  raw <- raw[raw$procedure != "TASK_ERROR" & is.finite(raw$estimate), ]
  grp <- interaction(raw$cell, raw$procedure, drop = TRUE)
  out <- do.call(rbind, lapply(split(raw, grp), function(d) {
    tr <- d$estimand_true[1]
    cov <- mean(d$ci_lo <= tr & d$ci_hi >= tr, na.rm = TRUE)
    data.frame(cell = d$cell[1], procedure = d$procedure[1], n_ok = nrow(d),
               true = tr, mean_est = mean(d$estimate),
               bias = mean(d$estimate) - tr,
               rel_bias = (mean(d$estimate) - tr) / tr,
               mcse_rel = (stats::sd(d$estimate) / sqrt(nrow(d))) / abs(tr),
               emp_se = stats::sd(d$estimate),
               coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
               secs = mean(d$secs), n_kept = mean(d$n_kept),
               nd_realised = mean(d$nd_realised),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out[order(match(out$cell, V36_CELLS),
            match(out$procedure, c("oracle", V36_ARMS))), ]
}

run_v36 <- function(n_rep = 100L, n_obs = 800L, nd_frac = 0.30, m = 30L,
                    base_seed = 20261109L, n_cores = NULL,
                    which_cells = V36_CELLS, arms = V36_ARMS,
                    sweeps = 3L, margin = "shash", project_root = NULL,
                    verbose = TRUE) {
  n_cores <- n_cores %||% max(1L, parallel::detectCores() - 2L)
  bad <- setdiff(arms, V36_ARMS)
  if (length(bad)) stop("Unknown arm(s): ", paste(bad, collapse = ", "))
  bad <- setdiff(which_cells, V36_CELLS)
  if (length(bad)) stop("Unknown cell(s): ", paste(bad, collapse = ", "))

  tasks <- list()
  for (ci in seq_along(which_cells)) {
    cl <- which_cells[ci]
    for (r in seq_len(n_rep))
      tasks[[length(tasks) + 1L]] <- list(
        sc = list(name = cl, sel = cl), rep = r,
        cell_index = match(cl, V36_CELLS))
  }
  if (verbose)
    cat(sprintf("V36: %d cells x %d reps x %d arms = %d tasks on %d cores\n",
                length(which_cells), n_rep, length(arms) + 1L,
                length(tasks), n_cores))

  ckpt_file <- getOption("v36.checkpoint", NULL)
  done <- list()
  if (identical(Sys.getenv("RESUME"), "1") && !is.null(ckpt_file) &&
      file.exists(ckpt_file)) {
    prev <- tryCatch(readRDS(ckpt_file), error = function(e) NULL)
    if (is.data.frame(prev) && nrow(prev)) {
      have <- stats::aggregate(list(n_arms = prev$procedure),
        by = list(cell = prev$cell, rep = prev$rep),
        FUN = function(z) length(unique(z)))
      full <- have[have$n_arms >= length(arms) + 1L, ]
      key_task <- vapply(tasks, function(tk) paste(tk$sc$name, tk$rep, sep = "\r"), "")
      key_done <- intersect(paste(full$cell, full$rep, sep = "\r"), key_task)
      keep <- !(key_task %in% key_done)
      prev <- prev[paste(prev$cell, prev$rep, sep = "\r") %in% key_done &
                   prev$procedure %in% c("oracle", arms), ]
      if (verbose)
        cat(sprintf("RESUME: %d/%d tasks already complete; %d to run.\n",
                    sum(!keep), length(tasks), sum(keep)))
      done[[1L]] <- list(prev); tasks <- tasks[keep]
      if (!length(tasks)) return(prev)
    }
  }

  chunk_size <- max(n_cores, 2L * n_cores)
  chunks <- split(seq_along(tasks), ceiling(seq_along(tasks) / chunk_size))
  t0 <- Sys.time()
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    res_chunk <- parallel::mclapply(
      tasks[idx], .v36_run_task,
      n_obs = n_obs, nd_frac = nd_frac, m = m, base_seed = base_seed,
      sweeps = sweeps, margin = margin, project_root = project_root,
      arms = arms, mc.cores = n_cores, mc.preschedule = FALSE)
    done[[length(done) + 1L]] <- res_chunk
    if (!is.null(ckpt_file)) {
      dfs <- Filter(is.data.frame, unlist(done, recursive = FALSE))
      if (length(dfs)) saveRDS(do.call(rbind, dfs), ckpt_file)
    }
    if (verbose) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("  chunk %d/%d done (%d/%d tasks, %.1f min elapsed)\n",
                  ci, length(chunks), max(idx), length(tasks), el))
      utils::flush.console()
    }
  }
  flat <- unlist(done, recursive = FALSE)
  dfs <- Filter(is.data.frame, flat)
  if (length(flat) - length(dfs) > 0)
    cat(sprintf("WARNING: %d worker(s) failed\n", length(flat) - length(dfs)))
  raw <- do.call(rbind, dfs)
  if (verbose)
    cat(sprintf("\nDone in %.1f min\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  raw
}

if (sys.nframe() == 0L) {
  tag  <- getenv("OUT_TAG", "v36")
  nrep <- as.integer(getenv("N_REP", 100L))
  nobs <- as.integer(getenv("N_OBS", 800L))
  mimp <- as.integer(getenv("M_IMP", 30L))
  ndf  <- as.numeric(getenv("ND_FRAC", 0.30))
  cls  <- if (nzchar(Sys.getenv("CELLS"))) split_csv(Sys.getenv("CELLS")) else V36_CELLS
  arms <- if (nzchar(Sys.getenv("ARMS")))  split_csv(Sys.getenv("ARMS"))  else V36_ARMS
  ncor <- as.integer(getenv("NCORES", max(1L, parallel::detectCores() - 2L)))
  results_dir <- "results"
  dir.create(results_dir, showWarnings = FALSE)
  options(v36.checkpoint = file.path(results_dir, paste0(tag, "_checkpoint.rds")))

  raw <- run_v36(n_rep = nrep, n_obs = nobs, nd_frac = ndf, m = mimp,
                 n_cores = ncor, which_cells = cls, arms = arms)
  s <- v36_summarise(raw)
  print(s, row.names = FALSE, digits = 4)
  saveRDS(raw, file.path(results_dir, paste0(tag, "_latest.rds")))
  write.csv(s,   file.path(results_dir, paste0(tag, "_summary.csv")), row.names = FALSE)
  write.csv(raw, file.path(results_dir, paste0(tag, "_raw.csv")), row.names = FALSE)
  cat(sprintf("\nWrote results/%s_{latest.rds,summary.csv,raw.csv}\n", tag))
}
