#!/usr/bin/env Rscript
# =============================================================================
# Track V3 runner -- BKMR mixture estimands (roadmap 04, design Phase 2b, R11)
# -----------------------------------------------------------------------------
# THE QUESTION. Phase 1 §7.7 found that linear congenial imputation stays biased
# on a mixture surface -- +19.6% with coverage 0.62. But it measured that on a
# SCAFFOLD estimand: the `logX1` coefficient of a matched linear model, which is
# the local main effect at the origin and not something a BKMR analysis reports.
# The claims ledger accordingly records the mixture verdict as "mechanism
# demonstrated, verdict not manuscript-final".
#
# This run replaces the scaffold with the estimands people actually report --
# overall mixture effect, single-exposure effects, the X1:X2 interaction, and the
# curvature in X1 -- each with its truth derived ANALYTICALLY from the generator
# (see R/estimands_bkmr.R for why analytic and not an oracle fit).
#
# THE ARMS
#   oracle_bkmr     complete data                   -- a GATE on the harness, not the yardstick
#   cc_bkmr         complete case
#   sub_lod2_bkmr   LOD/sqrt(2) substitution        -- standard applied practice
#   pipeline_bkmr   shipped block-FCS, m imputations, BKMR each, Rubin pooled
#
# WHAT WOULD FALSIFY THE EXPECTED RESULT. The hypothesis is that `pipeline_bkmr`
# is materially biased on `int_X1X2` and `curv_X1`, because the X block's
# conditional for the censored exposure is LINEAR in (Y, other X, Z) and cannot
# represent a surface with an interaction and a quadratic. If it comes back
# unbiased on those two, the §7.7 verdict does not generalise to the reported
# estimands and the mixture restriction in the README is too strong.
#
# ACCEPTANCE (see ../PLAN_pipeline_validation.md §8 for the full text and for the
# record of why the gate was widened from 5% to 10% BEFORE the run).
#   GATE   oracle_bkmr: |rel. bias| <= 10% and coverage within [0.90, 0.98] on all
#          seven estimands. If the gate fails the run says nothing about the other
#          arms and must be diagnosed before they are read.
#   MAIN   pipeline_bkmr on `int_X1X2` and `curv_X1`: report bias and coverage
#          with Monte-Carlo error, AND the excess over `oracle_bkmr` on the same
#          estimand -- the oracle measures the estimator's own floor, so raw bias
#          alone would charge the imputation for error BKMR makes on complete
#          data. "Materially biased" is |rel. bias| > 10% or coverage < 0.90.
#   REF    sub_lod2_bkmr and cc_bkmr are reported for context; no criterion.
#   Nothing here is a pass/fail on the pipeline: V3 establishes the SCOPE of the
#   mixture restriction, it does not gate a release.
#
# WHY THE ORACLE FLOOR MATTERS. A 44-rep oracle-only sweep at iter 1000/3000/8000
# found a STABLE +6-7% bias on `curv_X1`, `overall_q75_q50` and `singvar_X1_q75`
# that does not shrink with chain length (6.61 -> 7.06 from 3000 to 8000, within
# MC error). That is BKMR's own finite-sample bias on this surface at n = 800,
# not an MCMC artifact, and it sets the floor every other arm is measured against.
# Coverage passed everywhere (0.909-0.977) at all three chain lengths.
#
# CONFIG (environment variables)
#   NCORES     fork workers                              (default: cores - 2)
#   N_REP      replications per scenario                 (default 200)
#   N_OBS      sample size                               (default 800, matching V1-V8)
#   ITER       BKMR MCMC iterations, second half kept    (default 3000)
#              MEASURED: 1000 is too short (curv_X1 bias 9.99% vs 6.61% at
#              3000); 3000 -> 8000 changes nothing (6.61 -> 7.06, within MC
#              error). 3000 is where it converges, at 2.6x the cost of 1000.
#   M          imputations for the pipeline arm          (default 10)
#   KNOTS      GPP knots; 0/unset = full GP              (default unset)
#   SEED       base seed                                 (default 20260828)
#   SCENARIOS  comma list                       (default nd20,nd40,nd40_all)
#   ARMS       comma list                                (default: all four)
#   QLO / QHI  low/high evaluation quantiles          (default 0.25 / 0.75)
#              A DIAGNOSTIC, not a tuning knob -- narrowing them tests whether
#              the oracle floor is a sparse-upper-tail artifact. Narrowing also
#              SHRINKS every truth, so compare absolute and relative bias
#              together. Changing these changes the estimands.
#   SWEEPS / MARGIN / PIPELINE_ROOT                      as V4
#   OUT_TAG    result/checkpoint prefix                  (default "v3")
#
# COST. `pipeline_bkmr` fits M BKMR models per replication and dominates
# everything else. Measure before scaling: see run_v3_bkmr.sh's header for the
# timed figures on this machine.
# =============================================================================

suppressWarnings(suppressMessages({
  library(leftcens)
  library(parallel)
  library(bkmr)
}))

.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "validation/phase1"
}
for (f in c("dgp.R", "censoring.R", "procedures.R", "procedures_pipeline.R",
            "metrics.R", "robustness.R", "proper_impute.R", "mice_impute.R",
            "bart_impute.R", "estimands_bkmr.R", "smc_impute.R",
            "procedures_bkmr.R")) {
  source(file.path(.here, "R", f))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
getenv <- function(key, default) {
  v <- Sys.getenv(key, unset = NA); if (is.na(v) || !nzchar(v)) default else v
}
split_csv <- function(x) trimws(strsplit(x, ",")[[1]])

V3_ARMS <- c("oracle_bkmr", "cc_bkmr", "sub_lod2_bkmr", "pipeline_bkmr",
             # V9: the SMC arms. Not in the V3 default set -- V3 is closed and its
             # numbers must stay reproducible. Select them explicitly via ARMS.
             "smc_oracle_bkmr", "smc_plugin_bkmr")

# The arm set V3 was registered and run with. Kept as the default so
# `./run_v3_bkmr.sh` reproduces the published V3 result unchanged.
V3_DEFAULT_ARMS <- c("oracle_bkmr", "cc_bkmr", "sub_lod2_bkmr", "pipeline_bkmr")

# ---- scenarios ---------------------------------------------------------------

# The mixture surface throughout -- an additive generator would make every
# interaction and curvature estimand exactly zero and the track vacuous.
# `z_form` stays "linear", which makes BKMR's own model CORRECTLY SPECIFIED for
# this generator (Y = h(logX) + gamma'Z + eps). That is deliberate: any bias
# measured here is imputation-induced, not model misspecification leaking in.
.v3_scenario <- function(name, nd_frac, censor_all = FALSE) {
  list(name = name, nd_frac = nd_frac, censor_all = censor_all,
       erf_form = "mixture", z_form = "linear", rho = 0.4, skew = 0.0,
       mcar_frac = 0.0, y_frac = 0.0)
}

v3_scenarios <- function() list(
  .v3_scenario("nd20",       0.20),
  .v3_scenario("nd40",       0.40),
  # Every exposure censored. The mixture estimands are functions of ALL THREE
  # exposures, so this is the cell where an imputation failure has the most
  # surface to act on -- `int_X1X2` and `curv_X1` in particular depend on X1 and
  # X2 jointly, and here both are imputed rather than one. It is also the most
  # expensive: the X block loops over three exposures per sweep instead of one.
  .v3_scenario("nd40_all",   0.40, censor_all = TRUE),

  # ---- V15: the f-sweep, for the attenuation SHAPE -------------------------
  # THEORY -> PREDICTION -> TEST. Curvature attenuation grows FASTER than the
  # censored fraction (V3: 20% -> 11.7%, 40% -> 56.7%, a 4.85x rise for 2x the
  # censoring) and nothing explains why. No mechanistically motivated law
  # reproduces that ratio; `f^2` is the least-bad at 4.00 and is frank
  # curve-fitting to two points. These four cells are the UNMEASURED points that
  # let the shape be pinned rather than guessed.
  #
  # APPENDED, never inserted: nd20/nd40/nd40_all keep canonical indices 1/2/3, so
  # V3's published results still reproduce bit-for-bit from the same seeds. Those
  # two cells double as calibration checks on this run.
  #
  # Focal exposure only -- `curv_X1` depends on X1 alone, which V3 confirmed
  # (nd40 -56.7% against nd40_all -58.4%: censoring the others adds ~2 pp).
  .v3_scenario("mixf10",     0.10),
  .v3_scenario("mixf30",     0.30),
  .v3_scenario("mixf50",     0.50),
  .v3_scenario("mixf60",     0.60)
)

V3_SCENARIOS <- c("nd20", "nd40", "nd40_all")

# The V15 sweep: the two measured points plus the four unmeasured ones.
V15_SCENARIOS <- c("mixf10", "nd20", "mixf30", "nd40", "mixf50", "mixf60")

.v3_make_bundle <- function(sc, truth, n_obs) {
  s <- simulate_complete(n = n_obs, truth = truth, rho = sc$rho, sd_x = 1,
                         mu_x = 0, skew = sc$skew, sigma_y = 1,
                         z_form = sc$z_form)
  complete <- s$data
  p <- length(truth$b)
  cens_which <- if (isTRUE(sc$censor_all)) seq_len(p) else 1L
  censored <- inject_left_censoring(complete, nd_frac = sc$nd_frac,
                                    censor_which = cens_which)
  list(complete = complete, censored = censored, truth = truth)
}

# ---- one task ----------------------------------------------------------------

.v3_run_task <- function(task, grid, n_obs, m, iter, n_knots, base_seed,
                         sweeps, margin, project_root, arms) {
  sc <- task$sc
  r  <- task$rep
  set.seed(base_seed + task$canon_index * 100000L + r)

  t0 <- Sys.time()
  res <- tryCatch({
    truth  <- make_truth(p = 3L, q = 2L, erf_form = "mixture")
    bundle <- .v3_make_bundle(sc, truth, n_obs)

    parts <- list()
    if ("oracle_bkmr" %in% arms)
      parts[[length(parts) + 1L]] <- proc_bkmr_oracle(
        bundle, grid, iter = iter, n_knots = n_knots, seed = base_seed + r)
    if ("cc_bkmr" %in% arms)
      parts[[length(parts) + 1L]] <- proc_bkmr_complete_case(
        bundle, grid, iter = iter, n_knots = n_knots, seed = base_seed + r)
    if ("sub_lod2_bkmr" %in% arms)
      parts[[length(parts) + 1L]] <- proc_bkmr_sub_lod2(
        bundle, grid, iter = iter, n_knots = n_knots, seed = base_seed + r)
    if ("pipeline_bkmr" %in% arms)
      parts[[length(parts) + 1L]] <- proc_bkmr_pipeline(
        bundle, grid, m = m, iter = iter, n_knots = n_knots,
        seed = base_seed + r, n_cores = 1L, outer_sweeps = sweeps,
        margin = margin, project_root = project_root)
    for (md in c("oracle", "plugin")) {
      if (paste0("smc_", md, "_bkmr") %in% arms)
        parts[[length(parts) + 1L]] <- proc_bkmr_smc(
          bundle, grid, m = m, iter = iter, n_knots = n_knots,
          seed = base_seed + r, mode = md, sweeps = sweeps,
          rho = sc$rho, sd_x = 1, mu_x = 0, sigma_y = 1)
    }

    do.call(rbind, parts)
  }, error = function(e) {
    data.frame(procedure = "TASK_ERROR", estimand = NA_character_,
               estimate = NA_real_, se = NA_real_, ci_lo = NA_real_,
               ci_hi = NA_real_, note = substr(conditionMessage(e), 1, 120),
               ubar = NA_real_, b = NA_real_, fmi = NA_real_,
               stringsAsFactors = FALSE)
  })

  res$estimand_true <- unname(grid$true[res$estimand])
  res$secs     <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res$rep      <- r
  res$scenario <- sc$name
  res
}

# ---- summarising -------------------------------------------------------------

#' Summarise by scenario x estimand x arm.
#'
#' `rel_bias` is relative to the ANALYTIC truth. Estimands whose truth is exactly
#' zero would make a relative bias undefined; none of the seven is zero on the
#' mixture generator, but `rel_bias` is guarded anyway and `bias` is always
#' reported alongside so the absolute scale is never lost.
summarise_v3 <- function(raw) {
  raw <- raw[!is.na(raw$estimand), , drop = FALSE]
  grp <- interaction(raw$scenario, raw$estimand, raw$procedure, drop = TRUE)
  parts <- split(raw, grp)

  rows <- lapply(parts, function(df) {
    truth <- df$estimand_true[1]
    est <- df$estimate
    ok <- is.finite(est)
    emp_se <- stats::sd(est, na.rm = TRUE)
    ci_w <- mean(df$ci_hi - df$ci_lo, na.rm = TRUE)
    n_ok <- sum(ok)
    data.frame(
      scenario  = df$scenario[1], estimand = df$estimand[1],
      procedure = df$procedure[1],
      n_rep = nrow(df), n_ok = n_ok, true = truth,
      mean_est = mean(est, na.rm = TRUE),
      bias     = mean(est - truth, na.rm = TRUE),
      rel_bias = if (isTRUE(abs(truth) > 1e-8)) mean(est - truth, na.rm = TRUE) / truth else NA_real_,
      # Monte-Carlo SE of the bias -- reported so a difference is never read as
      # real without its own noise floor beside it (the V7 lesson).
      mcse_bias = emp_se / sqrt(max(n_ok, 1L)),
      rmse   = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      emp_se = emp_se,
      mean_ci_w = ci_w,
      width_se_ratio = if (isTRUE(emp_se > 0)) (ci_w / 3.919928) / emp_se else NA_real_,
      coverage = mean(df$ci_lo <= truth & truth <= df$ci_hi, na.rm = TRUE),
      # Coverage's own MC error, at the nominal 0.95.
      mcse_cov = sqrt(0.95 * 0.05 / max(n_ok, 1L)),
      ubar = mean(df$ubar, na.rm = TRUE),
      b    = mean(df$b, na.rm = TRUE),
      fmi  = mean(df$fmi, na.rm = TRUE),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$procedure <- factor(out$procedure, levels = c(V3_ARMS, "TASK_ERROR"))
  out <- out[order(out$scenario, out$estimand, out$procedure), ]
  rownames(out) <- NULL
  out
}

print_v3 <- function(s) {
  for (sc in unique(s$scenario)) {
    cat(sprintf("\n================ scenario: %s ================\n", sc))
    x <- s[s$scenario == sc, ]
    for (es in unique(x$estimand)) {
      y <- x[x$estimand == es, ]
      cat(sprintf("\n  -- %s  (analytic truth = %+.4f) --\n", es, y$true[1]))
      cat(sprintf("     %-16s %5s %9s %8s %8s %8s %8s\n",
                  "arm", "n_ok", "mean_est", "rel%", "mcse%", "width/se", "cov"))
      for (i in seq_len(nrow(y))) {
        cat(sprintf("     %-16s %5d %9.4f %8.2f %8.2f %8.3f %8.3f\n",
                    as.character(y$procedure[i]), y$n_ok[i], y$mean_est[i],
                    100 * y$rel_bias[i],
                    100 * y$mcse_bias[i] / abs(y$true[i]),
                    y$width_se_ratio[i], y$coverage[i]))
      }
    }
  }
  invisible(s)
}

# ---- runner ------------------------------------------------------------------

run_v3 <- function(n_rep = 200L, n_obs = 800L, m = 10L, iter = 3000L,
                   n_knots = NULL, base_seed = 20260828L, n_cores = NULL,
                   which_scenarios = V3_SCENARIOS, arms = V3_DEFAULT_ARMS,
                   sweeps = 3L, margin = "shash", project_root = NULL,
                   q_lo = 0.25, q_hi = 0.75, verbose = TRUE) {

  n_cores <- n_cores %||% max(1L, parallel::detectCores() - 2L)
  if (.Platform$OS.type == "windows") n_cores <- 1L

  bad <- setdiff(arms, V3_ARMS)
  if (length(bad)) {
    stop("Unknown arm(s): ", paste(bad, collapse = ", "),
         ". Valid: ", paste(V3_ARMS, collapse = ", "), call. = FALSE)
  }

  truth <- make_truth(p = 3L, q = 2L, erf_form = "mixture")
  grid  <- bkmr_estimand_grid(truth, mu_x = 0, sd_x = 1,
                              q_lo = q_lo, q_hi = q_hi)

  # The truth check runs EVERY time, not just in the smoke test: if the grid and
  # the closed forms ever disagree, every number downstream is measured against
  # the wrong target, and that must stop the run rather than be discovered later.
  chk <- check_bkmr_truth(truth, mu_x = 0, sd_x = 1, q_lo = q_lo, q_hi = q_hi)
  if (verbose) {
    cat(sprintf("Evaluation quantiles: lo = %.2f, mid = 0.50, hi = %.2f",
                q_lo, q_hi))
    if (!isTRUE(all.equal(c(q_lo, q_hi), c(0.25, 0.75))))
      cat("   <-- NON-DEFAULT; truths differ from the registered run")
    cat("\n")
    cat("Analytic estimand truth (verified against independent closed forms):\n")
    print(data.frame(estimand = chk$estimand, true = round(chk$grid, 5)),
          row.names = FALSE)
  }

  scs_all <- v3_scenarios()
  canon_names <- vapply(scs_all, `[[`, "", "name")
  keep <- canon_names %in% which_scenarios
  if (!any(keep)) stop("No scenarios matched.", call. = FALSE)
  scs <- scs_all[keep]

  if (verbose) cat("Loading pipeline environment once in parent ...\n")
  invisible(v1_load_pipeline(project_root, quiet = TRUE))

  tasks <- list(); k <- 0L
  for (sc in scs) for (r in seq_len(n_rep)) {
    k <- k + 1L
    tasks[[k]] <- list(sc = sc, rep = r, canon_index = match(sc$name, canon_names))
  }

  if (verbose) {
    cat(sprintf("Scenarios: %s\n", paste(vapply(scs, `[[`, "", "name"), collapse = ", ")))
    cat(sprintf("Arms: %s\n", paste(arms, collapse = ", ")))
    cat(sprintf("Tasks: %d | workers: %d | n = %d | m = %d | iter = %d | knots = %s\n",
                length(tasks), n_cores, n_obs, m, iter,
                if (is.null(n_knots)) "none (full GP)" else n_knots))
  }

  # Chunked + checkpointed, as V4: a kill costs one chunk, not the run.
  ckpt_file <- getOption("v3.checkpoint", NULL)
  done <- list()
  chunk_size <- max(n_cores, 2L * n_cores)
  chunks <- split(seq_along(tasks), ceiling(seq_along(tasks) / chunk_size))

  t0 <- Sys.time()
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    res_chunk <- parallel::mclapply(
      tasks[idx], .v3_run_task,
      grid = grid, n_obs = n_obs, m = m, iter = iter, n_knots = n_knots,
      base_seed = base_seed, sweeps = sweeps, margin = margin,
      project_root = project_root, arms = arms,
      mc.cores = n_cores, mc.preschedule = FALSE)
    done[[length(done) + 1L]] <- res_chunk

    if (!is.null(ckpt_file)) {
      flat <- unlist(done, recursive = FALSE)
      dfs <- Filter(is.data.frame, flat)
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
  n_fail <- length(flat) - length(dfs)
  if (n_fail > 0) cat(sprintf("WARNING: %d worker(s) failed\n", n_fail))
  raw <- do.call(rbind, dfs)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (verbose) cat(sprintf("\nDone in %.1f min\n", elapsed / 60))

  list(raw = raw, summary = summarise_v3(raw),
       meta = list(n_rep = n_rep, n_obs = n_obs, m = m, iter = iter,
                   n_knots = n_knots, base_seed = base_seed, n_cores = n_cores,
                   scenarios = which_scenarios, arms = arms, sweeps = sweeps,
                   margin = margin, elapsed_secs = elapsed,
                   q_lo = q_lo, q_hi = q_hi, quantiles = grid$quantiles,
                   truth = grid$true, truth_check = chk))
}

# --- script entry point ------------------------------------------------------
if (sys.nframe() == 0L) {
  n_rep  <- as.integer(getenv("N_REP", 200L))
  n_obs  <- as.integer(getenv("N_OBS", 800L))   # 800 = the calibrated setting; see the .sh header
  m      <- as.integer(getenv("M", 10L))
  iter   <- as.integer(getenv("ITER", 3000L))
  # 50 knots is the calibrated default -- the gate, the 10% bar and the cost
  # figures were all established with it. KNOTS=0 selects the full GP (~10x
  # slower) and is only for re-checking the approximation.
  knots  <- as.integer(getenv("KNOTS", 50L)); if (knots <= 0L) knots <- NULL
  seed   <- as.integer(getenv("SEED", 20260828L))
  ncores <- as.integer(getenv("NCORES", max(1L, parallel::detectCores() - 2L)))
  sweeps <- as.integer(getenv("SWEEPS", 3L))
  margin <- getenv("MARGIN", "shash")
  scs    <- if (nzchar(Sys.getenv("SCENARIOS"))) split_csv(Sys.getenv("SCENARIOS")) else V3_SCENARIOS
  arms   <- if (nzchar(Sys.getenv("ARMS")))      split_csv(Sys.getenv("ARMS"))      else V3_DEFAULT_ARMS
  proot  <- if (nzchar(Sys.getenv("PIPELINE_ROOT"))) Sys.getenv("PIPELINE_ROOT") else NULL
  q_lo   <- as.numeric(getenv("QLO", 0.25))
  q_hi   <- as.numeric(getenv("QHI", 0.75))

  tag <- getenv("OUT_TAG", "v3")
  results_dir <- file.path(.here, "results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  options(v3.checkpoint = file.path(results_dir, paste0(tag, "_checkpoint.rds")))

  res <- run_v3(n_rep = n_rep, n_obs = n_obs, m = m, iter = iter,
                n_knots = knots, base_seed = seed, n_cores = ncores,
                which_scenarios = scs, arms = arms, sweeps = sweeps,
                margin = margin, project_root = proot,
                q_lo = q_lo, q_hi = q_hi)

  cat("\n================ V3 BKMR mixture estimands ================\n")
  print_v3(res$summary)

  saveRDS(res, file.path(results_dir, paste0(tag, "_latest.rds")))
  utils::write.csv(res$summary, file.path(results_dir, paste0(tag, "_summary.csv")),
                   row.names = FALSE)
  utils::write.csv(res$raw, file.path(results_dir, paste0(tag, "_raw.csv")),
                   row.names = FALSE)
  cat(sprintf("\nWrote results to %s (%s_latest.rds, %s_summary.csv, %s_raw.csv)\n",
              results_dir, tag, tag, tag))
}
