#!/usr/bin/env Rscript
# =============================================================================
# Phase 2 -- the TRAP demonstration on the mixture benchmark   (DESIGN §8 step 4)
# -----------------------------------------------------------------------------
# Shows that half-limit SUBSTITUTION and non-surface-congenial IMPUTATION bias the
# recovered mixture-HBEB *benchmark*, and that the bias is worse on the non-linear
# surface (the §7.7 trap, now on the quantity that matters for policy). Every
# procedure feeds the SAME analysis (bkmr surface) + inversion (benchmark.R), so
# differences are attributable to censored-exposure handling alone.
#
# Procedures (surface-congenial procedure 5 is the FIX, built later -- not here):
#   oracle        complete uncensored exposures            (recovery ceiling)
#   substitution  ND -> MDL/2, DNQ -> geometric midpoint   (reference-paper method)
#   impute_blind  leftcens Y-blind draw (no outcome)        (uncongenial)
#   impute_linear leftcens Y-aware LINEAR draw (shash)      (congenial but linear -> the trap)
#
# Uses full `bkmr` at small n (the published fidelity engine; A-BKMR comes later).
# Config via env: N, ITER, M, NKNOTS, SWEEPS, ERF (mixture|additive), LIMIT_SCALE, SEED, REPS.
# =============================================================================

suppressWarnings(suppressMessages({library(MASS); library(leftcens); library(bkmr)}))

.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "R"
}
for (fn in c("benchmark.R", "dgp.R", "censoring.R")) source(file.path(.here, fn))

getenv <- function(k, d) { v <- Sys.getenv(k, NA); if (is.na(v) || !nzchar(v)) d else v }

# ---- adjustment covariates (confounders R,parity + BW precision predictors) ---
# R is the latent renal confounder; we adjust for it so every arm targets the
# CAUSAL surface and the arm differences isolate censored-exposure handling.
covariate_df <- function(d) {
  data.frame(R = d$R, parity = d$parity, ga = d$ga, male = d$male,
             smk_current = as.integer(d$smoking == "current"),
             smk_past    = as.integer(d$smoking == "past"),
             height_z = d$height_z, bmi_z = d$bmi_z, age = d$age)
}

# ---- procedures: each returns a list of completed data.frames -----------------

proc_oracle <- function(dcomplete, ...) list(dcomplete)

proc_substitution <- function(dcens, truth, ...) {
  dd <- dcens
  for (a in truth$analytes) {
    tier <- dcens[[paste0(a, "_tier")]]
    lo <- dcens[[paste0(a, "_lo")]]; hi <- dcens[[paste0(a, "_hi")]]
    x <- dcens[[a]]
    x[tier == "nd"]  <- (hi[tier == "nd"]) / 2                 # MDL/2
    x[tier == "dnq"] <- sqrt(lo[tier == "dnq"] * hi[tier == "dnq"])  # geometric mid
    dd[[a]] <- x
  }
  list(dd)
}

#' Block-FCS leftcens imputation. use_y toggles the outcome in the predictor set
#' (congenial vs blind). Linear conditioning (shash margin) -> the §7.7 trap.
block_fcs <- function(dcens, truth, use_y, m = 5L, sweeps = 3L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  an <- truth$analytes; Zc <- covariate_df(dcens); ybw <- dcens$bw; n <- nrow(dcens)
  # per-analyte log bounds + censored flag
  tier <- lapply(an, function(a) dcens[[paste0(a, "_tier")]]); names(tier) <- an
  lwr <- sapply(an, function(a) {
    t <- tier[[a]]; ifelse(t == "quant", NA, ifelse(t == "nd", -Inf,
                            log(dcens[[paste0(a, "_lo")]]))) })
  upr <- sapply(an, function(a) {
    t <- tier[[a]]; ifelse(t == "quant", NA, log(dcens[[paste0(a, "_hi")]])) })
  out <- vector("list", m)
  for (im in seq_len(m)) {
    L <- sapply(an, function(a) {
      t <- tier[[a]]; lo <- dcens[[paste0(a,"_lo")]]; hi <- dcens[[paste0(a,"_hi")]]
      ifelse(t == "quant", log(dcens[[a]]),
             ifelse(t == "nd", log(hi/2), (log(lo)+log(hi))/2)) })   # midpoint init
    for (s in seq_len(sweeps)) for (j in seq_along(an)) {
      cens <- tier[[an[j]]] != "quant"; if (!any(cens)) next
      preds <- cbind(Zc, L[, -j, drop = FALSE])
      if (use_y) preds <- cbind(bw = ybw, preds)
      yv <- L[, j]; yv[cens] <- NA
      imp <- tryCatch(leftcens::impute_censored_conditional(
               y = yv, x = preds,
               lower = ifelse(cens, lwr[, j], NA), upper = ifelse(cens, upr[, j], NA),
               m = 1L, margin = "shash"), error = function(e) NULL)
      if (!is.null(imp)) L[cens, j] <- imp[cens, 1]
    }
    dd <- dcens; for (j in seq_along(an)) dd[[an[j]]] <- exp(L[, j]); out[[im]] <- dd
  }
  out
}

# ---- bkmr surface + benchmark for one completed dataset ----------------------

fit_bkmr_surface <- function(dfilled, truth, iter = 1000L, nknots = 50L) {
  an <- truth$analytes
  Z <- log(as.matrix(dfilled[, an]))
  ctr <- colMeans(Z); scl <- apply(Z, 2, stats::sd); Zs <- scale(Z, ctr, scl)
  Xc <- as.matrix(covariate_df(dfilled))
  y  <- dfilled$bw; yc <- mean(y); ys <- stats::sd(y); yst <- (y - yc) / ys
  kn <- Zs[sample.int(nrow(Zs), min(nknots, nrow(Zs))), , drop = FALSE]
  invisible(utils::capture.output(         # kmbayes prints progress even at verbose=FALSE
    fit <- bkmr::kmbayes(y = yst, Z = Zs, X = Xc, iter = iter, verbose = FALSE,
                         varsel = FALSE, knots = kn)))
  sel <- seq(floor(iter/2), iter, by = max(1L, floor(iter/200)))   # thinned post-burnin
  predict_h <- function(Xmat) {
    Zn <- scale(log(Xmat[, an, drop = FALSE]), center = ctr, scale = scl)
    bkmr::ComputePostmeanHnew(fit, Znew = Zn, sel = sel, method = "approx")$postmean * ys
  }
  make_surface(predict_h, an)
}

#' Benchmark (proportional direction, 25% reduction from the arm's OWN current)
#' for one completed dataset.
bench_of <- function(dfilled, truth, iter, nknots, rel = 0.25) {
  an <- truth$analytes; Xmat <- as.matrix(dfilled[, an])
  tr <- default_translation(truth$bw_mean, truth$bw_sd, truth$lbw_cut, truth$N)
  surf <- fit_bkmr_surface(dfilled, truth, iter, nknots)
  dv <- dir_proportional(an)
  curr <- residual_cases(Xmat, surf, 0, dvec = dv, tr = tr)
  b <- benchmark_mixture(Xmat, surf, target_cases = curr * (1 - rel), dvec = dv, tr = tr)
  c(theta = b$theta, b$conc)
}

# ---- driver ------------------------------------------------------------------

#' One replication: returns a data.frame of the PFOS/PFHxS benchmark + theta per
#' method (proportional direction, `rel`-reduction from each arm's own current).
one_rep <- function(truth, n, iter, m, nknots, sweeps, limit_scale, seed, rel) {
  dcomplete <- simulate_complete(n, truth, seed = seed)$data
  # MAR covariates are a separately-solved pipeline step (miceRanger); this demo
  # isolates the CENSORED-EXPOSURE trap, so covariates stay complete.
  dcens <- inject_censoring(dcomplete, limit_scale = limit_scale)
  procs <- list(
    oracle        = proc_oracle(dcomplete),
    substitution  = proc_substitution(dcens, truth),
    impute_blind  = block_fcs(dcens, truth, use_y = FALSE, m = m, sweeps = sweeps, seed = seed),
    impute_linear = block_fcs(dcens, truth, use_y = TRUE,  m = m, sweeps = sweeps, seed = seed))
  rows <- lapply(names(procs), function(nm) {
    bs <- lapply(procs[[nm]], function(df)
      tryCatch(bench_of(df, truth, iter, nknots, rel), error = function(e) NULL))
    bs <- bs[!vapply(bs, is.null, logical(1))]
    est <- if (length(bs)) rowMeans(do.call(cbind, bs))
           else stats::setNames(rep(NA, length(truth$analytes)+1), c("theta", truth$analytes))
    cbind(data.frame(method = nm, stringsAsFactors = FALSE), as.data.frame(as.list(est)))
  })
  do.call(rbind, rows)
}

run_trap_demo <- function(n = 800L, iter = 2500L, m = 5L, nknots = 50L, sweeps = 3L,
                          erf_form = "mixture", limit_scale = 1.0, seed = 1L,
                          rel = 0.25, reps = 3L) {
  truth <- make_truth(erf_form, "real")
  tb <- true_benchmark(truth, rel_targets = rel, n_pop = 1e5, seed = 99L)
  tb <- tb[tb$direction == "prop", ]
  truth_pfos <- tb[[ "pfos" ]]; truth_theta <- tb$theta

  cat(sprintf("\n=== TRAP DEMO  n=%d erf=%s limit_scale=%.1f iter=%d m=%d reps=%d (%d%%-reduction) ===\n",
              n, erf_form, limit_scale, iter, m, reps, round(rel*100)))
  cat(sprintf("TRUTH: theta*=%.4f  PFOS benchmark=%.3f ng/mL\n", truth_theta, truth_pfos))

  all <- list()
  for (rp in seq_len(reps)) {
    t0 <- Sys.time()
    r <- one_rep(truth, n, iter, m, nknots, sweeps, limit_scale, seed + rp, rel)
    r$rep <- rp; all[[rp]] <- r
    cat(sprintf("  rep %d/%d done (%.0fs)\n", rp, reps, as.numeric(difftime(Sys.time(), t0, units="secs"))))
  }
  res <- do.call(rbind, all)

  # aggregate: mean PFOS benchmark per method, bias vs truth and vs ORACLE
  agg <- aggregate(pfos ~ method, res, function(x) c(mean = mean(x), sd = stats::sd(x)))
  agg <- do.call(data.frame, agg); names(agg) <- c("method","pfos_mean","pfos_sd")
  oracle_mean <- agg$pfos_mean[agg$method == "oracle"]
  agg$bias_truth  <- 100 * (agg$pfos_mean - truth_pfos)  / truth_pfos
  agg$bias_oracle <- 100 * (agg$pfos_mean - oracle_mean) / oracle_mean
  ord <- c("oracle","substitution","impute_blind","impute_linear")
  agg <- agg[match(ord, agg$method), ]

  cat(sprintf("\n%-14s %10s %8s %12s %13s\n", "method","PFOS","sd","bias_truth%","bias_oracle%"))
  for (i in seq_len(nrow(agg))) cat(sprintf("%-14s %10.3f %8.3f %+11.1f %+12.1f\n",
      agg$method[i], agg$pfos_mean[i], agg$pfos_sd[i], agg$bias_truth[i], agg$bias_oracle[i]))
  cat("\n(bias_oracle isolates censored-exposure handling from bkmr's common attenuation.)\n")

  outdir <- file.path(.here, "..", "results"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  tag <- sprintf("trap_%s_ls%.0f_n%d", erf_form, limit_scale*10, n)
  saveRDS(list(res = res, agg = agg, truth_pfos = truth_pfos, truth_theta = truth_theta),
          file.path(outdir, paste0(tag, ".rds")))
  cat(sprintf("saved: results/%s.rds\n", tag))
  invisible(list(res = res, agg = agg))
}

if (sys.nframe() == 0) {
  run_trap_demo(
    n = as.integer(getenv("N", 800)), iter = as.integer(getenv("ITER", 2500)),
    m = as.integer(getenv("M", 5)), nknots = as.integer(getenv("NKNOTS", 50)),
    sweeps = as.integer(getenv("SWEEPS", 3)), erf_form = getenv("ERF", "mixture"),
    limit_scale = as.numeric(getenv("LIMIT_SCALE", 1.0)),
    seed = as.integer(getenv("SEED", 1)), reps = as.integer(getenv("REPS", 3)))
}
