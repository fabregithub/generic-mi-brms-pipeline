# =============================================================================
# OLS analysis arms -- the same estimands, fitted on the TRUE basis instead of
# with BKMR.
# -----------------------------------------------------------------------------
# WHY THESE EXIST. P-B (2026-10-26, `DERIVATION_declared_basis.md`) built a
# parameter-free prediction of declared-basis bias and, with it, a replica that
# matches V33's imputation in every respect checked against the source and
# differs ONLY in using OLS on the true basis instead of BKMR. The replica
# agrees with the asymptotic fixed point to <= 0.025 in all twelve cells and
# still gives f-exponents of 0.83-1.33 against a measured 1.46-2.28. That
# localises the unexplained `f^2` scaling -- BOTH routes to it, V15's
# `excess% = -323.4*f^2` and V34's `f^1.85` -- to the ANALYSIS model's response
# to degraded curvature information, not to the imputation.
#
# These arms make that testable on V15's own ground: run the f-sweep with the
# BKMR arms and the OLS arms side by side, on the SAME simulated datasets and
# the SAME seeds, so the comparison is paired and the only thing that differs
# is the analysis model.
#
# WHAT AN OLS ARM IS AND IS NOT. It is not a competitor to BKMR and not a
# recommendation: it is given the generator's own basis, so it is correctly
# specified and efficient by construction, which BKMR is not asked to be. It is
# an INSTRUMENT -- the analysis model whose response to information loss is
# known exactly, so that whatever BKMR adds shows up as the difference.
#
# HOW THE ESTIMANDS ARE COMPUTED. Every estimand is a zero-sum contrast of h at
# fixed evaluation points (`R/estimands_bkmr.R`), so for a linear model it is a
# linear functional of beta:
#       est = (C %*% Xpts) %*% beta,     var = L V L',   L = C %*% Xpts
# with Xpts the design rows at the evaluation points. The covariates are set to
# zero when building Xpts: every contrast is zero-sum, so their contribution
# cancels exactly and the choice cannot affect any estimand. Intervals are
# normal-theory on the residual df -- exact here, where BKMR needed percentile
# intervals because its posteriors are not guaranteed symmetric.
# =============================================================================

#' Design rows for the estimand grid's evaluation points, under `dgp_formula()`.
#'
#' @return A matrix with one row per evaluation point and one column per model
#'   coefficient, in the order `lm()` produces them.
.lm_grid_design <- function(grid, truth) {
  pts <- as.data.frame(grid$Z)
  names(pts) <- paste0("logX", seq_len(grid$p))
  # Zero-sum contrasts: the covariate columns cancel, so any value does.
  pts$Z1 <- 0; pts$Z2 <- 0
  pts$Y  <- 0
  stats::model.matrix(dgp_formula(truth), data = pts)
}

#' One `lm` fit -> the estimand table, in `bkmr_rows()`'s contract.
.lm_tab <- function(fit, grid, truth) {
  Xp <- .lm_grid_design(grid, truth)
  b  <- stats::coef(fit)
  if (any(is.na(b)) || !identical(colnames(Xp), names(b))) return(NULL)
  L  <- grid$C %*% Xp
  est <- as.vector(L %*% b)
  V   <- stats::vcov(fit)
  vr  <- as.vector(rowSums((L %*% V) * L))
  df  <- fit$df.residual
  tq  <- stats::qt(0.975, df)
  data.frame(estimand = grid$names, est = est, se = sqrt(vr),
             ci_lo = est - tq*sqrt(vr), ci_hi = est + tq*sqrt(vr),
             ubar = NA_real_, b = NA_real_, fmi = NA_real_,
             stringsAsFactors = FALSE)
}

#' Per-imputation estimates and variances, for Rubin pooling.
.lm_moments <- function(d, grid, truth) {
  fit <- tryCatch(stats::lm(dgp_formula(truth), data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tb <- .lm_tab(fit, grid, truth)
  if (is.null(tb)) return(NULL)
  list(estimand = tb$estimand, est = tb$est, var = tb$se^2)
}

# ---- arm: complete data, OLS ------------------------------------------------

proc_lm_oracle <- function(bundle, grid, seed = NULL, label = "oracle_lm", ...) {
  fit <- tryCatch(stats::lm(dgp_formula(bundle$truth), data = bundle$complete),
                  error = function(e) NULL)
  if (is.null(fit)) return(bkmr_na_rows(label, grid, "lm fit failed"))
  tb <- .lm_tab(fit, grid, bundle$truth)
  if (is.null(tb)) return(bkmr_na_rows(label, grid, "rank-deficient lm fit"))
  bkmr_rows(label, tb, "ols on the true basis")
}

# ---- arm: shipped block-FCS imputation, OLS analysis ------------------------
# Deliberately mirrors `proc_bkmr_pipeline()` step for step -- same imputation
# call, same per-dataset loop, same `rubin_pool()` -- so the ONLY difference
# between `pipeline_lm` and `pipeline_bkmr` is the analysis model.

proc_lm_pipeline <- function(bundle, grid, m = 10L, seed = NULL, n_cores = 1L,
                             outer_sweeps = 3L, margin = "shash",
                             project_root = NULL, quiet = TRUE,
                             label = "pipeline_lm", ...) {
  if (!requireNamespace("leftcens", quietly = TRUE) ||
      !"impute_censored_conditional" %in% getNamespaceExports("leftcens")) {
    return(bkmr_na_rows(label, grid, "needs leftcens >= 0.9.0"))
  }
  env <- tryCatch(v1_load_pipeline(project_root, quiet = quiet),
                  error = function(e) NULL)
  if (is.null(env)) return(bkmr_na_rows(label, grid, "could not source pipeline"))
  if (!is.null(seed)) set.seed(seed)

  inp <- .v1_make_pipeline_inputs(bundle, outer_sweeps = outer_sweeps,
                                  margin = margin, n_cores = n_cores,
                                  mid = TRUE, z_imputer = "bart")
  imputed <- tryCatch(
    env$run_censored_exposure_block_fcs(
      data = inp$data, analysis_spec = inp$analysis_spec,
      var_dict = inp$var_dict, m = as.integer(m),
      seed = if (is.null(seed)) 1L else as.integer(seed)),
    error = function(e) structure(list(), err = conditionMessage(e)))
  if (!length(imputed))
    return(bkmr_na_rows(label, grid,
                        paste("failed:", attr(imputed, "err") %||% "no datasets")))

  ests <- vars <- matrix(NA_real_, nrow = length(imputed),
                         ncol = length(grid$names),
                         dimnames = list(NULL, grid$names))
  for (i in seq_along(imputed)) {
    di <- as.data.frame(imputed[[i]])
    for (k in seq_along(inp$cens_x)) di[[inp$cens_x[k]]] <- log(di[[inp$expo_names[k]]])
    mo <- .lm_moments(di, grid, bundle$truth)
    if (is.null(mo)) next
    ests[i, mo$estimand] <- mo$est
    vars[i, mo$estimand] <- mo$var
  }
  n_ok <- sum(is.finite(ests[, 1]))
  if (n_ok == 0L) return(bkmr_na_rows(label, grid, "all per-imputation lm fits failed"))

  tab <- do.call(rbind, lapply(grid$names, function(nm) {
    pl <- rubin_pool(ests[, nm], vars[, nm])
    data.frame(estimand = nm, est = unname(pl["est"]), se = unname(pl["se"]),
               ci_lo = unname(pl["ci_lo"]), ci_hi = unname(pl["ci_hi"]),
               ubar = unname(pl["ubar"]), b = unname(pl["b"]),
               fmi = unname(pl["fmi"]), stringsAsFactors = FALSE)
  }))
  note <- sprintf("m=%d; sweeps=%d; %s; cens=%d", m, outer_sweeps, margin,
                  length(inp$cens_x))
  if (n_ok < length(imputed))
    note <- paste0(note, sprintf("; %d/%d lm fits failed", length(imputed) - n_ok, length(imputed)))
  bkmr_rows(label, tab, substr(note, 1, 120))
}
