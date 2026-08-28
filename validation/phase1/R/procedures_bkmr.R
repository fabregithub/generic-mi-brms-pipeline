# =============================================================================
# Track V3 harness -- the BKMR arms
# -----------------------------------------------------------------------------
# Each arm returns ONE ROW PER ESTIMAND (not one row per replication, as the
# scalar-estimand tracks do), with the columns the rest of the harness uses:
#
#   procedure, estimand, estimate, se, ci_lo, ci_hi, note, ubar, b, fmi
#
#   oracle_bkmr    complete, uncensored data                (tests the HARNESS)
#   cc_bkmr        listwise-drop the censored rows          (the lazy baseline)
#   sub_lod2_bkmr  substitute LOD/sqrt(2)                    (what is done in practice)
#   pipeline_bkmr  the shipped block-FCS, m imputations, BKMR each, Rubin pooled
#                                                            (the ARM UNDER TEST)
#
# WHAT `oracle_bkmr` IS FOR. In every other track the oracle is the yardstick.
# Here it is not: the yardstick is the analytic truth from `estimands_bkmr.R`.
# The oracle instead tests the *harness* -- if BKMR on complete data does not
# recover the analytic truth with near-nominal coverage, then the estimand
# extraction, the quantile convention, or the MCMC length is wrong, and no
# statement about the other arms means anything. It is a gate, not a reference.
#
# WHY NO `leftcens_prestep` ARM. Phase 1 already established that imputing X
# without Y attenuates the ERF; repeating it here would spend BKMR fits to
# re-confirm a closed question. The open question is whether the *congenial*
# engine -- which does condition on Y -- survives a surface its linear
# conditional cannot represent.
# =============================================================================

# ---- fitting one BKMR model --------------------------------------------------

#' Fit one BKMR model and return posterior draws of every V3 estimand.
#'
#' @param d A data.frame with `Y`, `logX1..logXp`, `Z1`, `Z2` and no missing
#'   values in those columns.
#' @param grid A list from [bkmr_estimand_grid()].
#' @param iter MCMC iterations. `sel` is the second half.
#' @param n_knots If non-NULL, fit the Gaussian predictive process with this many
#'   knots (`fields::cover.design`) instead of the full GP. Cuts the per-iteration
#'   cost from O(n^3) to O(n * k^2) at the price of an approximation.
#' @param seed Optional RNG seed.
#' @return A draws x estimands matrix, or NULL on failure.
fit_bkmr_draws <- function(d, grid, iter = 1000L, n_knots = NULL, seed = NULL) {
  p <- grid$p
  xcols <- paste0("logX", seq_len(p))
  zcols <- intersect(c("Z1", "Z2"), names(d))
  need <- c("Y", xcols, zcols)
  if (!all(need %in% names(d))) return(NULL)

  d <- d[stats::complete.cases(d[, need, drop = FALSE]), , drop = FALSE]
  if (nrow(d) < 20L) return(NULL)

  Z <- as.matrix(d[, xcols, drop = FALSE])
  X <- as.matrix(d[, zcols, drop = FALSE])
  y <- d$Y

  if (!is.null(seed)) set.seed(seed)

  knots <- NULL
  if (!is.null(n_knots) && n_knots > 0 && n_knots < nrow(Z)) {
    knots <- tryCatch(fields::cover.design(Z, nd = as.integer(n_knots))$design,
                      error = function(e) NULL)
  }

  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      bkmr::kmbayes(y = y, Z = Z, X = X, iter = as.integer(iter),
                    verbose = FALSE, varsel = FALSE, knots = knots))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  bkmr_estimand_draws(fit, grid)
}

# ---- row assembly ------------------------------------------------------------

#' One result row per estimand, in the harness's column convention.
bkmr_rows <- function(procedure, tab, note = "") {
  data.frame(procedure = procedure, estimand = tab$estimand,
             estimate = tab$est, se = tab$se,
             ci_lo = tab$ci_lo, ci_hi = tab$ci_hi,
             note = note,
             ubar = tab$ubar %||% NA_real_,
             b = tab$b %||% NA_real_,
             fmi = tab$fmi %||% NA_real_,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Empty (all-NA) rows, so a failed arm still occupies its place in the grid
#' rather than silently vanishing from the summary.
bkmr_na_rows <- function(procedure, grid, note) {
  data.frame(procedure = procedure, estimand = grid$names,
             estimate = NA_real_, se = NA_real_, ci_lo = NA_real_,
             ci_hi = NA_real_, note = substr(note, 1, 120),
             ubar = NA_real_, b = NA_real_, fmi = NA_real_,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Single-fit arms: summarise the posterior directly.
#'
#' Uses the PERCENTILE credible interval rather than mean +/- 1.96 SD. BKMR
#' posteriors for interaction and curvature contrasts are not guaranteed
#' symmetric, and a normal-approximation interval would then mis-state coverage
#' for reasons that have nothing to do with the imputation under test.
.bkmr_single_tab <- function(draws) {
  mo <- bkmr_draws_moments(draws)
  data.frame(estimand = mo$estimand, est = mo$est, se = sqrt(mo$var),
             ci_lo = mo$cri_lo, ci_hi = mo$cri_hi,
             ubar = NA_real_, b = NA_real_, fmi = NA_real_,
             stringsAsFactors = FALSE)
}

# ---- arm 1: oracle -----------------------------------------------------------

proc_bkmr_oracle <- function(bundle, grid, iter = 1000L, n_knots = NULL,
                             seed = NULL) {
  dr <- fit_bkmr_draws(bundle$complete, grid, iter = iter, n_knots = n_knots,
                       seed = seed)
  if (is.null(dr)) return(bkmr_na_rows("oracle_bkmr", grid, "bkmr fit failed"))
  bkmr_rows("oracle_bkmr", .bkmr_single_tab(dr), sprintf("iter=%d", iter))
}

# ---- arm 2: complete case ----------------------------------------------------

proc_bkmr_complete_case <- function(bundle, grid, iter = 1000L, n_knots = NULL,
                                    seed = NULL) {
  # `logXj` is NA on censored rows, so complete.cases() inside fit_bkmr_draws()
  # performs the listwise drop.
  d <- bundle$censored
  dr <- fit_bkmr_draws(d, grid, iter = iter, n_knots = n_knots, seed = seed)
  if (is.null(dr)) return(bkmr_na_rows("cc_bkmr", grid, "bkmr fit failed"))
  n_used <- sum(stats::complete.cases(
    d[, c("Y", paste0("logX", seq_len(grid$p)), "Z1", "Z2")]))
  bkmr_rows("cc_bkmr", .bkmr_single_tab(dr),
            sprintf("iter=%d; n=%d", iter, n_used))
}

# ---- arm 3: LOD/sqrt(2) substitution ----------------------------------------

#' The default practice in the applied literature, included so the study says
#' something about what people actually do rather than only about the pipeline.
proc_bkmr_sub_lod2 <- function(bundle, grid, iter = 1000L, n_knots = NULL,
                               seed = NULL) {
  d <- bundle$censored
  cens_cols <- attr(d, "censored_cols") %||% "logX1"
  for (cc in cens_cols) {
    lod_col <- paste0(cc, "_lod")
    if (!lod_col %in% names(d)) next
    is_c <- is.na(d[[cc]])
    # LOD/sqrt(2) is defined on the CONCENTRATION scale; the columns here are
    # logs, so the substitution is log(LOD) - log(sqrt(2)).
    d[[cc]][is_c] <- d[[lod_col]][is_c] - log(sqrt(2))
  }
  dr <- fit_bkmr_draws(d, grid, iter = iter, n_knots = n_knots, seed = seed)
  if (is.null(dr)) return(bkmr_na_rows("sub_lod2_bkmr", grid, "bkmr fit failed"))
  bkmr_rows("sub_lod2_bkmr", .bkmr_single_tab(dr), sprintf("iter=%d", iter))
}

# ---- arm 4: the shipped block-FCS, then BKMR per imputation ------------------

#' The arm under test.
#'
#' Reuses `.v1_make_pipeline_inputs()` and the shipped
#' `run_censored_exposure_block_fcs()`, so this exercises pipeline code rather
#' than a harness copy of it -- the lesson V8 recorded about V7's BART arm.
#'
#' POOLING. Each imputed dataset yields a BKMR posterior per estimand; the
#' posterior mean and variance go into `rubin_pool()`. This deliberately matches
#' what the pipeline's Step 6 does in spirit (within + between decomposition)
#' rather than concatenating draws, so `ubar`, `b` and `fmi` are recorded per
#' estimand and an interval failure can be attributed the way V4 attributed one.
#'
#' @param m Completed datasets. NOTE the cost: this arm fits `m` BKMR models per
#'   replication, so it dominates the run.
proc_bkmr_pipeline <- function(bundle, grid, m = 10L, iter = 1000L,
                               n_knots = NULL, seed = NULL, n_cores = 1L,
                               outer_sweeps = 3L, margin = "shash",
                               project_root = NULL, quiet = TRUE,
                               label = "pipeline_bkmr") {
  if (!requireNamespace("leftcens", quietly = TRUE) ||
      !"impute_censored_conditional" %in% getNamespaceExports("leftcens")) {
    return(bkmr_na_rows(label, grid, "needs leftcens >= 0.9.0"))
  }
  env <- tryCatch(v1_load_pipeline(project_root, quiet = quiet),
                  error = function(e) NULL)
  if (is.null(env)) return(bkmr_na_rows(label, grid, "could not source pipeline"))

  if (!is.null(seed)) set.seed(seed)

  # `z_imputer = "bart"` is passed EXPLICITLY rather than left NULL. The V3
  # scenarios have complete covariates, so the Z block is idle and the choice is
  # currently inert -- but `.v1_make_pipeline_inputs()` defaults it to NULL,
  # which `resolve_z_imputer()` maps to the pre-v1.4.0 IMPROPER forest. Anyone
  # adding covariate missingness to a V3 scenario would then silently be testing
  # an imputer the pipeline stopped shipping two versions ago. Naming the
  # shipped default here makes the arm what it claims to be.
  inp <- .v1_make_pipeline_inputs(bundle, outer_sweeps = outer_sweeps,
                                  margin = margin, n_cores = n_cores, mid = TRUE,
                                  z_imputer = "bart")

  imputed <- tryCatch(
    env$run_censored_exposure_block_fcs(
      data = inp$data, analysis_spec = inp$analysis_spec,
      var_dict = inp$var_dict, m = as.integer(m),
      seed = if (is.null(seed)) 1L else as.integer(seed)),
    error = function(e) structure(list(), err = conditionMessage(e)))

  if (!length(imputed)) {
    return(bkmr_na_rows(label, grid,
                        paste("failed:", attr(imputed, "err") %||% "no datasets")))
  }

  # One BKMR fit per completed dataset; collect the per-estimand moments.
  ests <- vars <- matrix(NA_real_, nrow = length(imputed),
                         ncol = length(grid$names),
                         dimnames = list(NULL, grid$names))

  for (i in seq_along(imputed)) {
    di <- as.data.frame(imputed[[i]])
    # The module returns each censored exposure on the DATA scale; the estimand
    # grid and the generator both live on the log scale.
    for (k in seq_along(inp$cens_x)) {
      di[[inp$cens_x[k]]] <- log(di[[inp$expo_names[k]]])
    }
    dr <- fit_bkmr_draws(di, grid, iter = iter, n_knots = n_knots,
                         seed = if (is.null(seed)) NULL else seed + 7000L + i)
    if (is.null(dr)) next
    mo <- bkmr_draws_moments(dr)
    ests[i, mo$estimand] <- mo$est
    vars[i, mo$estimand] <- mo$var
  }

  n_ok_fits <- sum(is.finite(ests[, 1]))
  if (n_ok_fits == 0L) {
    return(bkmr_na_rows(label, grid, "all per-imputation bkmr fits failed"))
  }

  tab <- do.call(rbind, lapply(grid$names, function(nm) {
    pl <- rubin_pool(ests[, nm], vars[, nm])
    data.frame(estimand = nm, est = unname(pl["est"]), se = unname(pl["se"]),
               ci_lo = unname(pl["ci_lo"]), ci_hi = unname(pl["ci_hi"]),
               ubar = unname(pl["ubar"]), b = unname(pl["b"]),
               fmi = unname(pl["fmi"]), stringsAsFactors = FALSE)
  }))

  note <- sprintf("m=%d; iter=%d; sweeps=%d; %s; cens=%d", m, iter,
                  outer_sweeps, margin, length(inp$cens_x))
  if (n_ok_fits < length(imputed))
    note <- paste0(note, sprintf("; %d/%d bkmr fits failed",
                                 length(imputed) - n_ok_fits, length(imputed)))
  bkmr_rows(label, tab, substr(note, 1, 120))
}
