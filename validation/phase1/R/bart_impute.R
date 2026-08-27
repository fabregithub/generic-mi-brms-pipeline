# =============================================================================
# Phase 1 harness -- a BART Z-block, for requirement R8 / roadmap item 06
# -----------------------------------------------------------------------------
# WHY THIS EXISTS. R8 (proper multiple imputation) is only partially satisfied.
# Two arms have been tested and neither is both flexible and properly dispersed:
#
#   properBoot  bootstrap on a random forest -- flexible, but a forest is ALREADY
#               bagged, so an outer bootstrap shifts it only slightly. Under-
#               disperses: +2.0-2.6% bias in the `combined` cells, width/SE 0.921.
#   micePmm     parametric posterior draw -- properly dispersed, but LINEAR, so
#               misspecified under non-linear covariates: bias -4.70%.
#
# BART is the natural resolution: a sum-of-trees model with priors on both the
# tree structure and the leaf parameters, so it is non-parametric like a forest
# AND fully Bayesian. A posterior sample IS a proper draw -- no bootstrap
# approximation, no linearity assumption.
#
# THE DRAW. dbarts::bart() returns `yhat.test` as an (ndpost x n_missing) matrix
# of posterior draws of the conditional mean, and `sigma` as the matching
# residual-SD draws. A posterior PREDICTIVE draw is therefore
#
#     yhat.test[j, ] + rnorm(n_missing, 0, sigma[j])
#
# which is proper by construction. Binary targets are fitted as factors, giving
# the probit form: probabilities are pnorm(yhat.test[j, ]) and the draw is a
# Bernoulli at those probabilities.
#
# COST HONESTY. One BART fit yields m draws only when there is a SINGLE target
# whose predictors are complete -- true for the `missing_y20` cells (target Y),
# false for `mcar_z40` (Z1, Z2) and `combined` (Z1, Z2, Y), where fully
# conditional specification requires m separate chains because each chain's
# predictor values differ. The single-target fast path is taken where it is
# valid; otherwise m chains are run. No approximation is made to save time.
#
# This is a harness INSTRUMENT under test. Nothing is adopted unless the
# pre-registered criteria in ROADMAP.md item 06 are met.
# =============================================================================

# Tuning. 50 trees keeps BART flexible while staying fast; 100 burn-in
# iterations is dbarts' own default territory for problems this size. Exposed as
# options so a run can vary them without editing this file.
.v7_bart_ntree <- function() getOption("v7.bart.ntree", 50L)
.v7_bart_nskip <- function() getOption("v7.bart.nskip", 100L)

#' One posterior predictive draw for a continuous target.
#'
#' @param ndraw How many posterior draws to return (1 for an FCS step, m for the
#'   single-target fast path).
.v7_bart_draw <- function(y_obs, x_obs, x_mis, ndraw = 1L, binary = FALSE) {
  fit <- tryCatch(
    dbarts::bart(
      x.train = x_obs,
      y.train = if (binary) factor(y_obs, levels = c(0, 1)) else y_obs,
      x.test  = x_mis,
      ntree   = .v7_bart_ntree(),
      nskip   = .v7_bart_nskip(),
      ndpost  = as.integer(ndraw),
      keeptrees = FALSE,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$yhat.test)) return(NULL)

  yh <- fit$yhat.test
  if (!is.matrix(yh)) yh <- matrix(yh, nrow = 1L)
  n_mis <- ncol(yh)

  if (binary) {
    # Probit scale: convert to probabilities, then draw.
    lapply(seq_len(nrow(yh)), function(j) {
      stats::rbinom(n_mis, 1L, stats::pnorm(yh[j, ]))
    })
  } else {
    # `sigma` length varies across dbarts versions (ndpost, or nskip + ndpost);
    # take the trailing ndpost entries so draw j always pairs with yhat row j.
    sg <- utils::tail(fit$sigma, nrow(yh))
    lapply(seq_len(nrow(yh)), function(j) {
      yh[j, ] + stats::rnorm(n_mis, 0, sg[j])
    })
  }
}

#' BART replacement for `run_row_level_imputation()`.
#'
#' Same signature and return shape (a list of `m` completed data frames), so it
#' can be assigned over the pipeline's helper inside the sourced environment.
.v7_bart_row_imputation <- function(data, imputation_spec, analysis_spec,
                                    inner_iter = 3L) {
  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop("The BART arm requires the 'dbarts' package.", call. = FALSE)
  }

  vars <- imputation_spec$vars
  m <- as.integer(imputation_spec$m %||% 1L)

  vars <- vars[names(vars) %in% names(data)]
  vars <- purrr::map(vars, ~ intersect(.x, names(data)))
  vars <- vars[vapply(vars, length, integer(1)) > 0]

  if (length(vars) == 0L) {
    return(rep(list(tibble::as_tibble(data)), m))
  }

  targets <- names(vars)
  work0 <- as.data.frame(data)
  na_map <- lapply(work0[targets], is.na)
  names(na_map) <- targets
  is_bin <- vapply(targets, function(v) .v4_is_binary(work0[[v]]), logical(1))
  names(is_bin) <- targets

  design <- function(d, preds) {
    preds <- preds[vapply(d[preds], is.numeric, logical(1))]
    as.data.frame(d[, preds, drop = FALSE])
  }

  # ---- fast path: one target whose predictors are already complete ----------
  # One fit gives all m draws, which is both cheaper and cleaner -- there is no
  # chain to iterate because nothing else is being imputed.
  if (length(targets) == 1L) {
    v <- targets[1L]
    mis <- na_map[[v]]
    preds <- setdiff(intersect(vars[[v]], names(work0)), v)
    X <- design(work0, preds)

    if (any(mis) && !anyNA(X)) {
      draws <- .v7_bart_draw(work0[[v]][!mis], X[!mis, , drop = FALSE],
                             X[mis, , drop = FALSE], ndraw = m,
                             binary = is_bin[[v]])
      if (!is.null(draws)) {
        return(lapply(draws, function(dj) {
          d <- work0; d[[v]][mis] <- dj; tibble::as_tibble(d)
        }))
      }
    }
  }

  # ---- general path: m independent FCS chains -------------------------------
  out <- vector("list", m)
  for (j in seq_len(m)) {
    if (!is.null(imputation_spec$seed)) set.seed(imputation_spec$seed + j)

    d <- work0
    for (v in targets) d[[v]] <- .v4_init_fill(d[[v]])

    for (it in seq_len(inner_iter)) {
      for (v in targets) {
        mis <- na_map[[v]]
        if (!any(mis)) next
        preds <- setdiff(intersect(vars[[v]], names(d)), v)
        X <- design(d, preds)
        if (!ncol(X) || anyNA(X)) next

        dr <- .v7_bart_draw(d[[v]][!mis], X[!mis, , drop = FALSE],
                            X[mis, , drop = FALSE], ndraw = 1L,
                            binary = is_bin[[v]])
        # A failed fit is fatal rather than silently falling back to something
        # improper -- the whole point of this arm is that every draw is proper.
        if (is.null(dr)) {
          stop("BART fit failed for target '", v, "' (chain ", j, ").",
               call. = FALSE)
        }
        d[[v]][mis] <- dr[[1L]]
      }
    }
    out[[j]] <- tibble::as_tibble(d)
  }

  out
}
