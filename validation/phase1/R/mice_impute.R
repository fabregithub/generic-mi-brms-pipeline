# =============================================================================
# Phase 1 harness -- a `mice` Z-block, for Track 05
# -----------------------------------------------------------------------------
# WHY THIS EXISTS. V5 shipped a proper-MI fix for the Z block (bootstrap the
# training data per imputation, keep the random forest). It closed 66% of the
# calibration shortfall at 40% covariate missingness but only 23% under heavy
# combined missingness, because a forest is already bagged and an outer bootstrap
# shifts it only slightly.
#
# The obvious alternative -- `mice`'s parametric proper draws (`pmm`, `logreg`),
# which sample model parameters from a posterior -- calibrated markedly better in
# V4/V5. But that evidence came from a DGP whose covariates are all linear and
# Gaussian, so the parametric model was CORRECTLY SPECIFIED and was never
# penalised for the one risk it carries. Choosing it on that evidence would be
# choosing on a rigged test.
#
# This file supplies the `mice` arm. Run it against `z_form = "nonlinear"`
# (see R/dgp.R) and the two families finally compete on equal terms: the forest
# can capture structure the linear model cannot, and the linear model has proper
# parameter draws the forest lacks. Whichever wins, wins on merit.
#
# This is a harness INSTRUMENT, not pipeline code -- nothing here is adopted
# unless the comparison says so.
# =============================================================================

#' `mice`-based replacement for `run_row_level_imputation()`.
#'
#' Same signature and return shape (a list of `m` completed data frames), so it
#' can be assigned over the pipeline's helper inside the sourced environment.
#'
#' Uses `pmm` for continuous targets and `logreg` for binary ones -- both draw
#' imputation-model parameters from their posterior, which is what makes them
#' proper. The predictor matrix is built from the same `vars` list the pipeline
#' hands to miceRanger, so the two arms condition on exactly the same variables
#' and only the imputation *model* differs.
.v5_mice_row_imputation <- function(data, imputation_spec, analysis_spec) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("The mice arm requires the 'mice' package.", call. = FALSE)
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
  work <- as.data.frame(data)

  # mice's logreg needs a factor; the harness carries binary covariates as
  # numeric 0/1. Convert those targets on the way in and back on the way out, so
  # the completed datasets match what every other arm returns.
  binary_targets <- targets[vapply(targets, function(v) .v4_is_binary(work[[v]]), logical(1))]
  for (v in binary_targets) work[[v]] <- factor(work[[v]], levels = c(0, 1))

  nm <- names(work)
  pred <- matrix(0L, length(nm), length(nm), dimnames = list(nm, nm))
  for (v in targets) pred[v, intersect(vars[[v]], nm)] <- 1L

  meth <- stats::setNames(rep("", length(nm)), nm)
  for (v in targets) {
    meth[v] <- if (v %in% binary_targets) "logreg" else "pmm"
  }

  imp <- tryCatch(
    mice::mice(
      work, m = m, method = meth, predictorMatrix = pred,
      maxit = imputation_spec$maxiter %||% 5L,
      seed = imputation_spec$seed %||% NA_integer_,
      printFlag = FALSE
    ),
    error = function(e) NULL,
    warning = function(w) NULL
  )

  if (is.null(imp)) {
    stop("mice() failed in the Z block; no fallback is attempted so the failure ",
         "is visible rather than silently improper.", call. = FALSE)
  }

  completed <- mice::complete(imp, action = "all")

  lapply(completed, function(d) {
    for (v in binary_targets) d[[v]] <- as.numeric(as.character(d[[v]]))
    tibble::as_tibble(d)
  })
}
