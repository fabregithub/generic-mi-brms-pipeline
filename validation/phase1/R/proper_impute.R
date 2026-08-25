# =============================================================================
# Phase 1 harness -- a PROPER-MI covariate imputer, for Track V4
# -----------------------------------------------------------------------------
# Track V4 of ../PLAN_pipeline_validation.md §9.
#
# V2 found that the pipeline's 95% intervals are ~17% too narrow when covariate
# and outcome missingness are both substantial, and that raising `m` does not fix
# it (so it is a systematic shortfall in the variance estimator, not MC noise).
# Two mechanisms remain in play; this file isolates the first one.
#
# THE TEST. The pipeline's Z block uses miceRanger (random forest + predictive
# mean matching), which does NOT draw imputation-model parameters from a
# posterior. That is "improper" MI in Rubin's sense and is known to understate
# between-imputation variance. This file provides a drop-in replacement that IS
# proper -- same call signature as `run_row_level_imputation()`, so it can be
# swapped into the pipeline's private environment while EVERYTHING else (the X
# block, the outer alternation, MID, the pooling) stays byte-identical.
#
# If swapping it in restores width/SE to ~1.0, the Z-block imputer is implicated.
# If it does not, the imputer is exonerated and the cause lies elsewhere.
#
# This is a validation instrument, NOT a proposed replacement for miceRanger: it
# is deliberately simple (linear / logistic, no interactions, no non-linearity),
# which is fine here because the DGP's covariates are genuinely linear-Gaussian
# and binary. It would be a poor general-purpose imputer.
# =============================================================================

#' Is this column effectively binary (0/1)?
.v4_is_binary <- function(x) {
  u <- unique(x[!is.na(x)])
  length(u) <= 2L && all(u %in% c(0, 1))
}

#' Mean/mode fill, used only to initialise the FCS loop.
.v4_init_fill <- function(x) {
  na <- is.na(x)
  if (!any(na)) return(x)
  if (.v4_is_binary(x)) {
    tb <- table(x[!na])
    x[na] <- as.numeric(names(tb)[which.max(tb)])
  } else {
    x[na] <- mean(x[!na])
  }
  x
}

#' Proper Bayesian draw for one continuous target: draw sigma^2 from its scaled
#' inverse-chi-square posterior, beta from its conditional normal posterior, then
#' add residual noise. This is the `norm` method of classical MI -- the textbook
#' proper draw that PMM/random-forest imputation approximates without the
#' parameter uncertainty.
.v4_draw_norm <- function(y, X, obs, mis) {
  Xo <- X[obs, , drop = FALSE]; yo <- y[obs]
  qrX <- qr(Xo)
  if (qrX$rank < ncol(Xo)) return(NULL)                 # singular: caller falls back
  beta_hat <- qr.coef(qrX, yo)
  res <- yo - Xo %*% beta_hat
  dfree <- length(yo) - ncol(Xo)
  if (dfree < 1L) return(NULL)
  s2 <- sum(res^2) / dfree
  # sigma^2 | data  ~  s2 * dfree / chisq(dfree)
  sigma2 <- s2 * dfree / stats::rchisq(1L, dfree)
  V <- chol2inv(qr.R(qrX))                              # (X'X)^-1
  Rv <- tryCatch(chol(sigma2 * V), error = function(e) NULL)
  if (is.null(Rv)) return(NULL)
  beta <- as.vector(beta_hat + t(Rv) %*% stats::rnorm(ncol(Xo)))
  Xm <- X[mis, , drop = FALSE]
  as.vector(Xm %*% beta) + stats::rnorm(sum(mis), 0, sqrt(sigma2))
}

#' Proper draw for one binary target: logistic fit, beta drawn from its
#' asymptotic normal posterior, then a Bernoulli draw at the drawn probability.
.v4_draw_logit <- function(y, X, obs, mis) {
  Xo <- X[obs, , drop = FALSE]; yo <- y[obs]
  fit <- tryCatch(
    stats::glm.fit(Xo, yo, family = stats::binomial()),
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(fit) || any(!is.finite(fit$coefficients))) return(NULL)
  beta_hat <- fit$coefficients
  W <- fit$weights
  XtWX <- crossprod(Xo * sqrt(W))
  V <- tryCatch(chol2inv(chol(XtWX)), error = function(e) NULL)
  if (is.null(V)) return(NULL)
  Rv <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(Rv)) return(NULL)
  beta <- as.vector(beta_hat + t(Rv) %*% stats::rnorm(length(beta_hat)))
  pm <- stats::plogis(as.vector(X[mis, , drop = FALSE] %*% beta))
  stats::rbinom(sum(mis), 1L, pm)
}

#' Proper-MI replacement for `run_row_level_imputation()`.
#'
#' Same signature and same return shape (a list of `m` completed data frames), so
#' it can be assigned over the pipeline's own helper inside the sourced
#' environment. Runs a short fully-conditional-specification loop: mean/mode
#' initialise, then repeatedly re-impute each target from a proper Bayesian draw
#' given the current values of the others.
#'
#' @param data,imputation_spec,analysis_spec As the pipeline helper.
#' @param inner_iter FCS sweeps within one completed dataset.
.v4_proper_row_imputation <- function(data, imputation_spec, analysis_spec,
                                      inner_iter = 5L) {
  vars <- imputation_spec$vars
  m <- as.integer(imputation_spec$m %||% 1L)

  if (length(vars) == 0L) {
    return(rep(list(tibble::as_tibble(data)), m))
  }

  targets <- names(vars)
  na_map <- lapply(data[targets], is.na)
  names(na_map) <- targets

  out <- vector("list", m)
  for (k in seq_len(m)) {
    d <- data
    for (v in targets) d[[v]] <- .v4_init_fill(d[[v]])   # initialise

    for (it in seq_len(inner_iter)) {
      for (v in targets) {
        mis <- na_map[[v]]
        if (!any(mis)) next
        preds <- intersect(vars[[v]], names(d))
        preds <- preds[vapply(d[preds], is.numeric, logical(1))]
        if (!length(preds)) next

        X <- cbind(1, as.matrix(d[, preds, drop = FALSE]))
        if (anyNA(X)) next                               # predictor still missing
        y <- d[[v]]
        obs <- !mis

        imp <- if (.v4_is_binary(data[[v]])) {
          .v4_draw_logit(y, X, obs, mis)
        } else {
          .v4_draw_norm(y, X, obs, mis)
        }
        # A singular/failed draw leaves the current values in place rather than
        # silently producing garbage; that is visible as reduced between-variance
        # only if it happens often, which the caller can check.
        if (!is.null(imp) && all(is.finite(imp))) d[[v]][mis] <- imp
      }
    }
    out[[k]] <- tibble::as_tibble(d)
  }
  out
}
