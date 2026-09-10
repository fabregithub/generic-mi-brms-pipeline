# =============================================================================
# The DAG-factorised censored-exposure draw (THEORY.md §6b), for all four cases
# -----------------------------------------------------------------------------
# TARGET, from Bayes plus the DAG's factorisation:
#
#   p(x | rest, x<=L)  ∝  p(x | Pa(X))                          parents
#                       * p(Y | x, Pa(Y))                         the outcome child
#                       * PROD_{C in Ch(X)} p(C | x, ...)        other children
#                       * 1{x <= L}
#
# THREE THINGS THAT ARE EASY TO GET WRONG, and are encoded here deliberately:
#
# 1. THE OUTCOME FACTOR IS THE TRUE CONDITIONAL, NOT THE ANALYSIS MODEL. For a
#    total-effect estimand the analyst fits Y ~ X, but p(Y | x) is not the truth
#    -- the truth conditions on Z. Using the analyst's reduced model here would
#    put a misspecified factor in the target. So the algorithm fits the RICHEST
#    outcome model for the imputation and leaves the estimand to the analysis
#    step. §1b(i) measured that this is safe (+0.05% for the total effect) and
#    costs only mildly conservative intervals.
#
# 2. `Z` IS NEEDED HERE EVEN WHEN IT IS ABSENT FROM THE ANALYSIS MODEL. Z enters
#    the target because it is observed and informative about x -- as a parent in
#    a fork, as a child in a pipe or collider. Dropping it for a total-effect
#    analysis does not remove the non-linearity, it moves it from p(Z|x) into
#    p(Y|x) (§6b). Measured: +16.9% with Z in, +14.3% with it out, -0.3% with the
#    child factor on.
#
# 3. THE OUTCOME FACTOR CONDITIONS ON Pa(Y), NOT ON "ALL COVARIATES". Under a
#    collider Z1 is a CHILD of Y, so putting it in the outcome factor would be
#    conditioning on a collider inside the imputation model -- the fitted
#    coefficient is badly biased even though the structural one is zero, and
#    Z1's information about x would be double-counted (it already enters through
#    p(Z1 | x, Y)). So Z1 goes in the outcome factor exactly when it is a parent
#    of Y (fork, pipe, mixed) and never when it is a child of Y (collider).
#
# THE CHILD MODEL IS AN ASSUMPTION, NOT AN ESTIMATE. p(C | x) must be evaluated
# BELOW the LOD, where x is unobserved, and THEORY.md §0b measures that this is
# not identifiable -- a natural spline returns exactly zero curvature by
# construction, a quadratic is 4.6x low when the truth is not quadratic, a cubic
# gives false positives. So `child_form` is a declared functional form and the
# honest use of this code is a sensitivity analysis across it.
# =============================================================================

#' Which role each covariate plays, read from the DGP rather than passed in.
dag_spec <- function(truth) {
  role <- truth$z_role %||% "precision"
  list(role = role,
       # Z1's relationship to X1: "parent", "child", "child_xy", or "none".
       z1 = switch(role,
                   fork = "parent", mixed = "parent",
                   pipe = "child", pipe_nl = "child",
                   collider = "child_xy",
                   "none"),
       # Is Z1 a PARENT of Y? It is in a fork/pipe/mixed; in a collider it is a
       # child of Y, so it must stay out of the outcome factor (see (3) above).
       z1_in_y = !identical(role, "collider"),
       nl = identical(role, "pipe_nl"))
}

#' Fit the child model `p(Z1 | X1, .)` on the OBSERVED rows and return an
#' evaluator usable BELOW the LOD.
#'
#' Returns `mu(x, y)` and `sd`. `y` is used only when Z1 is a child of Y as well
#' (the collider), where the child factor is `p(Z1 | x, Y)`.
#'
#' @param child_form "true" uses the generator's own arrow -- an ORACLE, for
#'   testing whether the FACTORISATION is right, separately from whether the
#'   child model can be estimated. "linear" is the shippable control: it is the
#'   form the current conditional already implies, so a cell where "linear" is
#'   enough is a cell that needs no new assumption. "quad"/"cubic" are the
#'   sensitivity axis, and THEORY.md 0b measures why they can only ever be that:
#'   below the LOD the curvature of the arrow is NOT identified.
.dag_child_fit <- function(w, obs, truth, spec, child_form) {
  if (identical(child_form, "true")) {
    # Role-aware, because the generator's arrow differs by role -- using
    # ef_nl_g() for a linear `pipe` would make the "oracle" arm wrong and it
    # would still return believable numbers.
    if (spec$nl)
      return(list(mu = function(x, y = 0) ef_nl_g(x, truth), sd = truth$nl_sd))
    return(list(mu = function(x, y = 0)
                  truth$delta_xz * x +
                  (if (spec$z1 == "child_xy") truth$delta_yz * y else 0),
                sd = 0.80))
  }
  k <- switch(child_form, linear = 1L, quad = 2L, cubic = 3L,
              stop("unknown child_form: ", child_form, call. = FALSE))
  extra <- if (spec$z1 == "child_xy") "Y" else character(0)
  rhs <- paste0("poly(logX1, ", k, ", raw = TRUE)")
  if (length(extra)) rhs <- paste(c(rhs, extra), collapse = " + ")
  f <- stats::lm(stats::as.formula(paste("Z1 ~", rhs)),
                 data = w[obs, , drop = FALSE])
  cf <- stats::coef(f); cf[is.na(cf)] <- 0
  bx <- cf[grepl("^poly", names(cf))]
  by <- if (length(extra)) cf[["Y"]] else 0
  # `x` arrives as a MATRIX (rows = censored cells, cols = grid points), so the
  # polynomial is built by accumulation rather than through model.matrix().
  list(mu = function(x, y = 0) {
         out <- cf[[1L]] + by * y                 # y recycles down the columns
         for (j in seq_along(bx)) out <- out + bx[[j]] * x^j
         out
       },
       sd = stats::summary.lm(f)$sigma)
}

#' One DAG-factorised draw of the censored logX1 values.
dag_draw_x1 <- function(w, truth, idx, lod, spec, child, coefs, sig, obs_na,
                        use_y = TRUE, use_child = TRUE, ngrid = 400L) {
  if (!length(idx)) return(w$logX1)
  p  <- length(truth$b); xc <- paste0("logX", seq_len(p))
  zc <- intersect(c("Z1", "Z2"), names(w))
  sub <- w[idx, , drop = FALSE]

  # --- p(x | Pa(X)) ----------------------------------------------------------
  # The other exposures always; Z1 too when it is a PARENT (fork/mixed) and
  # never when it is a child, which is the whole point of the factorisation.
  #
  # FIT INTERVAL-CENSORED, NOT ON THE OBSERVED ROWS. An lm() on the above-LOD
  # rows is a regression on a response-truncated sample, which biases the slope
  # and shrinks sigma -- it would inject a defect into the PARENT factor and the
  # test would then be scoring the wrong thing. survreg's interval-censored
  # Gaussian uses every row and is the same fit `leftcens` performs internally.
  pa <- c(xc[-1], if (spec$z1 == "parent") "Z1")
  lo <- w$logX1; hi <- w$logX1
  lo[obs_na] <- -Inf; hi[obs_na] <- lod[obs_na]
  f_pa <- survival::survreg(
    stats::as.formula(paste("survival::Surv(lo, hi, type = \"interval2\") ~",
                            paste(pa, collapse = "+"))),
    data = cbind(w, lo = lo, hi = hi), dist = "gaussian")
  m_pa <- stats::predict(f_pa, newdata = sub, type = "response")
  s_pa <- f_pa$scale

  # --- where to put the grid ------------------------------------------------
  # A grid of L - 6 s_pa .. L is NOT enough. Both the parent and the outcome
  # factor are Gaussian in x, so their product is Gaussian with a known mean and
  # SD -- and for a cell with an extreme Y that mean can sit near the bottom of a
  # fixed window, which clips the target and biases the draw UPWARD while still
  # looking like a sensible truncated distribution. So the window is anchored on
  # that combined Gaussian, per cell, and only the child factor (the one
  # non-Gaussian piece, and a bounded one) is left outside the anchoring.
  zy <- if (use_y) intersect(zc, names(coefs)) else character(0)
  if (use_y) {
    eta0 <- coefs[["(Intercept)"]] +
      as.vector(as.matrix(sub[, xc[-1], drop = FALSE]) %*% coefs[xc[-1]]) +
      (if (length(zy)) as.vector(as.matrix(sub[, zy, drop = FALSE]) %*% coefs[zy]) else 0)
    b1  <- coefs[["logX1"]]
    prec <- 1 / s_pa^2 + b1^2 / sig^2
    m_c  <- (m_pa / s_pa^2 + b1 * (sub$Y - eta0) / sig^2) / prec
    s_c  <- sqrt(1 / prec)
  } else {
    m_c <- m_pa; s_c <- s_pa
  }
  hi_g <- lod[idx]
  lo_g <- pmin(m_c - 8 * s_c, hi_g - 8 * s_c)
  G <- lo_g + outer(hi_g - lo_g, seq(0, 1, length.out = ngrid))
  lw <- stats::dnorm(G, m_pa, s_pa, log = TRUE)

  # --- p(Y | x, Pa(Y); theta) : the TRUE conditional, not the analysis model.
  # `coefs` comes from a fit of Y on (all X, Pa(Y) covariates) -- see (1), (3).
  # `eta0` and `b1` were built above, where they set the grid.
  if (use_y) lw <- lw + stats::dnorm(sub$Y, eta0 + b1 * G, sig, log = TRUE)

  # --- the other children: p(Z1 | x, ...) -----------------------------------
  if (use_child && spec$z1 %in% c("child", "child_xy")) {
    # For a collider Z1 is a child of BOTH x and Y, so Y enters its mean; mu()
    # handles that itself and ignores `y` in every other role.
    lw <- lw + stats::dnorm(sub$Z1, child$mu(G, sub$Y), child$sd, log = TRUE)
  }

  lw <- lw - apply(lw, 1, max); W <- exp(lw)
  rs <- rowSums(W)
  cdf <- t(apply(W, 1, cumsum)); cdf <- cdf / cdf[, ngrid]
  u <- stats::runif(length(idx))
  out <- w$logX1
  out[idx] <- vapply(seq_along(idx), function(i)
    stats::approx(cdf[i, ], G[i, ], xout = u[i], rule = 2, ties = "ordered")$y, 0)
  # The grid spans (L - 6 s_pa, L]. Once Y and the children enter, the target is
  # NARROWER than the parent factor, so 6 s_pa is generous -- but "generous" is
  # an argument, not a measurement, and a grid that clipped the target would bias
  # the draw upward while still looking like a sensible truncated distribution.
  # Reported so test_v24_dag.R can assert on it instead of re-deriving the
  # weights, which is how two earlier test constructions went wrong.
  attr(out, "edge_weight") <- max(W[, 1] / rs)
  out
}

#' m completed datasets from the DAG-factorised draw.
dag_impute_datasets <- function(d, truth, m = 30L, sweeps = 3L,
                                child_form = "true", use_y = TRUE,
                                use_child = TRUE, ngrid = 400L) {
  spec <- dag_spec(truth)
  p <- length(truth$b); xc <- paste0("logX", seq_len(p))
  zc <- intersect(c("Z1", "Z2"), names(d))
  isc <- is.na(d$logX1); lod <- d$logX1_lod
  # MISSING COVARIATES ARE OUT OF SCOPE, AND MUST SAY SO. Every factor in the
  # target conditions on the covariates, so a missing Z1 makes the parent fit,
  # the outcome factor and the child factor all NA -- which surfaced as
  # "need at least two non-NA values to interpolate" from approx(), three calls
  # away from the cause. V24 isolates the EXPOSURE draw; the covariate block was
  # settled by V20/V21 and combining the two is a separate cell (ROADMAP).
  zmiss <- vapply(intersect(c("Z1", "Z2"), names(d)),
                  function(v) anyNA(d[[v]]), TRUE)
  if (any(zmiss))
    stop("dag_impute_datasets() needs fully observed covariates; missing: ",
         paste(names(zmiss)[zmiss], collapse = ", "),
         ". V24 isolates the exposure draw -- run these cells with mcar_frac = 0.",
         call. = FALSE)
  # The richest CORRECT outcome model -- all exposures plus the covariates that
  # are parents of Y -- used for the target regardless of the estimand.
  zy <- if (spec$z1_in_y) zc else setdiff(zc, "Z1")
  fo_y <- stats::as.formula(paste("Y ~", paste(c(xc, zy), collapse = "+")))

  out <- vector("list", m); edge <- 0
  for (i in seq_len(m)) {
    w <- d; w$logX1[isc] <- lod[isc] - 0.5
    for (t in seq_len(sweeps)) {
      fy <- stats::lm(fo_y, data = w)
      df <- fy$df.residual
      s2 <- sum(stats::residuals(fy)^2) / stats::rchisq(1, df)
      V  <- stats::summary.lm(fy)$cov.unscaled
      cf <- as.vector(stats::coef(fy) +
                        t(chol(V)) %*% stats::rnorm(length(stats::coef(fy))) * sqrt(s2))
      names(cf) <- names(stats::coef(fy))
      ch <- .dag_child_fit(w, !isc, truth, spec, child_form)
      dr <- dag_draw_x1(w, truth, which(isc), lod, spec, ch, cf, sqrt(s2),
                        obs_na = isc, use_y = use_y,
                        use_child = use_child, ngrid = ngrid)
      edge <- max(edge, attr(dr, "edge_weight") %||% 0)
      w$logX1 <- as.vector(dr)
    }
    out[[i]] <- w
  }
  attr(out, "edge_weight") <- edge
  out
}

# =============================================================================
# V24 procedure arms
# -----------------------------------------------------------------------------
# These are HARNESS INSTRUMENTS, not pipeline code -- exactly as exact_fork.R's
# arms are. They isolate the exposure draw: the covariate block is left at the
# exact conditional so that the only thing varying across arms is which factors
# enter p(x | rest, x <= L). `dag_noY` and `dag_Yonly` are the two draws that
# already exist in the literature and in this pipeline respectively, so the
# ladder is: no Y -> Y -> Y + children.
# =============================================================================

#' One V24 arm: run the DAG-factorised draw and pool.
#'
#' Pools through `fit_lm_estimand()` / `rubin_pool()` / `one_row()` -- the same
#' helpers every other arm uses -- so a V24 row is comparable to a V20/V21 row
#' and the pooling itself is not a new thing that could differ. `fit_lm_estimand`
#' fits `dgp_formula(truth)`, which under `target = "total"` has already dropped
#' Z1, so the estimand follows the scenario with nothing extra here.
proc_dag_arm <- function(bundle, m = 30L, seed = NULL, sweeps = 3L,
                         child_form = "true", use_y = TRUE, use_child = TRUE,
                         label = NULL, ...) {
  label <- label %||% paste0("dag_", if (!use_y) "noY" else
                             if (!use_child) "Yonly" else child_form)
  if (!is.null(seed)) set.seed(seed)
  truth <- bundle$truth
  imp <- tryCatch(
    dag_impute_datasets(bundle$censored, truth, m = as.integer(m),
                        sweeps = as.integer(sweeps), child_form = child_form,
                        use_y = use_y, use_child = use_child),
    error = function(e) structure(list(), err = conditionMessage(e)))
  if (!length(imp))
    return(one_row(label, NA, NA, NA, NA,
                   substr(paste("dag failed:", attr(imp, "err") %||% "none"), 1, 120)))
  es <- vs <- rep(NA_real_, length(imp))
  for (i in seq_along(imp)) {
    e <- fit_lm_estimand(as.data.frame(imp[[i]]), truth)
    es[i] <- e["est"]; vs[i] <- e["se"]^2
  }
  pl <- rubin_pool(es, vs)
  # `edge` is the grid-coverage diagnostic (see dag_draw_x1). It belongs in the
  # note rather than in a silent assumption: a run whose edge weight is not
  # negligible has a grid-truncated target, and its bias is an artefact.
  one_row(label, pl["est"], pl["se"], pl["ci_lo"], pl["ci_hi"],
          sprintf("child=%s; y=%d; ch=%d; m=%d; sweeps=%d; edge=%.1e",
                  child_form, as.integer(use_y), as.integer(use_child), m,
                  sweeps, attr(imp, "edge_weight") %||% NA_real_),
          ubar = pl["ubar"], b = pl["b"], fmi = pl["fmi"])
}

proc_dag_noY    <- function(b, ...) proc_dag_arm(b, use_y = FALSE, use_child = FALSE, ...)
proc_dag_Yonly  <- function(b, ...) proc_dag_arm(b, use_y = TRUE,  use_child = FALSE, ...)
# The labels must match the arm names registered in run_v4_variance.R's V4_ARMS
# exactly -- summarise_v4() factors `procedure` against that list, so a mismatch
# does not error, it silently produces an all-NA row named <NA>.
proc_dag_lin    <- function(b, ...) proc_dag_arm(b, child_form = "linear", label = "dag_lin",   ...)
proc_dag_quad   <- function(b, ...) proc_dag_arm(b, child_form = "quad",   label = "dag_quad",  ...)
proc_dag_cubic  <- function(b, ...) proc_dag_arm(b, child_form = "cubic",  label = "dag_cubic", ...)
proc_dag_true   <- function(b, ...) proc_dag_arm(b, child_form = "true",   label = "dag_true",  ...)
