# =============================================================================
# Track V36 -- outcome-dependent selection
# -----------------------------------------------------------------------------
# Every track before this one varies which COVARIATES are conditioned on
# (`z_role`). None varies which ROWS survive. This file adds that, and nothing
# else: selection is applied to the bundle before any arm sees it, so no
# procedure needs to know selection exists.
#
# WHAT IS AND IS NOT MODELLED HERE. The selection variable is a DESCENDANT OF
# THE EXPOSURE that shares an unmeasured cause `U` with the outcome:
#
#     logX1 -> M <- U -> Y          (select on M)
#
# Selecting on `M` opens the collider path logX1 -> M <- U -> Y, which is the
# mechanism a low-birth-weight exclusion triggers. `M` is deliberately NOT given
# an arrow into `Y`: if it had one it would be a mediator proper, the analysis
# model's `logX1` coefficient would become a total effect, and `truth$b[1]`
# would stop being the estimand. **So this prices collider stratification from
# selecting on a descendant of the exposure. It does NOT additionally price the
# blocked-mediation component** -- that changes the estimand and needs its own
# design. Say so in the findings; do not let the shorthand "mediator exclusion"
# imply both.
#
# `U` is added to `Y` after generation. That is legitimate because `U` is drawn
# independently of `(logX, Z)`, so it enters as extra residual variance and
# moves no coefficient -- `truth$estimand_true` remains the unselected
# population's `b[1]`, which DESIGN_v36_selection.md §2 requires as the target.
# =============================================================================

#' Selection cell specifications. `drop` is the fraction of rows removed;
#' 0.092 matches Itoh et al.'s 1,534 of 16,709.
v36_sel_spec <- function(name) {
  base <- list(name = name, drop = 0.092, alpha = 0.60, delta = 0.80,
               u_on_y = 0.60, sd_noise = 0.50)
  switch(name,
    none      = utils::modifyList(base, list(kind = "none")),
    # selection on a descendant of logX1 that shares U with Y
    desc      = utils::modifyList(base, list(kind = "desc")),
    # STRENGTH LEVER IS THE DROP FRACTION, NOT `delta`. Collider bias here is
    # NON-MONOTONE in `delta`: at delta = 0 selection is on a pure descendant of
    # the exposure (no bias -- see `descnull`), and as delta grows selection
    # tends to selection on `U` alone, which is independent of the exposure, so
    # the bias vanishes again. It peaks in between. Measured 2026-11-09 at
    # n = 4000, 400 reps: delta 0.80 -> -17.5%, delta 1.60 -> -12.8%. Raising
    # `delta` to "strengthen" a cell makes the bias SMALLER and reads as a
    # refutation. The drop fraction is monotone and is also the quantity an
    # epidemiologist observes.
    desc2     = utils::modifyList(base, list(kind = "desc", drop = 0.20)),
    # FALSIFICATION CONTROL: same structure, but M has no link to Y at all.
    # Selecting on a pure descendant of the exposure must NOT bias b[1].
    descnull  = utils::modifyList(base, list(kind = "desc", delta = 0.00,
                                             u_on_y = 0.60)),
    # selection on the outcome itself
    outcome   = utils::modifyList(base, list(kind = "outcome", sd_noise = 1.00)),
    outcome2  = utils::modifyList(base, list(kind = "outcome", sd_noise = 1.00,
                                             drop = 0.20)),
    # live-birth style: a collider of the exposure and the outcome
    collider  = utils::modifyList(base, list(kind = "collider", alpha = 0.60,
                                             delta = 0.60)),
    stop("unknown selection cell: ", name, call. = FALSE))
}

#' Build a V36 bundle: generate, add `U` to `Y`, censor, then select.
#'
#' Order matters. The LOD is a quantile of the FULL pre-selection sample,
#' because a detection limit is a property of the assay, not of who was kept --
#' so the realised non-detect fraction inside the selected sample can differ
#' from `nd_frac`, and it is recorded as an attribute rather than forced back.
v36_make_bundle <- function(sc, truth) {
  sel <- v36_sel_spec(sc$sel %||% "none")
  n   <- sc$n

  comp <- simulate_complete(n, truth, rho = sc$rho %||% 0.4,
                            skew = sc$skew %||% 0.0,
                            sd_x = sc$sd_x %||% 1.0,
                            z_form = sc$z_form %||% "linear")
  d <- comp$data

  u <- stats::rnorm(n)
  if (!identical(sel$kind, "none")) d$Y <- d$Y + sel$u_on_y * u

  keep <- rep(TRUE, n)
  if (!identical(sel$kind, "none")) {
    s_var <- switch(sel$kind,
      desc     = sel$alpha * d$logX1 + sel$delta * u +
                 stats::rnorm(n, 0, sel$sd_noise),
      outcome  = d$Y + stats::rnorm(n, 0, sel$sd_noise),
      collider = sel$alpha * d$logX1 + sel$delta * d$Y +
                 stats::rnorm(n, 0, sel$sd_noise))
    keep <- s_var >= stats::quantile(s_var, probs = sel$drop, names = FALSE)
  }

  cens <- inject_left_censoring(d, nd_frac = sc$nd_frac,
                                censor_which = 1L)
  lods <- attr(cens, "lods")
  cens <- inject_mcar_outcome(cens, sc$y_frac %||% 0.0)

  # Subset both frames identically, and carry the censoring attributes through
  # -- `[.data.frame` drops them, and the pipeline arms read `lods`.
  comp_s <- d[keep, , drop = FALSE]
  cens_s <- cens[keep, , drop = FALSE]
  attr(cens_s, "lods") <- lods
  attr(cens_s, "censored_cols") <- attr(cens, "censored_cols")
  attr(cens_s, "nd_frac_realised") <- mean(is.na(cens_s$logX1))
  rownames(comp_s) <- NULL; rownames(cens_s) <- NULL

  list(complete = comp_s, censored = cens_s, truth = truth,
       sel = sel, n_kept = sum(keep), nd_realised = mean(is.na(cens_s$logX1)))
}

#' LOD/sqrt(2) substitution on the additive path -- Itoh et al.'s actual method.
#'
#' The BKMR path has `proc_bkmr_sub_lod2`; the additive path never needed one.
#' Kept here rather than added to `run_procedures()` so no shared dispatch
#' changes and no existing track's RNG stream or arm set moves.
proc_sub_lod2_lm <- function(bundle) {
  d <- bundle$censored
  xc <- "logX1"; lc <- "logX1_lod"
  if (!lc %in% names(d)) return(one_row("sub_lod2", NA_real_, NA_real_,
                                        NA_real_, NA_real_, "no lod column"))
  is_c <- is.na(d[[xc]])
  # The substitution is defined on the CONCENTRATION scale; these columns are
  # logs, so LOD/sqrt(2) is log(LOD) - log(sqrt(2)).
  d[[xc]][is_c] <- d[[lc]][is_c] - log(sqrt(2))
  fit <- fit_lm_estimand(d, bundle$truth)
  if (is.na(fit["est"])) return(one_row("sub_lod2", NA_real_, NA_real_,
                                        NA_real_, NA_real_, "fit failed"))
  ci <- fit["est"] + c(-1, 1) * stats::qnorm(0.975) * fit["se"]
  one_row("sub_lod2", unname(fit["est"]), unname(fit["se"]), ci[1], ci[2])
}
