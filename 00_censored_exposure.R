# =============================================================================
# Censored-exposure block-FCS imputation  (PLAN §5 / Phase 4)
# -----------------------------------------------------------------------------
# A congenial imputation path for the special case where a focal EXPOSURE
# (predictor of the outcome) is left-/interval-censored below a reporting limit.
# Imputing such an exposure without the outcome Y biases the exposure-response
# coefficient (PLAN §2); this path handles it with a two-engine block-FCS:
#
#   * Z block (MAR covariates)      -> miceRanger, via the EXISTING pipeline
#                                       helper run_row_level_imputation() (reused,
#                                       not reimplemented), conditioning on the
#                                       current X and Y.
#   * X block (censored exposure)   -> leftcens::impute_censored_conditional()
#                                       (>= 0.9.0): skew-aware, bound-respecting,
#                                       Y-aware, proper-MI conditional draw.
#   * alternate for a few outer sweeps, per completed dataset.
#
# It is opt-in (strategy = "censored_exposure_block_fcs") and returns an
# `imputed_list` in the SAME shape as the other strategies, so 03_impute.R writes
# the standard imputed_###.rds + manifest and Steps 4-12 are untouched.
#
# Input convention (leftcens interval columns): each censored exposure `X` carries
# per-row bounds `X_lo` / `X_hi` (suffixes configurable) of the interval holding
# the true value: X_lo == X_hi where exactly observed; X_lo = -Inf (or 0) with
# X_hi = LOD for a left-censored non-detect; both finite for an interval.
# =============================================================================

#' Build (y, lower, upper) for one exposure from its lo/hi columns, on the
#' modelling scale (log if the variable's dictionary `scale` is "log").
.ce_exposure_bounds <- function(data, exposure, lo_col, hi_col, log_scale) {
  lo <- data[[lo_col]]; hi <- data[[hi_col]]
  if (is.null(lo) || is.null(hi)) {
    stop("Censored exposure '", exposure, "' needs columns '", lo_col,
         "' and '", hi_col, "' in the data.", call. = FALSE)
  }
  if (log_scale) {                       # impute on the log scale (where shash lives)
    lo <- ifelse(lo <= 0, -Inf, log(lo)) # 0/negative lower bound -> left-censored
    hi <- ifelse(is.finite(hi), log(hi), Inf)
  }
  exact <- is.finite(lo) & is.finite(hi) & abs(hi - lo) < 1e-9
  list(
    y     = ifelse(exact, lo, NA_real_),               # observed value where exact
    lower = ifelse(exact, NA_real_, lo),
    upper = ifelse(exact, NA_real_, hi),
    log_scale = log_scale
  )
}

#' Initialise a censored exposure column: observed value where exact, else a
#' point inside the interval (bounded away from an infinite end).
.ce_init_exposure <- function(b) {
  y <- b$y
  cens <- is.na(y)
  lo <- b$lower[cens]; hi <- b$upper[cens]
  fill <- ifelse(is.finite(lo) & is.finite(hi), (lo + hi) / 2,
          ifelse(is.finite(hi), hi - 0.5,                 # left-censored: just below LOD
          ifelse(is.finite(lo), lo + 0.5, 0)))            # right-censored
  y[cens] <- fill
  y
}

#' Default number of inner FCS iterations for the Z block.
#'
#' The OUTER block-FCS loop already alternates the Z and X blocks, so it was
#' assumed the Z-block imputer needed no inner iterations of its own: 3 outer x
#' 3 inner ran 9 alternations where 3 were designed, and tripled the BART fits.
#' V8 tested that assumption at 1000 reps per cell and it is only half true --
#' the outer loop alternates BETWEEN blocks, not among the Z block's own
#' targets, so with several targets one inner pass never lets them condition on
#' each other:
#'
#'   two targets, neither the outcome  paired bias diff  +0.04 pp (CI +/-0.13) -- free
#'   three targets incl. imputed Y     paired bias diff  -0.72 pp             -- breaches
#'
#' So inner = 1 is used only in the configuration measured to be free, and the
#' pre-fix inner = 3 everywhere else. Y is excluded from the cheap path
#' deliberately: the X block conditions on Y, so an under-conditioned Y feeds
#' straight into the exposure draw. V8's cells confound "3 targets" with
#' "Y is a target" (see FINDINGS_v8.md, Caveats), and this predicate is the
#' reading that is safe under either explanation.
#'
#' An explicit `analysis_spec$imputation$bart_inner_iter` always wins.
#'
#' @param targets Character vector of Z-block imputation targets this sweep.
#' @param y_var The outcome variable name.
#' @return 1L on the validated cheap path, else 3L.
#' @seealso validation/phase1/FINDINGS_v8.md
.ce_default_inner_iter <- function(targets, y_var) {
  if (length(targets) <= 2L && !(y_var %in% targets)) 1L else 3L
}

#' The automatic X-block predictor set.
#'
#' Outcome + in-model covariates + AUXILIARY covariates + the exposures. The
#' auxiliaries are the 2026-09-09 change: see the note at the call site and
#' validation/phase1/FINDINGS_v24.md. Both dictionary columns are optional,
#' hence the defensive accessors.
#'
#' Its own function because no bundled example reaches it -- the censored-exposure
#' example passes `predictors` explicitly, which overrides this entirely -- so a
#' regression here would otherwise be invisible to the test suite.
#' See test/test_censored_exposure.R.
ce_auto_predictors <- function(var_dict, y_var, exposures) {
  flag <- function(col) {
    v <- if (col %in% names(var_dict)) var_dict[[col]] else rep(FALSE, nrow(var_dict))
    v <- as.logical(v); v[is.na(v)] <- FALSE; v
  }
  role <- if ("role" %in% names(var_dict)) var_dict$role else rep("", nrow(var_dict))
  unique(c(y_var,
           var_dict$var[flag("use_in_model")],
           var_dict$var[flag("use_as_auxiliary") | role %in% "auxiliary"],
           exposures))
}

# =============================================================================
# WHY THERE IS NO DAG CHILD FACTOR HERE  --  removed 2026-09-09 after measuring
# -----------------------------------------------------------------------------
# The correct distribution for a below-LOD exposure factorises along the DAG:
#
#   p(x | rest, x<=L) ~ p(x | Pa(X)) . p(Y | x, Pa(Y)) . PROD p(C | x, .) . 1{x<=L}
#
# leftcens covers the first two as a SINGLE regression of x on its predictor
# set, which is exact wherever the joint is Gaussian -- confounder +0.20%,
# linear mediator -0.02% (direct) / +0.19% (total), collider +0.16%. It is
# badly biased in one case: a mediator whose relationship with the exposure is
# CURVED, at +22.0% (direct) and +17.9% (total), and that does not shrink with n.
#
# The third factor fixes it: a grid sampler that builds the product explicitly
# lands at -0.09% to -0.63% across every case, at three sample sizes
# (validation/phase1/FINDINGS_v24.md). So the ALGORITHM is validated.
#
# WHAT WAS TRIED HERE AND DOES NOT WORK. An implementation that asks leftcens
# for K candidate draws and reweights them by the child factor -- attractive
# because it leaves leftcens's shash margin untouched, and reimplementing that
# margin caused four separate defects (FINDINGS_v13.md). Measured against the
# grid sampler on identical data:
#
#                                   linear pipe    pipe_nl (the case that needs it)
#   weight p(C | x)                    + 7.61%        +27.37%
#   weight p(C | x, Y)   <- correct    + 0.46%        +32.53%
#   grid sampler                       - 0.60%        - 0.41%
#
# Two findings, in order. FIRST, the weight is p(C | x, Y), not p(C | x): the
# proposal already conditions on Y, and since
# p(Y|x,C).p(C|x) = p(Y|x).p(C|Y,x), the residual factor carries Y. That
# correction is real -- it fixes the linear pipe outright.
#
# SECOND, and fatally, it does not help where the feature is needed. With the
# child removed from the predictor set (required, or its information is counted
# twice), leftcens regresses x on Y LINEARLY -- but marginalising the child out
# makes Y = b1*x + gamma1*g(x) + ..., genuinely non-linear in x. So the PROPOSAL
# is misspecified, and importance reweighting can only ADD a missing factor, not
# repair a wrong proposal; that would need the full target/proposal ratio, and
# leftcens does not expose its density. The clean contrast is the table above:
# same code, same weight, linear proposal works and non-linear proposal does not.
#
# WHAT A CORRECT IMPLEMENTATION NEEDS. The grid sampler, whose outcome factor
# conditions on the child (where Y IS linear in x) -- which means evaluating the
# ANALYSIS model's linear predictor as a function of a candidate exposure value,
# per row. Feasible, but not a small addition: the analysis model here is an
# arbitrary brms formula (splines, interactions, mo() terms), and the margin
# would have to be preserved by working on the shash-z scale via leftcens's
# exported fit_shash_margin() / x_to_z() / z_to_x(). Tracked in ROADMAP.md.
#
# UNTIL THEN, and this is what ce_preflight_curvature() below tells the user:
# the curved-mediator case has a zero-assumption mitigation. Removing such a
# covariate from censored_exposure$predictors roughly HALVES the bias
# (+22.03% -> +9.07%), because forcing a linear conditional onto a curved
# relationship does more damage than discarding the information. Not a fix, but
# free and strictly better.
#
# Scripts: validation/phase1/derive_sir_weight.R, derive_sir_proposal.R,
#          derive_shipped_form.R, derive_dag_basis.R
# =============================================================================

#' Preflight: is this analysis in the cell where the exposure draw is unsafe?
#'
#' WHY THIS IS WORTH SHIPPING SEPARATELY FROM THE FIX. The exposure draw is
#' exact for a confounder, and exact for a mediator whose relationship with the
#' exposure is LINEAR (measured -0.02%). It is badly biased only when that
#' relationship is CURVED (+22.0% direct, +17.9% total, at 40% non-detects).
#' Whether a relationship is curved is visible ABOVE the LOD, so unlike the
#' child factor itself this check needs no assumption -- it reports what the
#' observed data shows, plus the non-detect rate that decides how far the draw
#' must extrapolate it.
#'
#' WHAT IT DOES NOT DO. It cannot measure the curvature BELOW the LOD; validation
#' measured that as not identifiable (THEORY.md 0b -- a natural spline returns
#' exactly zero there by construction). So a clean report is not a guarantee, it
#' is the absence of a visible warning sign.
#'
#' @return A data.frame, one row per (exposure, covariate), invisibly.
ce_preflight_curvature <- function(data, analysis_spec, var_dict,
                                   covariates = NULL, verbose = TRUE) {
  ce <- analysis_spec$imputation$censored_exposure
  if (is.null(ce) || length(ce$exposure_vars) == 0L) return(invisible(NULL))
  lo_suffix <- ce$lo_suffix %||% "_lo"; hi_suffix <- ce$hi_suffix %||% "_hi"
  y_var <- analysis_spec$outcome$y_var

  if (is.null(covariates)) {
    # Restrict to variables the DICTIONARY declares, not every numeric column --
    # otherwise row ids, weights and bookkeeping columns get screened and the
    # report fills with noise the user has to learn to ignore.
    dict_ok <- if ("role" %in% names(var_dict))
      var_dict$var[var_dict$role %in% c("covariate", "auxiliary", "exposure")]
      else var_dict$var
    num <- names(data)[vapply(data, is.numeric, logical(1))]
    covariates <- setdiff(intersect(num, dict_ok),
                          c(y_var, ce$exposure_vars,
                            paste0(ce$exposure_vars, lo_suffix),
                            paste0(ce$exposure_vars, hi_suffix)))
  }
  out <- list()
  for (x in ce$exposure_vars) {
    lo <- data[[paste0(x, lo_suffix)]]; hi <- data[[paste0(x, hi_suffix)]]
    if (is.null(lo) || is.null(hi)) next
    obs <- is.finite(lo) & is.finite(hi) & abs(hi - lo) < 1e-9
    ndr <- mean(!obs, na.rm = TRUE)
    xv <- ifelse(obs, lo, NA_real_)
    lg <- if (!is.null(ce$log_scale)) isTRUE(ce$log_scale) else
      identical(var_dict$scale[match(x, var_dict$var)], "log")
    if (lg) xv <- ifelse(!is.na(xv) & xv > 0, log(xv), NA_real_)
    for (z in covariates) {
      zv <- data[[z]]
      ok <- !is.na(xv) & !is.na(zv)
      # Need enough observed exposure values, and a covariate that actually
      # varies, before a curvature comparison means anything.
      if (sum(ok) < 50L || stats::sd(zv[ok]) == 0) next
      f_lin <- stats::lm(zv[ok] ~ xv[ok])
      f_crv <- tryCatch(stats::lm(zv[ok] ~ splines::ns(xv[ok], df = 4)),
                        error = function(e) NULL)
      if (is.null(f_crv)) next
      s_lin <- stats::sigma(f_lin); s_crv <- stats::sigma(f_crv)
      if (!is.finite(s_crv) || s_lin <= 0 || s_crv <= 0) next
      # DEGREES-OF-FREEDOM CORRECTION. A 4-df spline fits noise slightly better
      # than a straight line even when the truth IS a straight line, so the raw
      # excess variance is positive in expectation under the null and the bias
      # grows as n shrinks. Under the null E[RSS_lin - RSS_crv] = dfe * sigma^2,
      # so that much is subtracted before anything is reported. Without this the
      # screen false-positived at 0.084 on a genuinely linear covariate at
      # n = 800 (against 0.0024 for the same DGP at n = 20,000) -- the same
      # flexible-basis false positive THEORY.md 0b measured for cubic bases.
      rss_l <- sum(stats::residuals(f_lin)^2)
      rss_c <- sum(stats::residuals(f_crv)^2)
      dfe   <- f_lin$df.residual - f_crv$df.residual
      excess <- max(0, (rss_l - rss_c) - dfe * s_crv^2)
      # The correction removes the EXPECTED excess but not its sampling noise,
      # so at small n a linear truth still reads a non-zero effect size (0.068
      # at n = 800 against 0.012 at n = 20,000 on the same linear DGP). The F
      # test is what controls that: under a linear truth its p-value is uniform,
      # whatever n. Both are reported, and the screen requires BOTH -- the p
      # value stops small-n noise, the effect size stops a large-n study
      # flagging a curvature too small to matter.
      p_curv <- tryCatch(
        stats::anova(f_lin, f_crv)[["Pr(>F)"]][2L], error = function(e) NA_real_)
      out[[length(out) + 1L]] <- data.frame(
        exposure = x, covariate = z, nd_rate = ndr, n_obs = sum(ok),
        # `abs_cor` is how much the covariate knows about the exposure at all.
        abs_cor = abs(stats::cor(xv[ok], zv[ok])),
        # `curvature` is the NON-LINEAR part of the relationship in units of the
        # covariate's own residual noise, corrected for the spline's extra
        # degrees of freedom (see above).
        #
        # WHY NOT THE RELATIVE SD REDUCTION (1 - s_crv/s_lin), which is the
        # obvious choice and was the first version. It is diluted by the
        # covariate's own noise, so it reads NEAR ZERO on a genuinely dangerous
        # relationship. Calibrated against the two DGPs whose bias is known:
        #
        #   DGP                          measured bias   1 - s_crv/s_lin   this
        #   --------------------------   -------------   ---------------   ------
        #   validation linear pipe             -0.02%            0.0000    ~0.00
        #   validation pipe_nl                +22.03%            0.0165     0.18
        #
        # The relative form puts the KNOWN-BAD case at 0.0165 -- under any
        # sensible threshold, so the screen would have missed the one cell it
        # exists for. Caught by test/test_censored_exposure.R, not by
        # inspection.
        curvature = sqrt(excess / sum(ok)) / s_crv,
        p_curv = p_curv,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(invisible(NULL))
  res <- do.call(rbind, out)
  res <- res[order(-res$curvature), , drop = FALSE]
  rownames(res) <- NULL

  # SCREENING values, not gates, and three of them because no one of them is
  # enough: |cor| > 0.10 says the covariate knows something about the exposure
  # at all; p < 0.01 says the curvature is not sampling noise (this is what
  # holds the false-positive rate down at small n); curvature > 0.05 says it is
  # large enough to matter (this is what stops a huge study flagging a
  # negligible one). The non-detect rate is the multiplier -- validation
  # measured the effect as negligible at 5% non-detects and dominant at 40%+.
  flag <- res$abs_cor > 0.10 & res$curvature > 0.05 & res$nd_rate > 0.15 &
    !is.na(res$p_curv) & res$p_curv < 0.01
  if (verbose && any(flag)) {
    message("\n=== Censored-exposure preflight: covariate-exposure curvature ===")
    for (i in which(flag)) with(res[i, ], message(sprintf(
      "  %-18s ~ %-18s  |cor| %.2f  curvature %.2f (p %.1e)  non-detects %.0f%%",
      exposure, covariate, abs_cor, curvature, p_curv, 100 * nd_rate)))
    message(
      "  These covariates have a CURVED relationship with a censored exposure.\n",
      "  If the exposure CAUSES any of them (a mediator), the exposure draw is\n",
      "  misspecified below the LOD: measured +22% bias (direct effect) and +18%\n",
      "  (total effect) at 40% non-detects, and it does NOT shrink with n.\n",
      "  If the covariate is a CAUSE of the exposure (a confounder), the draw is\n",
      "  fine -- only the mediator direction is affected, and only you can say\n",
      "  which it is.\n",
      "  Options, in order of preference:\n",
      "    1. Transform the covariate so the relationship is straight, and pass\n",
      "       the transformed version via censored_exposure$predictors. This is\n",
      "       a real fix and costs nothing.\n",
      "    2. Remove the mediator from censored_exposure$predictors -- measured\n",
      "       to roughly HALVE the bias (+22% -> +9%) at no assumption cost,\n",
      "       because forcing a linear conditional onto a curved relationship\n",
      "       does more damage than discarding the information. Not a fix.\n",
      "    3. Report the estimate as biased away from the null. There is no\n",
      "       option to switch on: see the header of this file for the method\n",
      "       that works and why it is not implemented yet.\n",
      "  See docs/covariate-roles.md and validation/phase1/FINDINGS_v24.md.")
  } else if (verbose) {
    message("Censored-exposure preflight: no covariate shows both a curved ",
            "relationship with a censored exposure and a non-detect rate above ",
            "15%. (This is the absence of a warning sign, not a guarantee -- ",
            "curvature below the LOD is not identifiable.)")
  }
  invisible(res)
}

#' One block-FCS completed dataset (a single imputation).
.ce_one_imputation <- function(data, analysis_spec, var_dict, ce, seed) {
  set.seed(seed)
  y_var <- analysis_spec$outcome$y_var
  exposures <- ce$exposure_vars
  lo_suffix <- ce$lo_suffix %||% "_lo"
  hi_suffix <- ce$hi_suffix %||% "_hi"
  sweeps <- as.integer(ce$outer_sweeps %||% 5L)
  margin <- ce$margin %||% "shash"
  zb_scale <- match.arg(as.character(ce$zblock_scale %||% "model"),
                        c("model", "data"))

  # Whether to impute on the log scale: an explicit config flag wins; otherwise
  # fall back to the dictionary `scale == "log"`. Keeping this separate from the
  # dictionary scale avoids double-transforming when the model formula already
  # applies log() to the exposure.
  log_for <- function(x) {
    if (!is.null(ce$log_scale)) return(isTRUE(ce$log_scale))
    identical(var_dict$scale[match(x, var_dict$var)], "log")
  }

  # Per-exposure bounds on the modelling scale, and the working column filled in.
  bounds <- list(); work <- data; drop_cols <- character(0)
  for (x in exposures) {
    lo_col <- paste0(x, lo_suffix); hi_col <- paste0(x, hi_suffix)
    b <- .ce_exposure_bounds(data, x, lo_col, hi_col, log_scale = log_for(x))
    bounds[[x]] <- b
    # On the MODELLING scale (log where configured), which is the scale the Z
    # block must see -- see the note at the end of the sweep loop.
    work[[x]] <- .ce_init_exposure(b)                 # complete predictor for the Z block
    drop_cols <- c(drop_cols, lo_col, hi_col)
  }
  # The lo/hi columns have served their purpose; remove them so they never enter
  # the Z-block predictors or the completed dataset (downstream uses the exposure).
  work <- work[, setdiff(names(work), drop_cols), drop = FALSE]

  # Predictor set for the X block: reduced set if given, else auto (outcome +
  # in-model covariates + AUXILIARY covariates + the OTHER exposures), matching
  # the plan's guidance to avoid conditioning on all ~300 Z.
  #
  # WHY AUXILIARIES ARE IN THIS SET (V24, 2026-09-09). They were not, and the
  # omission was documented as a sub-1-percentage-point curiosity. It is not:
  # the variable most likely to be marked "impute only, not in the model" is a
  # MEDIATOR in a total-effect analysis, and a mediator is exactly the covariate
  # that carries information about a censored exposure (it is an effect of it).
  # Dropping such a covariate from this set was measured at 6-7% bias, FLAT in n
  # -- so interval coverage falls as the study grows: 0.945 -> 0.830 -> 0.535 for
  # a direct effect and 0.830 -> 0.480 -> 0.040 for a total effect over
  # n = 800 / 3200 / 12800. Keeping it costs nothing and needs no assumption:
  # conditioning on it LINEARLY is enough whenever its relationship with the
  # exposure is linear. See validation/phase1/FINDINGS_v24.md.
  #
  # `role == "auxiliary"` and `use_as_auxiliary == TRUE` both count, matching
  # `auxiliary_vars` in 00_common_functions.R rather than inventing a second
  # rule. Both columns are optional, hence the defensive accessor.
  x_predictors <- ce$predictors %||% ce_auto_predictors(var_dict, y_var, exposures)

  # Impute Y in the Z block too (MID): the X-block conditions on Y, so Y must be
  # complete for that step; the imputed-Y rows are deleted afterwards. Track the
  # ORIGINAL missing cells so every outer sweep re-imputes them given the current
  # X (miceRanger only targets currently-NA cells — without this the block would
  # stop alternating after the first sweep).
  z_spec_as <- analysis_spec
  z_spec_as$imputation$impute_y <- TRUE
  # Force the Z-block (miceRanger) to run SERIALLY: this function already
  # parallelises over the m completed datasets (mclapply), and a PSOCK cluster
  # started inside a fork fails ("cannot open server socket"). Cross-dataset
  # parallelism (n_cores) is the right level; the per-dataset Z-block stays serial.
  z_spec_as$parallel$impute_workers <- 1L
  # Z-block inner FCS count. This cannot be resolved here: it depends on which
  # targets the Z block actually has, which is only known once
  # make_row_level_imputation_spec() has run inside the sweep loop below. Capture
  # any explicit user/validation-arm setting now (it always wins) and apply
  # .ce_default_inner_iter() per sweep. Read from `analysis_spec`, not the
  # `z_spec_as` copy, so the per-sweep assignment below cannot be mistaken for a
  # user setting on the next sweep.
  user_inner_iter <- analysis_spec$imputation$bart_inner_iter
  # Reset only the originally-missing Z/Y cells each sweep — NOT the exposures,
  # whose missing (censored) cells are handled by the X block via the fixed bounds.
  orig_na <- lapply(data[intersect(names(data), names(work))], is.na)
  reset_cols <- setdiff(names(orig_na)[vapply(orig_na, any, logical(1))], exposures)

  for (t in seq_len(sweeps)) {
    # Re-open the originally-missing Z/Y cells so the block re-imputes them.
    for (cc in reset_cols) work[[cc]][orig_na[[cc]]] <- NA
    # --- Z block: miceRanger (reused helper), m = 1, conditioning on current X, Y.
    z_spec <- make_row_level_imputation_spec(work, z_spec_as, var_dict)
    if (length(z_spec$vars) > 0) {
      z_spec_as$imputation$bart_inner_iter <-
        user_inner_iter %||% .ce_default_inner_iter(names(z_spec$vars), y_var)
      z_spec$m <- 1L; z_spec$seed <- seed + 1000L + t
      work <- run_row_level_imputation(work, z_spec, z_spec_as)[[1]]
    }
    # --- X block: one congenial censored draw per exposure, given current Z, Y.
    for (x in exposures) {
      b <- bounds[[x]]
      preds <- setdiff(intersect(x_predictors, names(work)), x)
      xmat <- work[, preds, drop = FALSE]
      # A predictor that still carries NA at this point would reach leftcens as a
      # silent hole rather than an error. It can only happen if the variable is
      # not a Z-block target -- e.g. an auxiliary marked impute_target = FALSE
      # that has missing values of its own -- so name the cause rather than the
      # symptom. Checked here and not once up front because the Z block fills
      # cells in as the sweeps proceed.
      na_preds <- names(xmat)[vapply(xmat, anyNA, logical(1))]
      if (length(na_preds)) {
        stop("X-block predictor(s) still missing after the Z block: ",
             paste(na_preds, collapse = ", "),
             ". Mark them impute_target = TRUE in 00_variable_dictionary.csv so ",
             "the Z block completes them, or drop them from ",
             "censored_exposure$predictors.", call. = FALSE)
      }
      imp <- leftcens::impute_censored_conditional(
        y = b$y, x = xmat, lower = b$lower, upper = b$upper,
        m = 1L, margin = margin)[, 1]
      # STAY ON THE MODELLING SCALE UNTIL THE SWEEPS ARE DONE. This used to
      # exponentiate here, inside the loop, which meant every LATER sweep's Z
      # block regressed the covariates on `exp(log X)` while the analysis model
      # -- and the truth -- are linear in `log X`. The Z-block imputation model
      # was therefore misspecified in SCALE, and the whole block-FCS was
      # incongenial with the analysis for exactly the documented configuration
      # (`scale = "log"`).
      #
      # Measured 2026-09-10: both pipeline Z-block arms sat 3-5 pp above their
      # harness equivalents with two entirely different Z engines, and changing
      # only this scale in the harness reproduced it (+0.38% -> +5.38% under a
      # confounder, -0.45% -> +5.44% under a mediator). See
      # validation/phase1/derive_zblock_scale.R and ROADMAP.md.
      # `zblock_scale` exists so the pre-v1.6.0 defect stays REPRODUCIBLE as a
      # configuration rather than as a patch: "model" (default) is correct,
      # "data" restores the old in-loop conversion. Two reasons to keep it:
      # measuring the defect needs both versions on identical data in ONE run,
      # and a permanent switch is a regression guard -- a future refactor that
      # reintroduces the conversion is caught by test/test_censored_exposure.R
      # rather than by someone re-deriving it. NOT a user-facing knob; it is
      # undocumented in 00_config.R and has no reason to be set.
      work[[x]] <- if (identical(zb_scale, "data") && isTRUE(b$log_scale))
        exp(imp) else imp
    }
  }
  # Back to the data scale ONCE, after the last sweep, so the returned dataset
  # keeps the user-facing convention (a data-scale exposure column, which the
  # analysis formula logs itself).
  if (!identical(zb_scale, "data")) {
    for (x in exposures) {
      if (isTRUE(bounds[[x]]$log_scale)) work[[x]] <- exp(work[[x]])
    }
  }
  work
}

#' Run the censored-exposure block-FCS and return an `imputed_list` (m completed
#' datasets), the same shape 03_impute.R expects from the other strategies.
#'
#' @param data,analysis_spec,var_dict As in the other Step-3 strategies.
#' @param m Number of completed datasets to produce (the new batch size).
#' @param seed Base seed for this batch.
run_censored_exposure_block_fcs <- function(data, analysis_spec, var_dict, m, seed) {
  if (!requireNamespace("leftcens", quietly = TRUE) ||
      !"impute_censored_conditional" %in% getNamespaceExports("leftcens")) {
    stop("strategy 'censored_exposure_block_fcs' requires leftcens >= 0.9.0 ",
         "(exporting impute_censored_conditional()).", call. = FALSE)
  }
  ce <- analysis_spec$imputation$censored_exposure
  if (is.null(ce) || length(ce$exposure_vars) == 0) {
    stop("strategy 'censored_exposure_block_fcs' needs ",
         "analysis_spec$imputation$censored_exposure$exposure_vars.", call. = FALSE)
  }
  y_var <- analysis_spec$outcome$y_var
  mid <- isTRUE(ce$mid_delete_imputed_y %||% TRUE)
  y_missing <- if (!is.null(data[[y_var]])) is.na(data[[y_var]]) else logical(nrow(data))

  # The m completed datasets are independent -> parallelise over them (fork).
  # n_cores from censored_exposure$n_cores, else parallel$impute_workers, else 1.
  n_cores <- as.integer(ce$n_cores %||% analysis_spec$parallel$impute_workers %||% 1L)
  if (.Platform$OS.type == "windows") n_cores <- 1L        # mclapply forks: unix only
  n_cores <- max(1L, min(n_cores, m))

  log_msg("Censored-exposure block-FCS | exposures:",
          paste(ce$exposure_vars, collapse = ", "),
          "| outer_sweeps:", ce$outer_sweeps %||% 5L,
          "| margin:", ce$margin %||% "shash",
          "| m:", m, "| n_cores:", n_cores)

  build_one <- function(d) {
    work <- .ce_one_imputation(data, analysis_spec, var_dict, ce, seed = seed + d)
    # MID: the imputed-Y rows carry no information for beta; drop them before the
    # brms fit so pooling is over information-bearing rows only (von Hippel 2007).
    if (mid && any(y_missing)) work <- work[!y_missing, , drop = FALSE]
    # Invariants (asserted, not just commented): every censored exposure is fully
    # imputed, and MID has removed all imputed-Y rows (no Y missing remains).
    for (x in ce$exposure_vars) {
      if (anyNA(work[[x]]))
        stop("Censored exposure '", x, "' still has NA after imputation (dataset ", d, ").", call. = FALSE)
    }
    if (mid && any(y_missing) && !is.null(work[[y_var]]) && anyNA(work[[y_var]]))
      stop("MID invariant violated: imputed-Y rows were not removed (dataset ", d, ").", call. = FALSE)
    log_msg("  completed dataset", d, "of", m,
            if (mid && any(y_missing)) paste0("(MID dropped ", sum(y_missing), " imputed-Y rows)") else "")
    tibble::as_tibble(work)
  }

  if (n_cores > 1L) {
    imputed_list <- parallel::mclapply(seq_len(m), build_one, mc.cores = n_cores, mc.preschedule = FALSE)
    errs <- vapply(imputed_list, inherits, logical(1), what = "try-error")
    if (any(errs)) stop("Censored-exposure imputation failed on dataset(s) ",
                        paste(which(errs), collapse = ", "), ": ",
                        as.character(imputed_list[[which(errs)[1]]]), call. = FALSE)
  } else {
    imputed_list <- lapply(seq_len(m), build_one)
  }

  imputed_list
}
