# =============================================================================
# Phase 1 harness -- procedure 5: the SHIPPED pipeline engine
# -----------------------------------------------------------------------------
# Track V1 of ../PLAN_pipeline_validation.md; generalised for Track V2.
#
# Every other congenial procedure in this harness is a *prototype*
# (`cens_mi_y_shash`, ~40 lines written for the study). This one calls the code
# users actually run:
#
#     run_censored_exposure_block_fcs()   in <root>/00_censored_exposure.R
#
# It is imported by sourcing that file -- never copy-pasted -- so a change to the
# pipeline is picked up here automatically, and a divergence between this row and
# the `cens_mi_y_shash` row is evidence about the shipped code, not about a
# re-implementation of it.
#
# Everything downstream of imputation (the analysis model via fit_lm_estimand(),
# Rubin pooling via rubin_pool(), the returned one_row() shape) is the harness's
# existing machinery, so the comparison against oracle / complete_case /
# leftcens_prestep / cens_mi_y_shash is like-for-like: the ONLY thing that differs
# across procedures is how missing data were handled.
#
# V2 generalisation: any subset of exposures censored, MCAR covariates (which
# makes the miceRanger Z block real rather than a no-op), and a missing outcome
# (which fires MID). V1 handled only "X1 censored, everything else complete".
# =============================================================================

# Defined by the runner; guarded so this file is usable standalone too.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---- locating and loading the pipeline --------------------------------------

#' Find the pipeline project root (the directory holding 00_censored_exposure.R).
#'
#' Order: explicit argument, then $PIPELINE_ROOT, then walk up from this file.
.v1_project_root <- function(project_root = NULL) {
  if (!is.null(project_root)) return(normalizePath(project_root, mustWork = TRUE))

  env_root <- Sys.getenv("PIPELINE_ROOT", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = TRUE))

  # This file lives at <root>/validation/phase1/R/. `.here` is set by the runner
  # to <root>/validation/phase1 when it sources the component files.
  start <- if (exists(".here", inherits = TRUE)) {
    get(".here", inherits = TRUE)
  } else {
    getwd()
  }

  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "00_censored_exposure.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate the pipeline project root (no 00_censored_exposure.R ",
       "found walking up from '", start, "'). Set PIPELINE_ROOT.", call. = FALSE)
}

# Cache: sourcing the pipeline is done once per session, not once per replication.
.v1_cache <- new.env(parent = emptyenv())

#' Source the pipeline into a private environment and return it.
#'
#' Sourced into a dedicated env (not globalenv) so pipeline definitions cannot
#' shadow harness functions. `run_censored_exposure_block_fcs()` resolves its
#' helpers (`make_row_level_imputation_spec`, `run_row_level_imputation`,
#' `log_msg`, `%||%`) inside this same env, so calling it from here works.
#'
#' GOTCHA: 00_common_functions.R is NOT standalone-sourceable -- it reads
#' `paths$objects` at top level (the mi_runtime_override mechanism) and mutates
#' `analysis_spec`. Both are stubbed before sourcing.
#'
#' PERFORMANCE: call this ONCE in the parent process before forking. The forked
#' workers then inherit the loaded environment copy-on-write instead of each
#' re-sourcing ~57 KB of R. `run_v2_robustness.R` does exactly that.
#'
#' @param quiet Replace `log_msg()` with a no-op. The module logs several lines
#'   per completed dataset; across thousands of tasks that is pure noise.
#'   Errors still propagate normally.
v1_load_pipeline <- function(project_root = NULL, quiet = TRUE,
                             proper_z = FALSE, mice_z = FALSE, bart_z = FALSE) {
  key <- paste0("env_", if (isTRUE(quiet)) "quiet" else "loud",
                if (isTRUE(proper_z)) "_properz" else "",
                if (isTRUE(mice_z)) "_micez" else "",
                if (isTRUE(bart_z)) "_bartz" else "")
  if (!is.null(.v1_cache[[key]])) return(.v1_cache[[key]])

  root <- .v1_project_root(project_root)
  e <- new.env(parent = globalenv())

  # Stubs required before sourcing (see gotcha above).
  assign("paths", list(objects = tempdir(), cache = tempdir()), envir = e)
  assign("analysis_spec", list(), envir = e)

  sys.source(file.path(root, "00_common_functions.R"), envir = e)
  sys.source(file.path(root, "00_censored_exposure.R"), envir = e)

  if (!exists("run_censored_exposure_block_fcs", envir = e, inherits = FALSE)) {
    stop("run_censored_exposure_block_fcs() not found after sourcing the ",
         "pipeline from '", root, "'.", call. = FALSE)
  }

  if (isTRUE(quiet)) assign("log_msg", function(...) invisible(NULL), envir = e)

  # Track V4: swap ONLY the Z-block imputer for a proper Bayesian draw, leaving
  # the X block, the outer alternation, MID and the pooling untouched. If this
  # restores interval calibration, miceRanger's improper draw is the cause.
  if (isTRUE(proper_z)) {
    if (!exists(".v4_proper_row_imputation")) {
      stop("proper_z = TRUE needs R/proper_impute.R sourced.", call. = FALSE)
    }
    assign("run_row_level_imputation", .v4_proper_row_imputation, envir = e)
  }

  # Track 05: swap the Z block for mice's parametric proper draws. Replacing the
  # whole function also removes the pipeline's internal proper_draw dispatch, so
  # this arm is unambiguously "mice and nothing else".
  if (isTRUE(mice_z)) {
    if (!exists(".v5_mice_row_imputation")) {
      stop("mice_z = TRUE needs R/mice_impute.R sourced.", call. = FALSE)
    }
    assign("run_row_level_imputation", .v5_mice_row_imputation, envir = e)
  }

  # R8 / item 06: swap the Z block for BART -- flexible AND properly dispersed.
  if (isTRUE(bart_z)) {
    if (!exists(".v7_bart_row_imputation")) {
      stop("bart_z = TRUE needs R/bart_impute.R sourced.", call. = FALSE)
    }
    assign("run_row_level_imputation", .v7_bart_row_imputation, envir = e)
  }

  assign("project_root", root, envir = e)
  .v1_cache[[key]] <- e
  e
}

# ---- adapting the harness bundle to the pipeline's input convention ---------

#' Build the pipeline's expected inputs from a harness bundle.
#'
#' SCALE GOTCHA. The harness works in log space: `logXj` holds the true value and
#' is NA on censored rows, with `logXj_lod` carrying the log LOD. The pipeline's
#' user-facing convention is a **data-scale** exposure column plus `X_lo`/`X_hi`
#' bounds, with `censored_exposure$log_scale = TRUE` telling the X block to impute
#' on the log scale. We deliberately use that convention rather than feeding log
#' values directly, because exercising `.ce_exposure_bounds()`'s log handling is
#' part of validating the shipped path.
#'
#'   row        Xj            Xj_lo            Xj_hi
#'   observed   exp(logXj)    exp(logXj)       exp(logXj)
#'   censored   NA            0                exp(logXj_lod)
#'
#' `Xj_lo = 0` is the left-censored marker: `.ce_exposure_bounds()` maps a
#' non-positive lower bound to -Inf on the log scale.
#'
#' Censored exposures are discovered from the bundle's `censored_cols` attribute,
#' so this handles one, some, or all of them. Exposures that were NOT censored
#' stay as plain observed `logXj` predictor columns.
#'
#' Covariates carrying MCAR missingness are marked `impute_target = TRUE` so the
#' miceRanger Z block actually imputes them -- that is what turns the block-FCS
#' alternation on. A missing outcome is left as-is; the module sets `impute_y`
#' internally and MID deletes those rows before the analysis fit.
#'
#' @param bundle A harness bundle: list(complete, censored, truth).
#' @param outer_sweeps,margin,n_cores Pipeline knobs (see 00_config.R).
#' @return list(data, analysis_spec, var_dict, cens_x, expo_names, ...)
.v1_make_pipeline_inputs <- function(bundle, outer_sweeps = 3L,
                                     margin = "shash", n_cores = 1L,
                                     mid = TRUE, proper_draw = FALSE) {
  d <- bundle$censored
  truth <- bundle$truth
  p <- length(truth$b)
  q <- length(truth$gamma)

  all_x  <- paste0("logX", seq_len(p))
  cens_x <- attr(d, "censored_cols") %||% "logX1"   # e.g. c("logX1","logX2")
  cens_x <- intersect(cens_x, all_x)
  open_x <- setdiff(all_x, cens_x)                  # exposures never censored
  covars <- c("Z1", "Z2")[seq_len(q)]

  expo_names <- sub("^log", "", cens_x)             # logX1 -> X1

  out <- data.frame(Y = d$Y, row_id = seq_len(nrow(d)), check.names = FALSE)

  for (k in seq_along(cens_x)) {
    cc <- cens_x[k]
    nm <- expo_names[k]
    is_c <- d[[paste0(cc, "_cens")]] == "left"
    val  <- exp(d[[cc]])                            # NA where censored
    out[[nm]]                  <- val
    out[[paste0(nm, "_lo")]]   <- ifelse(is_c, 0, val)
    out[[paste0(nm, "_hi")]]   <- ifelse(is_c, exp(d[[paste0(cc, "_lod")]]), val)
  }
  for (cc in open_x) out[[cc]] <- d[[cc]]
  for (zz in covars) out[[zz]] <- d[[zz]]

  # Dictionary. Only `var`, `impute_target`, `use_in_model`, `use_as_auxiliary`,
  # `role`, `type`, `scale` are read (via get_col defaults). Covariates with
  # missingness become imputation targets so the Z block engages.
  z_has_na <- vapply(covars, function(z) anyNA(out[[z]]), logical(1))
  n_expo   <- length(expo_names) + length(open_x)

  var_dict <- data.frame(
    var              = c("Y", expo_names, open_x, covars),
    role             = c("outcome", rep("exposure", n_expo),
                         rep("covariate", length(covars))),
    type             = c("continuous", rep("continuous", n_expo),
                         ifelse(covars == "Z2", "binary", "continuous")),
    timing           = "single",
    scale            = "none",
    impute_target    = c(FALSE, rep(FALSE, n_expo), unname(z_has_na)),
    use_in_model     = TRUE,
    use_as_auxiliary = FALSE,
    stringsAsFactors = FALSE
  )

  # Predictor set for the X block: outcome + the other exposures + covariates.
  # In the V1 configuration this reproduces `cens_mi_y_shash`'s set exactly, so
  # any difference between those two rows is the engine, not the conditioning
  # set. It stays linear on the mixture surface -- that is the §7.7 gap, and
  # deliberately not corrected for here.
  ce <- list(
    exposure_vars        = expo_names,
    predictors           = c("Y", expo_names, open_x, covars),
    outer_sweeps         = as.integer(outer_sweeps),
    margin               = margin,
    log_scale            = TRUE,
    lo_suffix            = "_lo",
    hi_suffix            = "_hi",
    mid_delete_imputed_y = isTRUE(mid),
    n_cores              = as.integer(n_cores)
  )

  analysis_spec <- list(
    data     = list(row_id_var = "row_id"),
    outcome  = list(y_var = "Y"),
    imputation = list(
      m = 1L, maxiter = 5L, verbose = FALSE,
      mean_match_k = NULL, seed = NULL,
      impute_y = FALSE,
      # Track V4 candidate fix: the pipeline's own opt-in proper-MI Z block.
      proper_draw = isTRUE(proper_draw),
      censored_exposure = ce
    ),
    parallel = list(impute_workers = as.integer(n_cores))
  )

  list(data = out, analysis_spec = analysis_spec, var_dict = var_dict,
       cens_x = cens_x, expo_names = expo_names, open_x = open_x,
       covars = covars, n_y_missing = sum(is.na(d$Y)))
}

# ---- procedure 5 -------------------------------------------------------------

#' Procedure 5: the shipped censored-exposure block-FCS.
#'
#' Calls `run_censored_exposure_block_fcs()` to produce `m` completed datasets,
#' fits the matched ERF model to each with the harness's `fit_lm_estimand()`, and
#' pools with `rubin_pool()`.
#'
#' PARALLELISM. The module forks over the `m` completed datasets via
#' `mclapply(mc.cores = ce$n_cores)`. `run_phase1()`'s replication loop is serial,
#' so passing cores through is safe there. `run_v2_robustness.R` instead forks at
#' the *replication* level and passes `n_cores = 1` here -- forking inside a fork
#' is what must be avoided, not forking as such.
#'
#' MID. `mid_delete_imputed_y` is TRUE. When the outcome has missing values the
#' module imputes Y (so the X block can condition on it) and then deletes those
#' rows before the fit. The expected row count is checked below, so MID is
#' *tested* rather than assumed -- in V1, where Y was complete, it never fired.
#'
#' @param bundle A harness bundle: list(complete, censored, truth).
#' @param m Number of completed datasets.
#' @param seed Base seed.
#' @param n_cores Forked workers for the `m` datasets (use 1 under rep-level forking).
#' @param outer_sweeps,margin Pipeline knobs.
#' @param project_root Pipeline root; defaults to $PIPELINE_ROOT or an upward walk.
#' @param quiet Silence the module's per-dataset logging.
proc_pipeline_block_fcs <- function(bundle, m = 20L, seed = NULL, n_cores = 1L,
                                    outer_sweeps = 3L, margin = "shash",
                                    project_root = NULL, quiet = TRUE,
                                    mid = TRUE, proper_z = FALSE,
                                    proper_draw = FALSE, mice_z = FALSE,
                                    bart_z = FALSE, label = NULL) {
  label <- label %||% "pipeline_block_fcs"

  if (!requireNamespace("leftcens", quietly = TRUE) ||
      !"impute_censored_conditional" %in% getNamespaceExports("leftcens")) {
    return(one_row(label, NA, NA, NA, NA, "needs leftcens >= 0.9.0"))
  }

  env <- tryCatch(v1_load_pipeline(project_root, quiet = quiet,
                                   proper_z = proper_z, mice_z = mice_z,
                                   bart_z = bart_z),
                  error = function(e) NULL)
  if (is.null(env)) {
    return(one_row(label, NA, NA, NA, NA, "could not source pipeline"))
  }

  if (!is.null(seed)) set.seed(seed)
  truth <- bundle$truth

  inp <- .v1_make_pipeline_inputs(bundle, outer_sweeps = outer_sweeps,
                                  margin = margin, n_cores = n_cores, mid = mid,
                                  proper_draw = proper_draw)

  imputed <- tryCatch(
    env$run_censored_exposure_block_fcs(
      data          = inp$data,
      analysis_spec = inp$analysis_spec,
      var_dict      = inp$var_dict,
      m             = as.integer(m),
      seed          = if (is.null(seed)) 1L else as.integer(seed)
    ),
    error = function(e) structure(list(), err = conditionMessage(e))
  )

  if (!length(imputed)) {
    msg <- attr(imputed, "err") %||% "block-FCS returned nothing"
    return(one_row(label, NA, NA, NA, NA, substr(paste("failed:", msg), 1, 120)))
  }

  # MID expectation: with MID on, exactly the imputed-Y rows are deleted and
  # nothing else; with MID off, no rows are deleted at all.
  n_expected <- if (isTRUE(mid)) nrow(inp$data) - inp$n_y_missing else nrow(inp$data)
  mid_ok <- TRUE

  ests <- vars <- rep(NA_real_, length(imputed))

  for (i in seq_along(imputed)) {
    di <- as.data.frame(imputed[[i]])

    if (nrow(di) != n_expected) mid_ok <- FALSE

    # Back to the harness's analysis scale: the module returns each censored
    # exposure on the DATA scale; the matched model is specified in logs.
    for (k in seq_along(inp$cens_x)) {
      di[[inp$cens_x[k]]] <- log(di[[inp$expo_names[k]]])
    }

    e <- fit_lm_estimand(di, truth)
    ests[i] <- e["est"]
    vars[i] <- e["se"]^2
  }

  pooled <- rubin_pool(ests, vars)

  note <- sprintf("m=%d; sweeps=%d; %s; cens=%d", m, outer_sweeps, margin,
                  length(inp$cens_x))
  if (isTRUE(proper_z))     note <- paste0(note, "; properZ")
  if (isTRUE(proper_draw))  note <- paste0(note, "; properBoot")
  if (isTRUE(mice_z))       note <- paste0(note, "; micePmm")
  if (isTRUE(bart_z))       note <- paste0(note, "; bartMI")
  if (!isTRUE(mid))     note <- paste0(note, "; MID off")
  if (isTRUE(mid) && inp$n_y_missing > 0)
    note <- paste0(note, sprintf("; MID dropped %d", inp$n_y_missing))
  if (!mid_ok)             note <- paste0(note, "; WARN MID row count wrong")
  if (truth$erf_form == "mixture") note <- paste0(note, "; linear-approx (§7.7)")
  n_bad <- sum(!is.finite(ests))
  if (n_bad > 0) note <- paste0(note, sprintf("; %d/%d fits failed", n_bad, length(ests)))

  one_row(label, pooled["est"], pooled["se"], pooled["ci_lo"], pooled["ci_hi"],
          substr(note, 1, 120),
          ubar = pooled["ubar"], b = pooled["b"], fmi = pooled["fmi"])
}

# ---- Track V4 isolating variants --------------------------------------------
# Each changes exactly ONE thing relative to `pipeline_block_fcs`, so a shift in
# width/SE attributes the variance shortfall to that component and nothing else.

#' Same engine, Z block swapped for a PROPER Bayesian draw.
#' Isolates: is miceRanger's improper (RF/PMM) draw understating between-variance?
proc_pipeline_properZ <- function(bundle, m = 20L, seed = NULL, n_cores = 1L,
                                  outer_sweeps = 3L, margin = "shash",
                                  project_root = NULL, quiet = TRUE) {
  proc_pipeline_block_fcs(bundle, m = m, seed = seed, n_cores = n_cores,
                          outer_sweeps = outer_sweeps, margin = margin,
                          project_root = project_root, quiet = quiet,
                          mid = TRUE, proper_z = TRUE,
                          label = "pipeline_properZ")
}

#' THE CANDIDATE FIX: the pipeline's own proper-MI Z block, switched on via
#' `analysis_spec$imputation$proper_draw`. Unlike `pipeline_properZ` -- a harness
#' instrument that swaps in a simple linear/logistic draw -- this arm exercises
#' real shipped pipeline code (`run_row_level_imputation_proper()`), which
#' bootstraps the training data per imputation and keeps the random forest.
#' This is the arm whose result decides whether the fix is adopted.
proc_pipeline_properBoot <- function(bundle, m = 20L, seed = NULL, n_cores = 1L,
                                     outer_sweeps = 3L, margin = "shash",
                                     project_root = NULL, quiet = TRUE) {
  proc_pipeline_block_fcs(bundle, m = m, seed = seed, n_cores = n_cores,
                          outer_sweeps = outer_sweeps, margin = margin,
                          project_root = project_root, quiet = quiet,
                          mid = TRUE, proper_z = FALSE, proper_draw = TRUE,
                          label = "pipeline_properBoot")
}

#' Track 05 arm: Z block replaced by `mice`'s parametric proper draws (pmm /
#' logreg). The counterpart to `properBoot` -- proper by construction but linear,
#' so it is the arm that a non-linear covariate DGP can finally penalise.
proc_pipeline_micePmm <- function(bundle, m = 20L, seed = NULL, n_cores = 1L,
                                  outer_sweeps = 3L, margin = "shash",
                                  project_root = NULL, quiet = TRUE) {
  proc_pipeline_block_fcs(bundle, m = m, seed = seed, n_cores = n_cores,
                          outer_sweeps = outer_sweeps, margin = margin,
                          project_root = project_root, quiet = quiet,
                          mid = TRUE, proper_z = FALSE, proper_draw = FALSE,
                          mice_z = TRUE, label = "pipeline_micePmm")
}

#' R8 / item 06 arm: Z block replaced by BART -- non-parametric like a forest,
#' fully Bayesian like a parametric draw, so posterior samples are proper by
#' construction. The arm that has to beat BOTH incumbents to close R8.
proc_pipeline_bartMI <- function(bundle, m = 20L, seed = NULL, n_cores = 1L,
                                 outer_sweeps = 3L, margin = "shash",
                                 project_root = NULL, quiet = TRUE) {
  proc_pipeline_block_fcs(bundle, m = m, seed = seed, n_cores = n_cores,
                          outer_sweeps = outer_sweeps, margin = margin,
                          project_root = project_root, quiet = quiet,
                          mid = TRUE, proper_z = FALSE, proper_draw = FALSE,
                          mice_z = FALSE, bart_z = TRUE,
                          label = "pipeline_bartMI")
}

#' Same engine, MID disabled (imputed-Y rows are KEPT and pooled normally).
#' Isolates: is naive Rubin pooling on multiple-imputation-then-deletion data
#' understating variance (the correction von Hippel 2007 derives)?
proc_pipeline_noMID <- function(bundle, m = 20L, seed = NULL, n_cores = 1L,
                                outer_sweeps = 3L, margin = "shash",
                                project_root = NULL, quiet = TRUE) {
  proc_pipeline_block_fcs(bundle, m = m, seed = seed, n_cores = n_cores,
                          outer_sweeps = outer_sweeps, margin = margin,
                          project_root = project_root, quiet = quiet,
                          mid = FALSE, proper_z = FALSE,
                          label = "pipeline_noMID")
}

#' Both changes at once -- if the two mechanisms are additive this should be the
#' best-calibrated arm.
proc_pipeline_properZ_noMID <- function(bundle, m = 20L, seed = NULL, n_cores = 1L,
                                        outer_sweeps = 3L, margin = "shash",
                                        project_root = NULL, quiet = TRUE) {
  proc_pipeline_block_fcs(bundle, m = m, seed = seed, n_cores = n_cores,
                          outer_sweeps = outer_sweeps, margin = margin,
                          project_root = project_root, quiet = quiet,
                          mid = FALSE, proper_z = TRUE,
                          label = "pipeline_properZ_noMID")
}
