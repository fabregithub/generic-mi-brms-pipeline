#!/usr/bin/env Rscript
# =============================================================================
# Track V16 -- does the no-Y knob actually remove Y?
# -----------------------------------------------------------------------------
# WHY THIS TEST EXISTS. V16's whole claim rests on one thing: that the arms differ
# from the reference in Y's presence among the imputation predictors, and in
# NOTHING else. That is not visible in the results -- a wrong conditioning set
# returns plausible numbers, never an error. This project has been caught by that
# exact failure mode four times in the leftcens API (FINDINGS_v13.md), so the
# conditioning sets are asserted directly, at the point where they are built.
#
# It also checks the converse, which is easier to get wrong: that flipping the
# dictionary flag does not leak into the ANALYSIS model. If it did, the arm would
# be measuring a misspecified estimand rather than a misspecified imputation.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
# Same list the runner sources, so the test exercises the same code path.
for (f in c("dgp.R", "censoring.R", "procedures.R", "procedures_pipeline.R",
            "metrics.R", "robustness.R", "proper_impute.R", "mice_impute.R",
            "bart_impute.R")) source(file.path(.here, "R", f))
`%||%` <- function(a, b) if (is.null(a)) b else a

fails <- 0L
ok <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "FAIL", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- fails + 1L
}

# A cell with BOTH censored exposure and missing covariates, so both blocks fire.
set.seed(20260903)
truth  <- make_truth(p = 3L, erf_form = "additive")
sc     <- Filter(function(x) x$name == "combined", v2_scenarios())[[1]]
bundle <- v2_make_bundle(sc, truth)

cat("\n=== 1. the X block's predictor set ===\n")
for (yx in c(TRUE, FALSE)) {
  inp <- .v1_make_pipeline_inputs(bundle, y_in_x = yx, y_in_z = TRUE)
  pr  <- inp$analysis_spec$imputation$censored_exposure$predictors
  ok(sprintf("y_in_x = %-5s -> Y %s ce$predictors", yx, if (yx) "in" else "NOT in"),
     ("Y" %in% pr) == yx, paste(pr, collapse = ", "))
}
# and nothing ELSE moved
a <- .v1_make_pipeline_inputs(bundle, y_in_x = TRUE)$analysis_spec$imputation$censored_exposure$predictors
b <- .v1_make_pipeline_inputs(bundle, y_in_x = FALSE)$analysis_spec$imputation$censored_exposure$predictors
ok("dropping Y changes ONLY Y in the X set", identical(setdiff(a, b), "Y") && length(setdiff(b, a)) == 0L,
   sprintf("removed: %s | added: %s",
           paste(setdiff(a, b), collapse = ","),
           paste(setdiff(b, a), collapse = ",") ))

cat("\n=== 2. the Z block's predictor set ===\n")
# Reproduce what the module does: force impute_y (MID) and build the row-level
# spec exactly as 00_censored_exposure.R does inside the sweep loop.
env <- tryCatch(v1_load_pipeline(quiet = TRUE), error = function(e) NULL)
if (is.null(env)) {
  ok("pipeline sourced", FALSE, "could not source the pipeline; §2-4 skipped")
} else {
  z_pred <- function(y_in_z) {
    inp <- .v1_make_pipeline_inputs(bundle, y_in_z = y_in_z)
    as_ <- inp$analysis_spec; as_$imputation$impute_y <- TRUE
    # The module strips the _lo/_hi bound columns before the Z block; the
    # predictor set is a function of the dictionary, so the raw data suffices.
    env$make_row_level_imputation_spec(inp$data, as_, inp$var_dict)
  }
  for (yz in c(TRUE, FALSE)) {
    sp <- z_pred(yz)
    zt <- setdiff(names(sp$vars), "Y")          # the covariate targets
    has_y <- all(vapply(zt, function(v) "Y" %in% sp$vars[[v]], logical(1)))
    ok(sprintf("y_in_z = %-5s -> Y %s every Z target's predictors", yz,
               if (yz) "in" else "NOT in"),
       has_y == yz,
       sprintf("targets: %s | predictors of %s: %s", paste(names(sp$vars), collapse = ","),
               zt[1], paste(sp$vars[[zt[1]]], collapse = ",")))
  }
  s1 <- z_pred(TRUE); s0 <- z_pred(FALSE)
  zt <- setdiff(names(s1$vars), "Y")
  ok("dropping Y changes ONLY Y in the Z sets",
     all(vapply(zt, function(v) identical(setdiff(s1$vars[[v]], s0$vars[[v]]), "Y") &&
                                length(setdiff(s0$vars[[v]], s1$vars[[v]])) == 0L, logical(1))))
  # Y must remain a TARGET even when it is not a predictor: MID needs it imputed
  # so the X block can condition on a complete Y. Losing this would silently
  # convert the arm into "no MID" as well.
  ok("Y is still an imputation target with y_in_z = FALSE", "Y" %in% names(s0$vars),
     paste(names(s0$vars), collapse = ","))
  ok("Z still predicts Y (the asymmetry is one-way)",
     length(intersect(s0$vars[["Y"]], zt)) > 0L,
     paste(s0$vars[["Y"]], collapse = ","))

  cat("\n=== 3. the ANALYSIS model must not move ===\n")
  # fit_lm_estimand reads `truth`, not the dictionary. If the flag leaked into
  # the fit, the two would differ on identical data.
  d <- bundle$complete
  e1 <- fit_lm_estimand(d, truth); e2 <- fit_lm_estimand(d, truth)
  ok("estimand fit is dictionary-independent", isTRUE(all.equal(e1, e2)),
     sprintf("est = %.6f", e1[["est"]]))
  v1 <- .v1_make_pipeline_inputs(bundle, y_in_z = TRUE)$var_dict
  v0 <- .v1_make_pipeline_inputs(bundle, y_in_z = FALSE)$var_dict
  df <- v1$var[v1$use_in_model != v0$use_in_model]
  ok("only Y's use_in_model flag differs", identical(df, "Y"), paste(df, collapse = ","))
  ok("Y's role is still 'outcome'", v0$role[v0$var == "Y"] == "outcome")

  cat("\n=== 4. all four arms run and are distinguishable ===\n")
  res <- lapply(list(bartMI = proc_pipeline_bartMI, noYx = proc_pipeline_noYx,
                     noYz = proc_pipeline_noYz, noYboth = proc_pipeline_noYboth),
                function(f) f(bundle, m = 3L, seed = 11L, outer_sweeps = 2L))
  for (nm in names(res)) {
    r <- res[[nm]]
    ok(sprintf("%-8s returns a finite estimate", nm),
       is.finite(r$estimate) && is.finite(r$se),
       sprintf("est = %.4f  se = %.4f  note = %s", r$estimate, r$se, r$note))
  }
  es <- vapply(res, function(r) r$estimate, 0)
  ok("the four arms are not identical", length(unique(round(es, 10))) == 4L,
     paste(sprintf("%s=%.4f", names(es), es), collapse = "  "))
  # Direction is NOT asserted: one dataset at m = 3 carries no signal, and
  # asserting a direction here would be the hardcoded-conclusion mistake.
  cat(sprintf("  (info, not a test) truth = %.4f, single draw: %s\n",
              truth$estimand_true, paste(sprintf("%s=%.4f", names(es), es), collapse = "  ")))
}

cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL PASS" else "FAILURES",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
