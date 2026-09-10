#!/usr/bin/env Rscript
# =============================================================================
# Unit tests for 00_censored_exposure.R
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS, when the project's other tests are end-to-end bash
# scripts. `ce_auto_predictors()` is the 2026-09-09 fix that puts auxiliary and
# out-of-model covariates back into the exposure draw, and NO BUNDLED EXAMPLE
# REACHES IT -- the censored-exposure example passes `predictors` explicitly,
# which overrides it entirely. So without this file the shipped behaviour has no
# coverage and a regression would be invisible to the suite.
#
# An end-to-end smoke test catches none of those. Run with:
#   Rscript test/test_dag_child_factor.R
# =============================================================================
setwd(dirname(dirname(normalizePath(
  sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))))
`%||%` <- function(a, b) if (is.null(a)) b else a
log_msg <- function(...) invisible(NULL)          # quiet
source("00_censored_exposure.R")

fails <- 0L
ok <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "FAIL", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- fails + 1L
}

set.seed(20260909)
n <- 1200
xm <- rnorm(n)                                     # exposure, modelling scale
g  <- function(x) 1.2 * tanh(1.8 * x) + 0.35 * (x^2 - 1)
Z  <- g(xm) + rnorm(n, 0, 0.8)                     # a CHILD of the exposure
Y  <- 0.4 * xm + 0.5 * Z + rnorm(n)
lod <- as.numeric(quantile(xm, 0.4))
cens <- xm < lod
xobs <- ifelse(cens, NA_real_, xm)

cat("\n=== the automatic predictor set (the 2026-09-09 fix) ===\n")
vd <- data.frame(
  var              = c("Y", "Xexp", "Zmodel", "Zaux", "Zauxrole", "Zunused"),
  use_in_model     = c(FALSE, TRUE,  TRUE,     FALSE,  FALSE,      FALSE),
  use_as_auxiliary = c(FALSE, FALSE, FALSE,    TRUE,   FALSE,      FALSE),
  role             = c("outcome", "exposure", "covariate", "covariate",
                       "auxiliary", "covariate"),
  stringsAsFactors = FALSE)
preds <- ce_auto_predictors(vd, "Y", "Xexp")
ok("the outcome is always in (congeniality -- the point of the strategy)",
   "Y" %in% preds)
ok("in-model covariates are in", "Zmodel" %in% preds)
ok("use_as_auxiliary = TRUE is in (this is the fix)", "Zaux" %in% preds)
ok("role = 'auxiliary' is in (matches auxiliary_vars in 00_common_functions.R)",
   "Zauxrole" %in% preds)
ok("a covariate that is neither in-model nor auxiliary stays OUT",
   !"Zunused" %in% preds)
ok("the exposures are in", "Xexp" %in% preds)
ok("missing dictionary columns do not error (both are optional)",
   {v2 <- vd[, c("var", "use_in_model")]
    identical(sort(ce_auto_predictors(v2, "Y", "Xexp")),
              sort(c("Y", "Xexp", "Zmodel")))})
ok("NA flags are treated as FALSE, not dropped",
   {v3 <- vd; v3$use_as_auxiliary[4] <- NA
    !"Zaux" %in% ce_auto_predictors(v3, "Y", "Xexp")})

cat("\n=== the curvature preflight ===\n")
# Built from a DGP where the answer is known: Zcurve is a curved child of the
# exposure, Zlin a straight one, Znoise unrelated. The screen must separate
# them, and must not fire on the two that are fine.
set.seed(4)
n <- 800
x  <- rnorm(n)
lod <- as.numeric(quantile(x, 0.35))
cens <- x < lod
dat <- data.frame(
  y      = 0.4 * x + rnorm(n),
  X      = exp(x),
  X_lo   = ifelse(cens, 0, exp(x)),
  X_hi   = ifelse(cens, exp(lod), exp(x)),
  Zcurve = 1.2 * tanh(1.8 * x) + 0.35 * (x^2 - 1) + rnorm(n, 0, 0.5),
  Zlin   = 0.6 * x + rnorm(n, 0, 0.5),
  Znoise = rnorm(n),
  rowid  = seq_len(n))
vdp <- data.frame(
  var  = c("y", "X", "Zcurve", "Zlin", "Znoise", "rowid"),
  role = c("outcome", "exposure", "covariate", "covariate", "covariate", "id"),
  scale = c("no", "log", "no", "no", "no", "no"),
  stringsAsFactors = FALSE)
as_p <- list(outcome = list(y_var = "y"),
             imputation = list(censored_exposure = list(
               exposure_vars = "X", log_scale = TRUE,
               lo_suffix = "_lo", hi_suffix = "_hi")))
pf <- suppressMessages(ce_preflight_curvature(dat, as_p, vdp))
gv <- function(z, f) pf[[f]][pf$covariate == z]
ok("the preflight returns one row per screened covariate", nrow(pf) == 3L,
   paste(pf$covariate, collapse = ", "))
ok("a variable the dictionary does not call a covariate is NOT screened",
   !"rowid" %in% pf$covariate,
   "row ids and bookkeeping columns would otherwise fill the report")
ok("the CURVED child is flagged: real association and real curvature",
   gv("Zcurve", "abs_cor") > 0.10 && gv("Zcurve", "curvature") > 0.05,
   sprintf("|cor| %.2f, curvature %.3f", gv("Zcurve", "abs_cor"),
           gv("Zcurve", "curvature")))
ok("the STRAIGHT child is not flagged (it is the safe case)",
   !(gv("Zlin", "curvature") > 0.05 && gv("Zlin", "p_curv") < 0.01),
   sprintf("|cor| %.2f, curvature %.3f, p %.2g", gv("Zlin", "abs_cor"),
           gv("Zlin", "curvature"), gv("Zlin", "p_curv")))
# THE STATISTIC ITSELF, against the DGP whose bias is known. The first version
# (relative residual-SD reduction) read 0.0165 on this case -- below any
# sensible threshold, so the screen would have missed the one cell it exists
# for. Asserted so a future simplification cannot quietly reintroduce that.
ok("the statistic clears the threshold on a curved arrow by a wide margin",
   gv("Zcurve", "curvature") > 5 * 0.05,
   sprintf("%.3f against a 0.05 threshold", gv("Zcurve", "curvature")))
# The EFFECT SIZE alone is not enough at small n -- a 4-df spline fits noise
# better than a line, and the df correction removes only the expected excess,
# not its sampling noise. The F test is what holds the false-positive rate
# down, so it is asserted separately. Measured over 60 reps of the known-good
# linear DGP: 0% flagged at n = 800, 3% at n = 3200; detection on the
# known-bad DGP 75% and 100%.
ok("the F test separates a straight arrow from a curved one",
   gv("Zlin", "p_curv") > 0.01 && gv("Zcurve", "p_curv") < 1e-6,
   sprintf("straight p %.2g, curved p %.1e", gv("Zlin", "p_curv"),
           gv("Zcurve", "p_curv")))
ok("both the effect size and the p value are reported, so neither gate is hidden",
   all(c("curvature", "p_curv") %in% names(pf)) && !anyNA(pf$p_curv))
ok("an unrelated covariate is not flagged",
   gv("Znoise", "abs_cor") <= 0.10)
ok("the non-detect rate is reported and correct",
   all(abs(pf$nd_rate - mean(cens)) < 1e-9),
   sprintf("%.0f%%", 100 * pf$nd_rate[1]))
ok("the curved child sorts to the top (the report is read from the top)",
   pf$covariate[1] == "Zcurve")
ok("a config with no censored exposures returns NULL, not an error",
   is.null(ce_preflight_curvature(dat, list(outcome = list(y_var = "y"),
     imputation = list()), vdp)))

cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL PASS" else "FAILURES",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
