#!/usr/bin/env Rscript
# =============================================================================
# V17 -- do the covariate roles have the structure they claim?
# -----------------------------------------------------------------------------
# WHY THIS EXISTS. V16 registered a prediction that assumed Z was a confounder in
# the `nl_*` cells. It was not -- Z was a DESCENDANT of the exposures, and the
# partial correlation that mattered was +0.005. Nobody checked, because a DGP
# returns plausible numbers whatever its arrows are. Every role below is
# therefore verified by its CONSEQUENCES: which analysis is unbiased, and which
# is biased. A role that cannot be caught getting this wrong is not a role.
#
# The regression guard at the end matters just as much: `simulate_complete()`
# gained an argument, and every earlier track's results must still reproduce from
# their recorded seeds.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.here, "R", "dgp.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

fails <- 0L
ok <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "FAIL", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- fails + 1L
}
N <- 400000L
beta <- function(d, zterms) {
  rhs <- c(paste0("logX", 1:3), zterms)
  unname(coef(lm(as.formula(paste("Y ~", paste(rhs, collapse = "+"))), data = d))[["logX1"]])
}
near <- function(a, b, tol) abs(a - b) < tol

cat("\n=== 1. every role keeps b1 = 0.40 recoverable by ITS OWN correct model ===\n")
for (role in c("precision", "fork", "pipe", "collider", "mixed")) {
  set.seed(101)
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role)
  d  <- simulate_complete(N, tr)$data
  b  <- beta(d, tr$z_in_model)
  ok(sprintf("%-10s correct model recovers 0.40", role), near(b, 0.40, 0.02),
     sprintf("b1 = %.4f  (adjusts for %s)", b,
             paste(tr$z_in_model, collapse = "+") ))
}

cat("\n=== 2. each role is BROKEN by the wrong adjustment (the actual test) ===\n")
set.seed(202)
tr <- make_truth(p = 3L, erf_form = "additive", z_role = "fork")
d  <- simulate_complete(N, tr)$data
b_adj <- beta(d, c("Z1", "Z2")); b_un <- beta(d, "Z2")
ok("fork: OMITTING Z1 biases b1", !near(b_un, 0.40, 0.02) && near(b_adj, 0.40, 0.02),
   sprintf("adjusted %.4f  vs  unadjusted %.4f", b_adj, b_un))

set.seed(203)
tr <- make_truth(p = 3L, erf_form = "additive", z_role = "collider")
d  <- simulate_complete(N, tr)$data
b_adj <- beta(d, c("Z1", "Z2")); b_un <- beta(d, "Z2")
ok("collider: INCLUDING Z1 biases b1", near(b_un, 0.40, 0.02) && !near(b_adj, 0.40, 0.02),
   sprintf("without Z1 %.4f  vs  with Z1 %.4f", b_un, b_adj))
ok("collider: Z1 is not a cause of Y (gamma1 = 0)", tr$gamma[1] == 0)
ok("collider: Z1 is dropped from the analysis formula",
   !"Z1" %in% all.vars(dgp_formula(tr)), deparse(dgp_formula(tr)))

set.seed(204)
tr <- make_truth(p = 3L, erf_form = "additive", z_role = "pipe")
d  <- simulate_complete(N, tr)$data
b_dir <- beta(d, c("Z1", "Z2")); b_tot <- beta(d, "Z2")
ok("pipe: adjusted = DIRECT effect (0.40)", near(b_dir, 0.40, 0.02),
   sprintf("%.4f", b_dir))
ok("pipe: unadjusted = TOTAL effect (b1 + delta*gamma1)",
   near(b_tot, tr$total_effect, 0.03),
   sprintf("%.4f vs total_effect %.4f", b_tot, tr$total_effect))

set.seed(205)
tr <- make_truth(p = 3L, erf_form = "additive", z_role = "mixed")
d  <- simulate_complete(N, tr)$data
b_adj <- beta(d, c("Z1", "Z2")); b_noz1 <- beta(d, "Z2"); b_noz2 <- beta(d, "Z1")
ok("mixed: dropping the fork (Z1) biases b1", !near(b_noz1, 0.40, 0.02),
   sprintf("%.4f", b_noz1))
ok("mixed: dropping the pipe (Z2) shifts b1 toward the total effect",
   !near(b_noz2, b_adj, 0.02),
   sprintf("full %.4f  vs  no-Z2 %.4f", b_adj, b_noz2))

cat("\n=== 3. the arrows are where they are claimed to be ===\n")
pcor <- function(d, a, b, given) {
  ra <- resid(lm(as.formula(sprintf("%s ~ %s", a, paste(given, collapse = "+"))), data = d))
  rb <- resid(lm(as.formula(sprintf("%s ~ %s", b, paste(given, collapse = "+"))), data = d))
  cor(ra, rb)
}
set.seed(301)
for (role in c("precision", "fork", "pipe", "collider")) {
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role)
  d  <- simulate_complete(N, tr)$data
  pc <- pcor(d, "Z1", "logX1", c("logX2", "logX3", "Z2"))
  want_dep <- role != "precision"
  ok(sprintf("%-9s Z1 %s logX1 given the rest", role,
             if (want_dep) "IS associated with" else "is INDEPENDENT of"),
     (abs(pc) > 0.10) == want_dep, sprintf("partial cor = %+.4f", pc))
}

cat("\n=== 4. regression guard: the pre-V17 path is unchanged ===\n")
for (zf in c("linear", "nonlinear")) {
  set.seed(999)
  tr <- make_truth(p = 3L, erf_form = "additive")
  a <- simulate_complete(800L, tr, z_form = zf)$data
  set.seed(999)
  b <- simulate_complete(800L, tr, z_form = zf)$data
  ok(sprintf("z_form = %-9s is deterministic under a fixed seed", zf),
     identical(a, b))
}
set.seed(999); tr <- make_truth(p = 3L, erf_form = "additive")
d1 <- simulate_complete(800L, tr)$data
ok("default z_role is 'precision'", identical(tr$z_role, "precision"))
ok("default truth keeps both Z in the analysis model",
   identical(tr$z_in_model, c("Z1", "Z2")), paste(tr$z_in_model, collapse = ","))
ok("default formula is unchanged",
   identical(deparse(dgp_formula(tr)), "Y ~ logX1 + logX2 + logX3 + Z1 + Z2"),
   deparse(dgp_formula(tr)))
ok("estimand_true is still b[1] in every role",
   all(vapply(c("precision","fork","pipe","collider","mixed"),
              function(r) make_truth(p=3L, z_role=r)$estimand_true == 0.40, logical(1))))

cat("\n=== 5. the auxiliary knob, and the defect it measures ===\n")
for (f in c("censoring.R", "metrics.R", "procedures.R", "procedures_pipeline.R",
            "robustness.R", "proper_impute.R", "mice_impute.R", "bart_impute.R"))
  source(file.path(.here, "R", f))
env <- tryCatch(v1_load_pipeline(quiet = TRUE), error = function(e) NULL)
if (is.null(env)) {
  ok("pipeline sourced", FALSE, "could not source the pipeline; section 5 skipped")
} else {
  set.seed(501)
  trc <- make_truth(p = 3L, erf_form = "additive", z_role = "collider")
  scc <- Filter(function(x) x$name == "zr_collider", v2_scenarios())[[1]]
  bun <- v2_make_bundle(scc, trc)

  vd <- .v1_make_pipeline_inputs(bun, z_aux = TRUE)$var_dict
  ok("z_aux marks Z1 auxiliary and drops it from the model",
     !vd$use_in_model[vd$var == "Z1"] && vd$use_as_auxiliary[vd$var == "Z1"])
  ok("z_aux touches nothing but Z1",
     identical(vd$var[vd$use_in_model != .v1_make_pipeline_inputs(bun)$var_dict$use_in_model], "Z1"))

  # Z block: does the auxiliary reach it? (Documented behaviour -- it should.)
  as_ <- .v1_make_pipeline_inputs(bun, z_aux = TRUE)$analysis_spec
  as_$imputation$impute_y <- TRUE
  inp <- .v1_make_pipeline_inputs(bun, z_aux = TRUE)
  sp  <- env$make_row_level_imputation_spec(inp$data, as_, inp$var_dict)
  zt  <- setdiff(names(sp$vars), c("Y", "Z1"))
  ok("Z BLOCK uses the auxiliary Z1", length(zt) > 0 && "Z1" %in% sp$vars[[zt[1]]],
     sprintf("predictors of %s: %s", zt[1], paste(sp$vars[[zt[1]]], collapse = ",")))

  # X block: reproduce 00_censored_exposure.R's auto_preds exactly.
  vdx <- inp$var_dict
  auto <- unique(c("Y", vdx$var[vdx$use_in_model %in% TRUE], inp$expo_names))
  ok("X BLOCK auto predictors DROP the auxiliary Z1 (the defect)",
     !"Z1" %in% auto, paste(auto, collapse = ","))

  ok("auto_x_preds hands the module NULL so it computes its own set",
     is.null(.v1_make_pipeline_inputs(bun, z_aux = TRUE,
                                      auto_x_preds = TRUE)$analysis_spec$imputation$censored_exposure$predictors))

  a <- proc_pipeline_auxZ_shipped(bun, m = 3L, seed = 21L, outer_sweeps = 2L)
  b <- proc_pipeline_auxZ_asdoc(bun,   m = 3L, seed = 21L, outer_sweeps = 2L)
  ok("both auxiliary arms run", is.finite(a$estimate) && is.finite(b$estimate),
     sprintf("shipped %.4f  as-documented %.4f  (truth %.2f)",
             a$estimate, b$estimate, trc$estimand_true))
  ok("the two arms differ (the defect has an effect to measure)",
     !isTRUE(all.equal(a$estimate, b$estimate)))
}

cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL PASS" else "FAILURES",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
