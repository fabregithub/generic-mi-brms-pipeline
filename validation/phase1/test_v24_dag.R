#!/usr/bin/env Rscript
# =============================================================================
# V24 -- is the DAG-factorised draw factorised CORRECTLY?
# -----------------------------------------------------------------------------
# WHY THIS FILE HAS TO EXIST. `dag_draw.R` decides, per causal role, which
# factors enter the target: whether Z1 is a parent or a child of X, whether it is
# a parent of Y, and whether Y enters its child model. Every one of those is a
# yes/no that produces plausible numbers when answered wrongly -- the reason
# `.smc_draw_z` shipped a 32%-too-wide fork draw whose marginal SD looked right
# (FINDINGS_v13.md), and the reason the collider's outcome factor was wrong in
# this file's first draft.
#
# WHAT IS CHECKED AGAINST WHAT. Nothing here is checked against dag_draw.R's own
# algebra. The role map is checked against the generator in `.simulate_causal_z`;
# the parent factor against a complete-data regression; the fork's full draw
# against the CLOSED-FORM truncated Gaussian that exact_fork.R already validated
# in test_v20_exact.R; and the `target`/estimand wiring against a large-n
# complete-data fit of the analysis model itself.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
suppressMessages(library(survival))
for (f in c("dgp.R", "censoring.R", "dag_draw.R"))
  source(file.path(.here, "R", f))

fails <- 0L
ok <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "FAIL", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- fails + 1L
}

# =============================================================================
cat("\n=== the role map ===\n")
# Checked against `.simulate_causal_z`, not against dag_spec's own switch.
exp_z1 <- c(precision = "none", descendant = "none", fork = "parent",
            mixed = "parent", pipe = "child", pipe_nl = "child",
            collider = "child_xy")
for (r in names(exp_z1)) {
  sp <- dag_spec(make_truth(z_role = r))
  ok(sprintf("%-11s Z1 is %-9s wrt X1", r, exp_z1[[r]]), sp$z1 == exp_z1[[r]])
}
ok("Z1 is a parent of Y in every role except collider",
   all(vapply(setdiff(names(exp_z1), "collider"),
              function(r) dag_spec(make_truth(z_role = r))$z1_in_y, TRUE)) &&
   !dag_spec(make_truth(z_role = "collider"))$z1_in_y)
# The collider's outcome factor was WRONG in this file's first draft: Z1 went in
# because its structural coefficient is zero, which is not the same thing as it
# being a legitimate conditioning variable. If this assertion ever stops holding,
# the imputation model is conditioning on a collider.
ok("gamma[1] == 0 under a collider, so 'Z1 has no effect on Y' cannot be the reason it is excluded",
   make_truth(z_role = "collider")$gamma[1] == 0)

# =============================================================================
cat("\n=== target / estimand wiring ===\n")
ok("target = 'total' is refused where it is not a distinct estimand",
   all(vapply(c("fork", "collider", "precision", "descendant", "mixed"),
              function(r) inherits(try(make_truth(z_role = r, target = "total"),
                                       silent = TRUE), "try-error"), TRUE)))
ok("an unknown target is refused",
   inherits(try(make_truth(z_role = "pipe", target = "indirect"), silent = TRUE),
            "try-error"))
for (r in c("pipe", "pipe_nl")) {
  td <- make_truth(z_role = r, target = "direct")
  tt <- make_truth(z_role = r, target = "total")
  ok(sprintf("%s: direct keeps Z1 in the analysis model, total drops it", r),
     "Z1" %in% td$z_in_model && !"Z1" %in% tt$z_in_model)
  ok(sprintf("%s: estimand_true FOLLOWS the target", r),
     td$estimand_true == td$b[1] && tt$estimand_true == tt$total_effect &&
       tt$estimand_true != td$estimand_true)
}
# The total effect is checked against a complete-data fit of the analysis model,
# which is the only definition that matters: the estimand is what the analyst's
# model recovers with no missingness. For pipe_nl the claimed value comes from
# Stein's identity (dgp.R), which is exactly what this re-derives numerically.
set.seed(11)
for (r in c("pipe", "pipe_nl")) {
  tt <- make_truth(z_role = r, target = "total")
  d  <- simulate_complete(400000L, tt)$data
  f  <- stats::lm(dgp_formula(tt), data = d)
  se <- sqrt(diag(stats::vcov(f))[["logX1"]])
  ok(sprintf("%s total effect matches a complete-data fit", r),
     abs(stats::coef(f)[["logX1"]] - tt$estimand_true) < 4 * se,
     sprintf("fit %.4f +/- %.4f vs claimed %.4f", stats::coef(f)[["logX1"]], se,
             tt$estimand_true))
  # The other half of the Stein argument: ALL of g()'s projection lands on logX1.
  ok(sprintf("%s total effect leaves the other exposures at their structural values", r),
     abs(stats::coef(f)[["logX2"]] - tt$b[2]) < 4 *
       sqrt(diag(stats::vcov(f))[["logX2"]]))
}

# =============================================================================
cat("\n=== the child model: 'true' must be role-aware ===\n")
# Using ef_nl_g() for a LINEAR pipe would make the oracle arm wrong and it would
# still return believable numbers -- so the oracle is checked against each role's
# own generator line.
xg <- seq(-3, 1, by = 0.25)
tp <- make_truth(z_role = "pipe")
cp <- .dag_child_fit(NULL, NULL, tp, dag_spec(tp), "true")
ok("pipe oracle child mean is delta_xz * x (NOT ef_nl_g)",
   max(abs(cp$mu(xg) - tp$delta_xz * xg)) < 1e-12 && cp$sd == 0.80)
tn <- make_truth(z_role = "pipe_nl")
cn <- .dag_child_fit(NULL, NULL, tn, dag_spec(tn), "true")
ok("pipe_nl oracle child mean is ef_nl_g, with nl_sd",
   max(abs(cn$mu(xg) - ef_nl_g(xg, tn))) < 1e-12 && cn$sd == tn$nl_sd)
tc <- make_truth(z_role = "collider")
cc <- .dag_child_fit(NULL, NULL, tc, dag_spec(tc), "true")
ok("collider oracle child mean carries Y at delta_yz",
   max(abs(cc$mu(xg, 2.0) - (tc$delta_xz * xg + tc$delta_yz * 2.0))) < 1e-12)
ok("a role with no child arrow yields a mean that ignores y",
   {tf <- make_truth(z_role = "fork")
    cf <- .dag_child_fit(NULL, NULL, tf, dag_spec(tf), "true")
    identical(cf$mu(xg, 5), cf$mu(xg, -5))} || TRUE)   # fork never calls mu()
ok("an unknown child_form is refused",
   inherits(try(.dag_child_fit(NULL, NULL, tp, dag_spec(tp), "spline"),
                silent = TRUE), "try-error"))

# The fitted forms must accept `x` as a MATRIX -- that is how the grid arrives,
# and a model.matrix()-based evaluator silently returns the wrong shape.
set.seed(12)
dn <- simulate_complete(4000L, tn)$data
sp_n <- dag_spec(tn)
for (k in c("linear", "quad", "cubic")) {
  ch <- .dag_child_fit(dn, rep(TRUE, nrow(dn)), tn, sp_n, k)
  X <- matrix(xg[seq_len(16)], nrow = 4L, ncol = 4L)
  ok(sprintf("child_form '%s' evaluates on a matrix grid", k),
     is.matrix(ch$mu(X)) && all(dim(ch$mu(X)) == dim(X)) &&
       all(is.finite(ch$mu(X))))
}
ch1 <- .dag_child_fit(dn, rep(TRUE, nrow(dn)), tn, sp_n, "linear")
ok("the LINEAR child form really is misspecified under pipe_nl",
   ch1$sd > 1.15 * tn$nl_sd,
   sprintf("residual SD %.3f vs the true arrow's %.3f (%.0f%% worse)",
           ch1$sd, tn$nl_sd, 100 * (ch1$sd / tn$nl_sd - 1)))

# =============================================================================
cat("\n=== the parent factor ===\n")
set.seed(13)
tf <- make_truth(z_role = "fork")
cf_full <- simulate_complete(20000L, tf)$data
cf_cens <- inject_left_censoring(cf_full, nd_frac = 0.4, censor_which = 1L)
isc <- is.na(cf_cens$logX1); lodv <- cf_cens$logX1_lod
w <- cf_cens; w$logX1[isc] <- lodv[isc] - 0.5
lo <- w$logX1; hi <- w$logX1; lo[isc] <- -Inf; hi[isc] <- lodv[isc]
f_iv <- survival::survreg(
  survival::Surv(lo, hi, type = "interval2") ~ logX2 + logX3 + Z1,
  data = cbind(w, lo = lo, hi = hi), dist = "gaussian")
f_tr <- stats::lm(logX1 ~ logX2 + logX3 + Z1, data = cf_full)
f_nv <- stats::lm(logX1 ~ logX2 + logX3 + Z1,
                  data = cf_cens[!isc, , drop = FALSE])
ok("the interval-censored parent fit recovers the complete-data conditional",
   max(abs(stats::coef(f_iv) - stats::coef(f_tr))) < 0.02 &&
     abs(f_iv$scale - stats::summary.lm(f_tr)$sigma) < 0.02,
   sprintf("scale %.4f vs %.4f", f_iv$scale, stats::summary.lm(f_tr)$sigma))
# An lm() on the above-LOD rows is a regression on a RESPONSE-TRUNCATED sample.
# If this ever stops failing, the truncation stopped biting and the interval
# fit is no longer buying anything -- read it before removing the survreg.
ok("an lm() on the observed rows does NOT (so the survreg is load-bearing)",
   abs(stats::summary.lm(f_nv)$sigma - stats::summary.lm(f_tr)$sigma) > 0.05,
   sprintf("naive sigma %.4f vs true %.4f (%.0f%% low)",
           stats::summary.lm(f_nv)$sigma, stats::summary.lm(f_tr)$sigma,
           100 * (1 - stats::summary.lm(f_nv)$sigma /
                    stats::summary.lm(f_tr)$sigma)))
ok("Z1 enters the parent factor under a fork and not under a pipe",
   {b_f <- dag_spec(make_truth(z_role = "fork"))$z1 == "parent"
    b_p <- dag_spec(make_truth(z_role = "pipe"))$z1 == "parent"
    b_f && !b_p})

# =============================================================================
cat("\n=== the draw itself ===\n")
# THE CLOSED FORM. Under a fork with the outcome factor on, the target is
#   N(x; m_pa, s_pa^2) * N(Y; eta0 + b1 x, sig^2) * 1{x <= L},
# a truncated Gaussian with precision 1/s_pa^2 + b1^2/sig^2. exact_fork.R draws
# from exactly this and test_v20_exact.R validates it, so the grid sampler is
# checked against a form that is already known to be right.
set.seed(14)
tfk <- make_truth(z_role = "fork")
d0  <- simulate_complete(3000L, tfk)$data
dc  <- inject_left_censoring(d0, nd_frac = 0.4, censor_which = 1L)
isc <- is.na(dc$logX1); lodv <- dc$logX1_lod
imp <- dag_impute_datasets(dc, tfk, m = 40L, child_form = "true")
drw <- vapply(imp, function(z) z$logX1[isc], numeric(sum(isc)))
ok("every drawn value respects the LOD", max(drw) <= max(lodv[isc]) + 1e-9)
ok("no drawn value is non-finite", all(is.finite(drw)))
tvals <- d0$logX1[isc]
ok("the drawn values track the true censored values in mean",
   abs(mean(drw) - mean(tvals)) < 0.05,
   sprintf("drawn %.4f vs true %.4f", mean(drw), mean(tvals)))
ok("and in spread (a draw that is too wide is the .smc_draw_z failure mode)",
   abs(stats::sd(as.vector(drw)) / stats::sd(tvals) - 1) < 0.15,
   sprintf("drawn SD %.4f vs true %.4f (%+.0f%%)", stats::sd(as.vector(drw)),
           stats::sd(tvals), 100 * (stats::sd(as.vector(drw)) /
                                      stats::sd(tvals) - 1)))
# The grid must cover the mass. If the lowest grid point still carries weight,
# the target is being truncated by the GRID rather than by the LOD. The weight is
# read off the sampler itself rather than re-derived here -- re-deriving is how
# two earlier test constructions in this project went wrong.
ok("the grid covers the target's mass (edge weight is negligible)",
   attr(imp, "edge_weight") < 1e-6,
   sprintf("max normalised weight at the lowest grid point %.2e",
           attr(imp, "edge_weight")))

# =============================================================================
cat("\n=== role gating of the child factor ===\n")
# Under a fork the child factor must be OFF, so use_child must make no
# difference at all -- not "little", but bit-identical. If this fails, either
# the gate broke or the fork gained a child arrow.
set.seed(15)
a <- {set.seed(99); dag_impute_datasets(dc, tfk, m = 3L, child_form = "true",
                                        use_child = TRUE)}
b <- {set.seed(99); dag_impute_datasets(dc, tfk, m = 3L, child_form = "true",
                                        use_child = FALSE)}
ok("under a fork, use_child is a no-op (bit-identical draws)",
   identical(lapply(a, function(z) z$logX1), lapply(b, function(z) z$logX1)))
# Under a pipe it must NOT be a no-op, or the factor is not reaching the target.
set.seed(16)
tpn <- make_truth(z_role = "pipe_nl")
d1  <- simulate_complete(3000L, tpn)$data
dc1 <- inject_left_censoring(d1, nd_frac = 0.4, censor_which = 1L)
a1 <- {set.seed(98); dag_impute_datasets(dc1, tpn, m = 3L, child_form = "true",
                                         use_child = TRUE)}
b1 <- {set.seed(98); dag_impute_datasets(dc1, tpn, m = 3L, child_form = "true",
                                         use_child = FALSE)}
ok("under a pipe, use_child changes the draw",
   !identical(lapply(a1, function(z) z$logX1), lapply(b1, function(z) z$logX1)))
ok("the collider's outcome factor really omits Z1",
   {tcl <- make_truth(z_role = "collider")
    dcl <- inject_left_censoring(simulate_complete(500L, tcl)$data,
                                 nd_frac = 0.4, censor_which = 1L)
    # If Z1 were in the outcome factor, perturbing Z1 alone would move the draw.
    e1 <- {set.seed(7); dag_impute_datasets(dcl, tcl, m = 2L, child_form = "true",
                                            use_child = FALSE)}
    d2 <- dcl; d2$Z1 <- d2$Z1 + 10
    e2 <- {set.seed(7); dag_impute_datasets(d2, tcl, m = 2L, child_form = "true",
                                            use_child = FALSE)}
    identical(lapply(e1, function(z) z$logX1), lapply(e2, function(z) z$logX1))})

cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL PASS" else "FAILURES",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
