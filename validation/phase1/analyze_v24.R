#!/usr/bin/env Rscript
# =============================================================================
# Track V24 -- analysis against the criteria registered in run_v24_dagdraw.sh
# -----------------------------------------------------------------------------
# EVERY VERDICT BELOW IS COMPUTED. Nothing is pre-written, and no gate is
# reported as met unless the number that meets it is printed next to it. Two
# earlier tracks resolved a gate for the wrong reason -- V22's `bart >= +3%`
# missed at +2.84% while its pre-written interpretation stood, and V23's "CI
# contains 2" passed by being uninformative (the CI contained 1 too). PLAN §11
# came out of those, and this script implements it: each gate carries an
# INFORMATIVENESS check that can void a pass.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.here, "R", "dgp.R"))
source(file.path(.here, "R", "robustness.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

TAG <- Sys.getenv("TAG", "v24")
NS  <- as.integer(strsplit(Sys.getenv("NLEVELS", "800,3200,12800"), ",")[[1]])

# The truth for each cell comes from the SCENARIO DEFINITION, not from a literal
# typed here -- C4 exists because a hardcoded 0.400 would score a correct
# total-effect cell as a 43% bias, and a hardcoded table could drift the same way.
scs <- v2_scenarios()
scs <- scs[vapply(scs, function(s) identical(s$axis, "dag_draw"), TRUE)]
tv <- vapply(scs, function(s) make_truth(p = 3L, erf_form = s$erf_form,
                                         z_role = s$z_role,
                                         target = s$target %||% "direct",
                                         nl_a = s$nl_a, nl_c = s$nl_c
                                         )$estimand_true, 0)
names(tv) <- vapply(scs, `[[`, "", "name")
CELLS <- names(tv)
LADDER <- c("dag_noY", "dag_Yonly", "dag_lin", "dag_quad", "dag_cubic", "dag_true")

load_stage <- function(n) {
  f <- file.path(.here, "results", sprintf("%s_n%d_latest.rds", TAG, n))
  if (!file.exists(f)) return(NULL)
  s <- as.data.frame(readRDS(f)$summary)
  s$n <- n
  s$truth <- tv[as.character(s$scenario)]
  s$bias      <- s$mean_est - s$truth
  s$bias_pct  <- 100 * s$bias / s$truth
  s$bias_se   <- s$bias / s$claimed_se        # SIGNED: C5 reads the direction
  # Monte-Carlo error of the mean bias, in the same currency as the gate.
  s$mcse_bse  <- (s$emp_se / sqrt(s$n_ok)) / s$claimed_se
  s
}
all <- do.call(rbind, Filter(Negate(is.null), lapply(NS, load_stage)))
if (is.null(all) || !nrow(all)) stop("no V24 results found for TAG=", TAG)

g <- function(d, cell, arm, f) {
  v <- d[[f]][d$scenario == cell & d$procedure == arm]
  if (length(v) != 1L) NA_real_ else v
}

verdicts <- list()
for (n in sort(unique(all$n))) {
  d <- all[all$n == n, , drop = FALSE]
  cat(sprintf("\n=====================================================================\n"))
  cat(sprintf("  n = %d   (%d reps)\n", n, max(d$n_ok, na.rm = TRUE)))
  cat(sprintf("=====================================================================\n"))
  cat(sprintf("%-13s %-10s %7s %9s %9s %8s %7s %7s\n", "cell", "arm", "truth",
              "bias%", "bias/SE", "mcse", "cover", "w/SE"))
  for (cl in intersect(CELLS, unique(as.character(d$scenario)))) {
    for (ar in c("oracle", LADDER)) {
      r <- d[d$scenario == cl & d$procedure == ar, , drop = FALSE]
      if (!nrow(r)) next
      cat(sprintf("%-13s %-10s %7.4f %+9.2f %+9.2f %8.2f %7.3f %7.2f\n",
                  cl, ar, r$truth, r$bias_pct, r$bias_se, r$mcse_bse,
                  r$coverage, r$width_se_ratio))
    }
    cat("\n")
  }

  # ---- C1: the factorisation, with the ORACLE child model -------------------
  # Written so a vacuous pass is impossible: a cell counts as evidence only if
  # its dag_noY arm actually HAS bias to remove.
  c1 <- data.frame(cell = CELLS,
                   true_bse = vapply(CELLS, function(c0) g(d, c0, "dag_true", "bias_se"), 0),
                   noy_bse  = vapply(CELLS, function(c0) g(d, c0, "dag_noY",  "bias_se"), 0))
  c1$informative <- abs(c1$noy_bse) > 0.3
  c1$passes      <- abs(c1$true_bse) < 0.3
  cat("--- C1  dag_true: |bias/SE| < 0.3 in all six cells ------------------\n")
  for (i in seq_len(nrow(c1))) cat(sprintf("  %-13s dag_true %+.2f  dag_noY %+.2f  %s%s\n",
    c1$cell[i], c1$true_bse[i], c1$noy_bse[i],
    if (isTRUE(c1$passes[i])) "under 0.3" else "OVER 0.3",
    if (isTRUE(c1$informative[i])) "" else "  <- VACUOUS: dag_noY has no bias to remove"))
  c1_met <- all(c1$passes, na.rm = TRUE) && all(c1$informative, na.rm = TRUE)
  cat(sprintf("  C1: %s\n", if (isTRUE(c1_met)) "MET" else
    if (all(c1$passes, na.rm = TRUE)) "NOT MET -- passes but on uninformative cells" else "NOT MET"))

  # ---- C2: the child factor is load-bearing, and only where predicted ------
  c2_cells <- c("dag_pnl_dir", "dag_pnl_tot")
  gain <- vapply(c2_cells, function(c0)
    abs(g(d, c0, "dag_Yonly", "bias_se")) - abs(g(d, c0, "dag_true", "bias_se")), 0)
  fork_same <- isTRUE(all.equal(g(d, "dag_fork", "dag_Yonly", "bias_se"),
                                g(d, "dag_fork", "dag_true", "bias_se")))
  cat("\n--- C2  child factor gains >= 1.0 bias/SE in pipe_nl; no-op under fork ---\n")
  for (c0 in c2_cells) cat(sprintf("  %-13s Yonly %+.2f -> true %+.2f   gain %+.2f\n",
    c0, g(d, c0, "dag_Yonly", "bias_se"), g(d, c0, "dag_true", "bias_se"), gain[[c0]]))
  cat(sprintf("  %-13s Yonly %+.2f vs true %+.2f   %s\n", "dag_fork",
              g(d, "dag_fork", "dag_Yonly", "bias_se"),
              g(d, "dag_fork", "dag_true", "bias_se"),
              if (fork_same) "identical, as the role gate requires" else
                "DIFFERENT -- the fork gate leaked"))
  c2_met <- all(gain >= 1.0, na.rm = TRUE) && fork_same
  cat(sprintf("  C2: %s\n", if (isTRUE(c2_met)) "MET" else "NOT MET"))

  # ---- C3: is a misspecified child factor worse than none? ------------------
  cat("\n--- C3  dag_lin strictly worse than dag_Yonly in the pipe_nl cells ---\n")
  c3 <- vapply(c2_cells, function(c0)
    abs(g(d, c0, "dag_lin", "bias_se")) - abs(g(d, c0, "dag_Yonly", "bias_se")), 0)
  for (c0 in c2_cells) cat(sprintf("  %-13s Yonly %+.2f  lin %+.2f   lin is %s by %.2f\n",
    c0, g(d, c0, "dag_Yonly", "bias_se"), g(d, c0, "dag_lin", "bias_se"),
    if (c3[[c0]] > 0) "WORSE" else "better", abs(c3[[c0]])))
  # The linear cells are the control: there dag_lin is correctly specified, so it
  # must NOT be worse. A track where dag_lin is worse everywhere is measuring
  # something other than misspecification.
  ctrl <- vapply(c("dag_pipe_dir", "dag_pipe_tot"), function(c0)
    abs(g(d, c0, "dag_lin", "bias_se")) - abs(g(d, c0, "dag_Yonly", "bias_se")), 0)
  for (c0 in names(ctrl)) cat(sprintf("  %-13s (LINEAR control) lin is %s by %.2f\n",
    c0, if (ctrl[[c0]] > 0) "WORSE" else "better", abs(ctrl[[c0]])))
  c3_met <- all(c3 > 0, na.rm = TRUE)
  c3_clean <- all(ctrl <= 0, na.rm = TRUE)
  cat(sprintf("  C3: %s%s\n",
              if (isTRUE(c3_met)) "HOLDS -- a linear child model must NOT be the default" else
                "does NOT hold -- dag_lin degrades gracefully and is a defensible default",
              if (c3_clean) "" else
                "  <- but dag_lin is also worse in the LINEAR controls, so this is not misspecification"))

  # ---- C4: the estimand wiring ---------------------------------------------
  cat("\n--- C4  the *_tot cells recover their own estimand, not b1 --------\n")
  c4 <- TRUE
  for (c0 in c("dag_pipe_tot", "dag_pnl_tot")) {
    est <- d$mean_est[d$scenario == c0 & d$procedure == "oracle"]
    if (!length(est)) { c4 <- FALSE; next }
    ok_t <- abs(est - tv[[c0]]) / tv[[c0]] < 0.05
    ok_b <- abs(est - 0.40) / 0.40 < 0.05
    cat(sprintf("  %-13s oracle %.4f  vs own estimand %.4f (%s)  vs b1 0.4000 (%s)\n",
                c0, est, tv[[c0]], if (ok_t) "matches" else "NO",
                if (ok_b) "also matches -- WIRING SUSPECT" else "correctly does not"))
    c4 <- c4 && ok_t && !ok_b
  }
  cat(sprintf("  C4: %s\n", if (isTRUE(c4)) "MET" else "NOT MET"))

  # ---- C5: direction (descriptive) -----------------------------------------
  cat("\n--- C5  direction of the incongenial draw's bias, per cell (DESCRIPTIVE) ---\n")
  for (c0 in CELLS) {
    b <- g(d, c0, "dag_noY", "bias_se"); m <- g(d, c0, "dag_noY", "mcse_bse")
    lab <- if (is.na(b)) "unmeasured" else if (abs(b) < 2 * m) "indistinguishable from zero"
           else if (b < 0) "toward null (conservative)" else "away from null (ANTI-conservative)"
    cat(sprintf("  %-13s %+.2f +/- %.2f  %s\n", c0, b, m, lab))
  }
  cat("  C5 is descriptive: no pass/fail, and must not be reported as one.\n")

  verdicts[[as.character(n)]] <- c(C1 = c1_met, C2 = c2_met, C3 = c3_met, C4 = c4)
}

cat("\n=====================================================================\n")
cat("  criteria by n\n")
cat("=====================================================================\n")
cat(sprintf("%8s %6s %6s %6s %6s\n", "n", "C1", "C2", "C3", "C4"))
for (n in names(verdicts)) {
  v <- verdicts[[n]]
  cat(sprintf("%8s %6s %6s %6s %6s\n", n,
              ifelse(v["C1"], "MET", "no"), ifelse(v["C2"], "MET", "no"),
              ifelse(v["C3"], "HOLDS", "no"), ifelse(v["C4"], "MET", "no")))
}
cat("\nC3 HOLDING is the outcome that constrains the deliverable: the algorithm\n")
cat("would be correct (C1) while its one unidentifiable input decides whether it\n")
cat("helps or hurts -- making it a sensitivity procedure, not a drop-in fix.\n")

utils::write.csv(all, file.path(.here, "results", sprintf("%s_analysis.csv", TAG)),
                 row.names = FALSE)
cat(sprintf("\nwrote results/%s_analysis.csv\n", TAG))
