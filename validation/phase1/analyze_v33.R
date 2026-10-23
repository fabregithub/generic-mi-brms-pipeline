# =============================================================================
# V33 -- can substantive-model-compatible imputation be used when the surface
#        form is DECLARED rather than known?
# -----------------------------------------------------------------------------
# Gates registered in run_v33_mixture.sh BEFORE the run. Transferred from the
# covariate-arrow analogue already measured (SUMMARY.md), where the same
# question was settled for Z1 = g(logX1):
#     correct basis, coefficients fitted ... -0.51%
#     correct form, wrong constant ......... +3.22%
#     omit the factor entirely ............. +8.47%
#     correct family, ONE TERM DROPPED ..... +20.71%   <- worse than omitting
#     a straight line ...................... +25.43%
#
#   P1 RECOVERY  smc_dec_true beats pipeline_bkmr and sub_lod2_bkmr on
#                curvature bias, and lands near oracle_bkmr.
#   P2 INCOMPLETE BASIS IS WORSE THAN NONE  smc_dec_noint and smc_dec_noquad
#                are WORSE than smc_dec_lin. THE DECIDING GATE.
#   P3 OVER-DECLARING IS CHEAP  smc_dec_over is close to smc_dec_true.
#
#   FLOOR: the curvature estimand must be measurably non-zero in oracle_bkmr,
#   or every comparison below is a ratio of near-zero arms.
#
#   NOT A GATE: which arm "wins" overall. The contribution is the mechanism and
#   its limit, not a leaderboard.
#
#   Reads the checkpoint if the final summary is absent, so a partial run is
#   still analysable (V3 checkpoints every 44-task chunk).
#
#   ⚠️ EVERY GATE CARRIES MONTE-CARLO ERROR, and reports THREE states, not two:
#       MET        resolved (>2 sigma) and in the predicted direction
#       NOT MET    resolved and in the WRONG direction
#       UNDECIDED  not resolved -- the data cannot tell yet
#   The first draft of this script compared point estimates with no error and
#   printed "NOT MET" off a 1.1-sigma difference on the deciding gate. A gate
#   that cannot distinguish "refuted" from "underpowered" is the same failure
#   this programme has now hit eight times, just moved from design to analysis.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
TAG <- Sys.getenv("TAG", "v33")
EST <- Sys.getenv("ESTIMAND", "curv_X1")     # the surface feature the fix is about
# V34 runs TWO scenarios in one invocation. Filtering on estimand alone would
# silently average nd20 and nd40_all together, which is the whole comparison.
SCN <- Sys.getenv("SCENARIO", "")

res <- file.path(.here, "results", sprintf("%s_latest.rds", TAG))
ck  <- file.path(.here, "results", sprintf("%s_checkpoint.rds", TAG))
if (file.exists(res)) {
  s <- as.data.frame(readRDS(res)$summary); src <- "final"
} else if (file.exists(ck)) {
  raw <- readRDS(ck)
  if (!exists("summarise_v3")) source(file.path(.here, "run_v3_bkmr.R"), local = TRUE)
  s <- as.data.frame(summarise_v3(raw)); src <- "CHECKPOINT (run incomplete)"
} else stop("no V33 results at ", res)

if (nzchar(SCN)) s <- s[s$scenario == SCN, ]
if (length(unique(s$scenario)) > 1L)
  stop("results hold ", length(unique(s$scenario)), " scenarios (",
       paste(unique(s$scenario), collapse = ", "),
       "). Set SCENARIO= to pick one -- averaging them would destroy the comparison.")
d <- s[s$estimand == EST, ]
if (!nrow(d)) stop("estimand '", EST, "' not present; have: ",
                   paste(unique(s$estimand), collapse = ", "))
g <- function(arm, f) { v <- d[[f]][d$procedure == arm]; if (length(v) == 1L) v else NA_real_ }
# ⚠️ REFERENCE IS THE ORACLE, NOT THE TRUTH.
# BKMR at these settings is itself biased on curvature: the oracle arm -- complete,
# uncensored data -- reads +15.1% against truth at 200 reps, 6 sigma from zero. That is
# GPP/iteration/knot approximation error, and EVERY arm inherits it. Scoring arms against
# truth therefore measures BKMR's bias plus the imputation's, and the first dominates.
# The IMPUTATION-induced bias is (arm - oracle)/oracle. The first draft of this script used
# truth and reached the opposite verdict on P2.
.orc    <- g("oracle_bkmr", "mean_est")
.orc_se <- g("oracle_bkmr", "mcse_bias")
d$vs_oracle <- (d$mean_est - .orc) / .orc
d$vs_or_se  <- sqrt(d$mcse_bias^2 + .orc_se^2) / abs(.orc)
ab   <- function(a) abs(g(a, "vs_oracle"))
abse <- function(a) g(a, "vs_or_se")

#' Compare |rel_bias| of two arms with MC error.
#' @param better the arm the prediction says should be SMALLER
#' @return "MET" / "NOT MET" / "UNDECIDED", printing the comparison
cmp <- function(better, worse, label = "") {
  a <- ab(better); b <- ab(worse)
  se <- sqrt(abse(better)^2 + abse(worse)^2)
  if (!is.finite(a) || !is.finite(b) || !is.finite(se) || se <= 0) {
    cat(sprintf("  %-15s vs %-15s  VOID (missing)\n", better, worse)); return("VOID") }
  z <- (b - a) / se                      # positive = predicted direction
  v <- if (abs(z) <= 2) "UNDECIDED" else if (z > 0) "MET" else "NOT MET"
  cat(sprintf("  %-15s %.3f  vs %-15s %.3f   diff %+.3f +/- %.3f  %4.1f sigma  %s%s\n",
              better, a, worse, b, a - b, se, z, v, if (nzchar(label)) paste0("  [", label, "]") else ""))
  v
}
roll <- function(v) if (any(v == "NOT MET")) "NOT MET" else if (any(v == "UNDECIDED" | v == "VOID")) "UNDECIDED" else "MET"

cat(sprintf("\n=== V33  estimand = %s  [%s]  reps = %d ===\n",
            EST, src, suppressWarnings(max(d$n_ok, na.rm = TRUE))))

keep <- c("procedure","n_ok","mean_est","rel_bias","vs_oracle","vs_or_se","coverage")
cat("\n--- all arms  (rel_bias = vs truth; vs_oracle = IMPUTATION-induced) ---\n")
print(format(d[order(abs(d$vs_oracle)), intersect(keep, names(d))], digits = 3), row.names = FALSE)
cat(sprintf("\n  oracle vs truth: %+.3f -- BKMR's OWN curvature bias, shared by every arm\n",
            g("oracle_bkmr", "rel_bias")))

# ---- FLOOR ------------------------------------------------------------------
tru <- unique(d$true)[1]
orc <- g("oracle_bkmr", "rel_bias")
cat(sprintf("\n--- FLOOR: is %s measurably non-zero? true = %.4f, oracle rel_bias = %+.3f ---\n",
            EST, tru, orc))
floor_ok <- is.finite(tru) && abs(tru) > 1e-6 && is.finite(orc) && abs(orc) < 0.25
cat(sprintf("  %s\n", if (floor_ok) "ok -- comparisons below are defined"
                      else "⚠️ VOID -- the oracle itself does not recover this estimand"))

cat("\n--- P1  does declaring the TRUE basis recover the surface? ---\n")
cat(sprintf("  all figures are IMPUTATION-induced (vs oracle); MC error ~ %.3f at %d reps\n",
            mean(d$vs_or_se, na.rm = TRUE),
            suppressWarnings(max(d$n_ok, na.rm = TRUE))))
p1 <- roll(c(cmp("smc_dec_true", "pipeline_bkmr", "vs the shipped engine"),
             cmp("smc_dec_true", "sub_lod2_bkmr", "vs LOD/sqrt(2)")))
cat(sprintf("  P1: %s\n", p1))

cat("\n--- P2  is an INCOMPLETE basis worse than declaring NOTHING? (deciding) ---\n")
cat("  prediction: smc_dec_lin should be SMALLER than both incomplete arms\n")
p2 <- roll(c(cmp("smc_dec_lin", "smc_dec_noint",  "interaction dropped"),
             cmp("smc_dec_lin", "smc_dec_noquad", "quadratic dropped")))
cat(sprintf("  P2: %s%s\n", p2,
   if (p2 == "MET") "  -- declared-form SMC is UNSAFE where the basis is uncertain"
   else if (p2 == "NOT MET") "  -- an incomplete basis still helps" else "  -- more reps needed"))

cat("\n--- P3  is OVER-declaring cheap? ---\n")
cat("  prediction: smc_dec_over is NOT meaningfully worse than smc_dec_true\n")
p3v <- cmp("smc_dec_true", "smc_dec_over", "spurious terms")
p3 <- if (p3v == "UNDECIDED") "MET (indistinguishable -- which is the prediction)" else
      if (p3v == "MET") "MET (over-declaring is in fact better)" else
      "NOT MET -- spurious terms cost real accuracy"
cat(sprintf("  P3: %s\n", p3))

cat("\n--- reading ---\n")
cat("P1 + P2 -> the fix works but only with a basis you can defend; for a PFAS\n")
cat("mixture that is rarely, so report the limit rather than the method.\n")
cat("P1 + !P2 + P3 -> 'declare generously, never sparsely' is shippable advice.\n")
cat("!P1 -> declaring the form is not enough; the surface must be learned.\n")
cat("\nOther estimands: ESTIMAND=int_X1X2 Rscript analyze_v33.R\n")
