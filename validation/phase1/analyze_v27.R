#!/usr/bin/env Rscript
# =============================================================================
# Track V27 -- analysis against the criteria registered in run_v27_zimputer.sh
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED, and every gate carries an informativeness check
# that can void a pass. Four gates in this programme have resolved for the
# wrong reason (V22, V23, and both of V26's), and the lesson from the last two
# is that a guard needs a FLOOR ON THE QUANTITY BEING COMPARED, not on a spread
# or a difference. C2 below applies that directly.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.here, "R", "dgp.R")); source(file.path(.here, "R", "robustness.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
TAG <- Sys.getenv("TAG", "v27")
NS  <- as.integer(strsplit(Sys.getenv("NLEVELS", "800,3200,12800"), ",")[[1]])

INERT  <- c("mcar_z40", "mar_z40", "nl_mcar_z40", "nl_mar_z40")   # C1's ground
ACTIVE <- c("zr_fork", "zr_pipe", "zr_collider", "zr_mixed",
            "zrmar_fork", "zrmar_pipe")                           # C2's ground
NLARROW <- "zr_pipenl"
BART <- "pipeline_bartMI"; MICE <- "pipeline_miceShipped"; HARN <- "pipeline_micePmm"

scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
truth_of <- function(nm) { s <- scs[[nm]]
  make_truth(p = 3L, erf_form = s$erf_form, y_form = s$y_form %||% "linear",
             z_role = s$z_role %||% "precision", nl_a = s$nl_a, nl_c = s$nl_c,
             target = s$target %||% "direct", delta_xz = s$delta_xz,
             delta_zx = s$delta_zx, y_family = s$y_family %||% "gaussian")$estimand_true }
load_stage <- function(n) {
  f <- file.path(.here, "results", sprintf("%s_n%d_latest.rds", TAG, n))
  if (!file.exists(f)) return(NULL)
  d <- as.data.frame(readRDS(f)$summary); d$n <- n
  d$truth <- vapply(as.character(d$scenario), truth_of, 0)
  d$bias_pct <- 100 * (d$mean_est - d$truth) / d$truth
  d$bias_se  <- (d$mean_est - d$truth) / d$claimed_se
  d$mcse     <- (d$emp_se / sqrt(d$n_ok)) / d$claimed_se
  d
}
all <- do.call(rbind, Filter(Negate(is.null), lapply(NS, load_stage)))
if (is.null(all) || !nrow(all)) stop("no V27 results for TAG=", TAG)
g <- function(d, cl, ar, f) { v <- d[[f]][d$scenario == cl & d$procedure == ar]
                              if (length(v) != 1L) NA_real_ else v }

verdicts <- list()
for (n in sort(unique(all$n))) {
  d <- all[all$n == n, , drop = FALSE]
  cells <- intersect(c(INERT, ACTIVE, NLARROW), unique(as.character(d$scenario)))
  cat(sprintf("\n=====================================================================\n"))
  cat(sprintf("  n = %d   (%d reps)\n", n, max(d$n_ok, na.rm = TRUE)))
  cat(sprintf("=====================================================================\n"))
  cat(sprintf("%-16s %-9s | %18s %18s | %8s\n", "cell", "family",
              "bart bias/SE(cov)", "mice bias/SE(cov)", "gain"))
  for (cl in cells) {
    fam <- if (cl %in% INERT) "inert" else if (cl %in% ACTIVE) "active" else "nl-arrow"
    bb <- g(d, cl, BART, "bias_se"); bm <- g(d, cl, MICE, "bias_se")
    cat(sprintf("%-16s %-9s | %+11.2f (%.2f) %+11.2f (%.2f) | %+8.2f\n", cl, fam,
                bb, g(d, cl, BART, "coverage"), bm, g(d, cl, MICE, "coverage"),
                abs(bb) - abs(bm)))
  }

  # ---- C1: no loss on V6's ground. THE GATE THAT DECIDES THE DEFAULT -------
  cat("\n--- C1  mice within 0.3 bias/SE of bart on inert + non-linear cells ---\n")
  c1 <- logical(0)
  for (cl in intersect(INERT, cells)) {
    bb <- abs(g(d, cl, BART, "bias_se")); bm <- abs(g(d, cl, MICE, "bias_se"))
    loss <- bm - bb
    cat(sprintf("  %-14s bart %.2f  mice %.2f  loss %+.2f %s\n", cl, bb, bm, loss,
                if (loss <= 0.3) "ok" else "EXCEEDS 0.3"))
    c1 <- c(c1, loss <= 0.3)
  }
  c1_met <- all(c1, na.rm = TRUE)
  cat(sprintf("  C1: %s%s\n", if (isTRUE(c1_met)) "MET" else "NOT MET",
              if (isTRUE(c1_met)) "" else "  -> the default must NOT switch"))

  # ---- C2: real gain where it matters -------------------------------------
  cat("\n--- C2  mice improves |bias/SE| by >= 1.0 in causally-active cells ----\n")
  c2 <- logical(0)
  for (cl in intersect(ACTIVE, cells)) {
    bb <- abs(g(d, cl, BART, "bias_se")); bm <- abs(g(d, cl, MICE, "bias_se"))
    # VOID where bart has no bias to remove -- a floor on the quantity compared.
    if (!is.finite(bb) || bb < 0.5) {
      cat(sprintf("  %-14s VOID: bart's own |bias/SE| is %.2f, nothing to improve\n", cl, bb))
      next
    }
    cat(sprintf("  %-14s bart %.2f -> mice %.2f  gain %+.2f %s\n", cl, bb, bm, bb - bm,
                if (bb - bm >= 1.0) "" else "< 1.0"))
    c2 <- c(c2, (bb - bm) >= 1.0)
  }
  c2_met <- length(c2) > 0 && all(c2, na.rm = TRUE)
  cat(sprintf("  C2: %s%s\n", if (isTRUE(c2_met)) "MET" else "NOT MET",
              if (length(c2) == 0) "  (no informative cell -- VOID)" else ""))

  # ---- C3: coverage not traded away ---------------------------------------
  cat("\n--- C3  mice coverage >= bart coverage - 0.02, every cell -----------\n")
  worst <- Inf; worst_cl <- ""
  for (cl in cells) {
    dd <- g(d, cl, MICE, "coverage") - g(d, cl, BART, "coverage")
    if (is.finite(dd) && dd < worst) { worst <- dd; worst_cl <- cl }
  }
  cat(sprintf("  worst coverage change %+.3f (%s)\n", worst, worst_cl))
  c3_met <- worst >= -0.02
  cat(sprintf("  C3: %s\n", if (c3_met) "MET" else "NOT MET"))

  # ---- C4: the port is faithful -------------------------------------------
  cat("\n--- C4  shipped option == harness override, within MC error ---------\n")
  worst4 <- 0; w4 <- ""
  for (cl in cells) {
    a <- g(d, cl, MICE, "bias_se"); b <- g(d, cl, HARN, "bias_se")
    mc <- g(d, cl, MICE, "mcse")
    if (!is.finite(a) || !is.finite(b)) next
    z <- abs(a - b) / max(mc, 1e-9)
    if (z > worst4) { worst4 <- z; w4 <- cl }
  }
  cat(sprintf("  largest discrepancy %.2f MC errors (%s)\n", worst4, w4))
  c4_met <- worst4 < 3
  cat(sprintf("  C4: %s%s\n", if (c4_met) "MET" else "NOT MET",
              if (c4_met) "" else "  -> the PORT differs from the method; fix before shipping"))

  # ---- C5: the non-linear ARROW cell, descriptive -------------------------
  if (NLARROW %in% cells) {
    cat("\n--- C5  the non-linear arrow cell (V22's ground), DESCRIPTIVE -------\n")
    cat(sprintf("  %-14s bart %+.2f  mice %+.2f\n", NLARROW,
                g(d, NLARROW, BART, "bias_se"), g(d, NLARROW, MICE, "bias_se")))
    cat("  C5 is descriptive: no pass/fail.\n")
  }
  verdicts[[as.character(n)]] <- c(C1 = c1_met, C2 = c2_met, C3 = c3_met, C4 = c4_met)
}

cat("\n=====================================================================\n")
cat(sprintf("%8s %6s %6s %6s %6s\n", "n", "C1", "C2", "C3", "C4"))
for (n in names(verdicts)) { v <- verdicts[[n]]
  cat(sprintf("%8s %6s %6s %6s %6s\n", n, ifelse(v["C1"],"MET","no"),
              ifelse(v["C2"],"MET","no"), ifelse(v["C3"],"MET","no"),
              ifelse(v["C4"],"MET","no"))) }
cat("\nC1 AND C2 both MET at every n => recommend switching the default.\n")
cat("C1 failing anywhere  => do NOT switch; ship mice as an option only.\n")
utils::write.csv(all, file.path(.here, "results", sprintf("%s_analysis.csv", TAG)),
                 row.names = FALSE)
