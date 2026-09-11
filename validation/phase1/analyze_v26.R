#!/usr/bin/env Rscript
# =============================================================================
# Track V26 -- analysis against the criteria registered in run_v26_family.sh
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED. Nothing pre-written, and no gate reported as met
# unless the number meeting it is printed beside it. Two tracks have resolved a
# gate for the wrong reason (V22's `bart >= +3%` missed at +2.84%; V23's "CI
# contains 2" passed by being uninformative), which is why each gate here
# carries an INFORMATIVENESS check that can void a pass.
#
# THE CURRENCY MATTERS HERE MORE THAN USUAL. A logistic coefficient is a
# log-odds-ratio; a Gaussian one is a mean difference. **Relative bias is not
# comparable across the two.** Everything below is read in bias/SE and coverage.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.here, "R", "dgp.R")); source(file.path(.here, "R", "robustness.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

TAG <- Sys.getenv("TAG", "v26")
NS  <- as.integer(strsplit(Sys.getenv("NLEVELS", "800,3200,12800"), ",")[[1]])
# Each binomial cell and the Gaussian cell it is the twin of.
PAIRS <- list(c(bin = "bin_base", gau = "mcar_z40"),
              c(bin = "bin_fork", gau = "zr_fork"),
              c(bin = "bin_pipe", gau = "zr_pipe"))
ARMS <- c("complete_case", "leftcens_prestep", "pipeline_bartMI",
          "pipeline_noYx", "pipeline_noYz", "pipeline_properZ",
          "pipeline_properZ_ds")

scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
truth_of <- function(nm) {
  s <- scs[[nm]]
  make_truth(p = 3L, erf_form = s$erf_form, y_form = s$y_form %||% "linear",
             z_role = s$z_role %||% "precision", nl_a = s$nl_a, nl_c = s$nl_c,
             target = s$target %||% "direct", delta_xz = s$delta_xz,
             delta_zx = s$delta_zx,
             y_family = s$y_family %||% "gaussian")$estimand_true
}
load_stage <- function(n) {
  f <- file.path(.here, "results", sprintf("%s_n%d_latest.rds", TAG, n))
  if (!file.exists(f)) return(NULL)
  d <- as.data.frame(readRDS(f)$summary); d$n <- n
  d$truth   <- vapply(as.character(d$scenario), truth_of, 0)
  d$bias_se <- (d$mean_est - d$truth) / d$claimed_se
  d$mcse    <- (d$emp_se / sqrt(d$n_ok)) / d$claimed_se
  d
}
all <- do.call(rbind, Filter(Negate(is.null), lapply(NS, load_stage)))
if (is.null(all) || !nrow(all)) stop("no V26 results for TAG=", TAG)
g <- function(d, cl, ar, f) { v <- d[[f]][d$scenario == cl & d$procedure == ar]
                              if (length(v) != 1L) NA_real_ else v }

verdicts <- list()
for (n in sort(unique(all$n))) {
  d <- all[all$n == n, , drop = FALSE]
  cat(sprintf("\n=====================================================================\n"))
  cat(sprintf("  n = %d   (%d reps)\n", n, max(d$n_ok, na.rm = TRUE)))
  cat(sprintf("=====================================================================\n"))
  for (p in PAIRS) {
    cat(sprintf("\n--- %s (binomial)  vs  %s (gaussian) ---\n", p[["bin"]], p[["gau"]]))
    cat(sprintf("%-22s %18s %18s %8s\n", "arm",
                "binom bias/SE(cov)", "gauss bias/SE(cov)", "ratio"))
    for (ar in c("oracle", ARMS)) {
      bb <- g(d, p[["bin"]], ar, "bias_se"); bg <- g(d, p[["gau"]], ar, "bias_se")
      if (is.na(bb) && is.na(bg)) next
      cat(sprintf("%-22s %+11.2f (%.2f) %+11.2f (%.2f) %8s\n", ar,
                  bb, g(d, p[["bin"]], ar, "coverage"),
                  bg, g(d, p[["gau"]], ar, "coverage"),
                  if (is.na(bb)||is.na(bg)||abs(bg) < 1e-9) "--"
                  else sprintf("%.2fx", abs(bb)/abs(bg))))
    }
  }

  # ---- C1: does the arm ORDERING survive the family change? ----------------
  cat("\n--- C1  arm ordering by |bias/SE| survives ------------------------\n")
  c1 <- logical(0); c1_inf <- logical(0)
  for (p in PAIRS) {
    bb <- vapply(ARMS, function(a) abs(g(d, p[["bin"]], a, "bias_se")), 0)
    bg <- vapply(ARMS, function(a) abs(g(d, p[["gau"]], a, "bias_se")), 0)
    ok <- is.finite(bb) & is.finite(bg)
    # informative only if the GAUSSIAN cell actually separates the arms
    inf <- sum(ok) >= 3L && (max(bg[ok]) - min(bg[ok])) > 0.3
    rho <- if (sum(ok) >= 3L) stats::cor(rank(bb[ok]), rank(bg[ok])) else NA_real_
    same <- identical(order(bb[ok]), order(bg[ok]))
    cat(sprintf("  %-10s rank corr %+.2f | identical order: %-5s | gaussian spread %.2f%s\n",
                p[["bin"]], rho, same, max(bg[ok]) - min(bg[ok]),
                if (!inf) "  <- VOID: no ordering to preserve" else ""))
    c1 <- c(c1, isTRUE(rho >= 0.8)); c1_inf <- c(c1_inf, inf)
  }
  c1_met <- all(c1, na.rm = TRUE) && all(c1_inf)
  cat(sprintf("  C1: %s\n", if (isTRUE(c1_met)) "MET" else
    if (all(c1, na.rm = TRUE)) "NOT MET -- holds only on uninformative cells" else "NOT MET"))

  # ---- C2: does the SIZE survive, within a factor of 2? --------------------
  cat("\n--- C2  each arm's |bias/SE| within 2x of its Gaussian twin -------\n")
  worst <- 0; worst_lab <- ""
  for (p in PAIRS) for (ar in ARMS) {
    bb <- abs(g(d, p[["bin"]], ar, "bias_se")); bg <- abs(g(d, p[["gau"]], ar, "bias_se"))
    if (!is.finite(bb) || !is.finite(bg) || bg < 0.1) next   # ratio meaningless near 0
    r <- max(bb/bg, bg/bb)
    if (r > worst) { worst <- r; worst_lab <- sprintf("%s / %s", p[["bin"]], ar) }
  }
  cat(sprintf("  worst ratio %.2fx  (%s)\n", worst, worst_lab))
  c2_met <- worst < 2
  cat(sprintf("  C2: %s\n", if (c2_met) "MET" else "NOT MET"))

  # ---- C3: the root claim -- noYx still biased toward the null -------------
  cat("\n--- C3  pipeline_noYx still biased TOWARD THE NULL under binomial -\n")
  c3 <- logical(0)
  for (p in PAIRS) {
    bb <- g(d, p[["bin"]], "pipeline_noYx", "bias_se")
    mc <- g(d, p[["bin"]], "pipeline_noYx", "mcse")
    lab <- if (is.na(bb)) "unmeasured" else if (abs(bb) < 2*mc) "indistinguishable from zero"
           else if (bb < 0) "toward null" else "AWAY from null"
    cat(sprintf("  %-10s %+.2f +/- %.2f  %s\n", p[["bin"]], bb, mc, lab))
    c3 <- c(c3, isTRUE(bb < 0 && abs(bb) >= 2*mc))
  }
  c3_met <- all(c3, na.rm = TRUE)
  cat(sprintf("  C3: %s\n", if (c3_met) "MET" else "NOT MET"))

  # ---- C4: the scale defect by family (DESCRIPTIVE) ------------------------
  cat("\n--- C4  the pre-v1.6.0 scale defect, by family (DESCRIPTIVE) ------\n")
  for (p in PAIRS) for (k in c("bin","gau")) {
    a <- g(d, p[[k]], "pipeline_properZ_ds", "bias_se")
    b <- g(d, p[[k]], "pipeline_properZ", "bias_se")
    if (is.na(a) || is.na(b)) next
    cat(sprintf("  %-10s (%s) defect = %+.2f bias/SE\n", p[[k]],
                if (k == "bin") "binom" else "gauss", a - b))
  }
  cat("  C4 is descriptive: no pass/fail, and must not be reported as one.\n")
  verdicts[[as.character(n)]] <- c(C1 = c1_met, C2 = c2_met, C3 = c3_met)
}

cat("\n=====================================================================\n")
cat(sprintf("%8s %6s %6s %6s\n", "n", "C1", "C2", "C3"))
for (n in names(verdicts)) { v <- verdicts[[n]]
  cat(sprintf("%8s %6s %6s %6s\n", n, ifelse(v["C1"],"MET","no"),
              ifelse(v["C2"],"MET","no"), ifelse(v["C3"],"MET","no"))) }
cat("\nC1 FAILING is the consequential outcome: guidance from 25 Gaussian tracks\n")
cat("could not then be quoted for logistic analyses without re-derivation.\n")
cat("\nNOTE: this track does NOT exercise Step 6's pooling machinery -- it pools\n")
cat("with Rubin's rules in the harness. That gap stays open either way.\n")
utils::write.csv(all, file.path(.here, "results", sprintf("%s_analysis.csv", TAG)),
                 row.names = FALSE)
