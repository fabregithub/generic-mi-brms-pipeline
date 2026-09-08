#!/usr/bin/env Rscript
# =============================================================================
# V21 -- which property of a covariate draw has to be right?
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED. Criteria: the header of run_v21_zdraw_ladder.sh,
# committed before the run, and ../PLAN_pipeline_validation.md §8n.
#
# THE ANCHOR IS READ FIRST. z21_exact is the true conditional; V20 measured it at
# +-0.5%. If it has drifted here, the ladder has no zero point and nothing below
# is interpretable.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
setwd(.here)
NS <- c(800L, 3200L, 12800L); CELLS <- c("zr_fork", "zr_pipe"); TV <- 0.40
LADDER <- c("z21_exact", "z21_fit_proper", "z21_fit_improper", "z21_pmm",
            "z21_bart", "z21_bart_inner3")
fails <- character(0)
gate <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "MISS", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- c(fails, l)
}
S <- do.call(rbind, lapply(NS, function(n) {
  d <- read.csv(sprintf("results/v21_n%d_summary.csv", n), stringsAsFactors = FALSE)
  d$n <- n; d }))
g <- function(sc, pr, col) {
  v <- S[[col]][S$n %in% NS & S$scenario == sc & S$procedure == pr]
  v[order(S$n[S$scenario == sc & S$procedure == pr])]
}
rb <- function(sc, pr) g(sc, pr, "rel_bias") * 100

cat("=== ANCHOR: the true conditional (and the oracle) ===\n")
am <- 0; om <- 0
for (sc in CELLS) {
  ae <- rb(sc, "z21_exact"); oo <- rb(sc, "oracle")
  am <- max(am, max(abs(ae))); om <- max(om, max(abs(oo)))
  cat(sprintf("  %-9s exact %s | oracle %s\n", sc,
              paste(sprintf("%7.2f%%", ae), collapse = ""),
              paste(sprintf("%7.2f%%", oo), collapse = "")))
}
gate("oracle unbiased (<2%) at every n and cell", om < 2, sprintf("worst %.2f%%", om))
gate("z21_exact anchor within 1 pp of zero everywhere", am <= 1,
     sprintf("worst %.2f%%", am))
if (am > 1) cat("  !! No zero point. The ladder below is NOT interpretable.\n")

cat("\n=== THE LADDER: relative bias by arm and n ===\n")
for (sc in CELLS) {
  cat(sprintf("\n  --- %s ---\n  %-18s %8s %8s %8s %10s\n", sc, "arm",
              "n=800", "n=3200", "n=12800", "slope"))
  for (a in LADDER) {
    v <- rb(sc, a)
    sl <- if (min(abs(v)) > 0.8) sprintf("%+.3f", coef(lm(log(abs(v)) ~ log(NS)))[[2]])
          else "  (~0; n/a)"
    cat(sprintf("  %-18s %7.2f%% %7.2f%% %7.2f%% %10s\n",
                sub("^z21_", "", a), v[1], v[2], v[3], sl))
  }
}

cat("\n=== the registered gates ===\n")
fp <- vapply(CELLS, function(sc) rb(sc, "z21_fit_proper")[1], 0)
gate("fit_proper |bias| <= 2 pp at n = 800 (a correct ESTIMATED draw is safe)",
     all(abs(fp) <= 2), paste(sprintf("%s %+.2f%%", CELLS, fp), collapse = "  "))
bt <- vapply(CELLS, function(sc) rb(sc, "z21_bart")[1], 0)
gate("bart carries the reference bias (+3 to +7% at n = 800)",
     all(bt >= 3 & bt <= 7), paste(sprintf("%s %+.2f%%", CELLS, bt), collapse = "  "))
bsl <- vapply(CELLS, function(sc) coef(lm(log(abs(rb(sc, "z21_bart"))) ~ log(NS)))[[2]], 0)
gate("bart's bias decays near -1/3 (within 0.15)", all(abs(bsl + 1/3) <= 0.15),
     paste(sprintf("%s %+.3f", CELLS, bsl), collapse = "  "))

cat("\n--- properness: bias should NOT move, `b` and coverage should ---\n")
for (sc in CELLS) {
  d <- rb(sc, "z21_fit_improper") - rb(sc, "z21_fit_proper")
  bp <- g(sc, "z21_fit_proper", "b"); bi <- g(sc, "z21_fit_improper", "b")
  cp <- g(sc, "z21_fit_proper", "coverage"); ci <- g(sc, "z21_fit_improper", "coverage")
  cat(sprintf("  %-9s bias shift %s | b ratio (improper/proper) %s\n", sc,
              paste(sprintf("%+6.2f", d), collapse = ""),
              paste(sprintf("%6.2f", bi / bp), collapse = "")))
  cat(sprintf("  %-9s coverage proper %s | improper %s\n", "",
              paste(sprintf("%6.3f", cp), collapse = ""),
              paste(sprintf("%6.3f", ci), collapse = "")))
}
dmax <- max(abs(unlist(lapply(CELLS, function(sc)
  rb(sc, "z21_fit_improper") - rb(sc, "z21_fit_proper")))))
gate("dropping properness moves BIAS by <= 1.5 pp", dmax <= 1.5,
     sprintf("worst %+.2f pp", dmax))
bratio <- min(unlist(lapply(CELLS, function(sc)
  g(sc, "z21_fit_improper", "b") / g(sc, "z21_fit_proper", "b"))))
gate("dropping properness SHRINKS b (V4's improper-MI signature)", bratio < 1,
     sprintf("smallest b ratio %.2f", bratio))

cat("\n--- donor matching: does pmm's bias decay? ---\n")
for (sc in CELLS) {
  v <- rb(sc, "z21_pmm")
  f <- lm(log(abs(v)) ~ log(NS)); ci <- confint(f)["log(NS)", ]
  cat(sprintf("  %-9s %s  slope %+.3f [%+.3f, %+.3f]\n", sc,
              paste(sprintf("%7.2f%%", v), collapse = ""),
              coef(f)[[2]], ci[1], ci[2]))
}

cat("\n--- the inner-FCS loop (bart_inner3 vs bart) ---\n")
for (sc in CELLS)
  cat(sprintf("  %-9s %s pp\n", sc,
              paste(sprintf("%+7.2f", rb(sc, "z21_bart_inner3") - rb(sc, "z21_bart")),
                    collapse = "")))

cat("\n=== what the ladder attributes ===\n")
for (sc in CELLS) {
  ref <- rb(sc, "z21_bart")[1]; anc <- rb(sc, "z21_exact")[1]
  span <- ref - anc
  cat(sprintf("\n  %-9s reference (bart) %+.2f%%, anchor %+.2f%%  ->  span %.2f pp\n",
              sc, ref, anc, span))
  for (a in c("z21_fit_proper", "z21_fit_improper", "z21_pmm")) {
    cat(sprintf("    %-18s %+6.2f%%  = %5.1f%% of the span\n", sub("^z21_", "", a),
                rb(sc, a)[1], (rb(sc, a)[1] - anc) / span * 100))
  }
  cat(sprintf("    %-18s %+6.2f%%  = %5.1f%% of the span   <- everything the\n",
              "bart", ref, 100))
  cat("                                                        ladder must explain\n")
}

cat("\n=== verdict ===\n")
if (!length(fails)) cat("  Every registered gate passed.\n") else {
  cat(sprintf("  %d gate(s) missed:\n", length(fails)))
  for (f in fails) cat("    - ", f, "\n", sep = "")
}
