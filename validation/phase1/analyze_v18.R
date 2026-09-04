#!/usr/bin/env Rscript
# =============================================================================
# V18 -- does the shipped default's bias survive larger n?
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED. Criteria: ../PLAN_pipeline_validation.md §8k.
# CONTROLS ARE READ FIRST, on purpose: if the oracle drifts, the estimand or the
# DGP is n-dependent and nothing else in the run is readable, so there is no
# point looking at the headline before that is settled.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
setwd(.here)
NS <- c(800L, 1600L, 3200L, 6400L, 12800L)
BIASED <- c("zr_fork", "zr_pipe", "zrmar_fork", "zrmar_pipe")
CTRL   <- c("zr_collider", "mcar_z40")
fails <- character(0)
gate <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "MISS", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- c(fails, l)
}

S <- do.call(rbind, lapply(NS, function(n) {
  d <- read.csv(sprintf("results/v18_n%d_summary.csv", n), stringsAsFactors = FALSE)
  d$n <- n; d
}))
get <- function(sc, proc, col) vapply(NS, function(n)
  S[[col]][S$n == n & S$scenario == sc & S$procedure == proc], 0)

cat("=== CONTROL 1: the oracle, at every n ===\n")
cat(sprintf("  %-14s %s\n", "cell", paste(sprintf("%9d", NS), collapse = "")))
omax <- 0
for (sc in c(BIASED, CTRL)) {
  rb <- get(sc, "oracle", "rel_bias") * 100
  omax <- max(omax, max(abs(rb)))
  cat(sprintf("  %-14s %s\n", sc, paste(sprintf("%8.2f%%", rb), collapse = "")))
}
gate("oracle stays within 2% at every n and cell", omax < 2,
     sprintf("worst |rel bias| = %.2f%%", omax))

cat("\n=== HEADLINE: relative bias of the SHIPPED default (pipeline_bartMI) ===\n")
cat(sprintf("  %-14s %s %10s\n", "cell",
            paste(sprintf("%9d", NS), collapse = ""), "slope"))
slopes <- list()
for (sc in c(BIASED, CTRL)) {
  rb <- get(sc, "pipeline_bartMI", "rel_bias") * 100
  fit <- lm(log(abs(rb)) ~ log(NS)); ci <- confint(fit)["log(NS)", ]
  slopes[[sc]] <- list(b = coef(fit)[["log(NS)"]], ci = ci, rb = rb)
  cat(sprintf("  %-14s %s %7.3f\n", sc,
              paste(sprintf("%8.2f%%", rb), collapse = ""), slopes[[sc]]$b))
}

cat("\n--- PRIMARY: the slope test (constant = 0, finite-sample = -0.5) ---\n")
for (sc in BIASED) {
  x <- slopes[[sc]]
  excl0 <- x$ci[1] > 0 || x$ci[2] < 0
  exclH <- x$ci[1] > -0.5 || x$ci[2] < -0.5
  verdict <- if (!excl0 && exclH) "CONSTANT bias" else
             if (excl0 && !exclH) "FINITE-SAMPLE bias" else
             if (excl0 && exclH)  "NEITHER -- intermediate exponent" else
                                  "UNDETERMINED -- CI covers both"
  cat(sprintf("  %-14s slope %+.3f  95%% CI [%+.3f, %+.3f]  ->  %s\n",
              sc, x$b, x$ci[1], x$ci[2], verdict))
}
pooled <- do.call(rbind, lapply(BIASED, function(sc)
  data.frame(y = log(abs(slopes[[sc]]$rb)), x = log(NS), cell = sc)))
pf <- lm(y ~ x + cell, data = pooled); pci <- confint(pf)["x", ]
cat(sprintf("\n  POOLED across the four biased cells: slope %+.3f  95%% CI [%+.3f, %+.3f]\n",
            coef(pf)[["x"]], pci[1], pci[2]))
gate("pooled slope CI excludes 0 (constant-bias rejected)",
     pci[1] > 0 || pci[2] < 0, sprintf("[%+.3f, %+.3f]", pci[1], pci[2]))
gate("pooled slope CI excludes -0.5 (finite-sample rejected)",
     pci[1] > -0.5 || pci[2] < -0.5, sprintf("[%+.3f, %+.3f]", pci[1], pci[2]))

cat("\n--- per-cell bar at n = 12800: within 2 pp of the constant-bias value ---\n")
PRED <- c(zr_fork = 4.96, zr_pipe = 5.39, zrmar_fork = 8.52, zrmar_pipe = 10.40)
for (sc in BIASED) {
  got <- slopes[[sc]]$rb[length(NS)]
  gate(sprintf("%-14s n=12800 within 2 pp of %+.2f%%", sc, PRED[[sc]]),
       abs(got - PRED[[sc]]) <= 2,
       sprintf("%+.2f%%  (finite-sample would give %+.2f%%)", got, PRED[[sc]] / 4))
}

cat("\n=== CONTROL 2: the two cells that should not develop a large bias ===\n")
for (sc in CTRL) {
  x <- slopes[[sc]]
  gate(sprintf("%-14s stays below 3%% at every n", sc), max(abs(x$rb)) < 3,
       sprintf("worst %+.2f%%; slope %+.3f [%+.3f, %+.3f]",
               x$rb[which.max(abs(x$rb))], x$b, x$ci[1], x$ci[2]))
}

cat("\n=== CONSEQUENCE: coverage, measured against the projection ===\n")
PROJ <- list(zr_fork = c(.939,NA,.892,NA,.693), zr_pipe = c(.944,NA,.903,NA,.731),
             zrmar_fork = c(.927,NA,.806,NA,.358), zrmar_pipe = c(.936,NA,.807,NA,.321))
cat(sprintf("  %-14s %s\n", "cell", paste(sprintf("%9d", NS), collapse = "")))
for (sc in c(BIASED, CTRL)) {
  cv <- get(sc, "pipeline_bartMI", "coverage")
  cat(sprintf("  %-14s %s\n", sc, paste(sprintf("%9.3f", cv), collapse = "")))
  if (sc %in% names(PROJ))
    cat(sprintf("  %-14s %s   <- projected\n", "",
                paste(sprintf("%9s", ifelse(is.na(PROJ[[sc]]), "-",
                                            sprintf("%.3f", PROJ[[sc]]))), collapse = "")))
}
cat("\n  width/SE (the coverage projections assume this holds constant in n):\n")
wmin <- 9; wmax <- 0
for (sc in c(BIASED, CTRL)) {
  w <- get(sc, "pipeline_bartMI", "width_se_ratio")
  wmin <- min(wmin, min(w)); wmax <- max(wmax, max(w))
  cat(sprintf("  %-14s %s\n", sc, paste(sprintf("%9.3f", w), collapse = "")))
}
gate("width/SE stays in 0.9-1.2 (coverage projections remain valid)",
     wmin >= 0.9 && wmax <= 1.2, sprintf("range %.3f - %.3f", wmin, wmax))

cat("\n=== the corrected picture, and how far it extrapolates ===\n")
# The coverage model was not wrong -- its INPUT was. Fed the measured bias and
# width/SE it reproduces observed coverage; fed the assumed constant bias it did
# not. Re-run it here as a check, then extrapolate on the FITTED exponent.
mx <- 0
for (sc in BIASED) {
  q <- subset(S, procedure == "pipeline_bartMI" & scenario == sc)
  q <- q[order(q$n), ]
  r <- abs(q$rel_bias * 0.40) / q$emp_se
  mo <- stats::pnorm(1.96 * q$width_se_ratio - r) - stats::pnorm(-1.96 * q$width_se_ratio - r)
  mx <- max(mx, max(abs(mo - q$coverage)))
}
gate("coverage model reproduces observed coverage from MEASURED inputs (<0.03)",
     mx < 0.03, sprintf("worst |model - observed| = %.3f", mx))

cat(sprintf("\n  bias decays as n^%.2f while the SE decays as n^-0.50, so bias/SE grows\n",
            coef(pf)[["x"]]))
cat(sprintf("  as n^%.2f -- real, but %.1fx slower than the n^0.50 the V17 warning assumed.\n",
            coef(pf)[["x"]] + 0.5, 0.5 / (coef(pf)[["x"]] + 0.5)))
# "n for .90" is NA where coverage is already below 0.90 at n = 12800 -- the
# threshold is behind us, not ahead.
cat("  (NA in the .90 column = already below 0.90 by n = 12800.)\n")
cat("\n  EXTRAPOLATION from the fitted exponent (not measured -- the run stops at 12800):\n")
cat(sprintf("  %-14s %12s %12s %12s\n", "cell", "cov @ 1e5", "n for .90", "n for .80"))
for (sc in BIASED) {
  q <- subset(S, procedure == "pipeline_bartMI" & scenario == sc); q <- q[order(q$n), ]
  r0 <- abs(q$rel_bias[5] * 0.40) / q$emp_se[5]; w <- q$width_se_ratio[5]
  g  <- coef(pf)[["x"]] + 0.5
  cvf <- function(n) { r <- r0 * (n / 12800)^g
                       stats::pnorm(1.96 * w - r) - stats::pnorm(-1.96 * w - r) }
  nfor <- function(target) {
    f <- function(ln) cvf(exp(ln)) - target
    if (f(log(1.28e4)) < 0) return(NA_real_)
    tryCatch(exp(stats::uniroot(f, c(log(1.28e4), log(1e12)))$root), error = function(e) NA_real_)
  }
  cat(sprintf("  %-14s %12.3f %12s %12s\n", sc, cvf(1e5),
              formatC(nfor(0.90), format = "d", big.mark = ","),
              formatC(nfor(0.80), format = "d", big.mark = ",")))
}
cat("\n  Coverage does eventually fail -- bias/SE is unbounded -- but it takes a\n")
cat("  sample far beyond this pipeline's realistic range to get there.\n")

cat("\n=== verdict ===\n")
if (!length(fails)) cat("  Every registered gate passed.\n") else {
  cat(sprintf("  %d gate(s) missed:\n", length(fails)))
  for (f in fails) cat("    - ", f, "\n", sep = "")
}
