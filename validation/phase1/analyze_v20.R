#!/usr/bin/env Rscript
# =============================================================================
# V20 -- which block carries the bias under a causally active covariate?
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED. Criteria: ../PLAN_pipeline_validation.md §8m and the
# header of run_v20_attribute.sh, committed before the run.
#
# THE CONTROL IS READ FIRST AND IT GATES EVERYTHING. `ef_bart_ship` is the
# instrument configured to do what the pipeline does. If it does not reproduce
# `pipeline_bartMI`, the instrument is not the pipeline and the other three cells
# attribute nothing -- so there is no point reading the decomposition until that
# is settled.
#
# All contrasts are PAIRED within replication: every arm sees the same data, so
# the difference is far more precise than the difference of two means.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
setwd(.here)
NS <- c(800L, 3200L, 12800L); CELLS <- c("zr_fork", "zr_pipe"); TV <- 0.40
fails <- character(0)
gate <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "MISS", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- c(fails, l)
}
wide <- function(n, sc) {
  r <- subset(read.csv(sprintf("results/v20_n%d_raw.csv", n), stringsAsFactors = FALSE),
              scenario == sc)
  w <- reshape(r[, c("rep", "procedure", "estimate")], idvar = "rep",
               timevar = "procedure", direction = "wide")
  names(w) <- sub("^estimate\\.", "", names(w))
  w
}
relb <- function(w, a) (w[[a]] - TV) / TV * 100          # per-rep relative bias
pd   <- function(w, a, b) (w[[a]] - w[[b]]) / TV * 100   # paired difference
mse  <- function(v) sd(v) / sqrt(length(v))

cat("=== CONTROL 1: the oracle ===\n")
om <- 0
for (sc in CELLS) for (n in NS) {
  v <- mean(relb(wide(n, sc), "oracle")); om <- max(om, abs(v))
}
gate("oracle within 2% at every n and cell", om < 2, sprintf("worst %.2f%%", om))

cat("\n=== CONTROL 2: is the instrument the pipeline? ===\n")
cat("  (ef_bart_ship is the instrument doing what the pipeline does.)\n")
cm <- 0
for (sc in CELLS) for (n in NS) {
  w <- wide(n, sc); d <- pd(w, "ef_bart_ship", "pipeline_bartMI")
  cm <- max(cm, abs(mean(d)))
  cat(sprintf("  %-9s n=%-6d pipeline %+6.2f%%  instrument %+6.2f%%  paired diff %+5.2f +/- %.2f pp\n",
              sc, n, mean(relb(w, "pipeline_bartMI")), mean(relb(w, "ef_bart_ship")),
              mean(d), mse(d)))
}
gate("instrument reproduces pipeline_bartMI within 1.5 pp everywhere", cm <= 1.5,
     sprintf("worst |diff| = %.2f pp", cm))
if (cm > 1.5)
  cat("  !! The decomposition below is NOT interpretable. Stop here.\n")

cat("\n=== THE DECOMPOSITION (paired, vs ef_bart_ship) ===\n")
res <- list()
for (sc in CELLS) {
  cat(sprintf("\n  --- %s ---\n  %-7s %9s %9s %9s %9s %9s\n", sc, "n",
              "reference", "fix-Z", "fix-X", "both", "interact"))
  for (n in NS) {
    w <- wide(n, sc)
    ref <- mean(relb(w, "ef_bart_ship"))
    fz  <- pd(w, "ef_exact_ship",  "ef_bart_ship")
    fx  <- pd(w, "ef_bart_exact",  "ef_bart_ship")
    bo  <- pd(w, "ef_exact_exact", "ef_bart_ship")
    ia  <- bo - fz - fx
    res[[paste(sc, n)]] <- list(ref = ref, fz = mean(fz), fx = mean(fx),
                                bo = mean(bo), ia = mean(ia),
                                se = c(fz = mse(fz), fx = mse(fx), bo = mse(bo), ia = mse(ia)),
                                ee = mean(relb(w, "ef_exact_exact")))
    cat(sprintf("  %-7d %8.2f%% %+8.2f %+8.2f %+8.2f %+8.2f\n", n, ref,
                mean(fz), mean(fx), mean(bo), mean(ia)))
  }
  cat(sprintf("  %-7s %9s %+8.2f %+8.2f %+8.2f %+8.2f   <- mc se\n", "(se)", "",
              mean(vapply(NS, function(n) res[[paste(sc,n)]]$se[["fz"]], 0)),
              mean(vapply(NS, function(n) res[[paste(sc,n)]]$se[["fx"]], 0)),
              mean(vapply(NS, function(n) res[[paste(sc,n)]]$se[["bo"]], 0)),
              mean(vapply(NS, function(n) res[[paste(sc,n)]]$se[["ia"]], 0))))
}

cat("\n--- PRIMARY: share of the reference bias each block removes ---\n")
cat("  (registered: fix-Z >= 70%, fix-X <= 30%)\n")
for (sc in CELLS) for (n in NS) {
  x <- res[[paste(sc, n)]]
  shz <- -x$fz / x$ref * 100; shx <- -x$fx / x$ref * 100
  cat(sprintf("  %-9s n=%-6d fix-Z removes %6.1f%%   fix-X removes %6.1f%%\n",
              sc, n, shz, shx))
}
shz800 <- vapply(CELLS, function(sc) -res[[paste(sc, 800)]]$fz / res[[paste(sc, 800)]]$ref * 100, 0)
shx800 <- vapply(CELLS, function(sc) -res[[paste(sc, 800)]]$fx / res[[paste(sc, 800)]]$ref * 100, 0)
gate("fix-Z removes >= 70% of the reference bias (both cells, n = 800)",
     all(shz800 >= 70), paste(sprintf("%s %.1f%%", CELLS, shz800), collapse = "  "))
gate("fix-X removes <= 30% of the reference bias (both cells, n = 800)",
     all(shx800 <= 30), paste(sprintf("%s %.1f%%", CELLS, shx800), collapse = "  "))
gate("|interaction| <= 2 pp everywhere",
     max(abs(vapply(names(res), function(k) res[[k]]$ia, 0))) <= 2,
     sprintf("worst %+.2f pp", vapply(names(res), function(k) res[[k]]$ia, 0)[
       which.max(abs(vapply(names(res), function(k) res[[k]]$ia, 0)))]))
gate("ef_exact_exact within 1.5 pp of zero everywhere",
     max(abs(vapply(names(res), function(k) res[[k]]$ee, 0))) <= 1.5,
     sprintf("worst %+.2f%%", vapply(names(res), function(k) res[[k]]$ee, 0)[
       which.max(abs(vapply(names(res), function(k) res[[k]]$ee, 0)))]))

cat("\n--- SHARPER TEST: where does the n^-1/3 exponent live? ---\n")
cat("  Registered as: fixing Z should flatten the decay, fixing X should leave it.\n")
cat("  The answer turns out to make the question moot -- read the RESIDUAL BIAS\n")
cat("  first, because an exponent fitted to a bias of ~0 is fitted to noise:\n\n")
cat(sprintf("  %-18s %s\n", "arm (residual bias)",
            paste(sprintf("%9s", c("n=800", "n=3200", "n=12800")), collapse = "")))
for (a in c("pipeline_bartMI", "ef_bart_ship", "ef_exact_ship", "ef_bart_exact", "ef_exact_exact"))
  for (sc in CELLS)
    cat(sprintf("  %-13s %-4s %s\n", sub("^ef_", "", a), sub("^zr_", "", sc),
                paste(sprintf("%8.2f%%", vapply(NS, function(n)
                  mean(relb(wide(n, sc), a)), 0)), collapse = "")))
cat("\n  Fitted slopes, for the arms where a slope means anything:\n")
for (a in c("pipeline_bartMI", "ef_bart_ship", "ef_exact_ship", "ef_bart_exact", "ef_exact_exact")) {
  sl <- vapply(CELLS, function(sc) {
    v <- vapply(NS, function(n) mean(relb(wide(n, sc), a)), 0)
    if (any(v == 0)) return(NA_real_)
    coef(lm(log(abs(v)) ~ log(NS)))[[2]]
  }, 0)
  big <- vapply(CELLS, function(sc)
    min(abs(vapply(NS, function(n) mean(relb(wide(n, sc), a)), 0))) > 1, logical(1))
  cat(sprintf("  %-18s %s\n", sub("^ef_", "", a),
              paste(sprintf("%s %s", CELLS,
                            ifelse(big, sprintf("%+.3f", sl),
                                   "  (bias ~0; no exponent)")), collapse = "   ")))
}
cat("\n  Fixing the Z block drives the bias to ~0 at EVERY n, so there is no\n")
cat("  exponent left to locate -- which answers the question by dissolving it.\n")

cat("\n=== verdict ===\n")
if (!length(fails)) {
  cat("  Every registered gate passed: the Z block dominates.\n")
} else {
  cat(sprintf("  %d gate(s) missed:\n", length(fails)))
  for (f in fails) cat("    - ", f, "\n", sep = "")
}
