#!/usr/bin/env Rscript
# =============================================================================
# V19 -- does the bias exponent track the Z-block imputer?
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED. Criteria: ../PLAN_pipeline_validation.md §8l.
# The oracle control is read FIRST: if it drifts with n, the estimand or the DGP
# is n-dependent and no exponent below means anything.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
setwd(.here)
NS <- c(800L, 1600L, 3200L, 6400L, 12800L)
ONPATH <- c("zr_fork", "zr_pipe")
NONPAR <- c("pipeline_bartMI", "pipeline_properBoot")
PARAM  <- c("pipeline_micePmm", "pipeline_properZ")
ARMS   <- c(NONPAR, PARAM)
fails <- character(0)
gate <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "MISS", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- c(fails, l)
}
S <- do.call(rbind, lapply(NS, function(n) {
  d <- read.csv(sprintf("results/v19_n%d_summary.csv", n), stringsAsFactors = FALSE)
  d$n <- n; d
}))
rb <- function(sc, pr) vapply(NS, function(n)
  S$rel_bias[S$n == n & S$scenario == sc & S$procedure == pr] * 100, 0)

cat("=== CONTROL: the oracle ===\n")
om <- 0
for (sc in c(ONPATH, "mcar_z40")) {
  v <- rb(sc, "oracle"); om <- max(om, max(abs(v)))
  cat(sprintf("  %-10s %s\n", sc, paste(sprintf("%8.2f%%", v), collapse = "")))
}
gate("oracle within 2% at every n and cell", om < 2, sprintf("worst %.2f%%", om))

cat("\n=== relative bias by arm, on-path cells ===\n")
sl <- list()
for (sc in ONPATH) {
  cat(sprintf("\n  --- %s ---\n  %-22s %s %11s\n", sc, "arm",
              paste(sprintf("%9d", NS), collapse = ""), "slope [95% CI]"))
  for (a in ARMS) {
    v <- rb(sc, a); f <- lm(log(abs(v)) ~ log(NS)); ci <- confint(f)["log(NS)", ]
    sl[[paste(sc, a)]] <- list(b = coef(f)[["log(NS)"]], ci = ci, v = v)
    cat(sprintf("  %-22s %s   %+.3f [%+.3f,%+.3f]\n", sub("^pipeline_", "", a),
                paste(sprintf("%8.2f%%", v), collapse = ""),
                coef(f)[["log(NS)"]], ci[1], ci[2]))
  }
}

cat("\n=== PRIMARY: pooled slope per arm (both on-path cells) ===\n")
pool <- list()
for (a in ARMS) {
  d <- do.call(rbind, lapply(ONPATH, function(sc)
    data.frame(y = log(abs(rb(sc, a))), x = log(NS), cell = sc)))
  f <- lm(y ~ x + cell, data = d); ci <- confint(f)["x", ]
  pool[[a]] <- c(b = coef(f)[["x"]], lo = unname(ci[1]), hi = unname(ci[2]))
  fam <- if (a %in% NONPAR) "nonparametric" else "parametric"
  cat(sprintf("  %-22s %-14s slope %+.3f  95%% CI [%+.3f, %+.3f]\n",
              sub("^pipeline_", "", a), fam, coef(f)[["x"]], ci[1], ci[2]))
}

cat("\n--- each arm must sit NEAR its family's rate, not merely on the right side ---\n")
for (a in ARMS) {
  want <- if (a %in% NONPAR) -1/3 else -1/2
  other <- if (a %in% NONPAR) -1/2 else -1/3
  b <- pool[[a]][["b"]]
  gate(sprintf("%-22s within 0.08 of %+.3f, and >0.08 from %+.3f",
               sub("^pipeline_", "", a), want, other),
       abs(b - want) <= 0.08 && abs(b - other) > 0.08,
       sprintf("slope %+.3f", b))
}

cat("\n--- the family split: are the two groups separated? ---\n")
np <- vapply(NONPAR, function(a) pool[[a]][["b"]], 0)
pa <- vapply(PARAM,  function(a) pool[[a]][["b"]], 0)
gate("every parametric CI lies below every nonparametric CI",
     max(vapply(PARAM, function(a) pool[[a]][["hi"]], 0)) <
     min(vapply(NONPAR, function(a) pool[[a]][["lo"]], 0)),
     sprintf("nonpar %s  vs  par %s",
             paste(sprintf("%+.3f", np), collapse = "/"),
             paste(sprintf("%+.3f", pa), collapse = "/")))
gate("the two arms WITHIN each family agree (spread <= 0.08)",
     abs(diff(np)) <= 0.08 && abs(diff(pa)) <= 0.08,
     sprintf("nonpar spread %.3f, par spread %.3f", abs(diff(np)), abs(diff(pa))))

cat("\n--- what separates the arms instead: does the bias decay AT ALL? ---\n")
for (a in ARMS) {
  b <- pool[[a]]; dec <- b[["hi"]] < 0
  cat(sprintf("  %-22s slope %+.3f [%+.3f,%+.3f]  ->  %s\n", sub("^pipeline_", "", a),
              b[["b"]], b[["lo"]], b[["hi"]],
              if (dec) "DECAYS with n" else "ASYMPTOTIC -- CI includes 0"))
}
cat("\n  bias at the largest n tested (n = 12800), on-path cells:\n")
for (a in ARMS)
  cat(sprintf("  %-22s %s\n", sub("^pipeline_", "", a),
              paste(sprintf("%-12s %+.2f%%", ONPATH,
                            vapply(ONPATH, function(sc) rb(sc, a)[5], 0)), collapse = "   ")))

cat("\n=== CONTROL: mcar_z40, the causally inert cell ===\n")
cat("  (bias is small there, so the exponent is poorly determined -- this shows\n")
cat("   whether the family split appears WITHOUT a confounding path.)\n")
for (a in ARMS) {
  v <- rb("mcar_z40", a); f <- lm(log(abs(v)) ~ log(NS)); ci <- confint(f)["log(NS)", ]
  cat(sprintf("  %-22s %s   %+.3f [%+.3f,%+.3f]\n", sub("^pipeline_", "", a),
              paste(sprintf("%8.2f%%", v), collapse = ""),
              coef(f)[["log(NS)"]], ci[1], ci[2]))
}

cat("\n=== verdict ===\n")
if (!length(fails)) {
  cat("  Every registered gate passed: the exponent tracks the imputer family.\n")
} else {
  cat(sprintf("  %d gate(s) missed:\n", length(fails)))
  for (f in fails) cat("    - ", f, "\n", sep = "")
}
