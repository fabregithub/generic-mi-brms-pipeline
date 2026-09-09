#!/usr/bin/env Rscript
# =============================================================================
# V23 -- does the censored-exposure bias follow bias% = 300 * u^2?
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED. Criteria: run_v23_curvature.sh's header, committed
# before the run.
#
# THE NULL CELL IS READ FIRST. cv000 sets g == 0, so there is nothing for a
# linear draw to miss and the predictor is zero by construction. If cv000 is not
# near zero, `u` is not measuring what it claims and no other cell is readable.
#
# REPORTED IN bias/SE, NOT RELATIVE BIAS (PLAN §11, changed 2026-09-08), with the
# DIRECTION alongside -- away from the null is anti-conservative.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
setwd(.here)
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path(.here, "R", "dgp.R"))
source(file.path(.here, "R", "censoring.R"))
source(file.path(.here, "R", "robustness.R"))

NS <- c(800L, 3200L); TV <- 0.40; K <- 300
CELLS <- c("cv000","cv040","cv119","cv153","cv182","cv262","cv364")
CALIB <- "cv262"
fails <- character(0)
gate <- function(l, p, d = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(p)) "PASS" else "MISS", l,
              if (nzchar(d)) paste0("  --  ", d) else ""))
  if (!isTRUE(p)) fails <<- c(fails, l)
}
S <- do.call(rbind, lapply(NS, function(n) {
  d <- read.csv(sprintf("results/v23_n%d_summary.csv", n), stringsAsFactors = FALSE)
  d$n <- n; d }))
get <- function(sc, pr, col, n) S[[col]][S$n == n & S$scenario == sc & S$procedure == pr]

# u, recomputed from the scenario definitions rather than transcribed
scs <- v2_scenarios(); nm <- vapply(scs, `[[`, "", "name")
U <- vapply(CELLS, function(c) {
  sc <- scs[[match(c, nm)]]
  ef_unrep_curvature(make_truth(p = 3L, z_role = "pipe_nl",
                                nl_a = sc$nl_a, nl_c = sc$nl_c), nd_frac = sc$nd_frac)
}, 0)

cat("=== NULL CELL AND ORACLE (read first) ===\n")
for (n in NS)
  cat(sprintf("  n=%-5d cv000 exact %+6.2f%%  oracle %+6.2f%%\n", n,
              100*get("cv000","z21_exact","rel_bias",n),
              100*get("cv000","oracle","rel_bias",n)))
n0 <- max(abs(vapply(NS, function(n) 100*get("cv000","z21_exact","rel_bias",n), 0)))
om <- max(abs(unlist(lapply(NS, function(n) vapply(CELLS, function(c)
  100*get(c,"oracle","rel_bias",n), 0)))))
gate("oracle unbiased (<2%) in every cell and n", om < 2, sprintf("worst %.2f%%", om))
gate("cv000 null within 1.5 pp of zero", n0 <= 1.5, sprintf("worst %.2f%%", n0))

cat("\n=== THE LAW: X-block bias vs 300*u^2 ===\n")
cat(sprintf("  %-7s %8s %10s %10s %8s %9s %-15s\n",
            "cell","u","predicted","measured","miss pp","bias/SE","direction"))
meas <- numeric(0); bse <- numeric(0)
for (c in CELLS) {
  rb <- 100*get(c,"z21_exact","rel_bias",800L); se <- get(c,"z21_exact","emp_se",800L)
  pr <- K*U[[c]]^2
  meas[c] <- rb; bse[c] <- abs(rb/100*TV)/se
  cat(sprintf("  %-7s %8.4f %9.2f%% %9.2f%% %8.2f %9.2f %-15s\n", c, U[[c]], pr, rb,
              rb - pr, bse[c], if (rb > 0) "AWAY from null" else "toward null"))
}
pred <- K*U^2
off <- abs(meas - pred); off[CALIB] <- NA
gate("every predicted cell within 5 pp of 300*u^2", max(off, na.rm = TRUE) <= 5,
     sprintf("worst %s off by %+.2f pp", names(which.max(off)), max(off, na.rm = TRUE)))

cat("\n--- the quadratic form: log|bias| on log u ---\n")
ok <- U > 0 & abs(meas) > 0.5
f <- lm(log(abs(meas[ok])) ~ log(U[ok])); ci <- confint(f)["log(U[ok])", ]
cat(sprintf("  %d cells in fit | slope %+.3f  95%% CI [%+.3f, %+.3f]\n",
            sum(ok), coef(f)[[2]], ci[1], ci[2]))
gate("log-log slope CI CONTAINS 2 (the orthogonality argument survives)",
     ci[1] <= 2 && ci[2] >= 2, sprintf("[%+.3f, %+.3f]", ci[1], ci[2]))
# A "contains 2" gate can pass by being uninformative, so say whether it also
# contains the alternatives it was supposed to rule out.
alt <- c(`1 (linear)` = 1, `1.5` = 1.5, `2 (quadratic)` = 2, `3 (cubic)` = 3)
inside <- names(alt)[alt >= ci[1] & alt <= ci[2]]
cat(sprintf("  CI width %.2f; it also contains: %s\n", ci[2] - ci[1],
            if (length(inside) > 1) paste(inside, collapse = ", ") else "2 alone"))
if (1 >= ci[1] && 1 <= ci[2])
  cat("  -> the CI contains BOTH 1 and 2, so this test does not discriminate.\n     Read it as 'not refuted', not as 'confirmed'.\n")
kfit <- coef(lm(meas[U>0] ~ 0 + I(U[U>0]^2)))[[1]]
cat(sprintf("  freely fitted constant (no intercept, quadratic fixed): k = %.0f (registered 300)\n", kfit))

cat("\n--- direction stability, and the bar in bias/SE ---\n")
pos <- meas[U > 0.05] > 0
gate("sign positive (AWAY from null) in every cell with u > 0.05", all(pos),
     paste(sprintf("%s %+.1f", names(pos), meas[U > 0.05]), collapse = "  "))
gate("shipped-path bias/SE <= 0.3 (PLAN §11 default gate) in every cell",
     max(bse) <= 0.3, sprintf("worst %s at %.2f", names(which.max(bse)), max(bse)))

cat("\n--- flat in n? (the derivation says asymptotic) ---\n")
for (c in CELLS) {
  a <- 100*get(c,"z21_exact","rel_bias",800L); b <- 100*get(c,"z21_exact","rel_bias",3200L)
  cat(sprintf("  %-7s n=800 %+6.2f%%  n=3200 %+6.2f%%  shift %+5.2f pp\n", c, a, b, b - a))
}
sh <- max(abs(vapply(CELLS, function(c)
  100*(get(c,"z21_exact","rel_bias",3200L) - get(c,"z21_exact","rel_bias",800L)), 0)))
gate("bias flat between n = 800 and 3200 (within 3 pp)", sh <= 3,
     sprintf("worst shift %+.2f pp", sh))

cat("\n=== THE OTHER ARMS, and the warning about item 07 ===\n")
cat(sprintf("  %-7s %10s %10s %10s %10s\n", "cell","exact","bartMI","micePmm","smc_xgrid"))
for (c in CELLS)
  cat(sprintf("  %-7s %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n", c,
              100*get(c,"z21_exact","rel_bias",800L),
              100*get(c,"pipeline_bartMI","rel_bias",800L),
              100*get(c,"pipeline_micePmm","rel_bias",800L),
              100*get(c,"smc_xgrid","rel_bias",800L)))
xg <- vapply(CELLS, function(c) 100*get(c,"smc_xgrid","rel_bias",800L), 0)
xc <- vapply(CELLS, function(c) get(c,"smc_xgrid","coverage",800L), 0)
bm <- vapply(CELLS, function(c) 100*get(c,"pipeline_bartMI","rel_bias",800L), 0)
gate("smc_xgrid (item 07) is NOT worse than the shipped path in any cell",
     all(abs(xg) <= abs(bm) + 2),
     sprintf("worst: %s at %+.1f%% (coverage %.3f) vs bartMI %+.1f%%",
             names(which.max(abs(xg) - abs(bm))), xg[which.max(abs(xg)-abs(bm))],
             xc[which.max(abs(xg)-abs(bm))], bm[which.max(abs(xg)-abs(bm))]))

cat("\n=== verdict ===\n")
if (!length(fails)) cat("  Every registered gate passed.\n") else {
  cat(sprintf("  %d gate(s) missed:\n", length(fails)))
  for (f in fails) cat("    - ", f, "\n", sep = "")
}
