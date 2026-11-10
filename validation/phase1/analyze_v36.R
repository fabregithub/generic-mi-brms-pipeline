#!/usr/bin/env Rscript
# V36 gates. Three states (MET / NOT MET / UNDECIDED), because a gate that
# cannot tell "refuted" from "underpowered" has misled this programme eight
# times. Registered predictions: DESIGN_v36_selection.md §5.
#
# THE GATE IS THE INTERACTION, not selection bias and not imputation bias:
#   interaction(arm, cell) = (arm_cell - oracle_cell) - (arm_none - oracle_none)
# Within a cell the arm and its oracle are PAIRED on the same datasets, so the
# inner differences are taken per replication and only the outer contrast is
# across independent samples.
suppressPackageStartupMessages(library(stats))
tag <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "v36"
raw <- readRDS(sprintf("results/%s_latest.rds", tag))
raw <- raw[raw$procedure != "TASK_ERROR" & is.finite(raw$estimate), ]
tr  <- raw$estimand_true[1]
ARMS <- c("complete_case","sub_lod2","pipeline_block_fcs","smc_xgrid")

w <- reshape(raw[, c("cell","rep","procedure","estimate")],
             idvar = c("cell","rep"), timevar = "procedure", direction = "wide")
names(w) <- sub("^estimate[.]", "", names(w))
cells <- intersect(c("none","descnull","desc","desc2","outcome","outcome2","collider"),
                   unique(w$cell))

pd <- function(cell, arm) {            # paired per-rep difference, % of truth
  d <- w[w$cell == cell, ]
  100 * (d[[arm]] - d[["oracle"]]) / tr
}
state <- function(est, se, thresh = 0) {
  if (!is.finite(se) || se == 0) return("UNDECIDED")
  z <- (est - thresh) / se
  if (abs(z) > 2) (if (z > 0) "MET+" else "MET-") else "UNDECIDED"
}

cat("================ V36: selection x imputation ================\n")
cat(sprintf("target b1 = %.4f   cells: %s\n\n", tr, paste(cells, collapse=", ")))

cat("---- pure SELECTION bias (oracle arm, vs target) ----\n")
on <- 100 * (w[w$cell=="none","oracle"] - tr) / tr
for (cl in cells) {
  v <- 100 * (w[w$cell==cl,"oracle"] - tr) / tr
  cat(sprintf("  %-10s %+7.2f%% +/- %.2f   (n=%d)\n", cl, mean(v),
              sd(v)/sqrt(length(v)), length(v)))
}

cat("\n---- P1  instrument check ----\n")
p1ok <- TRUE
if ("descnull" %in% cells) {
  v <- 100 * (w[w$cell=="descnull","oracle"] - tr) / tr
  s <- state(mean(v) - mean(on), sqrt(var(v)/length(v) + var(on)/length(on)))
  cat(sprintf("  descnull vs none: %+.2f pp -> %s %s\n", mean(v)-mean(on), s,
              if (s == "UNDECIDED") "(required: no bias) OK" else "** LEAKING **"))
  if (s != "UNDECIDED") p1ok <- FALSE
}
for (cl in intersect(c("desc","outcome","collider"), cells)) {
  v <- 100 * (w[w$cell==cl,"oracle"] - tr) / tr
  s <- state(mean(v) - mean(on), sqrt(var(v)/length(v) + var(on)/length(on)))
  cat(sprintf("  %-9s vs none: %+.2f pp -> %s %s\n", cl, mean(v)-mean(on), s,
              if (s == "UNDECIDED") "** NO SELECTION BIAS -- mechanism not wired **" else "OK"))
  if (s == "UNDECIDED") p1ok <- FALSE
}
cat(sprintf("  P1: %s\n", if (p1ok) "MET" else "NOT MET -- read nothing below until fixed"))

cat("\n---- P2  THE GATE: does imputation make selection bias WORSE? ----\n")
# SIGN TRAP, fixed 2026-11-11. The first version of this scorer reported the
# SIGNED interaction and called a positive value "amplification". Selection bias
# in every live cell is NEGATIVE, so a positive interaction means the arm is
# nearer the truth, not further from it -- the opposite reading. The gate is
# therefore stated on |bias|, paired per replication against the SAME cell's
# oracle, which is sign-safe whatever direction the selection bias takes.
cat("  paired per rep: |arm - truth| - |oracle - truth|, pp of truth.\n")
cat("  POSITIVE = imputation makes it WORSE (P2's prediction).\n")
cat("  sub_lod2 borrows nothing from Y and is the mechanism control.\n")
cat(sprintf("\n  %-20s", "arm"))
for (cl in setdiff(cells, "none")) cat(sprintf(" %17s", cl))
cat("\n")
p2 <- list()
for (a in ARMS) {
  if (!a %in% names(w)) next
  cat(sprintf("  %-20s", a))
  for (cl in setdiff(cells, "none")) {
    d  <- w[w$cell == cl, ]
    dd <- 100 * (abs(d[[a]] - tr) - abs(d$oracle - tr)) / tr
    est <- mean(dd); se <- sd(dd)/sqrt(nrow(d))
    st  <- state(est, se)
    p2[[paste(a, cl)]] <- c(est, se)
    cat(sprintf(" %+6.2f+/-%-5.2f%-4s", est, se,
                switch(st, "MET+"=" W", "MET-"=" B", "UNDECIDED"="  ")))
  }
  cat("\n")
}
cat("\n  W = significantly WORSE than its own oracle (P2 direction)\n")
cat("  B = significantly BETTER -- P2 REFUTED in that cell: the imputation\n")
cat("      partially OFFSETS the selection bias rather than adding to it.\n")
key <- grep("^pipeline_block_fcs (desc|outcome|collider)$", names(p2), value = TRUE)
nW  <- sum(vapply(p2[key], function(v) v[1]/v[2] >  2, TRUE))
nB  <- sum(vapply(p2[key], function(v) v[1]/v[2] < -2, TRUE))
cat(sprintf("\n  P2 on pipeline_block_fcs: %d cell(s) WORSE, %d BETTER of %d -> %s\n",
            nW, nB, length(key),
            if (nW == length(key)) "MET" else if (nB > 0) "NOT MET (and reversed)" else "UNDECIDED"))

cat("\n---- P3  does the fix rescue selection bias? (it should NOT) ----\n")
cat("  smc_xgrid should land on oracle_sel, not on the target.\n")
for (cl in setdiff(cells, c("none","descnull"))) {
  if (!"smc_xgrid" %in% names(w)) break
  d  <- w[w$cell == cl, ]
  fx <- 100 * (d$smc_xgrid - tr) / tr          # fix vs TARGET
  oc <- 100 * (d$oracle    - tr) / tr          # selection bias that remains
  dif <- mean(fx) - mean(oc)
  se  <- sd(fx - oc) / sqrt(nrow(d))           # paired
  s   <- state(dif, se)
  cat(sprintf("  %-9s fix %+7.2f%%  oracle_sel %+7.2f%%  diff %+.2f +/- %.2f -> %s\n",
              cl, mean(fx), mean(oc), dif, se,
              if (s == "UNDECIDED") "P3 MET (fix does NOT rescue)" else
                "P3 NOT MET -- the fix moved it, which is a real finding"))
}
