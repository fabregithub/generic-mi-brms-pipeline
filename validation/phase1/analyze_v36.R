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

cat("\n---- P2  THE GATE: does imputation AMPLIFY selection bias? ----\n")
cat("  positive = imputation makes it worse.  sub_lod2 borrows nothing from Y\n")
cat("  and is the mechanism control: it must show NO interaction.\n")
cat(sprintf("\n  %-20s", "arm"))
for (cl in setdiff(cells, "none")) cat(sprintf(" %18s", cl))
cat("\n")
for (a in ARMS) {
  if (!a %in% names(w)) next
  dn <- pd("none", a)
  cat(sprintf("  %-20s", a))
  for (cl in setdiff(cells, "none")) {
    dd <- pd(cl, a)
    est <- mean(dd) - mean(dn)
    se  <- sqrt(var(dd)/length(dd) + var(dn)/length(dn))
    cat(sprintf(" %+7.2f+/-%-5.2f%-4s", est, se,
                switch(state(est, se), "MET+"="  *", "MET-"="  -", "UNDECIDED"="   ")))
  }
  cat("\n")
}
cat("\n  * = significantly positive (amplification)   - = significantly negative\n")

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
