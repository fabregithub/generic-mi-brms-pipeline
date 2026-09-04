#!/usr/bin/env Rscript
# =============================================================================
# V16 + V17 -- evaluate the registered predictions
# -----------------------------------------------------------------------------
# EVERY VERDICT IS COMPUTED. Nothing here is written as a literal: an earlier
# round of these scripts printed conclusions its own numbers contradicted, so
# each bar is compared in code and the pass/fail text is generated.
#
# The estimand is the PAIRED shift of each arm against pipeline_bartMI, within
# replication, in percentage points of the true b1 -- the same quantity
# predict_roles.R derived before the run. Pairing is what makes it a clean
# contrast: both arms see byte-identical data, so the reference's own bias
# cancels.
#
# Criteria: ../PLAN_pipeline_validation.md §8i (V16) and §8j (V17).
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
setwd(.here)

fails <- character(0)
gate <- function(label, pass, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(pass)) "PASS" else "MISS", label,
              if (nzchar(detail)) paste0("  --  ", detail) else ""))
  if (!isTRUE(pass)) fails <<- c(fails, label)
}

# ---- paired shifts, per cell ------------------------------------------------
shift <- function(raw, sc) {
  d <- subset(raw, scenario == sc)
  w <- reshape(d[, c("rep", "procedure", "estimate", "estimand_true")],
               idvar = c("rep", "estimand_true"), timevar = "procedure",
               direction = "wide")
  tv  <- w$estimand_true[1]
  ref <- w$`estimate.pipeline_bartMI`
  per <- function(arm) (w[[paste0("estimate.pipeline_", arm)]] - ref) / tv * 100
  arms <- c("noYx", "noYz", "noYboth")
  out <- lapply(arms, per); names(out) <- arms
  out$leak_v <- out$noYboth - out$noYx - out$noYz
  out$n <- nrow(w)
  out$ref_bias <- (mean(ref) - tv) / tv * 100
  out
}
mn <- function(v) mean(v); se <- function(v) sd(v) / sqrt(length(v))

report <- function(tag, cells) {
  raw <- read.csv(sprintf("results/%s_raw.csv", tag), stringsAsFactors = FALSE)
  st <- lapply(cells, function(sc) shift(raw, sc)); names(st) <- cells
  cat(sprintf("\n%-16s %8s %8s %8s %8s %8s\n", "cell", "noYx", "noYz", "noYboth", "leak", "ref"))
  for (sc in cells) {
    x <- st[[sc]]
    cat(sprintf("%-16s %8.2f %8.2f %8.2f %8.2f %7.2f%%\n", sc,
                mn(x$noYx), mn(x$noYz), mn(x$noYboth), mn(x$leak_v), x$ref_bias))
  }
  cat(sprintf("%-16s %8.2f %8.2f %8.2f %8.2f\n", "  (mc se)",
              mean(vapply(st, function(x) se(x$noYx), 0)),
              mean(vapply(st, function(x) se(x$noYz), 0)),
              mean(vapply(st, function(x) se(x$noYboth), 0)),
              mean(vapply(st, function(x) se(x$leak_v), 0))))
  st
}

cat("=== V16 (stage 1): the root claim ===\n")
V16 <- c("mcar_z40", "combined", "nl_mcar_z40", "nl_combined", "mar_z40", "nl_mar_z40")
s16 <- report("v16", V16)

cat("\n--- V16 gates (PLAN §8i) ---\n")
LIN <- c("mcar_z40", "combined")           # the linear cells the prediction was derived for
for (sc in LIN) {
  x <- s16[[sc]]
  gate(sprintf("%s: noYboth within 3 pp of -14.6", sc),
       abs(mn(x$noYboth) - (-14.6)) <= 3, sprintf("%.2f", mn(x$noYboth)))
  gate(sprintf("%s: noYx within 3 pp of -15.8", sc),
       abs(mn(x$noYx) - (-15.8)) <= 3, sprintf("%.2f", mn(x$noYx)))
  gate(sprintf("%s: |noYz| <= 1.5 pp", sc),
       abs(mn(x$noYz)) <= 1.5, sprintf("%.2f", mn(x$noYz)))
  gate(sprintf("%s: leak positive at 2x mc se", sc),
       mn(x$leak_v) - 2 * se(x$leak_v) > 0,
       sprintf("%.2f +/- %.2f", mn(x$leak_v), se(x$leak_v)))
}

cat("\n=== V17a (stage 2): the covariate's causal role ===\n")
V17 <- c("zr_fork", "zr_pipe", "zr_collider", "zr_mixed",
         "zrmar_fork", "zrmar_pipe", "zrmar_collider", "zrmar_mixed")
s17 <- report("v17a", V17)

# The precision row comes from stage 1 -- V17a has no precision cell, because
# mcar_z40 / mar_z40 already ARE that structure.
prec <- list(mcar = s16$mcar_z40, mar = s16$mar_z40)

cat("\n--- THE HEADLINE: the structural 2x2 (PLAN §8j) ---\n")
onpath <- c("fork", "pipe", "mixed")
for (mech in c("mcar", "mar")) {
  pre <- if (mech == "mcar") "zr_" else "zrmar_"
  big <- vapply(onpath, function(r) abs(mn(s17[[paste0(pre, r)]]$noYz)) > 10, logical(1))
  gate(sprintf("%s: |noYz| > 10 pp in fork, pipe and mixed", toupper(mech)),
       all(big),
       paste(sprintf("%s=%.1f", onpath,
                     vapply(onpath, function(r) mn(s17[[paste0(pre, r)]]$noYz), 0)),
             collapse = "  "))
  p <- prec[[mech]]
  gate(sprintf("%s: |noYz| < 3 pp in precision", toupper(mech)),
       abs(mn(p$noYz)) < 3, sprintf("%.2f", mn(p$noYz)))
}

cat("\n--- the mechanism term: MAR - MCAR on noYz ---\n")
# Paired per replication: seed_as gives each MCAR/MAR twin byte-identical
# complete data and identical exposure censoring.
for (r in c("fork", "pipe", "mixed", "collider")) {
  a <- s17[[paste0("zr_", r)]]$noYz; b <- s17[[paste0("zrmar_", r)]]$noYz
  d <- b - a; m <- mean(d); s <- sd(d) / sqrt(length(d))
  gate(sprintf("%-8s MAR - MCAR positive and <= 12 pp", r),
       (m - 2 * s > 0) && m <= 12, sprintf("%+.2f +/- %.2f", m, s))
}

cat("\n--- per-cell magnitudes (PLAN §8j) ---\n")
BAR <- list(
  list("zr_fork",         26.7, 5),  list("zrmar_fork",     32.4, 5),
  list("zr_pipe",         30.4, 5),  list("zrmar_pipe",     36.1, 5),
  list("zr_mixed",        19.5, 5),  list("zrmar_mixed",    23.8, 5),
  list("zrmar_collider",   3.1, 3))
for (b in BAR) {
  x <- s17[[b[[1]]]]; got <- mn(x$noYz)
  gate(sprintf("%-16s noYz within %g pp of %+.1f", b[[1]], b[[3]], b[[2]]),
       abs(got - b[[2]]) <= b[[3]], sprintf("%+.2f", got))
}
gate("zr_collider      |noYz| <= 1.5 pp",
     abs(mn(s17$zr_collider$noYz)) <= 1.5, sprintf("%+.2f", mn(s17$zr_collider$noYz)))
for (sc in c("zr_collider", "zrmar_collider"))
  gate(sprintf("%-16s leak negative at 2x mc se", sc),
       mn(s17[[sc]]$leak_v) + 2 * se(s17[[sc]]$leak_v) < 0,
       sprintf("%.2f +/- %.2f", mn(s17[[sc]]$leak_v), se(s17[[sc]]$leak_v)))
for (sc in c("zr_pipe", "zrmar_pipe"))
  gate(sprintf("%-16s noYz below the +75 pp total-effect bound", sc),
       mn(s17[[sc]]$noYz) < 75, sprintf("%+.2f", mn(s17[[sc]]$noYz)))

cat("\n=== V17b (stage 3): does use_as_auxiliary reach the X block? ===\n")
raw3 <- read.csv("results/v17b_raw.csv", stringsAsFactors = FALSE)
for (sc in c("zr_collider", "zrmar_collider")) {
  d <- subset(raw3, scenario == sc)
  w <- reshape(d[, c("rep", "procedure", "estimate", "estimand_true")],
               idvar = c("rep", "estimand_true"), timevar = "procedure", direction = "wide")
  tv <- w$estimand_true[1]
  sh <- (w$`estimate.pipeline_auxZ_asdoc` - w$`estimate.pipeline_auxZ_shipped`) / tv * 100
  cat(sprintf("  %-16s shipped %+.2f%%  as-documented %+.2f%%  |  paired difference %+.2f +/- %.2f pp\n",
              sc,
              (mean(w$`estimate.pipeline_auxZ_shipped`) - tv) / tv * 100,
              (mean(w$`estimate.pipeline_auxZ_asdoc`)   - tv) / tv * 100,
              mean(sh), sd(sh) / sqrt(length(sh))))
  # Significance and materiality are different questions, and at 500 reps the
  # first is easy to reach. Report both rather than letting a p-value stand in
  # for a decision.
  sig <- abs(mean(sh)) - 2 * sd(sh) / sqrt(length(sh)) > 0
  mat <- abs(mean(sh)) >= 1.0
  cat(sprintf("     -> %s at 2x mc se; %s (|shift| %s 1 pp)\n",
              if (sig) "DETECTABLE" else "not detectable",
              if (mat) "MATERIAL" else "not material",
              if (mat) ">=" else "<"))
}

cat("\n=== the shipped default's OWN bias, by covariate role ===\n")
# Every gate above is a paired contrast, in which the reference's bias cancels.
# That is what makes those contrasts clean -- and it also hides this. The
# reference arm IS the shipped configuration (Y in both blocks, bartMI), so its
# absolute bias is what a user actually gets.
cat(sprintf("  %-16s %10s %10s %10s %8s\n", "cell", "rel bias", "coverage", "width/se", "fmi"))
for (tag in c("v16", "v17a")) {
  sm <- read.csv(sprintf("results/%s_summary.csv", tag), stringsAsFactors = FALSE)
  sm <- subset(sm, procedure == "pipeline_bartMI")
  for (i in order(sm$scenario)) {
    r <- sm[i, ]
    cat(sprintf("  %-16s %9.2f%% %10.3f %10.3f %8.3f\n",
                r$scenario, 100 * r$rel_bias, r$coverage, r$width_se_ratio, r$fmi))
  }
}
sm <- read.csv("results/v17a_summary.csv", stringsAsFactors = FALSE)
sm <- subset(sm, procedure == "pipeline_bartMI")
worst <- sm[which.max(abs(sm$rel_bias)), ]
gate("shipped default stays within 10% bias in every V17a cell",
     max(abs(sm$rel_bias)) < 0.10,
     sprintf("worst: %s at %+.2f%% with coverage %.3f",
             worst$scenario, 100 * worst$rel_bias, worst$coverage))

cat("\n=== verdict ===\n")
if (length(fails) == 0L) {
  cat("  Every registered gate passed.\n")
} else {
  cat(sprintf("  %d gate(s) missed:\n", length(fails)))
  for (f in fails) cat("    - ", f, "\n", sep = "")
}
