#!/usr/bin/env Rscript
# V35 gates. Scored with MC error, in three states (MET / NOT MET / UNDECIDED),
# because this project has spoiled seven gates by scoring a trend against noise.
# The GATE is a comparison of SHAPES (log-log slopes), not of levels.
suppressPackageStartupMessages(library(stats))
tag <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "v35"
s <- read.csv(sprintf("results/%s_summary.csv", tag), stringsAsFactors = FALSE)
fmap <- c(mixf10=0.10, nd20=0.20, mixf30=0.30, nd40=0.40, mixf50=0.50, mixf60=0.60)
s$f <- fmap[s$scenario]
s <- s[!is.na(s$f), ]

excess <- function(arm, orc, est = "curv_X1") {
  a <- s[s$procedure==arm & s$estimand==est, ]
  o <- s[s$procedure==orc & s$estimand==est, ]
  if (!nrow(a) || !nrow(o)) return(NULL)
  m <- merge(a, o, by = "f", suffixes = c("_a","_o"))
  data.frame(f = m$f,
             excess = 100*((1+m$rel_bias_a)/(1+m$rel_bias_o) - 1),
             mcse = 100*sqrt((m$mcse_bias_a/abs(m$true_a))^2 +
                             (m$mcse_bias_o/abs(m$true_o))^2))
}
slope <- function(d) {
  d <- d[d$f >= 0.2 & is.finite(d$excess) & abs(d$excess) > 0, ]
  if (nrow(d) < 3) return(c(NA, NA, NA))
  fit <- lm(log(abs(excess)) ~ log(f), data = d)
  c(coef(fit)[2], confint(fit)[2, ])
}

cat("=========== V35: f-scaling by analysis model (curv_X1) ===========\n")
tab <- NULL
for (p in c("pipeline_bkmr","pipeline_lm")) {
  orc <- if (grepl("_lm$", p)) "oracle_lm" else "oracle_bkmr"
  d <- excess(p, orc); if (is.null(d)) next
  cat(sprintf("\n-- %s (vs %s) --\n", p, orc))
  print(d, row.names = FALSE, digits = 3)
  sl <- slope(d)
  cat(sprintf("   log-log slope (f>=0.2): %.2f  [%.2f, %.2f]\n", sl[1], sl[2], sl[3]))
  tab <- rbind(tab, data.frame(arm = p, slope = sl[1], lo = sl[2], hi = sl[3]))
}

cat("\n---------------- P1 (THE GATE) ----------------\n")
if (!is.null(tab) && nrow(tab) == 2) {
  b <- tab[tab$arm=="pipeline_bkmr", ]; l <- tab[tab$arm=="pipeline_lm", ]
  overlap <- !(b$lo > l$hi || l$lo > b$hi)
  cat(sprintf("bkmr %.2f [%.2f, %.2f]   lm %.2f [%.2f, %.2f]\n",
              b$slope, b$lo, b$hi, l$slope, l$lo, l$hi))
  if (overlap) {
    cat("P1 NOT MET -- the CIs OVERLAP. P-B's localisation is WRONG:\n")
    cat("   the exponent is not the analysis model's, and the imputation is back in play.\n")
  } else if (b$slope > l$slope) {
    cat("P1 MET -- disjoint CIs, BKMR steeper. The f^2 exponent is a property of\n")
    cat("   the ANALYSIS MODEL's response to degraded curvature information.\n")
  } else {
    cat("P1 NOT MET -- disjoint but in the WRONG DIRECTION (OLS steeper than BKMR).\n")
  }
} else cat("P1 UNDECIDED -- both arms needed.\n")

cat("\n---------------- P2 (instrument check) ----------------\n")
o <- s[s$procedure=="oracle_lm" & s$estimand %in% c("curv_X1","int_X1X2"), ]
if (nrow(o)) {
  o$rel_pct <- 100*o$rel_bias; o$mcse_pct <- 100*o$mcse_bias/abs(o$true)
  print(o[order(o$estimand, o$f), c("estimand","f","rel_pct","mcse_pct","coverage")],
        row.names = FALSE, digits = 3)
  for (e in unique(o$estimand)) {
    g <- o[o$estimand==e, ]
    tr <- summary(lm(rel_bias ~ f, data = g))$coefficients
    z  <- max(abs(g$rel_pct)/g$mcse_pct)
    cat(sprintf("   %-9s trend-in-f p = %.3f | max |z| from 0 = %.2f -> %s\n", e,
                tr[2,4], z,
                if (tr[2,4] > 0.05 && z < 3) "MET" else "CHECK: instrument not clean"))
  }
} else cat("P2 UNDECIDED -- oracle_lm not in this run.\n")

cat("\n---------------- P3 (V15 reproduction) ----------------\n")
d <- excess("pipeline_bkmr","oracle_bkmr")
if (!is.null(d)) {
  ref <- c("0.2" = -13.93, "0.4" = -55.62)
  for (ff in names(ref)) {
    r <- d[abs(d$f - as.numeric(ff)) < 1e-9, ]
    if (!nrow(r)) next
    z <- (r$excess - ref[[ff]])/r$mcse
    cat(sprintf("   f = %s: now %.1f%% +/- %.1f vs V15 %.1f%%  -> %.1f sigma %s\n",
                ff, r$excess, r$mcse, ref[[ff]], z,
                if (abs(z) < 3) "MET" else "NOT MET (V15 not comparable)"))
  }
} else cat("P3 UNDECIDED -- pipeline_bkmr not in this run.\n")
