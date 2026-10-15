# =============================================================================
# V31 -- is the bias expansion merely TOO LOW-ORDER at realistic |d|?
# -----------------------------------------------------------------------------
# The `d` instrument (R/instrument_d.R) is exact on contamination of KNOWN size
# (self-test 0.99-1.10 at |d| <= 0.15, all three families, both channels) and
# not on the real arms (pilot: median |pred/meas - 1| = 0.59). Shrinking the
# missing fraction shrinks |d| without touching the mechanism, so:
#
#   P1  (the deciding gate) pred/meas RISES toward 1 as the fraction falls, in
#       ALL THREE families -> the expansion is right and needs higher order.
#       Flat, or inconsistent across families -> the structure is wrong.
#   P1b at the 0.10 rung the ratio is within 25% of 1.
#   P0  pred_b matches the measured bias within 25% wherever it is measurable.
#
# FLOOR: a ratio is only formed where |measured| > 2 x its MC error. Seven gates
# in this programme have been spoiled by a denominator that can reach zero.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
suppressMessages({source(file.path(.here, "R", "dgp.R"))
                  source(file.path(.here, "R", "robustness.R"))})
`%||%` <- function(a, b) if (is.null(a)) b else a

TAG <- Sys.getenv("TAG", "v31")
N   <- as.integer(Sys.getenv("NOBS", "3200"))
LADDER <- list(
  gaussian = c("fr10_gau", "fr20_gau", "zr_fork"),
  binomial = c("fr10_bin", "fr20_bin", "bin_fork"),
  ordinal3 = c("fr10_or3", "fr20_or3", "ord3_fork"))
FRAC <- c(0.10, 0.20, 0.40)
ARMS <- c("pipeline_properZ", "pipeline_properZ_ds", "pipeline_bartMI", "pipeline_noYx")

scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
truth_of <- function(nm) {
  s <- scs[[nm]]
  make_truth(p = 3L, erf_form = s$erf_form, y_form = s$y_form %||% "linear",
             z_role = s$z_role %||% "precision", target = s$target %||% "direct",
             delta_xz = s$delta_xz, delta_zx = s$delta_zx,
             y_family = s$y_family %||% "gaussian", y_tau = s$y_tau)$estimand_true
}
f <- file.path(.here, "results", sprintf("%s_n%d_latest.rds", TAG, N))
if (!file.exists(f)) stop("no V31 results at ", f)
d <- as.data.frame(readRDS(f)$summary)
d$truth    <- vapply(as.character(d$scenario), truth_of, 0)
d$measured <- d$mean_est - d$truth
d$mcse     <- d$emp_se / sqrt(d$n_ok)
g <- function(cl, ar, fld) { v <- d[[fld]][d$scenario == cl & d$procedure == ar]
                             if (length(v) != 1L) NA_real_ else v }

cat(sprintf("\n=== V31  n = %d  (%d reps) ===\n", N, max(d$n_ok, na.rm = TRUE)))
cat("\n--- d_rms by rung (the lever) ---\n")
for (fam in names(LADDER)) cat(sprintf("  %-9s %s\n", fam,
  paste(sprintf("%.2f", vapply(LADDER[[fam]], function(c) mean(
    d$d_rms[d$scenario == c & d$procedure %in% ARMS], na.rm = TRUE), 0)), collapse = " -> ")))

tab <- NULL
for (fam in names(LADDER)) for (i in seq_along(FRAC)) {
  cl <- LADDER[[fam]][i]
  for (ar in ARMS) {
    m <- g(cl, ar, "measured"); s <- g(cl, ar, "mcse"); p <- g(cl, ar, "pred_b")
    tab <- rbind(tab, data.frame(fam = fam, frac = FRAC[i], cell = cl, arm = ar,
      measured = m, mcse = s, pred = p,
      ratio = if (is.finite(m) && is.finite(s) && abs(m) > 2 * s) p / m else NA_real_,
      borrow = g(cl, ar, "borrow"), atten = g(cl, ar, "atten"),
      d_rms = g(cl, ar, "d_rms")))
  }
}
cat("\n--- per cell x arm (ratio shown only above the floor) ---\n")
pr <- tab; pr$arm <- sub("pipeline_", "", pr$arm)
print(format(pr[, c("fam","frac","arm","measured","mcse","pred","ratio","d_rms")], digits = 3),
      row.names = FALSE)

cat("\n--- P1  does pred/meas rise toward 1 as the fraction falls? ---\n")
p1 <- TRUE
for (fam in names(LADDER)) {
  mr <- vapply(FRAC, function(fr) {
    v <- tab$ratio[tab$fam == fam & tab$frac == fr]; v <- v[is.finite(v)]
    if (length(v)) median(v) else NA_real_ }, 0)
  nn <- vapply(FRAC, function(fr) sum(is.finite(tab$ratio[tab$fam == fam & tab$frac == fr])), 0L)
  dist <- abs(mr - 1)
  closer <- all(diff(dist[!is.na(dist)]) > 0)   # distance from 1 grows with frac
  cat(sprintf("  %-9s frac 0.10/0.20/0.40 -> %s   (|ratio-1| = %s)  n_arms %s  %s\n",
      fam, paste(sprintf("%6.2f", mr), collapse = " "),
      paste(sprintf("%.2f", dist), collapse = " "), paste(nn, collapse = "/"),
      if (isTRUE(closer)) "rises toward 1" else "DOES NOT"))
  if (!isTRUE(closer)) p1 <- FALSE
}
cat(sprintf("  P1: %s\n", if (isTRUE(p1)) "MET -- E1, the expansion is right and too low-order"
                          else "NOT MET -- E2, the structure does not survive"))

cat("\n--- P1b  at the 0.10 rung, within 25% of 1? ---\n")
p1b <- TRUE
for (fam in names(LADDER)) {
  v <- tab$ratio[tab$fam == fam & tab$frac == 0.10]; v <- v[is.finite(v)]
  if (!length(v)) { cat(sprintf("  %-9s VOID (no arm above the floor)\n", fam)); p1b <- NA; next }
  ok <- abs(median(v) - 1) <= 0.25
  cat(sprintf("  %-9s median %.3f over %d arms  %s\n", fam, median(v), length(v),
              if (ok) "ok" else "MISS"))
  if (!ok) p1b <- FALSE
}
cat(sprintf("  P1b: %s\n", if (isTRUE(p1b)) "MET" else if (is.na(p1b)) "VOID" else "NOT MET"))

cat("\n--- P0  overall agreement where measurable ---\n")
v <- tab$ratio[is.finite(tab$ratio)]
cat(sprintf("  %d of %d cell x arm above the floor | median |ratio-1| = %.2f | within 25%%: %d\n",
            length(v), nrow(tab), median(abs(v - 1)), sum(abs(v - 1) <= 0.25)))
cat(sprintf("  P0: %s\n", if (median(abs(v - 1)) <= 0.25) "MET" else "NOT MET"))

cat("\n--- diagnostic: is the prediction a near-cancellation? ---\n")
# pred = borrow + atten, so the cancellation is |sum| / sum|.|, NOT |difference|.
tab$cancel <- abs(tab$borrow + tab$atten) / pmax(abs(tab$borrow) + abs(tab$atten), 1e-12)
for (fam in names(LADDER)) for (fr in FRAC) {
  z <- tab[tab$fam == fam & tab$frac == fr, ]
  cat(sprintf("  %-9s frac %.2f  |borrow|+|atten| = %.4f   |pred| = %.4f   pred/sum = %.2f\n",
      fam, fr, mean(abs(z$borrow) + abs(z$atten), na.rm = TRUE),
      mean(abs(z$pred), na.rm = TRUE), mean(z$cancel, na.rm = TRUE)))
}
