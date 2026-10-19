# =============================================================================
# V32 -- is the borrowing scale a FAMILY constant, or arm-dependent?
# -----------------------------------------------------------------------------
# Gates are on the CHANNEL SCALES in `measured = a*borrow + b*atten`, never on
# pred/meas: V31 showed that ratio is dominated by a fixed cancellation
# amplification (x1.4-x3.8) which says nothing about the channels.
#   P1 one (a,b) per family fits the 12 points at chi2/df <= 2
#   P2 an arm-specific scale does not significantly improve it (LRT p > 0.05)
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
suppressMessages({source(file.path(.here, "R", "dgp.R"))
                  source(file.path(.here, "R", "robustness.R"))})
`%||%` <- function(a, b) if (is.null(a)) b else a
TAG <- Sys.getenv("TAG", "v32"); N <- as.integer(Sys.getenv("NOBS", "3200"))
FAM <- list(gaussian = c("fr10_gau","fr20_gau"),
            binomial = c("fr10_bin","fr20_bin"),
            ordinal3 = c("fr10_or3","fr20_or3"))
ARMS <- c("pipeline_properZ","pipeline_properZ_ds",
          "pipeline_bartMI","pipeline_noYx","pipeline_noYz")

scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
tru <- function(nm) { s <- scs[[nm]]
  make_truth(p = 3L, erf_form = s$erf_form, z_role = s$z_role %||% "precision",
             target = s$target %||% "direct", delta_xz = s$delta_xz,
             delta_zx = s$delta_zx, y_family = s$y_family %||% "gaussian",
             y_tau = s$y_tau)$estimand_true }
f <- file.path(.here, "results", sprintf("%s_n%d_latest.rds", TAG, N))
if (!file.exists(f)) stop("no V32 results at ", f)
d <- as.data.frame(readRDS(f)$summary)
d$measured <- d$mean_est - vapply(as.character(d$scenario), tru, 0)
d$mcse <- d$emp_se / sqrt(d$n_ok)
d$surv <- abs(d$borrow + d$atten) / pmax(abs(d$borrow) + abs(d$atten), 1e-12)

cat(sprintf("\n=== V32  n = %d  (%d reps) ===\n", N, max(d$n_ok, na.rm = TRUE)))
cat("\n--- conditioning spread (|T1+T2| / (|T1|+|T2|)) by arm ---\n")
z <- d[d$procedure %in% ARMS, ]
print(round(sort(tapply(z$surv, sub("pipeline_", "", z$procedure), mean, na.rm = TRUE)), 3))

cat("\n--- P1  one (a,b) per family over 5 arms x 2 fractions ---\n")
p1 <- TRUE; fits <- list()
for (fam in names(FAM)) {
  w <- d[d$scenario %in% FAM[[fam]] & d$procedure %in% ARMS, ]
  w <- w[is.finite(w$borrow) & is.finite(w$atten) & is.finite(w$measured), ]
  m0 <- lm(measured ~ 0 + borrow + atten, data = w, weights = 1/w$mcse^2)
  co <- summary(m0)$coefficients
  chi <- sum((w$measured - predict(m0))^2 / w$mcse^2) / (nrow(w) - 2)
  fits[[fam]] <- list(m0 = m0, w = w)
  cat(sprintf("  %-9s a = %+.3f +/- %.3f | b = %+.3f +/- %.3f | chi2/df = %6.2f (n=%d)  %s\n",
      fam, co["borrow","Estimate"], co["borrow","Std. Error"],
      co["atten","Estimate"], co["atten","Std. Error"], chi, nrow(w),
      if (chi <= 2) "fits" else "DOES NOT FIT"))
  if (chi > 2) p1 <- FALSE
}
cat(sprintf("  P1: %s\n", if (isTRUE(p1))
  "MET -- `a` is a FAMILY constant; a transfer law exists"
  else "NOT MET -- `a` is not a single family constant"))

cat("\n--- P2  does an ARM-specific borrowing scale add anything? ---\n")
p2 <- TRUE
for (fam in names(FAM)) {
  w <- fits[[fam]]$w; m0 <- fits[[fam]]$m0
  if (length(unique(w$procedure)) < 3L) { cat(sprintf("  %-9s VOID\n", fam)); next }
  m1 <- lm(measured ~ 0 + borrow:procedure + atten, data = w, weights = 1/w$mcse^2)
  an <- anova(m0, m1)
  pv <- an[["Pr(>F)"]][2]
  cat(sprintf("  %-9s arm-specific a:  F = %5.2f  p = %.4f  -> %s\n", fam,
      an[["F"]][2], pv, if (is.na(pv) || pv > 0.05) "no improvement" else "ARM-DEPENDENT"))
  if (!is.na(pv) && pv <= 0.05) p2 <- FALSE
}
cat(sprintf("  P2: %s\n", if (isTRUE(p2)) "MET -- no arm dependence detected"
                          else "NOT MET -- the borrowing scale depends on the arm"))

cat("\n--- per cell x arm (for the record) ---\n")
z <- d[d$procedure %in% ARMS, c("scenario","procedure","measured","mcse","borrow","atten","pred_b","surv")]
z$procedure <- sub("pipeline_", "", z$procedure)
print(format(z[order(z$scenario, z$procedure), ], digits = 3), row.names = FALSE)
