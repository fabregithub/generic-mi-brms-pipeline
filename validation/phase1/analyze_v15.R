# =============================================================================
# Track V15 -- evaluate the registered prediction  excess% = 323.4 * f^2
# -----------------------------------------------------------------------------
# Persisted (unlike earlier tracks' ad-hoc analyses) because this is the first
# track whose headline is a VERDICT ON A PREDICTION, not a pass/fail against a
# threshold: the arithmetic that accepts or rejects the law has to be re-runnable.
#
# EVERY VERDICT BELOW IS COMPUTED. No conclusion is written as a literal -- an
# earlier round of these scripts printed narration that its own numbers
# contradicted, so the text is generated from the data or not at all.
#
# Estimator, reverse-engineered from V3 and checked against its published
# headline (-11.71 / -56.66 vs -11.7 / -56.7):
#     excess% = mean_reps(pipeline - oracle) / true * 100
# Paired within replication, so BKMR's own oracle floor cancels.
# =============================================================================
raw <- read.csv("results/v15_raw.csv", stringsAsFactors = FALSE)

F_OF <- c(mixf10 = 0.10, nd20 = 0.20, mixf30 = 0.30,
          nd40 = 0.40, mixf50 = 0.50, mixf60 = 0.60)
K_LAW    <- 323.4       # calibrated on V3's two points; curve-fitting, declared as such
TOL_PP   <- 10          # pre-declared per-cell tolerance, percentage points
V3_KNOWN <- c(nd20 = -11.71, nd40 = -56.66)

paired <- function(d, sc, est = "curv_X1") {
  d <- subset(d, scenario == sc & estimand == est &
                procedure %in% c("oracle_bkmr", "pipeline_bkmr"))
  w <- reshape(d[, c("rep", "procedure", "estimate", "estimand_true")],
               idvar = c("rep", "estimand_true"), timevar = "procedure",
               direction = "wide")
  w <- w[stats::complete.cases(w), ]
  diff <- (w$estimate.pipeline_bkmr - w$estimate.oracle_bkmr) / w$estimand_true[1] * 100
  list(n = length(diff), excess = mean(diff), se = sd(diff) / sqrt(length(diff)),
       per_rep = diff)
}

cells <- lapply(names(F_OF), function(sc) paired(raw, sc))
names(cells) <- names(F_OF)

cat("=== V15: measured vs predicted, curv_X1 paired excess over oracle ===\n\n")
cat(sprintf("%-8s %5s %4s %10s %10s %8s %8s %s\n",
            "cell", "f", "n", "measured", "predicted", "diff_pp", "mcse", "within 10pp?"))
tab <- do.call(rbind, lapply(names(F_OF), function(sc) {
  f <- F_OF[[sc]]; cc <- cells[[sc]]
  pred <- -K_LAW * f^2                     # sign: attenuation is negative
  dpp  <- cc$excess - pred
  ok   <- abs(dpp) <= TOL_PP
  cat(sprintf("%-8s %5.2f %4d %9.2f%% %9.2f%% %8.2f %8.2f %s\n",
              sc, f, cc$n, cc$excess, pred, dpp, cc$se,
              if (ok) "yes" else "NO"))
  data.frame(scenario = sc, f = f, n = cc$n, measured = cc$excess,
             mcse = cc$se, predicted = pred, diff_pp = dpp, within = ok)
}))

cat("\n--- calibration: do V3's cells reproduce? ---\n")
calib <- do.call(rbind, lapply(names(V3_KNOWN), function(sc) {
  m <- cells[[sc]]$excess; v3 <- V3_KNOWN[[sc]]
  # V3 ran 200 reps, V15 runs 100; compare against the pooled MC error.
  d <- m - v3
  data.frame(scenario = sc, v15 = m, v3 = v3, diff = d,
             mcse_v15 = cells[[sc]]$se,
             consistent = abs(d) <= 2 * cells[[sc]]$se)
}))
print(calib, row.names = FALSE, digits = 4)
cat(if (all(calib$consistent))
      "  -> V3 reproduced; the four new cells are comparable to the two old ones.\n"
    else
      "  -> DIVERGENCE from V3. Something changed; treat cross-track comparison as void.\n")

cat("\n--- falsification test 1: per-cell tolerance (", TOL_PP, "pp) ---\n", sep = "")
new <- subset(tab, !(scenario %in% names(V3_KNOWN)))
fail1 <- subset(new, !within)
cat(sprintf("  unmeasured cells tested: %d | outside tolerance: %d%s\n",
            nrow(new), nrow(fail1),
            if (nrow(fail1)) paste0(" (", paste(fail1$scenario, collapse = ", "), ")") else ""))

cat("\n--- falsification test 2: log-log slope, is the exponent 2? ---\n")
# Fitted across cells on |excess|; sign is uniform where attenuation holds.
fit_ok <- tab$measured < 0
if (!all(fit_ok))
  cat(sprintf("  NOTE: %d cell(s) have positive excess (%s); excluded from the log fit.\n",
              sum(!fit_ok), paste(tab$scenario[!fit_ok], collapse = ", ")))
lf <- tab[fit_ok, ]
m  <- lm(log(abs(measured)) ~ log(f), data = lf)
ci <- confint(m)["log(f)", ]
cat(sprintf("  cells in fit: %d | slope = %.3f  95%% CI [%.3f, %.3f]\n",
            nrow(lf), coef(m)[["log(f)"]], ci[1], ci[2]))
excl2 <- (ci[1] > 2 || ci[2] < 2)
cat(sprintf("  CI %s 2 -> test 2 %s\n",
            if (excl2) "EXCLUDES" else "contains",
            if (excl2) "REJECTS f^2" else "does not reject f^2"))

# Rep-level bootstrap: propagate Monte-Carlo error into the exponent.
set.seed(20260903)
bs <- replicate(2000, {
  ex <- sapply(lf$scenario, function(sc) mean(sample(cells[[sc]]$per_rep, replace = TRUE)))
  if (any(ex >= 0)) return(NA_real_)
  coef(lm(log(abs(ex)) ~ log(lf$f)))[[2]]
})
bci <- quantile(bs, c(.025, .975), na.rm = TRUE)
cat(sprintf("  bootstrap slope 95%% CI [%.3f, %.3f] (%d/%d draws usable)\n",
            bci[1], bci[2], sum(!is.na(bs)), length(bs)))

cat("\n--- verdict ---\n")
rejected <- nrow(fail1) > 0 || excl2
cat(sprintf("  registered law excess%% = %.1f * f^2 is %s\n", K_LAW,
            if (rejected) "REJECTED" else "NOT REJECTED"))
if (nrow(fail1))
  cat(sprintf("    by test 1: %s\n",
              paste(sprintf("%s off by %.1f pp", fail1$scenario, fail1$diff_pp),
                    collapse = "; ")))
if (excl2) cat("    by test 2: the exponent is not 2.\n")

cat("\n--- what the shape actually is (for the theory, not a new claim) ---\n")
# The log-log fit above minimises LOG residuals, so comparing it to the
# registered law on a PERCENTAGE-POINT scale would flatter the registered law
# for free. Refit k and p on the pp scale before any such comparison.
nl <- try(nls(measured ~ -k * f^p, data = lf,
              start = list(k = K_LAW, p = 2)), silent = TRUE)
alt <- data.frame(form = "k*f^2 (registered)", k = K_LAW, p = 2)
if (!inherits(nl, "try-error"))
  alt <- rbind(alt, data.frame(form = "k*f^p (fitted, pp scale)",
                               k = coef(nl)[["k"]], p = coef(nl)[["p"]]))
for (i in seq_len(nrow(alt))) {
  pr <- -alt$k[i] * lf$f ^ alt$p[i]
  cat(sprintf("  %-26s k=%6.1f p=%.2f  max|resid| = %5.2f pp  RMS = %5.2f pp\n",
              alt$form[i], alt$k[i], alt$p[i],
              max(abs(lf$measured - pr)), sqrt(mean((lf$measured - pr)^2))))
}
cat("  (a free exponent cannot fit worse than a fixed one on its own scale;\n")
cat("   the question is whether the improvement is worth the extra parameter.)\n")

cat(sprintf("\n  curvature destroyed, by f: %s\n",
            paste(sprintf("%.2f->%.0f%%", lf$f, lf$measured / -100 * 100),
                  collapse = "  ")))
# Does the estimate cross zero anywhere -- i.e. is curvature inverted, not just lost?
sg <- do.call(rbind, lapply(names(F_OF), function(sc) {
  d <- subset(raw, scenario == sc & estimand == "curv_X1" &
                procedure == "pipeline_bkmr")
  data.frame(scenario = sc, f = F_OF[[sc]], true = d$estimand_true[1],
             mean_est = mean(d$estimate), frac_negative = mean(d$estimate < 0))
}))
cat("\n  pipeline point estimate vs the true curvature (", sprintf("%.4f", sg$true[1]),
    "):\n", sep = "")
print(sg[, c("scenario", "f", "mean_est", "frac_negative")], row.names = FALSE, digits = 3)
inv <- subset(sg, mean_est < 0)
cat(if (nrow(inv))
      sprintf("  -> sign INVERTED at %s: the fitted surface curves the wrong way.\n",
              paste(inv$scenario, collapse = ", "))
    else "  -> no cell inverts the sign; curvature is attenuated but never reversed.\n")

write.csv(tab, "results/v15_prediction_test.csv", row.names = FALSE)
cat("\nwrote results/v15_prediction_test.csv\n")
