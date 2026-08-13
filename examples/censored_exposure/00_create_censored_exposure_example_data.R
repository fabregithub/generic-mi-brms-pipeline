# Create a simulated example for the censored-exposure block-FCS strategy.
#
# Two focal EXPOSURES are censored below reporting limits and represented with
# leftcens interval columns (X_lo / X_hi):
#   * expo1 — single left-censoring (non-detect below one LOD).
#   * expo2 — three-tier: non-detect (< MDL) and detected-not-quantified
#             (MDL <= x < LCMRL, an interval), plus quantified values.
# The outcome depends on log(expo1) + log(expo2); a covariate is MAR and some
# outcomes are missing (to exercise MID). Truth is known, so the exposure-response
# coefficients (0.6 and 0.4) can be checked after fitting.

dir.create("data", recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
n <- 500
z1 <- rnorm(n)
z2 <- rbinom(n, 1, 0.5)

expo1 <- exp(0.4 * z1 + rnorm(n))            # right-skewed, covariate-linked
expo2 <- exp(0.3 * z2 + rnorm(n))

y <- 0.6 * log(expo1) + 0.4 * log(expo2) + 0.5 * z1 - 0.3 * z2 + rnorm(n)

# expo1: single left-censoring at ~30% non-detects.
L1 <- as.numeric(quantile(expo1, 0.30)); nd1 <- expo1 < L1

# expo2: three-tier limits (MDL / LCMRL) -> ND (~20%) + DNQ interval (~20%).
MDL   <- as.numeric(quantile(expo2, 0.20))
LCMRL <- as.numeric(quantile(expo2, 0.40))
nd2  <- expo2 < MDL
dnq2 <- expo2 >= MDL & expo2 < LCMRL

dat <- data.frame(
  row_id = seq_len(n),
  y  = y,
  z1 = z1,
  z2 = z2,
  # exposure columns: observed value where quantified, NA where censored.
  expo1 = ifelse(nd1, NA, expo1),
  expo2 = ifelse(nd2 | dnq2, NA, expo2),
  # leftcens interval bounds (lo == hi where quantified; lo <= 0 => left-censored).
  expo1_lo = ifelse(nd1, 0, expo1),
  expo1_hi = ifelse(nd1, L1, expo1),
  expo2_lo = ifelse(nd2, 0, ifelse(dnq2, MDL, expo2)),
  expo2_hi = ifelse(nd2, MDL, ifelse(dnq2, LCMRL, expo2))
)

dat$y[sample(n, 30)]  <- NA                  # missing outcomes -> MID
dat$z1[sample(n, 70)] <- NA                  # MAR covariate -> miceRanger Z-block

saveRDS(dat, "data/censored_exposure_example.rds")

message("Saved: data/censored_exposure_example.rds")
message("Rows: ", nrow(dat), ", columns: ", ncol(dat))
message("expo1 non-detects: ", sum(nd1), " (", round(100 * mean(nd1)), "%)")
message("expo2 non-detects: ", sum(nd2), ", detected-not-quantified: ", sum(dnq2))
message("Missing outcomes: 30, MAR covariate cells: 70")
