# ============================================================
# 00_create_birthwt_spline_monotonic_example_data.R
# Public custom-formula demo dataset
#
# Dataset: MASS::birthwt
# Outcome: low, binary indicator of low birth weight
# Purpose: test brms custom_formula with
#   - s(age_z): spline smooth term
#   - mo(lwt_q): monotonic ordered quintile effect
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
})

dir.create("data", recursive = TRUE, showWarnings = FALSE)

make_ordered_quintile <- function(x) {
  q <- dplyr::ntile(x, 5)
  ordered(q, levels = 1:5)
}

add_missing <- function(x, prop = 0.08) {
  x <- x
  idx <- sample(seq_along(x), size = floor(length(x) * prop))
  x[idx] <- NA
  x
}

set.seed(20260521)

dat <- MASS::birthwt %>%
  as_tibble() %>%
  mutate(
    row_id = row_number(),

    # Outcome and categorical predictors
    low = factor(low, levels = c(0, 1)),
    race = factor(race),
    smoke = factor(smoke, levels = c(0, 1)),
    ht = factor(ht, levels = c(0, 1)),
    ui = factor(ui, levels = c(0, 1)),

    # Ordered quintile of maternal weight.
    # This is intentionally complete in this demo so mo(lwt_q) can be tested
    # without adding another derived-variable imputation dependency.
    lwt_q = make_ordered_quintile(lwt)
  )

# Add reproducible missingness to predictors only.
# Do not add missingness to the outcome or to lwt_q in this demo.
dat <- dat %>%
  mutate(
    age = add_missing(age, 0.06),
    lwt = add_missing(lwt, 0.08),
    race = add_missing(race, 0.05),
    smoke = add_missing(smoke, 0.05),
    ptl = add_missing(ptl, 0.08),
    ht = add_missing(ht, 0.04),
    ui = add_missing(ui, 0.04),
    ftv = add_missing(ftv, 0.08)
  )

saveRDS(
  dat,
  "data/birthwt_spline_monotonic_example.rds",
  compress = FALSE
)

message("Saved data/birthwt_spline_monotonic_example.rds")
message("Rows: ", nrow(dat))
message("Outcome distribution:")
print(table(dat$low, useNA = "ifany"))
message("lwt_q distribution:")
print(table(dat$lwt_q, useNA = "ifany"))
message("Missingness:")
print(colSums(is.na(dat)))
