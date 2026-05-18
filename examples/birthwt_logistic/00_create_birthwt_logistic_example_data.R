# ============================================================
# 00_create_birthwt_logistic_example_data.R
# Public logistic-regression demo dataset
#
# Dataset: MASS::birthwt
# Outcome: low, binary indicator of low birth weight
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
})

dir.create("data", recursive = TRUE, showWarnings = FALSE)

dat <- MASS::birthwt %>%
  as_tibble() %>%
  mutate(
    row_id = row_number(),

    # Convert categorical/binary variables to factors
    low = factor(low, levels = c(0, 1)),
    race = factor(race),
    smoke = factor(smoke, levels = c(0, 1)),
    ht = factor(ht, levels = c(0, 1)),
    ui = factor(ui, levels = c(0, 1))
  )

# Add reproducible missingness in predictors only.
# Do NOT add missingness to outcome low for the basic logistic demo.
set.seed(20260515)

add_missing <- function(x, prop = 0.08) {
  x <- x
  idx <- sample(seq_along(x), size = floor(length(x) * prop))
  x[idx] <- NA
  x
}

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
  "data/birthwt_logistic_example.rds",
  compress = FALSE
)

message("Saved data/birthwt_logistic_example.rds")
message("Rows: ", nrow(dat))
message("Outcome distribution:")
print(table(dat$low, useNA = "ifany"))

message("Missingness:")
print(colSums(is.na(dat)))
