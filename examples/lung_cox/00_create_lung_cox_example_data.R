# ============================================================
# 00_create_lung_cox_example_data.R
# Public survival-analysis demo dataset
#
# Dataset: survival::lung
# Outcome: time to death (days), right-censored
# Reference: Loprinzi et al. (1994) J Clin Oncol
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

dir.create("data", recursive = TRUE, showWarnings = FALSE)

dat <- survival::lung %>%
  as_tibble() %>%
  transmute(
    row_id  = row_number(),

    # Survival response:
    #   time     = days to death or last follow-up
    #   censored = 0 means event (death), 1 means right-censored
    #   brms cox() uses this convention.
    time     = time,
    censored = ifelse(status == 2, 0L, 1L),

    # Predictors
    age      = age,
    sex      = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    ph_ecog  = factor(ph.ecog, levels = 0:3),   # ECOG performance score
    ph_karno = ph.karno,                          # Karnofsky score (physician)
    wt_loss  = wt.loss,                           # Weight loss last 6 months (kg)
    meal_cal = meal.cal                            # Calories consumed at meals
  )

# Add reproducible extra missingness to predictors only.
# (ph_karno and wt_loss already have NAs in the source data.)
set.seed(20260515)

add_missing <- function(x, prop) {
  idx <- sample(seq_along(x), size = floor(length(x) * prop))
  x[idx] <- NA
  x
}

dat <- dat %>%
  mutate(
    age      = add_missing(age,      0.04),
    ph_karno = add_missing(ph_karno, 0.06),
    meal_cal = add_missing(meal_cal, 0.06)
  )

saveRDS(
  dat,
  "data/lung_cox_example.rds",
  compress = FALSE
)

message("Saved data/lung_cox_example.rds")
message("Rows: ", nrow(dat))
message("Events / censored:")
print(table(dat$censored, dnn = "censored (1=censored, 0=event)"))
message("Missingness:")
print(colSums(is.na(dat)))
