# Create a public-data example input file for the generic MI + brms pipeline.
# This uses the built-in datasets::airquality data.

dir.create("data", recursive = TRUE, showWarnings = FALSE)

dat <- datasets::airquality

# Make Month categorical for the example.
dat$Month <- factor(dat$Month)

saveRDS(dat, "data/airquality_example.rds")

message("Saved: data/airquality_example.rds")
message("Rows: ", nrow(dat), ", columns: ", ncol(dat))
