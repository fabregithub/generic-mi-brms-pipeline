# Run the full generic MI + brms pipeline.
# For the public airquality example, run 00_create_airquality_example_data.R once first.

steps <- c(
  "01_validate_config.R",
  "02_prepare_data.R",
  "03_impute.R",
  "04_fit_models.R",
  "05_diagnostics.R",
  "06_posterior_summary.R",
  "07_posterior_prediction.R",
  "08_publication_results.R"
)

for (s in steps) {
  cat("\n\n================ RUNNING", s, "================\n\n")
  source(s)
}

cat("\nPipeline completed successfully.\n")
