# ============================================================
# launch.R — interactive launcher for generic-mi-brms-pipeline
# ============================================================
# Open this file in RStudio and click Source (or press Ctrl+Shift+S).
# No terminal or command line needed.
# ============================================================

# Ensure we are in the project root (the folder containing this file)
if (file.exists("launch.R")) {
  # Already in the right place
} else if (file.exists(file.path(dirname(sys.frame(1)$ofile), "launch.R"))) {
  setwd(dirname(sys.frame(1)$ofile))
} else {
  stop(
    "Please open launch.R from the project root folder ",
    "(the folder that contains run_all.R and 00_config.R)."
  )
}

.run_clean <- function(script) {
  # Source in a fresh environment so successive menu choices don't bleed state
  e <- new.env(parent = globalenv())
  source(script, local = e)
  invisible(NULL)
}

.clean_outputs <- function(what) {
  targets <- switch(what,
    all = list(
      dirs  = c("objects", "fits", "results"),
      files = c("pipeline_error.flag", "pipeline_success.flag",
                "pipeline_progress.log", "pipeline_heartbeat.txt",
                "run_all_stdout.log")
    ),
    fits = list(
      dirs  = character(0),
      files = c(
        list.files("fits",    pattern = "fit_imp_.*\\.rds$", full.names = TRUE),
        list.files("objects", pattern = "fit_manifest|fit_status|fit_smoke",
                   full.names = TRUE),
        list.files("results", pattern = "fit_status|parameter_draws|parameter_summary|missing_y",
                   full.names = TRUE),
        "pipeline_error.flag", "pipeline_success.flag", "run_all_stdout.log"
      )
    )
  )
  for (d in targets$dirs)   if (dir.exists(d))   unlink(d, recursive = TRUE)
  for (f in targets$files)  if (file.exists(f))  file.remove(f)
  cat("Done.\n")
}

repeat {
  cat("
========================================================
  generic-mi-brms-pipeline — interactive launcher
========================================================
  1. Run full pipeline          (run_all.R)
  2. Validate config only       (01_validate_config.R)
  3. Prepare data               (02_prepare_data.R)
  4. Impute missing data        (03_impute.R)
  5. Fit models                 (04_fit_models.R)
  6. Posterior summary          (06_posterior_summary.R)
  7. Diagnostics                (05_diagnostics.R)
  8. Posterior prediction       (07_posterior_prediction.R)
  9. Publication results        (08_publication_results.R)
 10. Clean ALL outputs          (objects/, fits/, results/)
 11. Clean fits/posteriors only (keeps imputed data)
  q. Quit
========================================================
")

  choice <- trimws(readline("Enter choice: "))

  if (choice == "q" || choice == "Q") {
    cat("Goodbye.\n")
    break
  }

  if (choice == "1")  { cat("Running full pipeline...\n");       .run_clean("run_all.R") }
  else if (choice == "2")  { cat("Validating config...\n");            .run_clean("01_validate_config.R") }
  else if (choice == "3")  { cat("Preparing data...\n");               .run_clean("02_prepare_data.R") }
  else if (choice == "4")  { cat("Imputing...\n");                     .run_clean("03_impute.R") }
  else if (choice == "5")  { cat("Fitting models...\n");               .run_clean("04_fit_models.R") }
  else if (choice == "6")  { cat("Posterior summary...\n");            .run_clean("06_posterior_summary.R") }
  else if (choice == "7")  { cat("Running diagnostics...\n");          .run_clean("05_diagnostics.R") }
  else if (choice == "8")  { cat("Posterior prediction...\n");         .run_clean("07_posterior_prediction.R") }
  else if (choice == "9")  { cat("Publication results...\n");          .run_clean("08_publication_results.R") }
  else if (choice == "10") {
    confirm <- trimws(readline("This will delete objects/, fits/, results/. Type YES to confirm: "))
    if (confirm == "YES") { cat("Cleaning all outputs...\n"); .clean_outputs("all") }
    else cat("Cancelled.\n")
  }
  else if (choice == "11") {
    confirm <- trimws(readline("This will delete fits and posteriors (keeps imputed data). Type YES to confirm: "))
    if (confirm == "YES") { cat("Cleaning fits...\n"); .clean_outputs("fits") }
    else cat("Cancelled.\n")
  }
  else {
    cat("Unrecognised choice. Please enter a number from the menu or q to quit.\n")
  }
}
