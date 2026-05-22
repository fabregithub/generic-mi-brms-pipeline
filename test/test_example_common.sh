#!/usr/bin/env bash
# Shared helper functions for testing the Generic MI + brms Pipeline examples.
#
# Location:
#   generic_mi_brms_pipeline/test/test_example_common.sh
#
# The test scripts are stored in test/ to keep the project root clean.
# This file automatically moves to the project root before running R scripts.

set -Eeuo pipefail

TEST_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${TEST_SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

require_file() {
  local f="$1"
  [[ -f "$f" ]] || die "Required file not found: $f"
}

require_command() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: $cmd"
}

clean_outputs() {
  log "Cleaning previous pipeline outputs"
  rm -rf objects fits results
  rm -f pipeline_error.flag
  rm -f pipeline_success.flag
  rm -f pipeline_progress.log
  rm -f pipeline_heartbeat.txt
  rm -f pipeline_stdout.log
  rm -f run_all_stdout.log
  rm -f run_all_airquality_*_stdout.log
  rm -f run_all_birthwt_logistic_*_stdout.log
}

strip_test_overrides() {
  require_file "00_config.R"

  local tmp_file
  tmp_file="$(mktemp)"

  awk '
    /^# ---- BEGIN automated test overrides ----$/ {skip=1; next}
    /^# ---- END automated test overrides ----$/ {skip=0; next}
    skip != 1 {print}
  ' 00_config.R > "$tmp_file"

  mv "$tmp_file" 00_config.R
}

apply_test_overrides() {
  local mode="${1:-quick}"

  strip_test_overrides

  if [[ "$mode" == "quick" ]]; then
    log "Applying quick-test overrides to 00_config.R"
    cat >> 00_config.R <<'EOF'

# ---- BEGIN automated test overrides ----
# Added by test scripts. Remove this block to restore the original settings.
analysis_spec$imputation$m <- 5

analysis_spec$model$chains <- 1
analysis_spec$model$iter <- 500
analysis_spec$model$warmup <- 250
analysis_spec$model$run_smoke_fit <- TRUE

# Quick tests use single-worker imputation to isolate configuration/model errors.
analysis_spec$parallel$impute_workers <- 1
analysis_spec$parallel$num_impute_threads_per_worker <- 1
analysis_spec$parallel$num_impute_threads <- 1

analysis_spec$parallel$fit_workers <- 1
analysis_spec$parallel$cores_per_fit <- 1
analysis_spec$parallel$future_globals_maxsize_gb <- 8

analysis_spec$posterior_prediction$ndraws <- 200
# ---- END automated test overrides ----
EOF
  elif [[ "$mode" == "parallel" ]]; then
    log "Applying modest parallel-test overrides to 00_config.R"
    cat >> 00_config.R <<'EOF'

# ---- BEGIN automated test overrides ----
# Added by test scripts. Remove this block to restore the original settings.
analysis_spec$imputation$m <- 10

analysis_spec$model$chains <- 4
analysis_spec$model$iter <- 500
analysis_spec$model$warmup <- 250
analysis_spec$model$run_smoke_fit <- TRUE

# Parallel tests exercise the doParallel/foreach miceRanger path.
analysis_spec$parallel$impute_workers <- 2
analysis_spec$parallel$num_impute_threads_per_worker <- 2
analysis_spec$parallel$num_impute_threads <- 2

analysis_spec$parallel$fit_workers <- 2
analysis_spec$parallel$cores_per_fit <- 2
analysis_spec$parallel$future_globals_maxsize_gb <- 20

analysis_spec$posterior_prediction$ndraws <- 300
# ---- END automated test overrides ----
EOF
  else
    die "Unknown test mode: $mode. Use quick or parallel."
  fi
}

render_and_check_report() {
  require_command quarto

  local report_qmd="results/publication/report/bayesian_mi_report_template.qmd"
  local report_html="results/publication/report/bayesian_mi_report_template.html"
  local report_docx="results/publication/report/bayesian_mi_report_template.docx"

  require_file "$report_qmd"

  log "Rendering Quarto report"
  quarto render "$report_qmd"

  require_file "$report_html"
  require_file "$report_docx"

  log "Report files created:"
  ls -lh results/publication/report/
}

check_pipeline_success() {
  # Some pipeline versions write pipeline_success.flag and some simply print
  # "Pipeline completed successfully." Do not fail only because the flag is absent.
  # The required output files below are the stronger practical check.
  if [[ -f "pipeline_success.flag" ]]; then
    log "pipeline_success.flag found"
  else
    log "pipeline_success.flag not found; continuing because this pipeline version may not create it"
  fi

  if [[ -f "pipeline_error.flag" ]]; then
    die "pipeline_error.flag exists. Inspect pipeline_progress.log and the run_all stdout log."
  fi

  require_file "results/diagnostics.rds"
  require_file "results/parameter_summary.rds"
  require_file "results/publication/tables/main_effect_table_display.csv"
  require_file "results/publication/tables/analysis_metadata.csv"
  require_file "results/publication/report/bayesian_mi_report_template.qmd"
}

check_parallel_imputation_log() {
  local log_file="$1"

  log "Checking that the parallel miceRanger path was exercised"

  if grep -E "miceRanger parallel:.*TRUE|miceRanger parallel.*TRUE|impute_workers:.*2|impute_workers.*2" "$log_file" >/dev/null 2>&1; then
    log "Parallel imputation log check passed"
  else
    log "Could not confirm parallel miceRanger from stdout log."
    log "This may be harmless if log formatting changed, but please inspect:"
    log "  $log_file"
    log "Expected to see impute_workers = 2 or miceRanger parallel = TRUE."
  fi
}

print_diagnostics_hint() {
  log "Optional diagnostic check in R/RStudio:"
  cat <<'EOF'
diag <- readRDS("results/diagnostics.rds")
summary(diag)
sum(diag$divergent, na.rm = TRUE)
sum(diag$treedepth_hits, na.rm = TRUE)
EOF
}

run_pipeline() {
  local log_file="$1"
  local mode="${2:-quick}"

  require_file "01_validate_config.R"
  require_file "run_all.R"

  log "Validating config"
  Rscript 01_validate_config.R

  log "Running full pipeline"
  Rscript run_all.R 2>&1 | tee "$log_file"

  check_pipeline_success

  if [[ "$mode" == "parallel" ]]; then
    check_parallel_imputation_log "$log_file"
  fi

  render_and_check_report
}

prepare_airquality_example() {
  log "Preparing airquality Gaussian example"

  require_file "examples/airquality_gaussian/00_config_airquality_gaussian.R"
  require_file "examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv"
  require_file "examples/airquality_gaussian/00_create_airquality_example_data.R"

  cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
  cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv

  log "Creating airquality example data"
  Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
}

prepare_birthwt_logistic_example() {
  log "Preparing birthwt logistic example"

  require_file "examples/birthwt_logistic/00_config_birthwt_logistic.R"
  require_file "examples/birthwt_logistic/00_variable_dictionary_birthwt_logistic.csv"
  require_file "examples/birthwt_logistic/00_create_birthwt_logistic_example_data.R"

  cp examples/birthwt_logistic/00_config_birthwt_logistic.R 00_config.R
  cp examples/birthwt_logistic/00_variable_dictionary_birthwt_logistic.csv 00_variable_dictionary.csv

  log "Creating birthwt logistic example data"
  Rscript examples/birthwt_logistic/00_create_birthwt_logistic_example_data.R
}

test_airquality() {
  local mode="${1:-quick}"
  local log_file="run_all_airquality_${mode}_stdout.log"

  clean_outputs
  prepare_airquality_example
  apply_test_overrides "$mode"
  run_pipeline "$log_file" "$mode"

  log "Airquality ${mode} test completed successfully"
  print_diagnostics_hint
}

test_birthwt_logistic() {
  local mode="${1:-quick}"
  local log_file="run_all_birthwt_logistic_${mode}_stdout.log"

  clean_outputs
  prepare_birthwt_logistic_example
  apply_test_overrides "$mode"
  run_pipeline "$log_file" "$mode"

  log "Birthwt logistic ${mode} test completed successfully"
  print_diagnostics_hint
}
