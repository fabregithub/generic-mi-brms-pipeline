#!/usr/bin/env bash
# Shared helper functions for testing the Generic MI + brms Pipeline examples.
#
# Location:
#   generic_mi_brms_pipeline/test/test_example_common.sh
#
# The test scripts live in test/ to keep the project root clean.
#
# Important behaviour:
#   - Each example test runs in its own isolated project copy under test/runs/.
#   - Root-level objects/, fits/, and results/ are not deleted or overwritten.
#   - Previous test runs are kept for inspection.
#
# Example:
#   bash test/test_all_examples_quick.sh
#
# Outputs:
#   test/runs/<timestamp>_<example>_<mode>/

set -Eeuo pipefail

TEST_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${TEST_SCRIPT_DIR}/.." && pwd)"
TEST_RUNS_DIR="${PROJECT_ROOT}/test/runs"

log() {
  # Important: log to stderr.
  # Some helper functions return paths through stdout and are called through
  # command substitution. Logging to stdout would pollute those returned paths.
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

require_file() {
  local f="$1"
  [[ -f "$f" ]] || die "Required file not found: $f"
}

require_dir() {
  local d="$1"
  [[ -d "$d" ]] || die "Required directory not found: $d"
}

require_command() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: $cmd"
}

timestamp_id() {
  date '+%Y%m%d_%H%M%S'
}

make_run_dir() {
  local example="$1"
  local mode="$2"
  local ts
  ts="$(timestamp_id)"

  local run_dir="${TEST_RUNS_DIR}/${ts}_${example}_${mode}"

  # Avoid rare collision when two tests start in the same second.
  if [[ -e "$run_dir" ]]; then
    run_dir="${run_dir}_$$"
  fi

  mkdir -p "$run_dir"
  printf '%s\n' "$run_dir"
}

copy_if_exists() {
  local src="$1"
  local dest="$2"

  if [[ -e "$src" ]]; then
    cp -R "$src" "$dest"
  fi
}

prepare_runtime_project() {
  local run_dir="$1"
  local runtime_project="${run_dir}/project"
  # This function is called via command substitution, so stdout must contain
  # only the runtime project path. Use log(), which writes to stderr, for
  # human-readable messages.

  mkdir -p "$runtime_project"

  log "Creating isolated runtime project:"
  log "  $runtime_project"

  # Root R scripts
  local root_files=(
    "00_common_functions.R"
    "00_config.R"
    "00_create_airquality_example_data.R"
    "00_variable_dictionary.csv"
    "01_validate_config.R"
    "02_prepare_data.R"
    "03_impute.R"
    "04_fit_models.R"
    "05_diagnostics.R"
    "06_posterior_summary.R"
    "07_posterior_prediction.R"
    "08_publication_results.R"
    "09_check_mo_parameter_columns.R"
    "10_publication_mo_results.R"
    "fit_single_imputation.R"
    "run_all.R"
    "99_clean_fitting_results.sh"
    "99_cleanall.sh"
  )

  for f in "${root_files[@]}"; do
    if [[ -f "${PROJECT_ROOT}/${f}" ]]; then
      cp "${PROJECT_ROOT}/${f}" "${runtime_project}/${f}"
    fi
  done

  # Example configs/data scripts/dictionaries.
  require_dir "${PROJECT_ROOT}/examples"
  cp -R "${PROJECT_ROOT}/examples" "${runtime_project}/examples"

  # Create empty data directory. Example data scripts will write into this.
  mkdir -p "${runtime_project}/data"

  # Optional docs/license files. These are useful when inspecting a run.
  copy_if_exists "${PROJECT_ROOT}/README.md" "${runtime_project}/"
  copy_if_exists "${PROJECT_ROOT}/LICENSE" "${runtime_project}/"

  printf '%s\n' "$runtime_project"
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
    log "Applying quick-test overrides to runtime 00_config.R"
    cat >> 00_config.R <<'EOF'

# ---- BEGIN automated test overrides ----
# Added by test scripts inside an isolated runtime project.
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
    log "Applying modest parallel-test overrides to runtime 00_config.R"
    cat >> 00_config.R <<'EOF'

# ---- BEGIN automated test overrides ----
# Added by test scripts inside an isolated runtime project.
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
    die "pipeline_error.flag exists. Inspect pipeline_progress.log and the run_all stdout log in the run folder."
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

check_mo_outputs() {
  log "Checking mo()/s() example outputs"

  # The spline/monotonic example should generate supplementary summaries for
  # special brms parameters such as simo_ and smooth terms.
  if [[ -f "results/special_parameter_summary.rds" || -f "results/special_parameter_summary.csv" ]]; then
    log "special_parameter_summary output found"
  else
    die "Expected special_parameter_summary output was not found for mo()/s() example."
  fi

  # Conditional effect figures are optional but expected for this example.
  if [[ -d "results/publication/figures" ]]; then
    log "Publication figures found:"
    ls -lh results/publication/figures || true
  else
    log "No results/publication/figures directory found; continuing because figure generation may be disabled by config."
  fi

  # Show extracted special columns for easier debugging.
  if [[ -f "results/parameter_draws.rds" ]]; then
    Rscript -e 'x <- readRDS("results/parameter_draws.rds"); cat("Special columns:\n"); print(grep("^(simo_|bsp_|sds_|bs_)", names(x), value = TRUE))'
  fi
}

print_diagnostics_hint() {
  log "Optional diagnostic check in R/RStudio from this run folder:"
  cat <<'EOF'
diag <- readRDS("results/diagnostics.rds")
summary(diag)
sum(diag$divergent, na.rm = TRUE)
sum(diag$treedepth_hits, na.rm = TRUE)
EOF
}

write_run_metadata() {
  local run_dir="$1"
  local example="$2"
  local mode="$3"
  local runtime_project="$4"

  cat > "${run_dir}/RUN_INFO.txt" <<EOF
Generic MI + brms pipeline test run

Example: ${example}
Mode: ${mode}
Created: $(date '+%Y-%m-%d %H:%M:%S')
Project root: ${PROJECT_ROOT}
Runtime project: ${runtime_project}

Inspect outputs in:
  ${runtime_project}/objects
  ${runtime_project}/fits
  ${runtime_project}/results

Stdout log:
  ${run_dir}/logs/run_all_${example}_${mode}_stdout.log
EOF
}

run_pipeline() {
  local log_file="$1"
  local mode="${2:-quick}"
  local example="${3:-unknown}"

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

  if [[ "$example" == "birthwt_spline_monotonic" ]]; then
    check_mo_outputs
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

prepare_birthwt_spline_monotonic_example() {
  log "Preparing birthwt spline + monotonic example"

  require_file "examples/birthwt_spline_monotonic/00_config_birthwt_spline_monotonic.R"
  require_file "examples/birthwt_spline_monotonic/00_variable_dictionary_birthwt_spline_monotonic.csv"
  require_file "examples/birthwt_spline_monotonic/00_create_birthwt_spline_monotonic_example_data.R"

  cp examples/birthwt_spline_monotonic/00_config_birthwt_spline_monotonic.R 00_config.R
  cp examples/birthwt_spline_monotonic/00_variable_dictionary_birthwt_spline_monotonic.csv 00_variable_dictionary.csv

  log "Creating birthwt spline + monotonic example data"
  Rscript examples/birthwt_spline_monotonic/00_create_birthwt_spline_monotonic_example_data.R
}

run_example_in_isolated_project() {
  local example="$1"
  local mode="$2"

  local run_dir
  local runtime_project
  local log_file

  run_dir="$(make_run_dir "$example" "$mode")"
  runtime_project="$(prepare_runtime_project "$run_dir")"

  write_run_metadata "$run_dir" "$example" "$mode" "$runtime_project"

  cd "$runtime_project"

  log "Running ${example} ${mode} test in isolated folder:"
  log "  $runtime_project"

  if [[ "$example" == "airquality" ]]; then
    prepare_airquality_example
  elif [[ "$example" == "birthwt_logistic" ]]; then
    prepare_birthwt_logistic_example
  elif [[ "$example" == "birthwt_spline_monotonic" ]]; then
    prepare_birthwt_spline_monotonic_example
  else
    die "Unknown example: $example"
  fi

  apply_test_overrides "$mode"

  mkdir -p "${run_dir}/logs"
  log_file="${run_dir}/logs/run_all_${example}_${mode}_stdout.log"
  run_pipeline "$log_file" "$mode" "$example"

  # If an older script or manual command created a stdout log in the runtime
  # project root, move it into the run-level logs folder.
  if [[ -f "run_all_${example}_${mode}_stdout.log" ]]; then
    mkdir -p "${run_dir}/logs"
    mv "run_all_${example}_${mode}_stdout.log" "${run_dir}/logs/"
  fi

  log "${example} ${mode} test completed successfully"
  log "Run outputs preserved in:"
  log "  $runtime_project"
  print_diagnostics_hint

  # Return to the real project root for the next test.
  cd "$PROJECT_ROOT"
}

test_airquality() {
  local mode="${1:-quick}"
  run_example_in_isolated_project "airquality" "$mode"
}

test_birthwt_logistic() {
  local mode="${1:-quick}"
  run_example_in_isolated_project "birthwt_logistic" "$mode"
}

test_birthwt_spline_monotonic() {
  local mode="${1:-quick}"
  run_example_in_isolated_project "birthwt_spline_monotonic" "$mode"
}

list_test_runs() {
  if [[ -d "$TEST_RUNS_DIR" ]]; then
    log "Existing test runs:"
    find "$TEST_RUNS_DIR" -maxdepth 2 -name RUN_INFO.txt -print | sort
  else
    log "No test runs found yet."
  fi
}
