#!/usr/bin/env bash
# ============================================================
# launch.sh — interactive launcher for generic-mi-brms-pipeline
# ============================================================
# Run from the project root:
#   bash launch.sh
# Or make executable and double-click (Mac):
#   chmod +x launch.sh
# ============================================================
set -euo pipefail

# Ensure we are in the project root
cd "$(dirname "${BASH_SOURCE[0]}")"

if [[ ! -f "run_all.R" ]]; then
  echo "ERROR: run_all.R not found. Run launch.sh from the project root folder."
  exit 1
fi

require_rscript() {
  if ! command -v Rscript &>/dev/null; then
    echo "ERROR: Rscript not found. Please install R and ensure it is on your PATH."
    exit 1
  fi
}

run_step() {
  require_rscript
  echo ""
  echo "======== Running $1 ========"
  Rscript "$1"
}

clean_all() {
  read -rp "This will delete objects/, fits/, results/. Type YES to confirm: " confirm
  if [[ "$confirm" == "YES" ]]; then
    rm -rf objects fits results
    rm -f pipeline_error.flag pipeline_success.flag \
          pipeline_progress.log pipeline_heartbeat.txt run_all_stdout.log
    echo "Done."
  else
    echo "Cancelled."
  fi
}

clean_fits() {
  read -rp "Delete fits and posteriors (keeps imputed data)? Type YES to confirm: " confirm
  if [[ "$confirm" == "YES" ]]; then
    rm -f fits/fit_imp_*.rds
    rm -f objects/fit_manifest.rds objects/fit_status.rds objects/fit_smoke_status.rds
    rm -f results/fit_status.csv results/fit_smoke_status.csv
    rm -f results/parameter_draws.rds results/parameter_summary.rds results/parameter_summary.csv
    rm -f results/parameter_draws_imp_*.rds objects/parameter_manifest.rds
    rm -f results/missing_y_draws.rds results/missing_y_summary.rds results/missing_y_summary.csv
    rm -f pipeline_error.flag pipeline_success.flag run_all_stdout.log
    echo "Done."
  else
    echo "Cancelled."
  fi
}

while true; do
  echo ""
  echo "========================================================"
  echo "  generic-mi-brms-pipeline — interactive launcher"
  echo "========================================================"
  echo "  1.  Run full pipeline          (run_all.R)"
  echo "  2.  Validate config only       (01_validate_config.R)"
  echo "  3.  Prepare data               (02_prepare_data.R)"
  echo "  4.  Impute missing data        (03_impute.R)"
  echo "  5.  Fit models                 (04_fit_models.R)"
  echo "  6.  Posterior summary          (06_posterior_summary.R)"
  echo "  7.  Diagnostics                (05_diagnostics.R)"
  echo "  8.  Posterior prediction       (07_posterior_prediction.R)"
  echo "  9.  Publication results        (08_publication_results.R)"
  echo " 10.  Clean ALL outputs          (objects/, fits/, results/)"
  echo " 11.  Clean fits/posteriors only (keeps imputed data)"
  echo "  q.  Quit"
  echo "========================================================"
  read -rp "Enter choice: " choice

  case "$choice" in
    1)  run_step "run_all.R" ;;
    2)  run_step "01_validate_config.R" ;;
    3)  run_step "02_prepare_data.R" ;;
    4)  run_step "03_impute.R" ;;
    5)  run_step "04_fit_models.R" ;;
    6)  run_step "06_posterior_summary.R" ;;
    7)  run_step "05_diagnostics.R" ;;
    8)  run_step "07_posterior_prediction.R" ;;
    9)  run_step "08_publication_results.R" ;;
    10) clean_all ;;
    11) clean_fits ;;
    q|Q) echo "Goodbye."; exit 0 ;;
    *) echo "Unrecognised choice. Please enter a number from the menu or q to quit." ;;
  esac
done
