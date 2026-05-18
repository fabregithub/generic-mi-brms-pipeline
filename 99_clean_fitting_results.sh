#!/bin/bash

rm -f fits/fit_imp_*.rds
rm -f objects/fit_manifest.rds
rm -f objects/fit_status.rds
rm -f objects/fit_smoke_status.rds
rm -f results/fit_status.csv
rm -f results/fit_smoke_status.csv
rm -f results/worker_logs/fit_worker_imp_*.log

rm -f results/parameter_draws.rds
rm -f results/parameter_summary.rds
rm -f results/parameter_summary.csv
rm -f results/parameter_draws_imp_*.rds
rm -f objects/parameter_manifest.rds

rm -f results/missing_y_draws.rds
rm -f results/missing_y_summary.rds
rm -f results/missing_y_summary.csv
rm -f results/missing_y_draws_imp_*.rds

rm -f pipeline_error.flag
rm -f pipeline_success.flag
rm -f run_all_stdout.log

