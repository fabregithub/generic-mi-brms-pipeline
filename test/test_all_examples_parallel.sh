#!/usr/bin/env bash
# Modest parallel tests for all bundled examples, including s()/mo().
# These exercise parallel miceRanger imputation with impute_workers = 2.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "All parallel example tests failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_airquality parallel
test_birthwt_logistic parallel
test_birthwt_spline_monotonic parallel
test_lung_cox parallel

list_test_runs
log "All parallel example tests completed successfully"
