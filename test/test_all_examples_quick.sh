#!/usr/bin/env bash
# Quick tests for all bundled examples, including s()/mo().
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "All quick example tests failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_airquality quick
test_birthwt_logistic quick
test_birthwt_spline_monotonic quick
test_lung_cox quick

list_test_runs
log "All quick example tests completed successfully"
