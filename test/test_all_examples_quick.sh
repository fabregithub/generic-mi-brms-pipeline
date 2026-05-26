#!/usr/bin/env bash
# Quick tests for both standard examples.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "All quick example tests failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_airquality quick
test_birthwt_logistic quick

list_test_runs
log "All quick example tests completed successfully"
