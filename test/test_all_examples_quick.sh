#!/usr/bin/env bash
# Quick tests for both standard examples.
set -Eeuo pipefail
cd "$(dirname "$0")"

source ./test_example_common.sh
trap 'die "All quick example tests failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_airquality quick
test_birthwt_logistic quick

log "All quick example tests completed successfully"
