#!/usr/bin/env bash
# Modest parallel tests for both standard examples.
# These exercise parallel miceRanger imputation with impute_workers = 2.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "All parallel example tests failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_airquality parallel
test_birthwt_logistic parallel

log "All parallel example tests completed successfully"
