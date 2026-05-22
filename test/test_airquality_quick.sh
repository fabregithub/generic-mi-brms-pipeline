#!/usr/bin/env bash
# Quick test for the airquality Gaussian example.
set -Eeuo pipefail
cd "$(dirname "$0")"

source ./test_example_common.sh
trap 'die "Airquality quick test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_airquality quick
