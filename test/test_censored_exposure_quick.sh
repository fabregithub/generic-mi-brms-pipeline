#!/usr/bin/env bash
# Quick test for the censored-exposure block-FCS example.
# Two censored exposures (expo1 left-censored; expo2 three-tier ND + DNQ interval),
# imputed with the Y-aware leftcens X-block + miceRanger Z-block, then fit with brms.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "Censored-exposure quick test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_censored_exposure quick
