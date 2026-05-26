#!/usr/bin/env bash
# Quick test for the birthwt spline + monotonic brms example.
# This exercises custom brms formula support for s() and mo().
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "Birthwt spline/monotonic quick test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_birthwt_spline_monotonic quick
