#!/usr/bin/env bash
# Modest parallel test for the birthwt spline + monotonic brms example.
# This exercises parallel miceRanger plus custom brms s()/mo() support.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "Birthwt spline/monotonic parallel test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_birthwt_spline_monotonic parallel
