#!/usr/bin/env bash
# Parallel test for the lung Cox proportional hazards example.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "Lung Cox parallel test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_lung_cox parallel
