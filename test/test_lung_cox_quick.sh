#!/usr/bin/env bash
# Quick test for the lung Cox proportional hazards example.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "Lung Cox quick test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_lung_cox quick
