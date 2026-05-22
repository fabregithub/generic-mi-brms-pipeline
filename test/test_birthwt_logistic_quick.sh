#!/usr/bin/env bash
# Quick test for the birthwt Bernoulli/logit example.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "Birthwt logistic quick test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_birthwt_logistic quick
