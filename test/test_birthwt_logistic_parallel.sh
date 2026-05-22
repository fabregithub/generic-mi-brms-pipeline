#!/usr/bin/env bash
# Modest parallel test for the birthwt Bernoulli/logit example.
# This exercises parallel miceRanger imputation with impute_workers = 2.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/test_example_common.sh"

trap 'die "Birthwt logistic parallel test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_birthwt_logistic parallel
