#!/usr/bin/env bash
# Modest parallel test for the birthwt Bernoulli/logit example.
set -Eeuo pipefail
cd "$(dirname "$0")"

source ./test_example_common.sh
trap 'die "Birthwt logistic parallel test failed at line $LINENO"' ERR

require_command Rscript
require_command quarto

test_birthwt_logistic parallel
