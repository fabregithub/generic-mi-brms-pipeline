#!/usr/bin/env bash
# =============================================================================
# V16 alone -- now stage 1 of run_v16_v17.sh.
# -----------------------------------------------------------------------------
# This file used to carry its own scenario list. It no longer does, because V16
# gained the MAR cells (mar_z40, nl_mar_z40) after this script was written, and
# two drivers each claiming to "run V16" with different cell sets is exactly the
# drift that makes an old result impossible to interpret later. There is one
# definition of the track and it lives in run_v16_v17.sh.
#
# Registered predictions and falsification criteria: ../PLAN_pipeline_validation.md §8i.
# =============================================================================
SELF="$(cd "$(dirname "$0")" && pwd)"
echo "V16 is stage 1 of the sequence driver. Running it via:"
echo "  STAGES=1 $SELF/run_v16_v17.sh"
echo ""
exec env STAGES=1 bash "$SELF/run_v16_v17.sh" "$@"
