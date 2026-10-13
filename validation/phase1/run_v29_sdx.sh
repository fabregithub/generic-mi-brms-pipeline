#!/usr/bin/env bash
# =============================================================================
# V29 -- s^2 on its own, at a FIXED arrow.
# -----------------------------------------------------------------------------
# V28 left exactly one live hypothesis. At n = 12800 and delta >= 0.4 the
# chain's derived/measured ratio TRENDS with the arrow in the fork sweep
# (0.944, 0.992, 1.112) and is FLAT in the pipe sweep (0.942, 0.953, 0.952).
# The arrow moves the exposure conditional's Z1 coefficient in both roles but
# moves Var(logX1) only in the fork, so the flat pipe ACQUITS the coefficient
# and the trending fork implicates Var(logX1).
#
# `sd_x` moves s^2 in BOTH roles with the arrow fixed at 0.60 -- the one lever
# that tests that channel on its own. Cells and their Var(logX1):
#     sxf07/10/13  fork  0.85 / 1.36 / 2.05   (1.36 = canonical zr_fork)
#     sxp07/10/13  pipe  0.49 / 1.00 / 1.68
#
# REGISTERED PREDICTION (2026-10-13, before the run; see R/robustness.R):
#   P1 ratio rises with sd_x in BOTH roles  -> s^2 is the driver
#   P2 rises in the FORK only               -> something fork-specific survives
#   P3 flat in both                         -> s^2 is not the driver at all
# The pipe's s^2 range (0.49 -> 1.68) is WIDER than the fork arrow sweep's
# (1.16 -> 1.64), so a flat pipe cannot be excused as too small a lever.
#
# This re-opens a question THEORY.md records as refuted. That refutation used a
# two-block simulation which misses the pipeline's level by 31% -- a limit
# THEORY states itself. These cells use the real pipeline arms.
#
# Cost: measured full-wave pilot for the identical shape (6 cells x 200 reps,
# n = 12800, 22 workers) ran 52 min as v28b. Not extrapolated.
#
#   bash run_v29_sdx.sh        # detach; logs/v29_<stamp>.{log,pid}
#   bash run_status.sh v29     # is it alive?
#   Rscript analyze_v29.R      # after it lands
# =============================================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
TAG=v29
STAMP="$(date '+%Y%m%d-%H%M%S')"
mkdir -p logs results

if [[ -z "${V29_DETACHED:-}" ]]; then
  V29_DETACHED=1 V29_STAMP="$STAMP" nohup bash "$0" "$@" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V29 detached (pid $PID)"
  echo "  log:    logs/${TAG}_${STAMP}.log"
  echo "  status: bash run_status.sh ${TAG}"
  exit 0
fi

STAMP="${V29_STAMP:-$STAMP}"
LOG="logs/${TAG}_${STAMP}.log"
NC="${NCORES:-$(( $(sysctl -n hw.ncpu) - 2 ))}"

{ echo "[$(date '+%H:%M:%S')] V28b start -- arrow sweep at n=12800, ${NC} workers"
  echo "cells: fkd04,fkd06,fkd08,ppd04,ppd06,ppd08 (delta=0.2 cells EXCLUDED --"
  echo "  their measured effect is ~0.4 pp, so derived/measured is a ratio of"
  echo "  near-zero arms and carries no information. Same trap as V26 C1/C2.)"
} > "$LOG"

SCENARIOS=sxf07,sxf10,sxf13,sxp07,sxp10,sxp13 \
ARMS=pipeline_properZ,pipeline_properZ_ds \
N_REP=200 N_OBS=12800 OUT_TAG=${TAG}_n12800 NCORES="$NC" \
  Rscript run_v4_variance.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?

{ echo ""; echo "[$(date '+%H:%M:%S')] V29 finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v29.R"; } >> "$LOG"
exit $rc
