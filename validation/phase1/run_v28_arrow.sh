#!/usr/bin/env bash
# =============================================================================
# V28 -- the arrow sweep at a THIRD sample size.
# -----------------------------------------------------------------------------
# V28's first two waves (n=800, n=3200) were launched as a bare background
# command with no `.pid` file. run_status.sh reports only what has one, so a
# live 24-core run showed as "Nothing running" and was very nearly killed and
# restarted on top of itself. Every wave gets a runner from here on.
#
# WHY A THIRD n. The closed-form chain is an ASYMPTOTIC statement; the measured
# wrapper effect is not. Between n=800 and n=3200 the fork moved +4.55 -> +4.99
# (toward its derived +5.05) and the pipe +5.63 -> +5.92 (PAST its derived
# +5.44). The "fork residual" that THEORY records as a fork-specific 15% miss is
# therefore, at minimum, entangled with sample size -- at n=3200 the fork is the
# role that AGREES. One more n decides between:
#
#   (P1) both converge to the chain's asymptote  -> the chain is right and the
#        n=800 comparison was the error; the fork residual dissolves.
#   (P2) the pipe keeps rising past +5.44        -> the chain under-predicts the
#        pipe, and the role contrast is real but the OPPOSITE of what is on
#        record.
#
# These are distinguishable: P1 needs the pipe to turn around, P2 needs it to
# keep climbing. Registered before the run (2026-10-13).
#
# Cost: MEASURED full-wave pilot, 22 reps x 2 cells = 115 s on 22 workers
# (results/v28pilot_n12800_*). 6 cells x 200 reps ~ 52 min. Not extrapolated
# from a part-wave.
#
#   bash run_v28_arrow.sh            # detach, write logs/v28b_<stamp>.{log,pid}
#   bash run_status.sh v28b          # is it alive?
#   Rscript analyze_v28.R            # after it lands
# =============================================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
TAG=v28b
STAMP="$(date '+%Y%m%d-%H%M%S')"
mkdir -p logs results

if [[ -z "${V28_DETACHED:-}" ]]; then
  V28_DETACHED=1 V28_STAMP="$STAMP" nohup bash "$0" "$@" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V28b detached (pid $PID)"
  echo "  log:    logs/${TAG}_${STAMP}.log"
  echo "  status: bash run_status.sh ${TAG}"
  exit 0
fi

STAMP="${V28_STAMP:-$STAMP}"
LOG="logs/${TAG}_${STAMP}.log"
NC="${NCORES:-$(( $(sysctl -n hw.ncpu) - 2 ))}"

{ echo "[$(date '+%H:%M:%S')] V28b start -- arrow sweep at n=12800, ${NC} workers"
  echo "cells: fkd04,fkd06,fkd08,ppd04,ppd06,ppd08 (delta=0.2 cells EXCLUDED --"
  echo "  their measured effect is ~0.4 pp, so derived/measured is a ratio of"
  echo "  near-zero arms and carries no information. Same trap as V26 C1/C2.)"
} > "$LOG"

SCENARIOS=fkd04,fkd06,fkd08,ppd04,ppd06,ppd08 \
ARMS=pipeline_properZ,pipeline_properZ_ds \
N_REP=200 N_OBS=12800 OUT_TAG=${TAG}_n12800 NCORES="$NC" \
  Rscript run_v4_variance.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?

{ echo ""; echo "[$(date '+%H:%M:%S')] V28b finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v28.R"; } >> "$LOG"
exit $rc
