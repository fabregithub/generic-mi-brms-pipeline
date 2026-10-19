#!/usr/bin/env bash
# =============================================================================
# Phase 2 censoring sweep -- driver.  RESUME.md step 1.
# -----------------------------------------------------------------------------
# Maps how mixture-HBEB benchmark bias grows as detection degrades
# (limit_scale 1.0 -> 3.0), answering the open framing decision in RESUME.md:
# the August trap demo showed substitution only -3.9% vs oracle at realistic
# detection, so the honest message may be CONDITIONAL -- "robust when the
# benchmark-driving exposures are well detected; the trap grows as detection
# degrades" -- rather than a universal claim.
#
# RUN IT FROM ANYWHERE, including the repo root:
#   bash ~/dev/r-projects/generic-mi-brms-pipeline/validation/phase2/run_trap_sweep.sh
#
# It resolves its own absolute path before changing directory, so the working
# directory never matters. (A driver that re-invokes itself through a relative
# "$0" dies silently after the cd while still printing "detached" -- that bug
# cost a 21-hour run on 2026-10-15.)
#
# Config -- env vars, passed straight through to R/run_trap_sweep.R:
#   LIMITS (1.0,1.6,2.2,3.0)  N (800)   ITER (2500)  M (5)
#   REPS (3)  NKNOTS (50)  SWEEPS (3)   ERF (mixture)  SEED (1)  REL (0.25)
#   FG=1  run in the foreground instead of detaching
#
# COST -- the R script's own estimate, which is self-consistent:
#   per limit_scale = (1 + M + M) fits x REPS, plus REPS for the one-off oracle
#   default = (11 x 3 x 4) + 3 = 135 fits at ~1 bkmr fit/min = ~2.25 h
#   faster: fewer LIMITS, M=3, or N=600.  NOT independently re-timed here.
#
# ⚠️ It is SEQUENTIAL -- ~1 fit/min on one core, so 21 of 22 cores sit idle.
#    Parallelising REPS would cut it a lot, but that means editing a validated
#    script, so it is left alone. Worth doing deliberately, not in passing.
#
# RUNS AND MANUSCRIPT ARE SEPARATE (2026-10-19). Everything executable lives
# here under ~/dev; the Dropbox folder holds manuscript text, refs and the plan
# (NOTES.md / RESUME.md) and nothing that runs. Simulating inside a synced
# folder had already produced two failures in one morning: reads blocking on
# online-only placeholders, and a smoke test silently overwriting the 14-Aug
# n=200 results. ./pipeline-snapshot-2026-08-14 is frozen provenance only --
# nothing sources it; phase2 is self-contained.
#
# Watch:   tail -f <the log path printed below> | grep -v '^Iteration:'
# Alive?   kill -0 $(cat logs/trap_sweep_<stamp>.pid)
# Stop:    kill $(cat logs/trap_sweep_<stamp>.pid)
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
PROJ="$(dirname "$SELF")"
P2="$PROJ"
cd "$P2" || { echo "phase2 not found at $P2" >&2; exit 1; }
mkdir -p "$PROJ/logs" "$PROJ/results"

STAMP="$(date '+%Y%m%d-%H%M%S')"
LOG="$PROJ/logs/trap_sweep_${STAMP}.log"
PID="$PROJ/logs/trap_sweep_${STAMP}.pid"

export LIMITS="${LIMITS:-1.0,1.6,2.2,3.0}" N="${N:-800}" ITER="${ITER:-2500}" \
       M="${M:-5}" REPS="${REPS:-3}" NKNOTS="${NKNOTS:-50}" SWEEPS="${SWEEPS:-3}" \
       ERF="${ERF:-mixture}" SEED="${SEED:-1}" REL="${REL:-0.25}"

# GUARD: R/run_trap_sweep.R writes results/trap_sweep_<ERF>_n<N>.{rds,csv} and
# will silently overwrite an existing file of that name. On 2026-10-19 a smoke
# test at N=200 destroyed the 14-Aug n=200 sweep results exactly this way
# (recovered from Dropbox version history). Refuse rather than clobber.
OUT_RDS="$P2/results/trap_sweep_${ERF}_n${N}.rds"
if [[ -e "$OUT_RDS" && -z "${OVERWRITE:-}" ]]; then
  echo "REFUSING TO RUN: $OUT_RDS already exists ($(date -r "$OUT_RDS" '+%F %H:%M'))." >&2
  echo "  It would be overwritten. Either:" >&2
  echo "    - move it aside, or" >&2
  echo "    - run with a different N, or" >&2
  echo "    - re-run with OVERWRITE=1 if you genuinely want it replaced." >&2
  exit 2
fi

{ echo "[$(date '+%F %H:%M:%S')] phase2 censoring sweep"
  echo "LIMITS=$LIMITS N=$N ITER=$ITER M=$M REPS=$REPS NKNOTS=$NKNOTS ERF=$ERF SEED=$SEED REL=$REL"
  echo "cwd=$P2"
  echo "expected fits = (1 + $M + $M) x $REPS x $(echo "$LIMITS" | tr ',' '\n' | grep -c .) + $REPS"
  echo "---"; } > "$LOG"

if [[ -n "${FG:-}" ]]; then
  Rscript R/run_trap_sweep.R >> "$LOG" 2>&1; rc=$?
  echo "[$(date '+%F %H:%M:%S')] finished rc=$rc" >> "$LOG"; exit $rc
fi

nohup Rscript R/run_trap_sweep.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "$PID"
command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$RPID" >/dev/null 2>&1 &

echo "phase2 censoring sweep detached (pid $RPID)"
echo "  log:     $LOG"
echo "  pid:     $PID"
echo "  watch:   tail -f '$LOG' | grep -v '^Iteration:'"
echo "  alive?:  kill -0 \$(cat '$PID') && echo ALIVE || echo DONE"
echo "  results: $P2/results/trap_sweep_${ERF}_n${N}.{rds,csv}"
