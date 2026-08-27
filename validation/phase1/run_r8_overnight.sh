#!/usr/bin/env bash
# =============================================================================
# R8 overnight chain — three phases, run sequentially, one command.
# -----------------------------------------------------------------------------
# Phase 1 answers R8. Phases 2 and 3 are pure upside: each is worth running
# whichever way Phase 1 lands, so nothing here depends on a result we do not
# have yet.
#
#   P1  decision grid            ~7 h   Does BART close R8? 7 scenarios x 4 arms
#                                       x 300 reps at the tuning already measured.
#   P2  BART tuning check        ~4 h   Is the Phase-1 verdict fragile to BART's
#                                       tree count? Same non-linear cells at
#                                       ntree=200 (the dbarts/literature default)
#                                       against ntree=50. Matters BOTH ways: if
#                                       BART wins at 50 we need to know it also
#                                       wins at the standard setting; if it loses
#                                       at 50 it may win at 200, and we learn that
#                                       tonight instead of losing another day.
#   P3  pilot m / FMI            ~3 h   Closes requirement R12 (design-plan
#                                       Phase 5): where does interval width
#                                       stabilise as m grows? m in {10, 20, 50}
#                                       against the m=30 already on record.
#
# Total ~14 h with margin inside an 18 h window. Each phase writes its own
# tagged results and checkpoint, so nothing overwrites anything and a phase that
# dies costs only itself — the later phases still run.
#
#   ./run_r8_overnight.sh              # all three phases
#   PHASES=1 ./run_r8_overnight.sh     # phase 1 only
#   PHASES=1,2 ./run_r8_overnight.sh   # skip the pilot-m phase
#
# Watch:  tail -f validation/phase1/logs/r8chain_<stamp>.log
# Results: results/r8p1_latest.rds, r8p2_latest.rds, r8p3m{10,20,50}_latest.rds
#
# If a phase dies, do NOT resume — resume is not trusted (see
# run_v4_variance.R). Read its checkpoint instead:
#   summarise_v4(readRDS("results/r8p1_checkpoint.rds"))
# =============================================================================
set -uo pipefail          # deliberately NOT -e: one phase failing must not
                          # abort the phases after it

SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

# --- self-detach ------------------------------------------------------------
# Each phase runs in the FOREGROUND so the phases are strictly ordered, which
# means the chain as a whole blocks its terminal. For a 14 h run that is a
# liability, not an inconvenience: closing the terminal would SIGHUP the whole
# chain. So on first invocation the script re-execs itself under nohup, prints
# the pid, and returns the terminal. Set R8_NO_DETACH=1 to run in the
# foreground deliberately (useful when debugging the chain itself).
if [[ "${R8_DETACHED:-0}" != "1" && "${R8_NO_DETACH:-0}" != "1" ]]; then
  mkdir -p logs
  export R8_CHAIN_STAMP="$(date +%Y%m%d-%H%M%S)"
  CHAIN_LOG="logs/r8chain_${R8_CHAIN_STAMP}.log"
  R8_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  CHAIN_PID=$!
  echo "${CHAIN_PID}" > "logs/r8chain_${R8_CHAIN_STAMP}.pid"

  # Idle-sleep protection as a sibling, never wrapping the chain (a caffeinate
  # wrapper was implicated in earlier run deaths).
  if command -v caffeinate >/dev/null 2>&1; then
    nohup caffeinate -i -w "${CHAIN_PID}" >/dev/null 2>&1 &
  fi

  echo "R8 chain detached.  PID ${CHAIN_PID}"
  echo "Log:   validation/phase1/${CHAIN_LOG}"
  echo "Watch: tail -f validation/phase1/${CHAIN_LOG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/r8chain_${R8_CHAIN_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260825}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"

PHASES="${PHASES:-1,2,3}"
STAMP="${R8_CHAIN_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/r8chain_${STAMP}.log"
mkdir -p results logs

want() { [[ ",${PHASES}," == *",$1,"* ]]; }

say() { echo "[$(date '+%H:%M:%S')] $*" >> "$LOG"; }

# One phase = one Rscript invocation, foreground, so phases are strictly ordered.
run_phase() {
  local tag="$1"; shift
  say "=== phase ${tag} starting: $* ==="
  ( export OUT_TAG="${tag}" "$@"
    Rscript run_v4_variance.R ) >> "$LOG" 2>&1
  local rc=$?
  if [[ $rc -eq 0 ]]; then
    say "=== phase ${tag} done (results/${tag}_latest.rds) ==="
  else
    say "=== phase ${tag} FAILED rc=${rc}; checkpoint results/${tag}_checkpoint.rds may hold partial results; continuing ==="
  fi
}

say "R8 chain starting. phases=${PHASES} ncores=${NCORES} log=${LOG}"

# ---- P1: the decision grid -------------------------------------------------
if want 1; then
  run_phase r8p1 \
    SCENARIOS=base,mcar_z40,missing_y20,combined,nl_mcar_z40,nl_missing_y20,nl_combined \
    ARMS=pipeline_block_fcs,pipeline_properBoot,pipeline_micePmm,pipeline_bartMI \
    N_REP=300 M=30 BART_NTREE=50
fi

# ---- P2: BART tuning sensitivity -------------------------------------------
# Non-linear cells only, and only the two arms still in contention, so 4x the
# per-fit cost of ntree=200 stays affordable.
if want 2; then
  run_phase r8p2 \
    SCENARIOS=nl_mcar_z40,nl_missing_y20,nl_combined \
    ARMS=pipeline_properBoot,pipeline_bartMI \
    N_REP=300 M=30 BART_NTREE=200
fi

# ---- P3: pilot m / FMI (requirement R12) -----------------------------------
# m = 30 is already on record from P1, so this adds the points either side of it.
# properBoot only: this asks how width behaves in m for the SHIPPED default, not
# which imputer wins.
if want 3; then
  for mm in 10 20 50; do
    run_phase "r8p3m${mm}" \
      SCENARIOS=combined,nl_combined \
      ARMS=pipeline_properBoot \
      N_REP=200 "M=${mm}"
  done
fi

say "R8 chain finished. Results:"
ls -1 results/r8p*_latest.rds >> "$LOG" 2>/dev/null
