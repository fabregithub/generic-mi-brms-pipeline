#!/usr/bin/env bash
# =============================================================================
# Launch the Phase-1 study DETACHED (walk away / close the terminal).
# Confirms the ERF bias of the no-Y leftcens pre-step and (if brms added) the
# joint-model fix. See README.md and ../PLAN_leftcensored_exposure_integration.md.
#
#   ./run_phase1.sh                                  # full no-brms grid, NCORES=14 (~5-8 min)
#   NCORES=8 ./run_phase1.sh                          # fewer forks
#   N_REP=150 ./run_phase1.sh                         # quicker look
#   CONFIG=quick ./run_phase1.sh                      # small/fast sanity grid
#   PROCS=oracle,complete_case,leftcens_prestep,brms_joint ./run_phase1.sh   # + Stan
#
# Watch:  tail -f validation/phase1/logs/phase1_<stamp>.log
# Stop:   kill $(cat validation/phase1/logs/phase1_<stamp>.pid)
# Done:   validation/phase1/results/latest.rds (+ phase1_summary.csv, phase1_raw.csv)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# Keep BLAS single-threaded so it doesn't fight the R-level loop.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1

# Study config (all overridable via environment). Defaults = full grid with the
# congenial cens_mi_y gold standard (fast, no Stan). brms_joint is deprecated in
# the scaffold (collapses to LOD-substitution) -- add it only for reference.
export CONFIG="${CONFIG:-full}"
export PROCS="${PROCS:-oracle,complete_case,leftcens_prestep,cens_mi_y_shash}"
export ERF="${ERF:-additive,mixture}"
export ND="${ND:-0.2,0.4}"
export N_REP="${N_REP:-300}"
export M="${M:-30}"
export N="${N:-800}"
export NCORES="${NCORES:-14}"   # fork workers for gsimp_mi (leave headroom on 16 cores)
export SEED="${SEED:-20260813}"

mkdir -p results logs
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="logs/phase1_${STAMP}.log"
PIDF="logs/phase1_${STAMP}.pid"

echo "Launching Phase 1: CONFIG=${CONFIG} PROCS=${PROCS} ERF=${ERF} ND=${ND} N_REP=${N_REP} M=${M} N=${N} NCORES=${NCORES}"
nohup Rscript run_phase1.R > "${LOG}" 2>&1 &
echo $! > "${PIDF}"
echo "PID $(cat "${PIDF}")  |  log: ${LOG}"
echo "Watch: tail -f validation/phase1/${LOG}"
echo "Stop:  kill \$(cat validation/phase1/${PIDF})"
