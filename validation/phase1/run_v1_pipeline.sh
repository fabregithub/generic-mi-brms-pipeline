#!/usr/bin/env bash
# =============================================================================
# Track V1 -- validate the SHIPPED censored-exposure engine against known truth.
# Launches DETACHED (walk away / close the terminal).
#
# Adds procedure 5 (`pipeline_block_fcs`, calling the real
# run_censored_exposure_block_fcs()) alongside the existing four, so the shipped
# code is measured like-for-like against oracle / complete-case / the no-Y
# pre-step / the `cens_mi_y_shash` prototype.
#
# See ../PLAN_pipeline_validation.md §6 for the acceptance criteria.
#
#   ./run_v1_pipeline.sh                    # full grid, 300 reps  (LONG)
#   CONFIG=quick N_REP=20 ./run_v1_pipeline.sh   # smoke: check the adapter first
#   ERF=additive ./run_v1_pipeline.sh       # additive only (the pass/fail case)
#   SWEEPS=1 ./run_v1_pipeline.sh           # 1 outer sweep (see note below)
#
# Watch:  tail -f validation/phase1/logs/v1_<stamp>.log
# Stop:   kill $(cat validation/phase1/logs/v1_<stamp>.pid)
# Done:   validation/phase1/results/latest.rds (+ phase1_summary.csv)
#
# NOTE on SWEEPS. With the harness's default design the covariates and outcome
# are fully observed, so the miceRanger Z block is a no-op and each outer sweep
# is just another independent draw of the censored exposure -- only the last is
# kept. SWEEPS=3 (the shipped example's setting) is therefore ~3x the cost of
# SWEEPS=1 for a statistically identical result *in this configuration*. Keep 3
# for fidelity to the shipped default; drop to 1 when you want speed.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# Keep BLAS single-threaded so it doesn't fight the R-level loop.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1

# Where to find 00_censored_exposure.R / 00_common_functions.R.
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"

export CONFIG="${CONFIG:-full}"
export PROCS="${PROCS:-oracle,complete_case,leftcens_prestep,cens_mi_y_shash,pipeline_block_fcs}"
export ERF="${ERF:-additive,mixture}"
export ND="${ND:-0.2,0.4}"
export N_REP="${N_REP:-300}"
export M="${M:-30}"
export N="${N:-800}"
export NCORES="${NCORES:-20}"   # 24-core box; leaves headroom
export SEED="${SEED:-20260825}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"

if [[ ! -f "${PIPELINE_ROOT}/00_censored_exposure.R" ]]; then
  echo "ERROR: PIPELINE_ROOT='${PIPELINE_ROOT}' has no 00_censored_exposure.R" >&2
  exit 1
fi

mkdir -p results logs
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="logs/v1_${STAMP}.log"
PIDF="logs/v1_${STAMP}.pid"

echo "Launching V1: CONFIG=${CONFIG} PROCS=${PROCS} ERF=${ERF} ND=${ND}"
echo "              N_REP=${N_REP} M=${M} N=${N} NCORES=${NCORES} SWEEPS=${SWEEPS} MARGIN=${MARGIN}"
echo "              PIPELINE_ROOT=${PIPELINE_ROOT}"
nohup Rscript run_phase1.R > "${LOG}" 2>&1 &
echo $! > "${PIDF}"
echo "PID $(cat "${PIDF}")  |  log: ${LOG}"
echo "Watch: tail -f validation/phase1/${LOG}"
echo "Stop:  kill \$(cat validation/phase1/${PIDF})"
