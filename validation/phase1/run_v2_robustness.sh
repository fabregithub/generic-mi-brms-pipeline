#!/usr/bin/env bash
# =============================================================================
# Track V2 -- robustness sweep for the shipped censored-exposure engine.
# Launches DETACHED (walk away / close the terminal).
#
# Sweeps the axes V1 left untouched, including the two that were dead code
# there: MCAR covariates (turns the miceRanger Z block on, so the block-FCS
# alternation is finally exercised) and a missing outcome (fires MID).
#
# See ../PLAN_pipeline_validation.md §7 for the acceptance criteria.
#
#   ./run_v2_robustness.sh                          # all scenarios, 300 reps
#   N_REP=20 ./run_v2_robustness.sh                 # smoke: check it runs
#   SCENARIOS=base,missing_y20 ./run_v2_robustness.sh   # subset
#   BIG_N=80000 BIG_N_REP=30 ./run_v2_robustness.sh # push the scale cell
#
# Watch:  tail -f validation/phase1/logs/v2_<stamp>.log
# Stop:   kill $(cat validation/phase1/logs/v2_<stamp>.pid)
# Done:   validation/phase1/results/v2_latest.rds (+ v2_summary.csv, v2_raw.csv)
#
# RESOURCES. Parallelism is at the REPLICATION level (one fork per task, cores
# set below); every procedure runs single-threaded inside its fork, so there is
# no nested forking. Memory is dominated by the large_n cell: each worker holds
# m completed datasets of n rows, i.e. roughly
#     n x 10 cols x 8 B x M  bytes  per worker
# ~50 MB/worker at BIG_N=20000 (M=30) and ~190 MB/worker at BIG_N=80000 --
# about 1 GB and 4 GB across 22 workers respectively, both comfortable on this
# box. BLAS is pinned to one thread so it does not fight the fork pool.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# Single-threaded BLAS: the R-level fork pool owns the parallelism.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1

export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"

export N_REP="${N_REP:-300}"
export M="${M:-30}"
export SEED="${SEED:-20260825}"
export NCORES="${NCORES:-22}"        # 24-core box; leave 2 for the OS and parent
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export BIG_N="${BIG_N:-20000}"
export BIG_N_REP="${BIG_N_REP:-50}"

if [[ ! -f "${PIPELINE_ROOT}/00_censored_exposure.R" ]]; then
  echo "ERROR: PIPELINE_ROOT='${PIPELINE_ROOT}' has no 00_censored_exposure.R" >&2
  exit 1
fi

mkdir -p results logs
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="logs/v2_${STAMP}.log"
PIDF="logs/v2_${STAMP}.pid"

echo "Launching V2: N_REP=${N_REP} M=${M} NCORES=${NCORES} SWEEPS=${SWEEPS} MARGIN=${MARGIN}"
echo "              BIG_N=${BIG_N} BIG_N_REP=${BIG_N_REP} SCENARIOS=${SCENARIOS:-<all>}"
echo "              PIPELINE_ROOT=${PIPELINE_ROOT}"
nohup Rscript run_v2_robustness.R > "${LOG}" 2>&1 &
echo $! > "${PIDF}"
echo "PID $(cat "${PIDF}")  |  log: ${LOG}"
echo "Watch: tail -f validation/phase1/${LOG}"
echo "Stop:  kill \$(cat validation/phase1/${PIDF})"
