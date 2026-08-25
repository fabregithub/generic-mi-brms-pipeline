#!/usr/bin/env bash
# =============================================================================
# Track V4 -- attribute the V2 interval under-coverage. Launches DETACHED.
#
# V2 found 95% intervals ~17% too narrow when covariate AND outcome missingness
# are both substantial, and proved it is NOT a small-M artefact (m=30 vs m=100 on
# identical data: no change). This run attributes the shortfall by changing one
# component at a time, and reports the Rubin variance decomposition (ubar / b /
# fmi) that V2 did not store.
#
# See ../PLAN_pipeline_validation.md §9 for what each outcome would mean.
#
#   ./run_v4_variance.sh                                # 4 scenarios x 4 arms, 150 reps
#   N_REP=20 ./run_v4_variance.sh                       # smoke
#   SCENARIOS=combined ./run_v4_variance.sh             # the worst cell only
#   ARMS=pipeline_block_fcs,pipeline_properZ ./run_v4_variance.sh   # one contrast
#
# Watch:  tail -f validation/phase1/logs/v4_<stamp>.log
# Stop:   kill $(cat validation/phase1/logs/v4_<stamp>.pid)
# Done:   validation/phase1/results/v4_latest.rds (+ v4_summary.csv, v4_raw.csv)
#
# COST. Four pipeline arms per replication, so roughly 4x V2's per-scenario cost
# (V2 measured: combined 266 s/rep, mcar_z40 133, missing_y20 98, base 5.6). At
# 150 reps that is very roughly 4 x 150 x (266+133+98+6) / 3600 ~= 84 CPU-hours,
# about 4 hours on 22 workers. Trim with SCENARIOS= or ARMS= if that is too long;
# `combined` alone with all four arms is ~45 min.
#
# The proper-Z arms are CHEAPER than the miceRanger reference (a linear/logistic
# draw costs far less than 90 random-forest fits), so the estimate above is
# conservative.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1

export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"

export N_REP="${N_REP:-150}"
export M="${M:-30}"
export SEED="${SEED:-20260825}"      # matches V2 -> reference arm sees identical data
export NCORES="${NCORES:-22}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"

if [[ ! -f "${PIPELINE_ROOT}/00_censored_exposure.R" ]]; then
  echo "ERROR: PIPELINE_ROOT='${PIPELINE_ROOT}' has no 00_censored_exposure.R" >&2
  exit 1
fi

mkdir -p results logs
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="logs/v4_${STAMP}.log"
PIDF="logs/v4_${STAMP}.pid"

# Record the config INTO the log, not just the terminal: a detached run's
# settings must be recoverable from its artefacts afterwards.
{
  echo "=== V4 run configuration ==="
  echo "started:       $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES} SWEEPS=${SWEEPS} MARGIN=${MARGIN}"
  echo "SCENARIOS=${SCENARIOS:-<default: base,mcar_z40,missing_y20,combined>}"
  echo "ARMS=${ARMS:-<default: all four>}"
  echo "PIPELINE_ROOT=${PIPELINE_ROOT}"
  echo "==========================="
} > "${LOG}"

echo "Launching V4: N_REP=${N_REP} M=${M} NCORES=${NCORES}"
echo "              SCENARIOS=${SCENARIOS:-<all 4>}  ARMS=${ARMS:-<all 4>}"
nohup Rscript run_v4_variance.R >> "${LOG}" 2>&1 &
R_PID=$!
echo "${R_PID}" > "${PIDF}"

# Idle-sleep protection as a SIBLING, not a parent. Wrapping R in
# `caffeinate -i Rscript ...` was tried and both detached runs then died right
# after their first chunk, while an identical foreground run survived -- so the
# wrapper is implicated. `-w` holds the assertion while the given pid lives
# without sitting in R's process ancestry.
if command -v caffeinate >/dev/null 2>&1; then
  nohup caffeinate -i -w "${R_PID}" >/dev/null 2>&1 &
  echo "caffeinate: idle sleep suppressed while pid ${R_PID} lives (sibling)"
fi
echo "PID $(cat "${PIDF}")  |  log: ${LOG}"
echo "Watch: tail -f validation/phase1/${LOG}"
echo "Stop:  kill \$(cat validation/phase1/${PIDF})"
