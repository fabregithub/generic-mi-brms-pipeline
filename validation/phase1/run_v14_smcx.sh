#!/usr/bin/env bash
# =============================================================================
# Track V14 -- does a SHIPPABLE exposure draw work? (roadmap item 07)
# -----------------------------------------------------------------------------
# V13 showed the exposure draw carries most of the non-linear-outcome penalty
# (3.1-3.4 pp of ~3.9); V9 showed it carries the whole mixture failure. Both were
# measured against a sampler that KNEW the generator's surface. This run asks
# whether a draw knowing only the analysis FORMULA can match it.
#
#   arm            exposure draw
#   pipeline_bartMI  shipped: leftcens, conditions on Y LINEARLY
#   smc_zexact       exact: knows the generator's coefficients (the ceiling)
#   smc_xgrid        NEW: grid draw, knows only the FORMULA, estimates the rest
#                    and draws its parameters from their posterior each sweep
#
# WHY A GRID AND NOT IMPORTANCE SAMPLING -- tested and rejected. With a non-linear
# exposure-response, mu(x) = Y can have a second root far from the linear
# solution (+2.0 and -4.667 here; censoring admits only -4.667). Importance
# sampling from a linear-conditional proposal misses it, and ESS DOES NOT WARN
# YOU: at sigma_y = 0.1 the draw was off by 4.6 while ESS/K read 0.999, because
# uniformly-bad candidates give uniform weights. A grid has no proposal to
# misplace. See test_v14_smcx.R.
#
# ACCEPTANCE (full text in ../PLAN_pipeline_validation.md §8g)
#   CONTROL   in the four LINEAR cells smc_xgrid must agree with smc_zexact and
#             pipeline_bartMI within MC error -- a linear outcome makes the
#             shipped conditional correct, so there is nothing to improve there.
#   PRIMARY   paired smc_xgrid - smc_zexact in the non-linear cells.
#             Pre-declared adequate if within +/-1 pp.
#   SECONDARY how much of the ~3 pp it recovers, AND whether b / b_share stay
#             comparable to the shipped arm -- recovering bias by collapsing
#             between-imputation variance would be the V4 defect in new clothing.
#
# WHAT THE OUTCOMES MEAN. Within 1 pp: item 07 has a shippable method and the rest
# is engineering (wire .ce_smc_x_grid into 00_censored_exposure.R behind a flag).
# Materially worse: the gap between knowing the surface and estimating it IS the
# obstacle, and 07 needs re-scoping again.
#
# SCALAR PATH ONLY. Whether the same draw recovers the MIXTURE estimands needs
# the BKMR harness (~10 h, separate). This runs first because it is cheap and a
# failure here would make the expensive run pointless.
#
# BEFORE RUNNING: `Rscript test_v14_smcx.R` -- 24 checks, no model fitting, ~10 s.
#
# COST (MEASURED: 88 tasks, 3 arms + oracle, m=30, 22 workers -> 2.1 min)
#   8 cells x 1000 reps = 8000 tasks -> ~3.2 h   <- the registered design
#   4 cells x 1000 reps = 4000 tasks -> ~1.6 h   (non-linear cells only)
#
# USAGE
#   ./run_v14_smcx.sh                  # the registered run (~3.2 h)
#   N_REP=100 ./run_v14_smcx.sh        # short pilot
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v14}"

if [[ "${V14_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V14_STAMP="$(date +%Y%m%d-%H%M%S)"
  V14_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V14_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V14 shippable-exposure-draw run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V14_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V14_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V14_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260825}"
export N_REP="${N_REP:-1000}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export OUT_TAG="${TAG}"
export SCENARIOS="${SCENARIOS:-mcar_z40,ynl_mcar_z40,missing_y20,ynl_missing_y20,combined,ynl_combined,mar_z40,ynl_mar_z40}"
export ARMS="${ARMS:-pipeline_bartMI,smc_zexact,smc_xgrid}"

STAMP="${V14_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

{
  echo "=== Track V14: a shippable exposure draw (item 07) ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
  echo "CONTROL: in the LINEAR cells smc_xgrid must agree with smc_zexact"
  echo "         and pipeline_bartMI within MC error."
  echo "=================================================="
} > "$LOG"

# Record the PID of the worker, not this wrapper -- see run_v11_ynl.sh.
Rscript run_v4_variance.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
echo "[$(date '+%H:%M:%S')] done rc=$rc -> results/${TAG}_latest.rds" >> "$LOG"
