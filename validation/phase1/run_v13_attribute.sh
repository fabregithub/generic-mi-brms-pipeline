#!/usr/bin/env bash
# =============================================================================
# Track V13 -- attributing the REST of the non-linear-outcome penalty
# -----------------------------------------------------------------------------
# V12 confirmed the Z-draw SHAPE mechanism but could attribute only 32-49% of the
# ~3 pp penalty to it. The rest lies in the Z conditional MEAN, the EXPOSURE
# draw, or both -- V12 varied all of them at once. This run separates them with a
# 2x2 plus the shipped reference:
#
#   arm                  Z draw                        X draw
#   pipeline_bartMI      BART mean + Gaussian          shipped (leftcens)
#   smc_zgauss_xship     EXACT mean + Gaussian         shipped (leftcens)
#   smc_zexact_xship     EXACT draw (mean + shape)     shipped (leftcens)
#   smc_zgauss           exact mean + Gaussian         EXACT
#   smc_zexact           EXACT draw                    EXACT
#
# giving a complete additive decomposition, every step paired on identical data:
#   bartMI      -> zgauss_xship   the Z conditional MEAN
#   zgauss_xship-> zexact_xship   the Z draw SHAPE
#   zexact_xship-> zexact         the EXPOSURE draw
#
# ACCEPTANCE (full text in ../PLAN_pipeline_validation.md §8f)
#   CONTROL 1  in the LINEAR cells all four SMC arms must agree within MC error.
#   CONTROL 2  the _xship arms must reproduce bartMI's `b` and `b_share`,
#              confirming they really do share its exposure draw.
#   PRIMARY    the three-way decomposition, which must roughly sum to the total.
#   SECONDARY  is a Z-ONLY fix sufficient? If smc_zexact_xship still carries most
#              of the penalty, roadmap item 08 must be re-scoped to include the
#              exposure draw -- which drags in item 07 and stops being cheap.
#
# A BUG THIS DESIGN ALREADY CAUGHT. The first implementation passed leftcens
# filled values AND point bounds on every row. The pipeline's convention is the
# opposite: observed rows supply `y` with NA bounds, censored rows supply bounds
# with y = NA. The wrong version put the _xship arms 5 pp adrift IN THE LINEAR
# CELLS, where they should have matched -- caught only because the smoke test
# compared them against a cell whose answer was already known.
#
# COST (MEASURED: 88 tasks, 5 arms + oracle, m=30, 22 workers -> 1.4 min)
#   8 cells x 1000 reps = 8000 tasks -> ~3.9 h   <- the registered design
#
# USAGE
#   ./run_v13_attribute.sh              # the registered run (~3.9 h)
#   N_REP=100 ./run_v13_attribute.sh    # short pilot
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v13}"

if [[ "${V13_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V13_STAMP="$(date +%Y%m%d-%H%M%S)"
  V13_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V13_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V13 attribution run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V13_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V13_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V13_STAMP}.pid)"
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
export ARMS="${ARMS:-pipeline_bartMI,smc_zgauss_xship,smc_zexact_xship,smc_zgauss,smc_zexact}"

STAMP="${V13_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

{
  echo "=== Track V13: attributing the non-linear-outcome penalty ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
  echo "CONTROL: in the LINEAR cells all four SMC arms must agree, and the"
  echo "         _xship arms must reproduce bartMI's b / b_share."
  echo "=================================================="
} > "$LOG"

# Record the PID of the worker, not this wrapper -- see run_v11_ynl.sh.
Rscript run_v4_variance.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
echo "[$(date '+%H:%M:%S')] done rc=$rc -> results/${TAG}_latest.rds" >> "$LOG"
