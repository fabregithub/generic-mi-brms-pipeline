#!/usr/bin/env bash
# =============================================================================
# Track V15 -- the attenuation shape. THE FIRST PREDICTION-LED TRACK.
# -----------------------------------------------------------------------------
# THE METHOD, adopted 2026-09-01 (../PLAN_pipeline_validation.md §7b): state the
# theory, derive a NUMBER for a cell nobody has measured, run, and revise the
# theory against the result. Tracks V0-V14 registered acceptance thresholds --
# "is the pipeline good enough?" This registers a PREDICTION -- "is our
# understanding right?" -- which can be wrong, and that is the point.
#
# THE OPEN QUESTION. Curvature attenuation grows FASTER than the censored
# fraction. V3 measured 20% non-detects -> 11.7% excess and 40% -> 56.7%: a
# 4.85x rise for 2x the censoring. THEORY.md 4 records this as unexplained, and
# it matters -- it says the penalty for uncongeniality accelerates with the
# amount of missing data.
#
# THE REGISTERED PREDICTION.  excess% = 323.4 * f^2
#
#     f      predicted     status
#   0.10        3.2%       UNMEASURED
#   0.20       12.9%       measured 11.7%   (calibration point)
#   0.30       29.1%       UNMEASURED
#   0.40       51.8%       measured 56.7%   (calibration point)
#   0.50       80.9%       UNMEASURED
#   0.60      116.4%       UNMEASURED
#
# THIS LAW IS FRANK CURVE-FITTING and is registered as such. Seven candidates
# were tested against the observed 4.85x ratio; the mechanistically motivated
# ones fit WORSE (f*E[x^2|cens] -> 1.14, f*|q_f| -> 0.60) than the atheoretical
# f^2 -> 4.00. A weak-but-falsifiable prediction moves the theory; a post-hoc
# description does not. What is NOT allowed is presenting the fitted constant as
# a derivation.
#
# FALSIFICATION, PRE-DECLARED. f^2 is rejected if
#   * the measured excess at ANY of f = 0.10, 0.30, 0.50, 0.60 differs from the
#     prediction by more than 10 percentage points, OR
#   * a log-log regression of excess on f gives a slope whose 95% CI excludes 2.
# The f = 0.60 cell is the sharpest test: 116.4% requires the estimate to have
# flipped sign and overshot, which is easy to refute.
#
# EITHER OUTCOME IS USEFUL. If f^2 survives there is a quantitative law to
# explain and an exponent to derive. If it fails, four new points still pin the
# shape well enough to constrain what a derivation must produce -- more than two
# points can do. The prediction stays in the record either way, with whatever
# supersedes it beside it.
#
# CALIBRATION CHECK BUILT IN. nd20 and nd40 are the cells V3 already measured,
# and they keep canonical indices 1-3, so this run reproduces them from the same
# seeds. If they do not come back at -11.7% and -56.7%, something changed and the
# four new points are not comparable to the two old ones.
#
# COST (derived from V3's measured 5.80 s per BKMR fit at 22 workers)
#   6 cells x 100 reps x 11 fits = 6600 fits -> ~10.6 h
#   MC error ~3.5 pp at 100 reps, against 17-35 pp gaps between candidate laws.
#
# USAGE
#   ./run_v15_fsweep.sh                 # the registered run (~10.6 h)
#   N_REP=10 ./run_v15_fsweep.sh        # short pilot
#   SCENARIOS=nd20,nd40 ./run_v15_fsweep.sh   # calibration cells only
#
# Results: results/v15_latest.rds, v15_summary.csv, v15_raw.csv
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v15}"

if [[ "${V15_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V15_STAMP="$(date +%Y%m%d-%H%M%S)"
  V15_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V15_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V15 f-sweep run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V15_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V15_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V15_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
# SEED must match V3's so the nd20/nd40 calibration cells reproduce.
export SEED="${SEED:-20260828}"
export N_REP="${N_REP:-100}"
export N_OBS="${N_OBS:-800}"
export ITER="${ITER:-3000}"
export M="${M:-10}"
export KNOTS="${KNOTS:-50}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export OUT_TAG="${TAG}"
export SCENARIOS="${SCENARIOS:-mixf10,nd20,mixf30,nd40,mixf50,mixf60}"
export ARMS="${ARMS:-oracle_bkmr,pipeline_bkmr}"

STAMP="${V15_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

{
  echo "=== Track V15: the attenuation shape (first prediction-led track) ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} N_OBS=${N_OBS} M=${M} ITER=${ITER} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
  echo ""
  echo "REGISTERED PREDICTION: curv_X1 excess% = 323.4 * f^2"
  echo "   f=0.10 -> 3.2%   f=0.30 -> 29.1%   f=0.50 -> 80.9%   f=0.60 -> 116.4%"
  echo "REJECTED IF: any cell misses by >10 pp, or a log-log slope CI excludes 2."
  echo "CALIBRATION: nd20 and nd40 must reproduce V3 (-11.7%, -56.7%)."
  echo "This law is curve-fitting to two points and is expected to be refuted."
  echo "==================================================================="
} > "$LOG"

# Record the worker's PID, not this wrapper's -- see run_v11_ynl.sh.
Rscript run_v3_bkmr.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
echo "[$(date '+%H:%M:%S')] done rc=$rc -> results/${TAG}_latest.rds" >> "$LOG"
