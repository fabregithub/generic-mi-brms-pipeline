#!/usr/bin/env bash
# =============================================================================
# Track V11 -- does the imputer choice survive a NON-LINEAR OUTCOME?
# -----------------------------------------------------------------------------
# THE QUESTION. `z_imputer = "bart"` has been the shipped default since v1.5.0,
# chosen over `mice pmm` and `forest_boot` in V6/V7. Every cell in that
# comparison had an outcome LINEAR in (logX, Z) by construction -- and
# FINDINGS_v7.md says so itself:
#
#   "Y is linear in (logX, Z) by construction, which is why parametric arms do
#    well in the outcome-dominated cells. Untested under a non-linear outcome."
#
# The Z block imputes Z1 conditioning on Y. If Y is linear in Z1 then
# p(Z1 | Y, X) is linear-Gaussian, a parametric imputer is correctly specified
# against the outcome, and the comparison can only ever reward calibration --
# never punish getting the conditional mean wrong. That is a confound in a
# SHIPPED DEFAULT, and this run removes it.
#
# THE MANIPULATION. The outcome gains b_zq * (Z1^2 - 1), b_zq = 0.40 (comparable
# to gamma[1] = 0.50). `dgp_formula()` gains the matching I(Z1^2) term, so the
# ANALYSIS model stays correctly specified and b_logX1 = 0.4 remains exactly
# recoverable -- any bias is the imputation model's, never analysis
# misspecification. Verified before registering: the oracle recovers b_logX1 at
# +0.14%, and adding I(Y^2) to a linear Z1 model gives F = 0.0 (p = 0.83) under
# the old design against F = 521.6 (p < 2e-16) under the new one.
#
# DISTINCT FROM V6'S AXIS. `z_form = "nonlinear"` makes Z1 a non-linear function
# of the OTHER PREDICTORS. `y_form = "nonlinear"` makes the OUTCOME non-linear in
# Z1, which is what conditions the imputation draw. Different confounds; only the
# first was ever addressed.
#
# PAIRED BY CONSTRUCTION. Each ynl_* cell carries `seed_as`, so it draws its data
# with its matched linear cell's seed: the pair sees byte-identical exposures,
# covariates and outcome noise and differs ONLY by the curvature term. That makes
# the "how much does pmm degrade?" contrast paired rather than unpaired, cutting
# its Monte-Carlo error 4.3x (measured) -- the difference between a 2.90 pp and a
# 0.65 pp minimum detectable effect, without adding a single replication.
#
# ACCEPTANCE (registered before the run; full text in ../PLAN_pipeline_validation.md §8d)
#   CONTROL   mcar_z40, combined and mar_z40 with pipeline_bartMI must reproduce
#             V10 BIT-IDENTICALLY on reps 1-300. NOT V7: V7's raw numbers are not
#             reproducible from the current tree (V7 and V8 already disagree by
#             up to 0.0128 on the same cell/arm/reps, and V8 predates all V11
#             work -- most likely v1.5.0's z_imputer selector, introduced between
#             them). V8 and V10 ARE bit-identical across 1000 reps, so the
#             post-v1.5.0 tree is the stable reference. missing_y20 has no
#             reproducible baseline and is uncontrolled.
#             A control failure INVALIDATES the run.
#   PRIMARY   is bartMI still best or joint-best on bias in the ynl cells?
#             Paired arm-vs-arm MC error 0.13 pp; min detectable 0.37 pp.
#   SECONDARY the degradation of micePmm from its linear cell to its non-linear
#             one. Paired MC error 0.23 pp; min detectable 0.65 pp.
#   DECISION  if bartMI is no longer best, the v1.5.0 default was decided on a
#             biased comparison and must be revisited. If it remains best, the
#             default is confirmed on a design that could have refuted it.
#   Either answer is useful. This is not expected to change the default -- it is
#   expected to EARN it.
#
# COST (MEASURED: 44 tasks, 4 arms, m=30, 22 workers -> 14.3 min = 0.325 min/task)
#   8 cells x  300 reps = 2400 tasks -> ~13 h   <- the registered design
#   8 cells x  500 reps = 4000 tasks -> ~22 h
#   8 cells x 1000 reps = 8000 tasks -> ~43 h
#   300 is ample because the comparisons that matter are PAIRED (see above).
#   Do not cut M below 30: the CONTROL requires reproducing V7, which used 30.
#
# USAGE
#   ./run_v11_ynl.sh                        # the registered run (~13 h)
#   N_REP=30 ./run_v11_ynl.sh               # short pilot first -- recommended
#   SCENARIOS=mcar_z40,ynl_mcar_z40 ./run_v11_ynl.sh    # one matched pair only
#   ARMS=pipeline_bartMI,pipeline_micePmm ./run_v11_ynl.sh  # the head-to-head only
#
# Results: results/v11_latest.rds, v11_summary.csv, v11_raw.csv
# Checkpoint after every chunk: results/v11_checkpoint.rds -- readable mid-run
#   with summarise_v4(readRDS(...)). Resume is NOT trusted; see run_v4_variance.R.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

# Log and PID names follow OUT_TAG, so results and logs share one identity.
TAG="${OUT_TAG:-v11}"

if [[ "${V11_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V11_STAMP="$(date +%Y%m%d-%H%M%S)"
  V11_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V11_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V11 non-linear-outcome run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V11_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V11_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V11_STAMP}.pid)"
  exit 0
fi

# BLAS single-threaded or it fights the R-level fork loop.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
# SEED is explicit, not defaulted: the CONTROL criterion depends on it matching V7.
export SEED="${SEED:-20260825}"
export N_REP="${N_REP:-300}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export OUT_TAG="${TAG}"
export SCENARIOS="${SCENARIOS:-mcar_z40,ynl_mcar_z40,missing_y20,ynl_missing_y20,combined,ynl_combined,mar_z40,ynl_mar_z40}"
export ARMS="${ARMS:-pipeline_block_fcs,pipeline_properBoot,pipeline_micePmm,pipeline_bartMI}"

STAMP="${V11_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

{
  echo "=== Track V11: non-linear outcome, imputer comparison ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
  echo "CONTROL: mcar_z40 / combined / mar_z40 with pipeline_bartMI must"
  echo "         reproduce V10 bit-identically on reps 1-300 (NOT V7 -- see header)."
  echo "         A control failure invalidates the run."
  echo "=========================================================="
} > "$LOG"

# Record the PID of the process that does the WORK, not this wrapper. If the
# wrapper dies while Rscript continues, a wrapper PID refers to nothing:
# run_status.sh would report "done" for a live run and `kill` would silently
# kill nothing. That happened on the first V9 launch.
Rscript run_v4_variance.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
echo "[$(date '+%H:%M:%S')] done rc=$rc -> results/${TAG}_latest.rds" >> "$LOG"
