#!/usr/bin/env bash
# =============================================================================
# Track V12 -- is the non-linear-outcome penalty the SHAPE of the Z-block draw?
# -----------------------------------------------------------------------------
# THE QUESTION. V11 found a non-linear outcome costs ~3 pp of bias for EVERY
# imputer -- BART, mice pmm, forest, bootstrapped forest -- with a spread of only
# 0.44 pp. That is strange if the conditional MEAN were the problem, since those
# arms differ enormously in how flexibly they model it.
#
# THE DIAGNOSIS. They all draw a covariate as
#
#     (fitted conditional mean)  +  homoscedastic Gaussian noise
#
# Computed exactly from the generator, p(Z1 | Y, X) has constant SD 0.894 and
# zero skew at every Y under a LINEAR outcome -- that draw is exactly right. Under
# a NON-LINEAR outcome the same conditional runs SD 0.626 -> 1.216 with skew
# reaching -1.338. The imputers share the part that is wrong (the noise) and
# differ only in the part that is not (the mean), which is precisely the pattern
# V11 measured.
#
# TWO HYPOTHESES WERE TESTED AND REJECTED FIRST, and are recorded so nobody
# re-runs them: the conditional does NOT become bimodal (the N(0,1) prior keeps
# it unimodal at every Y), and the harm is NOT in the X block (the exposure's
# predictors are complete in these cells, and the oracle arm -- which imputes
# nothing -- degrades by 0.01 pp).
#
# THE MANIPULATION. Two Z-block draws sharing the SAME, CORRECT conditional mean:
#
#   smc_zgauss   exact mean + HOMOSCEDASTIC GAUSSIAN noise -- what every shipped
#                imputer effectively does, but handed the mean exactly, so no
#                mean-modelling error remains to confound the comparison
#   smc_zexact   a draw from the TRUE conditional, spread and skew included
#
# Both use the same exact exposure draw and the same outcome model, so their
# difference is the shape of the Z draw and nothing else. pipeline_bartMI is
# carried as the reference point.
#
# ACCEPTANCE (registered before the run; full text in ../PLAN_pipeline_validation.md §8e)
#   CONTROL   in the LINEAR cells the two draws must be indistinguishable -- the
#             true conditional is Gaussian there, so shape has nothing to add.
#             A difference means they differ for some other reason: run is void.
#   PRIMARY   the PAIRED smc_zgauss - smc_zexact difference in the four
#             non-linear cells: the cost of getting shape wrong, mean held right.
#   SECONDARY how much of V11's ~3 pp that accounts for.
#   If the contrast is large, the fix is a shape-aware Z-block sampler -- a
#   concrete change to run_row_level_imputation_bart(), not research. If small,
#   the mechanism is elsewhere again.
#
# NOT THE SAME BUG AS V3/V9. That one is the X block drawing the exposure from a
# conditional with the wrong functional FORM. This is the Z block drawing a
# covariate with the right mean and the wrong SHAPE. One fix will not close both.
#
# COST (MEASURED: 88 tasks, 3 arms + oracle, m=30, 22 workers -> 1.7 min)
#   8 cells x  300 reps = 2400 tasks -> ~46 min
#   8 cells x 1000 reps = 8000 tasks -> ~2.6 h   <- the registered design
#   Cheap because the analysis model is lm, not MCMC. The shape contrast is
#   paired on identical data, so its MC error should be well under 0.2 pp.
#
# BEFORE RUNNING: `Rscript test_v12_zdraw.R` -- 16 checks, no model fitting, ~2 s.
# It verifies that "exact" reproduces the true conditional's mean, SD AND skew,
# and that "gaussian" matches the mean while losing the shape. Without that the
# contrast would confound shape with mean and answer nothing.
#
# USAGE
#   ./run_v12_zdraw.sh                    # the registered run (~2.6 h)
#   N_REP=100 ./run_v12_zdraw.sh          # short pilot
#   SCENARIOS=mcar_z40,ynl_mcar_z40 ./run_v12_zdraw.sh   # one matched pair
#
# Results: results/v12_latest.rds, v12_summary.csv, v12_raw.csv
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v12}"

if [[ "${V12_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V12_STAMP="$(date +%Y%m%d-%H%M%S)"
  V12_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V12_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V12 Z-draw-shape run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V12_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V12_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V12_STAMP}.pid)"
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
export ARMS="${ARMS:-pipeline_bartMI,smc_zgauss,smc_zexact}"

STAMP="${V12_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

{
  echo "=== Track V12: Z-block draw, shape versus mean ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
  echo "CONTROL: in the LINEAR cells smc_zexact and smc_zgauss must be"
  echo "         indistinguishable. A difference there voids the run."
  echo "=================================================="
} > "$LOG"

# Record the PID of the worker, not this wrapper -- see run_v11_ynl.sh.
Rscript run_v4_variance.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
echo "[$(date '+%H:%M:%S')] done rc=$rc -> results/${TAG}_latest.rds" >> "$LOG"
