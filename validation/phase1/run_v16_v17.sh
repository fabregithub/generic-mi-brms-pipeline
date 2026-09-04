#!/usr/bin/env bash
# =============================================================================
# V16 + V17 -- the root claim, then the covariate's causal role. THREE STAGES,
#              run in sequence, each with its own registered prediction.
# -----------------------------------------------------------------------------
# WHY THREE STAGES AND NOT ONE RUN. run_v4_variance.R applies one arm set to
# every scenario it is given, and these questions do not share an arm set:
#
#   Stage 1  V16   Is the root claim true, in the structure the pipeline ships?
#                  4 cells x the noY 2x2. The covariate is causally inert in all
#                  four -- which is exactly why stage 2 exists.
#
#   Stage 2  V17a  Does the answer depend on what Z IS? Same 2x2, but Z is now a
#                  confounder / mediator / collider / both. The fork cell is the
#                  first time in this project that a covariate has HAD to be
#                  adjusted for, so it is the first time the Z block's imputation
#                  can bite the focal estimand at all -- and the derivation says
#                  it bites hard: noYz goes from -0.5 pp in the precision
#                  structure every earlier track used, to +26.7 pp in the fork
#                  and +30.4 pp in the pipe.
#
# MCAR *AND* MAR, IN EVERY STAGE. The covariate mechanism is not a detail here.
# Under MCAR, omitting Y from the Z block costs information -- the congeniality
# penalty this project has been measuring since V1. But
# `inject_mar_covariates()` drives missingness off **Y** (and logX2), so under
# MAR the same omission drops the variable the MECHANISM depends on: the
# imputation no longer conditions on what makes the data MAR, which is a
# VALIDITY failure, not an efficiency one. Those are different claims and they
# get separate cells, paired by `seed_as` so each MCAR/MAR contrast sees
# byte-identical complete data and byte-identical exposure censoring.
#
#   Stage 3  V17b  Does `use_as_auxiliary` reach the censored draw? Only
#                  meaningful in the collider cell, and only against the
#                  bartMI reference -- so it cannot share stage 2's arm set.
#                  This is the one stage that exercises 00_censored_exposure.R's
#                  own `auto_preds`, a code path NO track since V1 has run.
#
# Stages run strictly in order and each writes its own results/<tag>_latest.rds,
# so a failure in one leaves the others' output intact and re-runnable. A stage
# that exits non-zero STOPS the sequence -- a later stage's numbers are not worth
# having if an earlier one broke the harness.
#
# THE REGISTERED PREDICTIONS live in ../PLAN_pipeline_validation.md §8i (V16) and
# §8j (V17) and are echoed into each stage's log at launch. They are DERIVED
# (predict_roles.R mirrors the engine's alternation with correctly-specified
# conditionals), not fitted -- so a miss is informative rather than embarrassing.
#
# COST -- MEASURED, not extrapolated. A full-wave pilot (22 reps = one task per
# worker per cell, 352 tasks, 0 errors) at the run's own settings (m = 30,
# sweeps = 3, 22 workers):
#
#   stage 1  V16    103 min   stage 2  V17a  67 min   stage 3  V17b  12 min
#   TOTAL ~3.1 h
#
# Two cells are 38% of the whole run: `combined` and `nl_combined` cost ~93 s per
# task against ~22 s for every other cell, because they censor all THREE
# exposures rather than the focal one. Dropping them
#
#   SCENARIOS is per-stage, so use STAGES plus an explicit list -- e.g.
#   STAGES=1 SCENARIOS=mcar_z40,nl_mcar_z40,mar_z40,nl_mar_z40 ./run_v16_v17.sh
#
# turns stage 1 into ~33 min and the sequence into ~1.9 h. Earlier quotes of
# 1.8 h and 5.0 h in this file's history came from 3-rep pilots, where per-worker
# setup dominates; both were noise, and this is the figure to trust.
#
# USAGE
#   ./run_v16_v17.sh                  # all three stages, in order (~3.1 h)
#   STAGES=1 ./run_v16_v17.sh         # V16 only
#   STAGES=2,3 ./run_v16_v17.sh       # the V17 pair only
#   N_REP=10 ./run_v16_v17.sh         # pilot every stage (use this to time it)
#
# Watch:  tail -f validation/phase1/logs/v16v17_<stamp>.log
# Alive?: bash validation/phase1/run_status.sh v16v17
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v16v17}"

if [[ "${SEQ_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export SEQ_STAMP="$(date +%Y%m%d-%H%M%S)"
  SEQ_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${SEQ_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V16 + V17 sequence detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${SEQ_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${SEQ_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${SEQ_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260903}"
export N_REP="${N_REP:-500}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"

STAMP="${SEQ_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NOY_ARMS="pipeline_bartMI,pipeline_noYx,pipeline_noYz,pipeline_noYboth"

stage_tag=(  ""  "v16"                                          "v17a"                                    "v17b" )
stage_name=( ""  "V16 - the root claim"                         "V17a - the covariate's causal role"      "V17b - auxiliary reaches the X block?" )
stage_scen=( ""  "mcar_z40,combined,nl_mcar_z40,nl_combined,mar_z40,nl_mar_z40" \
                 "zr_fork,zr_pipe,zr_collider,zr_mixed,zrmar_fork,zrmar_pipe,zrmar_collider,zrmar_mixed" \
                 "zr_collider,zrmar_collider" )
stage_arms=( ""  "$NOY_ARMS"                                    "$NOY_ARMS"                               "pipeline_bartMI,pipeline_auxZ_shipped,pipeline_auxZ_asdoc" )

{
  echo "=== V16 + V17 sequence ==="
  echo "started: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SWEEPS=${SWEEPS} SEED=${SEED} NCORES=${NCORES}"
  echo ""
  echo "PREDICTIONS (paired shift vs pipeline_bartMI, percentage points)."
  echo "Derived by predict_roles.R; criteria in ../PLAN_pipeline_validation.md §8i, §8j."
  echo ""
  echo "  Stage 1  V16, linear cells:  noYboth -14.6 | noYx -15.8 | noYz -0.5"
  echo "           leak = noYboth-(noYx+noYz) = +1.6, POSITIVE: the mixed"
  echo "           configuration should be WORSE than omitting Y from both."
  echo "           (predict_roles.R supersedes predict_v16.R, which imputed only"
  echo "           Z1; the numbers moved <0.3 pp, so the registration stands.)"
  echo ""
  echo "  Stage 2  V17a -- noYz shift vs bartMI, in pp, by role and mechanism:"
  echo "                        MCAR      MAR     diff"
  echo "           precision    -0.5     +0.1     +0.5"
  echo "           fork        +26.7    +32.4     +5.7"
  echo "           pipe        +30.4    +36.1     +5.7"
  echo "           collider     -0.7     +3.1     +3.8"
  echo "           mixed       +19.5    +23.8     +4.4"
  echo ""
  echo "           TWO SEPARABLE RESULTS."
  echo "           (a) STRUCTURAL 2x2: the penalty is large exactly when Z is"
  echo "               adjusted for AND on an open X-Y path -- confounding"
  echo "               (fork) or mediating (pipe) alike -- and about zero when"
  echo "               Z is adjusted for but off-path (precision), or on a path"
  echo "               the analysis correctly omits (collider, MCAR)."
  echo "           (b) MECHANISM: MAR-on-Y adds +4 to +6 pp in every on-path"
  echo "               cell -- a validity failure on top of the congeniality"
  echo "               one -- but adds NOTHING off-path. The precision row is"
  echo "               the informative null: the Z imputation is genuinely"
  echo "               invalid there and b1 still does not move."
  echo ""
  echo "           Rejected if |noYz| is not >10 pp in fork/pipe/mixed or not"
  echo "           <3 pp in precision, under EITHER mechanism; or if MAR-MCAR"
  echo "           is not positive in all four on-path cells. Bars: PLAN 8j."
  echo ""
  echo "  Stage 3  V17b. auxZ_shipped drops the auxiliary from the X block;"
  echo "           auxZ_asdoc keeps it. Any difference IS the defect's size."
  echo "           A null here says the fix is not worth making."
  echo "==================================================================="
} > "$LOG"

WANT="${STAGES:-1,2,3}"
rc=0
for i in 1 2 3; do
  case ",$WANT," in *",$i,"*) ;; *) continue ;; esac
  t="${stage_tag[$i]}"
  {
    echo ""
    echo "-------------------------------------------------------------------"
    echo "[$(date '+%H:%M:%S')] STAGE $i : ${stage_name[$i]}"
    echo "  scenarios: ${stage_scen[$i]}"
    echo "  arms:      ${stage_arms[$i]}"
    echo "  results:   results/${t}_latest.rds"
    echo "-------------------------------------------------------------------"
  } >> "$LOG"

  OUT_TAG="$t" SCENARIOS="${stage_scen[$i]}" ARMS="${stage_arms[$i]}" \
    Rscript run_v4_variance.R >> "$LOG" 2>&1 &
  RPID=$!
  # Record the WORKER's pid, not this wrapper's -- see run_v11_ynl.sh. Rewritten
  # per stage so run_status.sh and `kill` always point at what is running now.
  printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
  wait "$RPID"; rc=$?

  if [[ $rc -ne 0 ]]; then
    echo "[$(date '+%H:%M:%S')] STAGE $i FAILED (rc=$rc) -- stopping the sequence." >> "$LOG"
    echo "  Earlier stages' results are complete and re-runnable with STAGES=." >> "$LOG"
    break
  fi
  echo "[$(date '+%H:%M:%S')] stage $i done -> results/${t}_latest.rds" >> "$LOG"
done

echo "" >> "$LOG"
echo "[$(date '+%H:%M:%S')] sequence finished rc=$rc" >> "$LOG"
