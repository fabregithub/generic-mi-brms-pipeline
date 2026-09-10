#!/usr/bin/env bash
# =============================================================================
# Track V24 -- the DAG-factorised exposure draw, all four cases in ONE sweep
# -----------------------------------------------------------------------------
# WHAT THIS TESTS. THEORY.md §6b derives the target for a left-censored exposure
# from Bayes plus the DAG's factorisation:
#
#   p(x | rest, x<=L)  ∝  p(x | Pa(X)) · p(Y | x, Pa(Y)) · ∏_{C∈Ch(X)} p(C | x, ·) · 1{x<=L}
#
# and concludes the third factor is needed in FOUR cases, which differ along two
# axes -- the covariate's causal role, and (for a pipe) which estimand:
#
#   dag_fork        fork      direct   Z1 adjusted (backdoor);  Z1 a PARENT of X
#   dag_pipe_dir    pipe      direct   Z1 adjusted;             Z1 a CHILD of X
#   dag_pipe_tot    pipe      total    Z1 NOT adjusted;         Z1 a CHILD of X
#   dag_pnl_dir     pipe_nl   direct   as pipe_dir, arrow non-linear
#   dag_pnl_tot     pipe_nl   total    as pipe_tot, arrow non-linear
#   dag_collider    collider  direct   Z1 NOT adjusted;   child of X *and* of Y
#
# WHY ALL SIX IN ONE SWEEP, rather than a track per case. The claim under test is
# not "the draw is unbiased here", it is "bias direction is DAG-aware and the
# factorisation fixes each case for its own reason". Cases run separately cannot
# rule out a sign being an artefact of one cell, and the two linear `pipe` cells
# are only meaningful AS controls for the two `pipe_nl` ones. `seed_as` pairs
# every pipe-family cell to `dag_pipe_dir`, so the contrasts are paired.
#
# THE ARM LADDER, in increasing order of what the target conditions on:
#
#   dag_noY     parents + truncation only            the incongenial draw
#   dag_Yonly   + p(Y | x, Pa(Y))                    what V0-V23 all used
#   dag_lin     + p(Z1 | x) fitted LINEAR            the only assumption-free form
#   dag_quad    + p(Z1 | x) fitted QUADRATIC         sensitivity
#   dag_cubic   + p(Z1 | x) fitted CUBIC             sensitivity
#   dag_true    + p(Z1 | x) at the TRUE arrow        ORACLE -- upper bound only
#
# `dag_true` IS NOT SHIPPABLE and is not proposed as a method. It is the ceiling:
# it answers "is the FACTORISATION right", separately from "can the child model
# be estimated". THEORY.md §0b measures that the second answer is no -- below the
# LOD the arrow's curvature is not identified (a natural spline returns exactly
# zero by construction, a quadratic is 4.6x low, a cubic false-positives). So
# dag_lin/quad/cubic are a SENSITIVITY AXIS, never a tuning choice.
#
# WHY THE COVARIATES ARE FULLY OBSERVED (mcar_frac = 0, against 0.4 in every
# zr_* cell). All three factors condition on Z, so a missing Z1 is not a harder
# version of this question -- it is a different one, and it would reintroduce
# V21's covariate-draw ladder as a confound. V20 and V21 already settled the
# covariate block. `dag_impute_datasets()` refuses missing covariates outright
# rather than returning NA. The combined case is open; see ROADMAP.md.
#
# REGISTERED CRITERIA (PLAN §11: each must be able to fail for the reason it
# names, and none may be satisfiable by an uninformative result).
#
#   C1  FACTORISATION. With `dag_true`, all six cells have |bias/SE| < 0.3.
#       Fails if the factorisation is wrong in any case, including the two that
#       V0-V23 never estimated (`*_tot`). Cannot pass vacuously: `dag_noY` in
#       the same cells must exceed 0.3, or the cell has no bias to remove.
#   C2  THE CHILD FACTOR IS LOAD-BEARING, and only where the theory says. In
#       `dag_pnl_dir` and `dag_pnl_tot`, `dag_true` improves on `dag_Yonly` by
#       at least 1.0 in |bias/SE|; in `dag_fork` the two are BIT-IDENTICAL
#       (Z1 is a parent, so `dag_spec()` gates the factor off).
#   C3  A MISSPECIFIED CHILD FACTOR IS WORSE THAN NONE. In the two `pipe_nl`
#       cells, `dag_lin` has strictly larger |bias/SE| than `dag_Yonly`.
#       If C3 HOLDS, a linear child model must not be the shipped default and
#       docs/covariate-roles.md gains a declared-assumption section. If C3
#       FAILS, `dag_lin` degrades gracefully and is a defensible default --
#       which is the more useful outcome, so this gate is written to be
#       falsifiable in the direction that would help.
#   C4  ESTIMAND. `dag_pipe_tot` and `dag_pnl_tot` recover their own
#       `estimand_true` (0.700 and 0.830), not b1 = 0.400. A cell scoring near
#       -43% is the estimand wiring failing, not the draw.
#   C5  DIRECTION. Across the six cells, `dag_noY`'s bias sign is recorded per
#       cell. THEORY.md §4b claims toward-null only when the outcome-borrowing
#       term T2 vanishes; this measures whether any cell is anti-conservative.
#       Descriptive: no pass/fail, and it must not be reported as one.
#
# MOST CONSEQUENTIAL OUTCOME. If C1 holds and C3 holds, the algorithm is correct
# but its one unidentifiable input decides whether it helps or hurts -- which
# makes the deliverable a sensitivity procedure, not a drop-in fix.
#
# COST -- measured by full-wave pilots (N_REP=4 -> 24 tasks ~ 22 workers), never
# extrapolated from a single rep. See CLAUDE.md.
#     n     pilot wave    predicted    ACTUAL (2026-09-09 run)
#             (24 tasks)     stage
#    800       0.4 min        20 min       13 min
#   3200       2.0 min       100 min       71 min
#  12800       7.5 min       375 min      237 min
#   TOTAL                    ~8.2 h      320 min (5.3 h) -- 35% under
#
# The full-wave pilots were CONSERVATIVE by 35%, which is the right direction to
# be wrong in and the first time in this project an estimate has not needed a
# caveat. A wave of 24 tasks on 22 workers is 2 waves, not 1, so 50x the wave
# time over-counts the tail; that is the whole of the discrepancy.
#
# Measured 2026-09-09 on 22 workers, with nothing else on the machine. An earlier
# set of pilots taken while a 20-worker job was running read 0.6 / 2.9 / 11.2 min
# -- 45-50% high. Contention is the single biggest source of error in these
# numbers, so a pilot run alongside other work is an upper bound, not a cost.
#
# n = 12800 is 77% of the run, and every gate here is a LEVEL rather than a
# slope, so STAGES=1,2 (2.2 h) settles C1-C4 and still spans 4x in n. The third
# level matters for one reason only: bias/SE grows as n^(1/6) (V18), so an arm
# that sits just inside 0.3 at n = 800 is at 0.3 x 16^(1/6) = 0.48 at n = 12800.
# C1 is therefore not decided until stage 3 for any arm that passes narrowly.
#
# USAGE
#   ./run_v24_dagdraw.sh                 # all three n levels
#   STAGES=1 ./run_v24_dagdraw.sh        # n = 800 only
#   N_REP=4 ./run_v24_dagdraw.sh         # full-wave timing pilot
#
# Results: results/v24_n<N>_latest.rds per stage.  Analyse with analyze_v24.R.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v24}"

if [[ "${V24_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V24_STAMP="$(date +%Y%m%d-%H%M%S)"
  V24_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V24_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V24 DAG-draw run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V24_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V24_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V24_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260909}"
export N_REP="${N_REP:-200}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export SCENARIOS="${SCENARIOS:-dag_fork,dag_pipe_dir,dag_pipe_tot,dag_pnl_dir,dag_pnl_tot,dag_collider}"
export ARMS="${ARMS:-dag_noY,dag_Yonly,dag_lin,dag_quad,dag_cubic,dag_true}"

# N_REP = 200 comes from the CURRENCY, not from habit. The gates are on bias/SE
# at 0.3, and the Monte-Carlo error of an R-rep mean bias/SE is about 1/sqrt(R):
# 0.16 at R = 40, which is half the gate and cannot resolve it; 0.07 at R = 200,
# which can. C2's 1.0-point improvement bar needs less, but it is read off the
# same rows.

STAMP="${V24_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 3200 12800)

{
  echo "=== Track V24: the DAG-factorised exposure draw, four cases at once ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "THE TARGET (THEORY.md 6b):"
  echo "  p(x | rest, x<=L) ~ p(x|Pa(X)) . p(Y|x,Pa(Y)) . PROD p(C|x,.) . 1{x<=L}"
  echo "The third factor is what V0-V23 never had. It is needed in FOUR cases:"
  echo "  fork (parent, no child factor) | pipe direct | pipe total | collider"
  echo ""
  echo "ARM LADDER: noY -> Yonly -> child factor at lin/quad/cubic/TRUE."
  echo "  dag_true is an ORACLE and is NOT a proposed method. It separates"
  echo "  'is the factorisation right' from 'can the child model be estimated',"
  echo "  and THEORY.md 0b already measured the second answer as NO below the LOD."
  echo ""
  echo "REGISTERED CRITERIA"
  echo "  C1 dag_true: |bias/SE| < 0.3 in all six cells (and dag_noY > 0.3, or"
  echo "     the cell had no bias to remove and C1 passed vacuously)"
  echo "  C2 child factor gains >= 1.0 in |bias/SE| in the pipe_nl cells, and is"
  echo "     BIT-IDENTICAL to dag_Yonly under the fork"
  echo "  C3 dag_lin is strictly WORSE than dag_Yonly in the pipe_nl cells"
  echo "     -> if so, a linear child model cannot be the shipped default"
  echo "  C4 the *_tot cells recover 0.700 / 0.830, not b1 = 0.400"
  echo "  C5 direction of dag_noY's bias per cell -- DESCRIPTIVE, not a gate"
  echo ""
  echo "MOST CONSEQUENTIAL OUTCOME: C1 and C3 both holding makes the deliverable"
  echo "  a SENSITIVITY PROCEDURE over a declared child model, not a drop-in fix."
  echo "==========================================================="
} > "$LOG"

WANT="${STAGES:-1,2,3}"
rc=0
for i in 1 2 3; do
  case ",$WANT," in *",$i,"*) ;; *) continue ;; esac
  N="${NLEVELS[$i]}"
  {
    echo ""
    echo "-------------------------------------------------------------------"
    echo "[$(date '+%H:%M:%S')] STAGE $i : n = ${N}"
    echo "  results: results/${TAG}_n${N}_latest.rds"
    echo "-------------------------------------------------------------------"
  } >> "$LOG"

  OUT_TAG="${TAG}_n${N}" N_OBS="$N" Rscript run_v4_variance.R >> "$LOG" 2>&1 &
  RPID=$!
  # Worker PID, not this wrapper's -- see run_v11_ynl.sh. Rewritten per stage.
  printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
  wait "$RPID"; rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[$(date '+%H:%M:%S')] STAGE $i (n=${N}) FAILED rc=$rc -- stopping." >> "$LOG"
    break
  fi
  echo "[$(date '+%H:%M:%S')] stage $i done -> results/${TAG}_n${N}_latest.rds" >> "$LOG"
done

{
  echo ""
  echo "[$(date '+%H:%M:%S')] V24 finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v24.R"
} >> "$LOG"
exit $rc
