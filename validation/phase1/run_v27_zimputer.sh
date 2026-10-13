#!/usr/bin/env bash
# =============================================================================
# Track V27 -- should `z_imputer` default to "mice" instead of "bart"?
# -----------------------------------------------------------------------------
# THE LIVE PROBLEM. The RECOMMENDED configuration (`z_imputer = "bart"`, Y in
# both blocks) carries +4.96% (fork), +5.39% (pipe), +8.52% and +10.40% (both
# under MAR) whenever a covariate is a cause of the exposure or an effect of
# it. The oracle is unbiased in all four, so it is the imputation's bias, not
# the design's -- and COVERAGE DOES NOT DETECT IT (0.93-0.97). It decays as
# n^-1/3, so 2.0-3.6% is still there at n = 12,800. This is the only remaining
# item where a user running the documented configuration is exposed.
#
# THE CANDIDATE, AND WHY IT IS CREDIBLE NOW. V21 attributed **100%** of that
# bias to BART's FLEXIBILITY: a correctly specified parametric draw is
# +0.26%/+0.58%, and improper and donor-matched variants are too, so properness
# moves the interval and not the bias. Measured post-v1.6.1 on matched cells:
# bart +4.83%/+6.03% against mice -0.08%/+1.45%.
#
# WHAT USED TO BLOCK IT. A parametric Z block is misspecified when the
# covariate's conditional is NON-LINEAR, and V6 measured `mice pmm` at -4.70%
# there -- which is exactly why BART was adopted. V22 built the deciding cell
# and the parametric draw won there too, but V22's cell has the non-linearity
# ON an X-Y path. The `nl_*` cells below are V6's original ground, and they are
# the cells this track exists to settle. **If mice loses there, the default
# stays `bart` and mice ships as an option for causally-active covariates.**
#
# CELLS, chosen so the decision cannot be made on friendly ground alone:
#   mcar_z40, mar_z40            inert covariate, linear     -- the control
#   nl_mcar_z40, nl_mar_z40      inert covariate, NON-LINEAR -- V6's ground
#   zr_fork, zr_pipe,
#   zr_collider, zr_mixed        causally active, linear     -- the gain
#   zrmar_fork, zrmar_pipe       causally active, MAR        -- the worst cells
#   zr_pipenl                    non-linear arrow ON-path    -- V22's cell
#
# ARMS: oracle (auto), pipeline_bartMI (today's default), pipeline_miceShipped
# (the real `z_imputer = "mice"` config field, not the harness override -- so a
# difference from `pipeline_micePmm` would be the PORT rather than the method).
#
# REGISTERED CRITERIA (PLAN 11: each must be able to fail for the reason it
# names; none may be satisfiable by an uninformative result).
#
#   C1  NO LOSS ON V6's GROUND. In the four `nl_*` and inert cells, mice's
#       |bias/SE| does not exceed bart's by more than 0.3. **This is the gate
#       that decides the default.** It can fail: V6 measured the parametric
#       draw at -4.70% on exactly this DGP.
#   C2  REAL GAIN WHERE IT MATTERS. In the six causally-active cells, mice
#       improves |bias/SE| by at least 1.0 against bart. Void in any cell where
#       bart's own |bias/SE| is under 0.5 -- there is no bias to remove there.
#   C3  COVERAGE IS NOT TRADED AWAY. In every cell, mice's coverage is at least
#       bart's minus 0.02. A point-estimate gain bought with intervals that
#       stop covering is not a gain.
#   C4  THE PORT IS FAITHFUL. `pipeline_miceShipped` and `pipeline_micePmm`
#       agree to within Monte-Carlo error in every cell. Fails if the shipped
#       implementation differs from the one the evidence was gathered on.
#   C5  n-SCALING, descriptive. Both arms at three n. No pass/fail.
#
# MOST CONSEQUENTIAL OUTCOME. If C1 and C2 both hold, the default should switch
# and a documented 5-10% bias leaves the recommended path. If C1 fails, it must
# not switch, and the honest result is a conditional recommendation keyed on
# whether the covariate is causally active -- which users cannot always know.
#
# COST -- measured by full-wave pilots, never extrapolated (CLAUDE.md).
#     n     pilot (24 tasks)   full stage (11 cells x 500 reps = 5500 tasks)
#    800          24 s              1.5 h
#   3200          52 s              3.3 h
#  12800         174 s             11.1 h
#   TOTAL                          ~15.9 h
#
# Measured 2026-10-11 on 22 workers, machine otherwise idle. n = 12800 is 70%
# of the run. STAGES=1,2 is 4.8 h and settles C1-C4 as levels; the third stage
# is for C5's n-scaling and for the fact that the bias this track targets
# decays as n^-1/3, so the GAIN shrinks with n and the largest-n cell is where
# the case for switching is weakest.
#
# USAGE
#   ./run_v27_zimputer.sh              # all stages
#   STAGES=1 ./run_v27_zimputer.sh     # n = 800 only
#   N_REP=4 ./run_v27_zimputer.sh      # full-wave timing pilot
#
# Results: results/v27_n<N>_latest.rds per stage. Analyse with analyze_v27.R.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"
TAG="${OUT_TAG:-v27}"

if [[ "${V27_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V27_STAMP="$(date +%Y%m%d-%H%M%S)"
  V27_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V27_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V27 z-imputer run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V27_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V27_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20261011}"
export N_REP="${N_REP:-500}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export SCENARIOS="${SCENARIOS:-mcar_z40,mar_z40,nl_mcar_z40,nl_mar_z40,zr_fork,zr_pipe,zr_collider,zr_mixed,zrmar_fork,zrmar_pipe,zr_pipenl}"
export ARMS="${ARMS:-pipeline_bartMI,pipeline_miceShipped,pipeline_micePmm}"

STAMP="${V27_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs
NLEVELS=(0 800 3200 12800)

{
  echo "=== Track V27: should z_imputer default to mice instead of bart? ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "THE LIVE PROBLEM: the RECOMMENDED config carries +5% to +10.4% under a"
  echo "  causally active covariate, and coverage does not detect it."
  echo "THE CANDIDATE: V21 attributes 100% of it to BART's flexibility; a"
  echo "  correctly specified parametric draw is +0.26%/+0.58%."
  echo "WHAT COULD STOP IT: V6 measured mice pmm at -4.70% where the covariate"
  echo "  conditional is NON-LINEAR -- which is why BART was adopted. The nl_*"
  echo "  cells are that ground, and C1 is that test."
  echo ""
  echo "REGISTERED CRITERIA"
  echo "  C1 no loss on V6's ground (nl_* and inert): mice within 0.3 bias/SE"
  echo "     of bart. THIS GATE DECIDES THE DEFAULT."
  echo "  C2 real gain in the causally-active cells: >= 1.0 bias/SE better"
  echo "     (void where bart's own |bias/SE| < 0.5)"
  echo "  C3 coverage not traded away: mice >= bart - 0.02 everywhere"
  echo "  C4 the port is faithful: miceShipped == micePmm within MC error"
  echo "  C5 n-scaling, descriptive"
  echo ""
  echo "IF C1 FAILS the default must NOT switch, and the honest outcome is a"
  echo "  conditional recommendation keyed on something users cannot always know."
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
    echo "-------------------------------------------------------------------"
  } >> "$LOG"
  OUT_TAG="${TAG}_n${N}" N_OBS="$N" Rscript run_v4_variance.R >> "$LOG" 2>&1 &
  RPID=$!
  printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
  wait "$RPID"; rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[$(date '+%H:%M:%S')] STAGE $i (n=${N}) FAILED rc=$rc -- stopping." >> "$LOG"
    break
  fi
  echo "[$(date '+%H:%M:%S')] stage $i done -> results/${TAG}_n${N}_latest.rds" >> "$LOG"
done
{ echo ""; echo "[$(date '+%H:%M:%S')] V27 finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v27.R"; } >> "$LOG"
exit $rc
