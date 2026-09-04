#!/usr/bin/env bash
# =============================================================================
# Track V19 -- does the bias exponent track the Z-block IMPUTER?
# -----------------------------------------------------------------------------
# THE CONJECTURE. V18 measured the shipped default's bias under a causally
# active covariate decaying as n^-0.335 (CI [-0.378, -0.291]) -- rejecting both
# the constant-bias account (slope 0) and the finite-sample one (-0.5).
# n^(-1/3) is the classical NONPARAMETRIC convergence rate for a Lipschitz
# function in one dimension, and the Z block is BART. So: the bias may simply be
# the Z-block imputer's own convergence rate, showing up in the estimand.
#
# WHY IT IS TESTABLE HERE, AND CHEAPLY. In the `fork` and `pipe` cells the true
# Z1 conditional is LINEAR-GAUSSIAN (Z1 ~ N(0,1); X1 and Y are linear in it), so
# a PARAMETRIC imputer is correctly specified there and should converge at the
# parametric rate. Different imputer families therefore predict visibly
# different curves on the SAME data:
#
#   arm                   imputer family                    predicted slope
#   pipeline_bartMI       BART, nonparametric               ~ -1/3   (measured)
#   pipeline_properBoot   bootstrapped random forest,       ~ -1/3
#                         also nonparametric
#   pipeline_micePmm      mice pmm / logreg, parametric     ~ -1/2
#                         and CORRECT in these cells
#   pipeline_properZ      linear/logistic posterior draw,   ~ -1/2
#                         parametric and correct
#
# TWO OF EACH FAMILY, on purpose. One of each would leave any difference
# attributable to the specific implementation; two of each makes the split a
# property of the FAMILY or not at all. V18 also measured bartMI's exponent
# independently (-0.335, CI [-0.378, -0.291]) under a different arm set, so the
# nonparametric side carries a replication the parametric side does not.
#
# The arms sit at very different bias LEVELS (a 4-rep smoke put them at +9.2,
# +4.9, +4.0 and -20.3%), which is the point: the conjecture is about the
# EXPONENT, not the magnitude, and an exponent shared by arms that disagree on
# level is a much stronger result than one fitted to a single curve.
#
# COST, AND WHY THE STAGE ORDER MATTERS. Measured per arm at n = 800, 22 reps,
# one cell, 22 workers: bartMI 9 s, micePmm 9 s, properZ 5 s, and
# **properBoot 156 s** -- 17x the next arm and ~80% of the four-arm total. It is
# kept because a weekend is available; on a weekday budget, dropping it cuts the
# run by ~8x and costs only the second nonparametric arm.
#
# Stages run in ASCENDING n and each writes its own results file, so an overrun
# is graceful: whatever levels have finished still support a slope fit, with a
# wider CI. There is no need to guess the total right.
#
# WHAT A NULL WOULD MEAN, and it is not nothing. The X block is `leftcens` in
# every arm, so its contribution is common to all four. If every arm shows the
# same exponent, the bias is NOT the Z imputer's rate -- it is coming from the X
# block or from the alternation, and the next step is a V13-style decomposition
# rather than more imputers. Either outcome moves the theory.
#
# FALSIFICATION, PRE-DECLARED
#   * PRIMARY: pooled slope per arm, log|rel bias| on log n.
#     The conjecture is REJECTED if the two parametric arms' slopes are not
#     separated from the two nonparametric arms' by more than their combined
#     95% CIs, i.e. if the family split is not visible.
#   * It is also rejected if micePmm or properZ lands within 0.08 of -1/3, or if
#     bartMI or properBoot lands within 0.08 of -1/2 -- each arm has to sit near
#     the rate its family predicts, not merely on the correct side.
#   * CONTROL: `oracle` unbiased (|rel bias| < 2%) at every n and cell. If it
#     drifts, nothing else is readable.
#   * CONTROL: `mcar_z40`, the causally inert cell, is carried at every n. Its
#     bias is small (-2.0% to -0.7% in V18), so its exponent is poorly
#     determined -- it is here to show whether the family split appears where
#     there is no confounding path, NOT as a second test of the conjecture.
#
# ARM SET NOTE (PLAN §10b): procedures consume the task RNG stream in order, so
# an arm's numbers are only comparable across runs with an identical arm set.
# V18's bartMI figures are therefore NOT the baseline -- bartMI is re-run here
# inside this arm set, and every comparison is within this run.
#
# USAGE
# COST -- from MEASURED per-arm times at n = 800 plus a MEASURED stage-time
# exponent of 1.03 over n = 800/1600/3200 (linear in n; V18's 1-arm sweep was
# 0.82). Extrapolated to the top two levels, which is an extrapolation and is
# flagged as one:
#
#   n =   800   ~2.0 h      n = 6400   ~16 h
#   n =  1600   ~4.1 h      n = 12800  ~33 h
#   n =  3200   ~8.1 h      TOTAL     ~63 h   (3 cells x 300 reps, 4 arms)
#
# That fits a Friday-morning start with little slack. STAGES=1,2,3,4 stops at
# n = 6400 for ~30 h and still spans 8x. N_REP=200 cuts the total by a third.
#
# USAGE
#   ./run_v19_imputer_rate.sh              # all five n levels (~63 h)
#   STAGES=1,2,3,4 ./run_v19_imputer_rate.sh   # stop at n = 6400
#   N_REP=22 ./run_v19_imputer_rate.sh     # full-wave timing pilot
#
# Results: results/v19_n<N>_latest.rds per stage.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v19}"

if [[ "${V19_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V19_STAMP="$(date +%Y%m%d-%H%M%S)"
  V19_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V19_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V19 imputer-rate run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V19_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V19_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V19_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260904}"
export N_REP="${N_REP:-300}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export SCENARIOS="${SCENARIOS:-zr_fork,zr_pipe,mcar_z40}"
export ARMS="${ARMS:-pipeline_bartMI,pipeline_properBoot,pipeline_micePmm,pipeline_properZ}"

STAMP="${V19_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 1600 3200 6400 12800)

{
  echo "=== Track V19: does the bias exponent track the Z-block imputer? ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "CONJECTURE: V18's n^-1/3 is the Z-block imputer's own convergence rate."
  echo "  In zr_fork and zr_pipe the true Z1 conditional is LINEAR-GAUSSIAN, so"
  echo "  a parametric imputer is CORRECTLY SPECIFIED and should beat -1/3."
  echo ""
  echo "  arm                  family            predicted slope"
  echo "  pipeline_bartMI      nonparametric        ~ -1/3  (V18 measured)"
  echo "  pipeline_properBoot  nonparametric        ~ -1/3"
  echo "  pipeline_micePmm     parametric, correct  ~ -1/2"
  echo "  pipeline_properZ     parametric, correct  ~ -1/2"
  echo ""
  echo "REJECTED IF: the two families' pooled slopes are not separated beyond"
  echo "  their combined 95% CIs; or a parametric arm lands within 0.08 of -1/3;"
  echo "  or a nonparametric arm lands within 0.08 of -1/2."
  echo "A NULL IS INFORMATIVE: the X block (leftcens) is common to all four arms,"
  echo "  so one shared exponent says the bias is NOT the Z imputer's rate and"
  echo "  the next step is a V13-style decomposition, not more imputers."
  echo "NOTE: V18's bartMI numbers are NOT the baseline -- different arm set,"
  echo "  different RNG stream (PLAN §10b). bartMI is re-run here."
  echo "==========================================================="
} > "$LOG"

WANT="${STAGES:-1,2,3,4,5}"
rc=0
for i in 1 2 3 4 5; do
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

echo "" >> "$LOG"
echo "[$(date '+%H:%M:%S')] sequence finished rc=$rc" >> "$LOG"
