#!/usr/bin/env bash
# =============================================================================
# Track V23 -- the censored-exposure defect, DERIVED FIRST
# -----------------------------------------------------------------------------
# WHY THIS ONE IS DIFFERENT. V22 found ~20% asymptotic bias in the censored
# exposure draw by SURPRISE -- it was built to answer a covariate-draw question
# and this fell out of a secondary cell. That is the pattern the cycle adopted at
# V15 was meant to replace, and it recurred. So before building anything, the
# mechanism was derived and probed:
#
#   1. WHY it must exist. `leftcens` draws a censored exposure from a conditional
#      LINEAR in its predictors. The true conditional for a censored logX1
#      contains p(Z1 | logX1) = N(Z1; g(logX1), s^2), whose log contributes
#      -(Z1 - g(logX1))^2 / 2s^2 -- non-linear in logX1 wherever g is curved.
#      A linear-Gaussian draw cannot represent it. Derivable from reading what
#      leftcens conditions on; no simulation needed.
#
#   2. WHY it is asymptotic. The fitted coefficients converge to the best LINEAR
#      approximation of a non-linear target. That approximation stays wrong at
#      every n -- which is exactly what V22 measured (+20.8/+20.3/+20.3% over a
#      16x range of n) and what a derivation would have predicted.
#
#   3. WHAT governs it. Not "how big is g" but how much of g a straight line
#      CANNOT express where the censored mass actually sits: the density-weighted
#      residual SD of g after its best linear fit below the LOD. Call it u --
#      `ef_unrep_curvature()` in dgp.R. A 5-setting, 40-rep probe (~4 min)
#      confirmed the direction: monotone-g settings and small-amplitude settings
#      both collapse toward zero, and only settings with substantial
#      UNREPRESENTABLE curvature are large.
#
#   4. WHY IT SHOULD ENTER SQUARED. u is by construction ORTHOGONAL, in the
#      density-weighted L2 sense, to the span of the linear predictors. A
#      first-order expansion of the bias functional pairs the misspecification
#      with the score, and orthogonality kills that term -- so the leading
#      contribution is second order. **This is the first time this project has
#      had a REASON for one of its quadratic laws.** V15's f^2 was frank curve
#      fitting; this one has an argument behind it, and if it holds the same
#      argument may explain V15's.
#
# THE REGISTERED LAW
#
#     X-block bias% = 300 * u^2
#
#   calibrated on ONE point (cv262, V22's own arrow, u = 0.2621, probe bias
#   +22.26% at 40 reps). Everything else is a prediction:
#
#     cell    a     c      u        predicted
#     cv000  0.0  0.00   0.0000       0.00%   <- the null: no arrow, nothing to miss
#     cv040  0.4  0.00   0.0399       0.48%
#     cv119  0.8  0.10   0.1190       4.25%
#     cv153  0.0  0.35   0.1531       7.04%
#     cv182  0.8  0.25   0.1818       9.91%
#     cv262  1.2  0.35   0.2621      20.60%   <- calibration
#     cv364  1.6  0.50   0.3636      39.66%
#
#   The 40-rep probe already disagrees with the law at two points (u = 0.1195
#   measured +0.09% against 4.3% predicted; u = 0.1531 measured +4.24% against
#   7.0%), each about 2 MC se away. **That is the main thing this run settles**:
#   whether a one-parameter u^2 law survives at 500 reps or whether the residuals
#   are real and the account needs a second term.
#
# TWO CONVENTIONS THIS RUN ADOPTS (PLAN §11, changed 2026-09-08)
#   1. THE BAR IS IN bias/SE, NOT RELATIVE BIAS. Relative bias decays as n^-1/3
#      while bias/SE GROWS as n^0.17 (V18), so a relative bar gets easier to
#      pass exactly as the bias becomes more consequential. At n = 12800 the old
#      10% relative bar permits a bias worth 95% of a 2.8-SE detectable effect.
#      Default gate: **bias/SE <= 0.3**, about 10% of a detectable effect. The
#      target effect size is a CHOICE -- 5% is an epidemiological convention, and
#      a no-threshold setting has no such floor -- so bias/SE is reported and the
#      conversion is the reader's.
#   2. DIRECTION IS REPORTED, NOT JUST MAGNITUDE. Away from the null inflates an
#      exposure-response (anti-conservative); toward the null risks missing one
#      (conservative); a sign that FLIPS with a design parameter is worse than
#      either, and V15 found exactly that. The probe says this defect is AWAY
#      from the null (+20.8% at u = 0.262), i.e. the dangerous direction, and
#      whether that sign is stable across the u sweep is a registered question.
#
# FALSIFICATION, PRE-DECLARED
#   * any predicted cell off by more than 5 pp -> the one-parameter law fails
#   * the SIGN is not positive in every cell with u > 0 -> the direction is not
#     stable in u, and no single caveat covers the defect
#   * log|bias| on log u: slope 95% CI excluding 2 -> the quadratic form fails,
#     and with it the orthogonality argument
#   * cv000 above 1.5 pp -> the null leaks and the whole predictor is suspect
#   * bias not flat in n (800 vs 3200) -> it is not the asymptotic
#     misspecification the derivation describes
#
# ARMS, and one of them is a warning
#   z21_exact         exact covariate draw + shipped leftcens -> ISOLATES the X
#                     block, which is what the law is about
#   pipeline_bartMI   the shipped configuration -- V22's secondary cell had NO
#                     shipped baseline, which this fixes
#   pipeline_micePmm  a 22-rep probe put it 10 pp better than the default here
#   smc_xgrid         ITEM 07's grid exposure draw. The same probe put it at
#                     +93.7% with coverage 0.000 -- INVERTED. Mechanism: it
#                     proposes from the exposure prior and reweights by
#                     p(Y | x, rest), so it never uses the Z1 = g(logX1)
#                     information at all, discarding the linear-but-partly-right
#                     covariate conditioning leftcens does use. **Item 07 is the
#                     only open build track. It must not ship without this cell
#                     in its acceptance set**, and confirming or clearing that
#                     inversion at 500 reps is a first-class goal of this run.
#
# COST -- MEASURED (full-wave pilot, 22 reps per cell per stage, 308 tasks, 0
# errors, at the run's own settings m = 30, sweeps = 3, 22 workers):
#
#   n =  800   29.3 s/task ->  78 min
#   n = 3200   79.0        -> 209 min
#   TOTAL ~4.8 h
#
# STAGES=1 costs 78 min and settles the LAW, which is a set of levels; the second
# stage only confirms flatness in n -- and V22 already measured that over a 16x
# range in one of these cells.
#
# USAGE
#   ./run_v23_curvature.sh                  # both n levels (~4.8 h)
#   STAGES=1 ./run_v23_curvature.sh         # n = 800 only
#   N_REP=22 ./run_v23_curvature.sh         # full-wave timing pilot
#
# Results: results/v23_n<N>_latest.rds per stage.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v23}"

if [[ "${V23_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V23_STAMP="$(date +%Y%m%d-%H%M%S)"
  V23_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V23_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V23 curvature-sweep run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V23_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V23_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V23_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260909}"
export N_REP="${N_REP:-500}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export SCENARIOS="${SCENARIOS:-cv000,cv040,cv119,cv153,cv182,cv262,cv364}"
export ARMS="${ARMS:-z21_exact,pipeline_bartMI,pipeline_micePmm,smc_xgrid}"

STAMP="${V23_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 3200)

{
  echo "=== Track V23: the censored-exposure defect, derived first ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "THE DERIVATION. leftcens draws a censored exposure from a conditional"
  echo "  LINEAR in its predictors. The true conditional contains"
  echo "  p(Z1 | logX1) = N(Z1; g(logX1), s^2), non-linear in logX1 wherever g"
  echo "  is curved. So what governs the damage is u = how much of g a straight"
  echo "  line cannot express where the censored mass sits -- and u should enter"
  echo "  SQUARED, because it is orthogonal to the linear predictors' span and"
  echo "  the first-order term in the bias expansion vanishes."
  echo ""
  echo "REGISTERED LAW:  X-block bias% = 300 * u^2"
  echo "  cv000 u=0.0000 ->  0.00% (the null)   cv182 u=0.1818 ->  9.91%"
  echo "  cv040 u=0.0399 ->  0.48%              cv262 u=0.2621 -> 20.60% (calib)"
  echo "  cv119 u=0.1190 ->  4.25%              cv364 u=0.3636 -> 39.66%"
  echo "  cv153 u=0.1531 ->  7.04%"
  echo ""
  echo "A 40-rep probe already disagrees at two points (u=0.1195 gave +0.09% vs"
  echo "  4.3% predicted; u=0.1531 gave +4.24% vs 7.0%), each ~2 MC se. Whether"
  echo "  a ONE-PARAMETER u^2 law survives at 500 reps is what this settles."
  echo ""
  echo "BAR IS IN bias/SE (PLAN 11, changed 2026-09-08): default gate 0.3, about"
  echo "  10% of a 2.8-SE detectable effect. Relative bias is the wrong currency --"
  echo "  it decays as n^-1/3 while bias/SE GROWS as n^0.17."
  echo "DIRECTION IS REPORTED: the probe puts this defect AWAY from the null"
  echo "  (+20.8% at u=0.262) -- the anti-conservative direction. Whether that"
  echo "  sign holds across the u sweep is a registered question."
  echo ""
  echo "REJECTED IF: any predicted cell off by >5 pp; the sign is not positive in"
  echo "  every cell with u > 0; or the log-log slope CI"
  echo "  excludes 2 (which would also kill the orthogonality argument); or"
  echo "  cv000 exceeds 1.5 pp; or the bias is not flat between n=800 and 3200."
  echo ""
  echo "WARNING ARM: smc_xgrid is ITEM 07's grid exposure draw. A 22-rep probe"
  echo "  put it at +93.7% with coverage 0.000 in this cell -- INVERTED, because"
  echo "  it reweights by p(Y|x) and never uses the Z1 = g(logX1) information."
  echo "  Item 07 is the only open build track; confirming or clearing that is a"
  echo "  first-class goal here."
  echo "==========================================================="
} > "$LOG"

WANT="${STAGES:-1,2}"
rc=0
for i in 1 2; do
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
