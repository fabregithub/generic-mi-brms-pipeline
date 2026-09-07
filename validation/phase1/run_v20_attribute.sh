#!/usr/bin/env bash
# =============================================================================
# Track V20 -- WHICH BLOCK carries the bias under a causally active covariate?
# -----------------------------------------------------------------------------
# THE QUESTION, and it is the main one left open. V17 found the shipped default
# biased +5% to +10.4% wherever the covariate is a confounder or a mediator. V18
# measured how that bias scales (n^-1/3, permanent). V19 ruled out the account
# that it was the imputer family's convergence rate -- and found BART is the only
# Z-block imputer whose bias decays at all. **Nothing yet says whether the bias
# comes from the X block, the Z block, or their alternation.** Every fix depends
# on that answer, so it is the step before any code change.
#
# THE DESIGN: a 2x2, swapping each block's draw for the TRUE conditional.
#
#                          X block: shipped (leftcens)   X block: exact
#   Z block: shipped BART      ef_bart_ship  (CONTROL)     ef_bart_exact
#   Z block: exact             ef_exact_ship               ef_exact_exact
#
#   fix-Z effect   = ef_exact_ship  - ef_bart_ship
#   fix-X effect   = ef_bart_exact  - ef_bart_ship
#   both           = ef_exact_exact - ef_bart_ship
#   interaction    = both - (fix-Z + fix-X)      <- the alternation's share
#
# `ef_bart_ship` IS THE CONTROL AND IT DECIDES WHETHER ANYTHING ELSE COUNTS. It
# is the instrument configured to do what the pipeline does, so it must reproduce
# `pipeline_bartMI`. If it does not, the instrument is not the pipeline and the
# other three cells attribute nothing. This is the V8 port check applied to a new
# instrument.
#
# WHY THE EXACT DRAWS NEEDED NEW CODE (R/exact_fork.R). V12/V13's exact Z draw
# computes p(Z1 | Y) ~ p(Y | Z1, rest) * N(Z1 | 0, 1) -- correct in the PRECISION
# structure every track before V17 used, where Z1 is exogenous. **Under a fork it
# is wrong**: Z1 -> logX1 with delta = 0.60, so logX1 is a child of Z1 and carries
# information about it that the missing p(logX1 | Z1) factor throws away.
# Measured: its draws are 32% too wide around the true conditional mean, while
# their marginal SD (0.998 against a true 0.997) looks perfectly fine. Reusing it
# would have produced a believable decomposition from a broken "exact" arm --
# the leftcens failure mode from FINDINGS_v13.md, again.
#
# The replacement conditions the JOINT GAUSSIAN of (Z1, logX, Y) in closed form,
# which is exact rather than approximate, and is checked against known answers in
# test_v20_exact.R (20 assertions: the analytic conditional mean must equal the
# population regression and the analytic SD its residual SD, because a Gaussian
# joint makes those identities).
#
# THE REGISTERED PREDICTION: the Z BLOCK dominates.
#   Under a fork the X block's true conditional IS linear-Gaussian, and leftcens
#   draws from a linear conditional with a flexible (shash) margin -- so it is
#   approximately correctly specified and should carry little. The Z block's true
#   conditional is also linear-Gaussian, but BART approximates it
#   NONPARAMETRICALLY and pays a smoothing error -- and n^-1/3 is exactly the rate
#   V18 measured for the estimand. So:
#
#     fix-Z removes >= 70% of the reference bias
#     fix-X removes <= 30%
#     |interaction| <= 2 pp
#     ef_exact_exact within 1.5 pp of ZERO (the true conditional throughout)
#
# THE SHARPER TEST, and why n is swept. If the bias is the Z draw's smoothing
# error, then **the n^-1/3 exponent should live in the Z block and vanish once it
# is fixed**: ef_exact_ship's bias should be flat in n (already ~0), while
# ef_bart_exact's should still decay at ~-1/3. Three levels make that visible.
#
# WHAT WOULD REFUTE IT. If fix-X carries most of it, the X block's linear
# conditional is at fault even where it is correctly specified, and item 07's
# grid draw becomes the fix. If neither block alone carries it but both together
# do, the alternation is the cause and no single-block fix will work -- which is
# V13's lesson repeating in a structure V13 could not see.
#
# NOTE (PLAN §10b): arms consume the task RNG stream in order, so V17/V18/V19
# figures are NOT the baseline. `pipeline_bartMI` is re-run inside this arm set.
#
# USAGE
#   ./run_v20_attribute.sh                 # all three n levels
#   STAGES=1 ./run_v20_attribute.sh        # n = 800 only
#   N_REP=22 ./run_v20_attribute.sh        # full-wave timing pilot
#
# Results: results/v20_n<N>_latest.rds per stage.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v20}"

if [[ "${V20_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V20_STAMP="$(date +%Y%m%d-%H%M%S)"
  V20_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V20_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V20 block-attribution run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V20_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V20_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V20_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260907}"
export N_REP="${N_REP:-500}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export SCENARIOS="${SCENARIOS:-zr_fork,zr_pipe}"
export ARMS="${ARMS:-pipeline_bartMI,ef_bart_ship,ef_exact_ship,ef_bart_exact,ef_exact_exact}"

STAMP="${V20_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 3200 12800)

{
  echo "=== Track V20: which block carries the bias? ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "2x2, swapping each block for the TRUE conditional:"
  echo "                      X shipped        X exact"
  echo "  Z shipped (BART)    ef_bart_ship     ef_bart_exact"
  echo "  Z exact             ef_exact_ship    ef_exact_exact"
  echo ""
  echo "CONTROL FIRST: ef_bart_ship must reproduce pipeline_bartMI (within"
  echo "  1.5 pp). If it does not, the instrument is not the pipeline and the"
  echo "  other three cells attribute nothing."
  echo ""
  echo "REGISTERED PREDICTION -- the Z BLOCK dominates:"
  echo "  fix-Z removes >= 70% of the reference bias"
  echo "  fix-X removes <= 30%"
  echo "  |interaction| <= 2 pp"
  echo "  ef_exact_exact within 1.5 pp of ZERO"
  echo "  Reasoning: under a fork the X block's true conditional IS"
  echo "  linear-Gaussian and leftcens is approximately correct there, while"
  echo "  BART approximates a linear truth nonparametrically and pays a"
  echo "  smoothing error at the n^-1/3 rate V18 measured for the estimand."
  echo ""
  echo "SHARPER TEST across n: if the bias is the Z draw's smoothing error, the"
  echo "  n^-1/3 exponent should VANISH once Z is fixed and survive when only X"
  echo "  is fixed."
  echo ""
  echo "REFUTED IF: fix-X carries most of it (then item 07's grid draw is the"
  echo "  fix), or neither block alone carries it but both together do (then the"
  echo "  alternation is the cause and no single-block fix works)."
  echo "NOTE: V17-V19 numbers are NOT the baseline -- different arm set, so a"
  echo "  different RNG stream (PLAN §10b). pipeline_bartMI is re-run here."
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

echo "" >> "$LOG"
echo "[$(date '+%H:%M:%S')] sequence finished rc=$rc" >> "$LOG"
