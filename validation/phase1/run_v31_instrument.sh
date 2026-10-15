#!/usr/bin/env bash
# =============================================================================
# V31 -- does the RECORDED imputation error `d` predict the bias it causes?
# -----------------------------------------------------------------------------
# V30 closed the standard-error half of outcome-family transfer exactly and left
# the bias half open: a defect's ABSOLUTE bias is not family-invariant
# (0.26x-2.2x), and both derivable reductions to a function of the information
# ratio were REJECTED (FINDINGS_v30 §4a, worst chi2/df 73.6 and 101.1). The
# failing step depends on the imputation error distribution `d`, which is not
# recoverable from summary statistics -- so `d` is now recorded directly
# (R/instrument_d.R, RECORD_D=1).
#
# THE PREDICTION, per (cell, arm):
#     [I_bb I_bt; I_tb I_tt] [b_b; b_t]
#         = [E[d g]; 0] - E[J (d'beta0)] - [E[w d d'] beta0; 0]
# and `pred_b` (its focal component) should equal the MEASURED bias
# `mean_est - estimand_true` in the same row.
#
# SELF-TESTED BEFORE WIRING, on contamination of KNOWN size, which caught a real
# omission: a FIRST-order expansion predicts the borrowing channel well
# (ratio 1.01-1.05) and the attenuation channel not at all (0.003-0.15, sign
# flipped under binomial), because for an error with E[d|x] = 0 attenuation is
# exactly what first order kills. With the second-order term `-E[w d d']beta0`
# restored, predicted/actual is 0.99-1.10 at |d| <= 0.15 and <=1.21 at |d| = 0.40
# across gaussian, binomial and ordinal, both channels.
#
# ⚠️ REAL |d| IS LARGE. The arms run d_rms ~ 0.53, at or past where the self-test
# showed the expansion straining. If the prediction misses, check `d_rms` before
# concluding the structure is wrong -- an inadequate expansion order and a wrong
# mechanism look the same in the ratio.
#
# ⚠️ SCOPE CHANGED BEFORE LAUNCH, on the pilot's evidence. The first design ran
# all 14 V30 cells to see whether `pred_b` matches. The pilot (4 cells x 120
# reps, results/dcheck_*) already answers that: it does NOT -- median
# |pred/meas - 1| = 0.59 -- and more reps cannot fix a biased predictor. The
# run instead sweeps the MISSING FRACTION (0.10 / 0.20 / 0.40) across gaussian,
# binomial and ordinal, which distinguishes "too low-order" from "wrong
# structure" and can fail either way.
#
# REGISTERED (2026-10-15, before the run):
#   P1 pred/meas RISES toward 1 as the missing fraction falls, in ALL THREE
#      families. That is E1 -- the expansion is right and needs higher order at
#      realistic |d|. If it stays off 1, or moves inconsistently across
#      families, that is E2 and the structure itself is wrong.
#   P1b at the 0.10 rung the ratio should be within 25% of 1, since d_rms there
#      is ~a quarter of the 0.40 rung's and the self-test holds to 1.10 at
#      |d| <= 0.15.
#   P0 pred_b matches the measured bias within 25% wherever the measured bias
#      is large enough to divide by (|measured| > 2 x its MC error -- a FLOOR,
#      because seven gates in this programme have been spoiled by a denominator
#      that can reach zero).
#   P2 the borrow/atten SPLIT explains the family dependence: an arm's
#      bias_fam/bias_gauss tracks its borrow share, since only the borrowing
#      term carries the outcome family's weights.
#   NOT A GATE: any single cell-arm ratio where the measured bias is near zero.
#
# ARMS: only the six that route through the instrumented core. `oracle`,
# `complete_case` and `leftcens_prestep` take other code paths and return NA for
# `pred_b` by design -- they are kept for context, not for the gates.
#
# Cost: MEASURED, 4 cells x 120 reps x 5 arms at n = 3200 ran 17.8 min
# (results/dcheck_*). 9 cells x 300 reps x 5 arms scales to ~1.7 h.
#
#   bash run_v31_instrument.sh        # detach; logs/v31_<stamp>.{log,pid}
#   bash run_status.sh v31
#   Rscript analyze_v31.R
# =============================================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
TAG=v31
STAMP="$(date '+%Y%m%d-%H%M%S')"
mkdir -p logs results

if [[ -z "${V31_DETACHED:-}" ]]; then
  V31_DETACHED=1 V31_STAMP="$STAMP" nohup bash "$0" "$@" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V31 detached (pid $PID)"
  echo "  log:    logs/${TAG}_${STAMP}.log"
  echo "  status: bash run_status.sh ${TAG}"
  echo "  then:   Rscript analyze_v31.R"
  exit 0
fi

STAMP="${V31_STAMP:-$STAMP}"
LOG="logs/${TAG}_${STAMP}.log"
export NCORES="${NCORES:-$(( $(sysctl -n hw.ncpu) - 2 ))}"
export N_REP="${N_REP:-300}"
export RECORD_D="${RECORD_D:-1}"
export M="${M:-30}"
export SCENARIOS="${SCENARIOS:-fr10_gau,fr20_gau,zr_fork,fr10_bin,fr20_bin,bin_fork,fr10_or3,fr20_or3,ord3_fork}"
export ARMS="${ARMS:-pipeline_properZ,pipeline_properZ_ds,pipeline_bartMI,pipeline_noYx}"
NLEVELS=(0 3200 12800)

{ echo "[$(date '+%H:%M:%S')] V31 start -- recording d (RECORD_D=${RECORD_D})"
  echo "N_REP=${N_REP} M=${M} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
} > "$LOG"

WANT="${STAGES:-1}"   # n = 3200 only; the instrument tests a mechanism, and
                      # reps buy more here than a second n does
rc=0
for i in 1 2; do
  case ",$WANT," in *",$i,"*) ;; *) continue ;; esac
  N="${NLEVELS[$i]}"
  { echo ""; echo "-----------------------------------------------------------"
    echo "[$(date '+%H:%M:%S')] STAGE $i : n = ${N}"
    echo "-----------------------------------------------------------"; } >> "$LOG"
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
{ echo ""; echo "[$(date '+%H:%M:%S')] V31 finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v31.R"; } >> "$LOG"
exit $rc
