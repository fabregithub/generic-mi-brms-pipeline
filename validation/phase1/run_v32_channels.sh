#!/usr/bin/env bash
# =============================================================================
# V32 -- is the borrowing scale a FAMILY constant, or arm-dependent?
# -----------------------------------------------------------------------------
# V31 localised outcome-family bias transfer to one channel: fitting
# `measured = a*borrow + b*atten`, the attenuation scale `b` is 0.94-1.03 under
# gaussian, binomial and ordinal alike, while the borrowing scale `a` falls to
# 0.63 and 0.41. Under a Gaussian outcome the whole instrument fits
# (chi2/df = 1.1; 0.1 at the cleanest rung).
#
# WHAT V31 COULD NOT SETTLE. Whether `a` is a property of the FAMILY -- in which
# case a transfer law exists and can be predicted -- or of the ARM, in which
# case there is no law to find. V31 had only 4 arms, and its pooled fit worked
# at the clean rung (chi2/df 0.5 / 1.2) while failing at the dirty one
# (115.8 / 79.4), which is consistent with either reading.
#
# THE DESIGN. Six instrumented arms instead of four, chosen to SPAN the
# borrow/atten mix rather than to repeat it: `noYx` drops Y from the X block and
# `noYz` from the Z block, which move the two channels in different directions
# (V31 measured surviving fractions |T1+T2|/(|T1|+|T2|) of 0.73 for noYx against
# 0.19 for properZ_ds). Six arms x 2 fractions = 12 points per family against 2
# parameters, so a single (a,b) per family is a real constraint rather than a
# fit.
#
# ONLY THE CLEAN RUNGS (0.10, 0.20). V31 showed the expansion itself degrades at
# 0.40 (chi2/df 115.8 binomial, 79.4 ordinal) -- `a` drifts with |d| there
# because g_eta_eta is non-zero for a GLM and identically zero for a Gaussian.
# Including 0.40 would re-measure that degradation, not the question asked.
#
# REGISTERED (2026-10-15, before the run):
#   P1 a single (a, b) per family fits the 10 points at chi2/df <= 2
#      -> `a` is a FAMILY constant and a transfer law exists.
#   P2 adding an arm-specific scale does not significantly improve the fit
#      (likelihood-ratio p > 0.05).
#   If P1 fails while the Gaussian cells still fit, `a` is arm-dependent and
#   there is no single transfer constant to find -- which is equally decisive
#   and would close this line rather than extend it.
#
#   NOT A GATE: `pred/meas` on any single cell x arm. V31 showed that statistic
#   is dominated by a fixed cancellation amplification (x1.4-x3.8) that says
#   nothing about the channels. The gates are on the CHANNEL SCALES.
#
# COST -- MEASURED by a timed FULL-WAVE pilot with the final arm set
# (6 cells x 22 reps x 6 arms = 132 tasks in 432 s = 3.27 s/task):
#     6 cells x 300 reps -> 1.6 h      6 cells x 500 reps -> 2.7 h
#     6 cells x 400 reps -> 2.2 h      6 cells x 600 reps -> 3.3 h
#
# ⚠️ HOW THE FIRST ATTEMPT WENT WRONG, because the cause is reusable.
# The 2026-10-15 launch quoted 1.6 h and was running at ~21 h when killed --
# out by 13x. `pipeline_block_fcs` costs 942 s per wave against 10-22 s for
# every other arm (45-94x), and alone was ~94% of the total: it is the one arm
# that runs the SHIPPED pipeline path (mice Z-block, m = 30, 3 sweeps, ~43 s
# per rep) rather than a harness arm. The estimate had scaled V31's per-task
# cost by the ARM COUNT, which silently assumes arms cost alike.
# The evidence was already in hand and was not read: an 18-task smoke test
# blew through a 600 s timeout, implying >0.5 s... >0.5 MIN per task. It was
# treated as pass/fail without looking at its clock.
# CLAUDE.md requires a timed full-wave pilot for exactly this reason.
#
# `block_fcs` IS DROPPED, and costs the design nothing: it was there to supply
# a well-conditioned arm, and `noYz` already sits at 1.00. The remaining five
# span 0.10-1.00 in |T1+T2|/(|T1|+|T2|), a 10x range in conditioning.
#
#   bash run_v32_channels.sh        # detach; logs/v32_<stamp>.{log,pid}
#   bash run_status.sh v32
#   Rscript analyze_v32.R
# =============================================================================
set -uo pipefail
# Resolve to an ABSOLUTE path BEFORE the cd, and re-invoke through it. Using
# "$0" here is a trap: the cd below changes the working directory, so a relative
# "$0" no longer resolves and the detached child dies instantly -- while the
# parent still prints "detached (pid N)" and writes a .pid, so the launch looks
# successful and nothing runs. Latent until someone launches from the repo root
# instead of from phase1/ (found 2026-10-15).
SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
cd "$(dirname "$SELF")"
TAG=v32
STAMP="$(date '+%Y%m%d-%H%M%S')"
mkdir -p logs results

if [[ -z "${V32_DETACHED:-}" ]]; then
  V32_DETACHED=1 V32_STAMP="$STAMP" nohup bash "$SELF" "$@" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V32 detached (pid $PID)"
  echo "  log:    logs/${TAG}_${STAMP}.log"
  echo "  status: bash run_status.sh ${TAG}"
  echo "  then:   Rscript analyze_v32.R"
  exit 0
fi

STAMP="${V32_STAMP:-$STAMP}"
LOG="logs/${TAG}_${STAMP}.log"
export NCORES="${NCORES:-$(( $(sysctl -n hw.ncpu) - 2 ))}"
export N_REP="${N_REP:-300}"
export RECORD_D="${RECORD_D:-1}"
export M="${M:-30}"
export SCENARIOS="${SCENARIOS:-fr10_gau,fr20_gau,fr10_bin,fr20_bin,fr10_or3,fr20_or3}"
export ARMS="${ARMS:-pipeline_properZ,pipeline_properZ_ds,pipeline_bartMI,pipeline_noYx,pipeline_noYz}"
NLEVELS=(0 3200 12800)

{ echo "[$(date '+%H:%M:%S')] V32 start -- channel scales, recording d (RECORD_D=${RECORD_D})"
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
{ echo ""; echo "[$(date '+%H:%M:%S')] V32 finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v32.R"; } >> "$LOG"
exit $rc
