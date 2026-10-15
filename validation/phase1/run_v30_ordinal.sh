#!/usr/bin/env bash
# =============================================================================
# V30 -- does the record transfer to an ORDINAL outcome, at the DERIVED size?
# -----------------------------------------------------------------------------
# V26 settled that the Gaussian record transfers to a BINARY outcome in
# direction and ordering but at ~0.474x in bias/SE, and THEORY.md derived that
# factor as Fisher information rather than anything about imputation. V30 asks
# whether the generalisation to a cumulative-logit outcome, which has NO free
# parameters, predicts the measured attenuation.
#
# THE PREDICTION IS ALREADY VALIDATED ON COMPLETE DATA. Before launching, the
# exact form was checked against the complete-data SE ratio in all six ordinal
# cells: misses -0.0000 .. +0.0005, mixed signs. That pre-flight caught two
# real errors in the naive sqrt(E[w]) form -- the Cov(w, xtilde^2) term and the
# non-orthogonality of the ESTIMATED cutpoints, the latter a flat ~1.3% effect
# that does NOT decay with n. So this run tests transfer THROUGH THE IMPUTATION
# machinery, not the information algebra, which is the point.
# See derive_ordinal_information.R.
#
# REGISTERED (2026-10-14, before the run; also in R/robustness.R):
#   P1 LEVEL  each cell's |bias/SE| ratio to its Gaussian twin is within 0.05 of
#             K=2 fork .4459 pipe .4488 | K=3 .5091/.5120
#             K=5 .5449/.5477           | K=7 .5559/.5586
#   P2 TREND  the ratio RISES with K, in both roles.
#   P3 ANCHOR bin_fork/bin_pipe reproduce V26's 0.4731 / 0.4527 within 0.06.
#   P4 UNBALANCED at K=3 the ratio FALLS as the split skews --
#             balanced .5091/.5120 -> 70/20/10 .4306/.4327 -> 85/10/05
#             .3358/.3382 -- and the 85/10/05 cell lands BELOW its own K=2
#             cell. This is the gate that makes the formula earn its keep:
#             it moves the prediction with K HELD FIXED, which "ratios just
#             rise with K" cannot do, over a span (0.173) wider than the
#             whole balanced K sweep (0.110). All four cells were validated
#             against the complete-data SE ratio to <=0.0033, with zero polr
#             convergence failures at the 84/10/05 split.
#
#   FLOOR, registered up front at 0.5: a ratio is only computed where the
#   GAUSSIAN TWIN HAS A DEFECT (|bias/SE| > 0.5). V26 registered 0.1, found it
#   far too weak -- two near-zero arms gave 23.56x -- and corrected mid-analysis.
#   Six gates in this programme have now been spoiled by a denominator that can
#   approach zero; this one carries the floor from the start.
#
#   NOT A GATE: the fork-vs-pipe contrast at fixed K. The derivation puts the
#   two roles within 0.003, far under achievable MC error, so such a criterion
#   could not fail for the reason it names. The roles are REPLICATION of P2.
#
# ARMS are V26's exactly, so the two tracks are directly comparable. (The
# shipped default moved to mice in v1.7.0, but these are explicit harness arms,
# and changing the set would break the comparison this track exists to make.)
#
# COST -- MEASURED by a full-wave pilot (N_REP = 22 = worker count), both
# stages, never extrapolated from a part-wave:
#     n =  3200 : 220 tasks in 12.9 min   (results/v30pilot_n3200_*)
#     n = 12800 : 9.2 min per 44-task chunk, i.e. 46 min for the same wave
# A task is one (scenario, rep) across all 8 arms; 10 scenarios.
#
#     N_REP      stage 1   stage 2    TOTAL
#       300        2.9 h    10.5 h    13.4 h   <- default; comfortable overnight
#       400        3.9 h    14.0 h    17.9 h
#       500        4.9 h    17.4 h    22.3 h   (V26 used 500)
#
# N_REP = 300 was chosen against the GATE, not for speed: bias/SE carries MC
# error ~1/sqrt(N_REP), so pooled over ~7 arms the ratio lands at ~0.025 where
# the Gaussian twin sits just above the 0.5 floor, and much tighter where it
# does not. Against a +/-0.05 gate that is ~2 sigma -- enough that P1 CAN fail,
# which is the point, without being so tight that noise fails it.
#
#   bash run_v30_ordinal.sh              # detach; logs/v30_<stamp>.{log,pid}
#   STAGES=1 bash run_v30_ordinal.sh     # n=3200 only
#   bash run_status.sh v30               # is it alive?
#   Rscript analyze_v30.R                # after it lands
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
TAG=v30
STAMP="$(date '+%Y%m%d-%H%M%S')"
mkdir -p logs results

if [[ -z "${V30_DETACHED:-}" ]]; then
  V30_DETACHED=1 V30_STAMP="$STAMP" nohup bash "$SELF" "$@" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V30 detached (pid $PID)"
  echo "  log:    logs/${TAG}_${STAMP}.log"
  echo "  status: bash run_status.sh ${TAG}"
  echo "  then:   Rscript analyze_v30.R"
  exit 0
fi

STAMP="${V30_STAMP:-$STAMP}"
LOG="logs/${TAG}_${STAMP}.log"
export NCORES="${NCORES:-$(( $(sysctl -n hw.ncpu) - 2 ))}"
export N_REP="${N_REP:-300}"
export M="${M:-30}"
export SCENARIOS="${SCENARIOS:-ord3_fork,ord5_fork,ord7_fork,ord3_pipe,ord5_pipe,ord7_pipe,ord3u_fork,ord3v_fork,ord3u_pipe,ord3v_pipe,bin_fork,bin_pipe,zr_fork,zr_pipe}"
export ARMS="${ARMS:-complete_case,leftcens_prestep,pipeline_bartMI,pipeline_noYx,pipeline_noYz,pipeline_properZ,pipeline_properZ_ds}"
NLEVELS=(0 3200 12800)

{ echo "[$(date '+%H:%M:%S')] V30 start -- ordinal transfer"
  echo "N_REP=${N_REP} M=${M} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
} > "$LOG"

WANT="${STAGES:-1,2}"
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
{ echo ""; echo "[$(date '+%H:%M:%S')] V30 finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v30.R"; } >> "$LOG"
exit $rc
