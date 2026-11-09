#!/usr/bin/env bash
# =============================================================================
# Track V36 -- outcome-dependent selection: does censored-exposure imputation
#              AMPLIFY it, or are the two defects merely additive?
# -----------------------------------------------------------------------------
# Design and REGISTERED PREDICTIONS: DESIGN_v36_selection.md (2026-11-09,
# written before any V36 code ran). Read it before interpreting output.
#
# WHY. Itoh et al. 2025 excluded 1,534 of 16,709 mothers on low birth weight,
# which Mezawa et al.'s DAG draws as a mediator. Nothing in this harness removes
# rows as a function of anything -- every track to date varies which COVARIATES
# are conditioned on, never which ROWS survive. This blocks the manuscript's
# §6.3 third row.
#
# THE GATE IS THE INTERACTION, not selection bias on its own:
#     interaction = (arm_sel - oracle_sel) - (arm_none - oracle_none)
# positive = imputation makes selection bias worse. Zero is a clean result too
# (the defects are additive) and is a NOT MET, not a failure.
#
# INSTRUMENT CHECK ALREADY PASSED (P1, 2026-11-09, oracle arm only, n = 4000,
# 400 reps): `descnull` -- selection on a pure descendant of the exposure with
# no link to Y -- reads -0.20% +/- 0.30, i.e. no bias, while `desc` reads
# -17.2%, `outcome` -13.4% and `collider` -32.8%.
#
# STRENGTH LEVER IS THE DROP FRACTION, NOT `delta`. Collider bias is
# NON-MONOTONE in delta (0.80 -> -17.5%, 1.60 -> -12.8%): as delta grows,
# selection tends to selection on U alone, which is independent of the
# exposure. Raising delta to "strengthen" a cell makes the bias smaller and
# reads as a refutation. See R/selection.R.
#
# PILOT FIRST. Do not quote a duration from anything smaller than a full wave;
# this project's extrapolations have been wrong by 3x, 8x and 13x.
#   PILOT=1 ./run_v36_selection.sh      # 2 cells x 22 reps = 44 tasks = 1 chunk
# =============================================================================
set -uo pipefail
: "${STAGES:=}" "${NCORES:=22}" "${N_OBS:=800}" "${M_IMP:=30}" "${ND_FRAC:=0.30}"
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

if [[ "${PILOT:-0}" == "1" ]]; then
  TAG="${OUT_TAG:-v36pilot}"; N_REP="${N_REP:-22}"; CELLS="${CELLS:-none,desc}"
else
  TAG="${OUT_TAG:-v36}";      N_REP="${N_REP:-100}"; CELLS="${CELLS:-none,descnull,desc,desc2,outcome,outcome2,collider}"
fi

if [[ "${V36_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V36_STAMP="$(date +%Y%m%d-%H%M%S)"
  export PILOT OUT_TAG="$TAG" N_REP CELLS STAGES NCORES N_OBS M_IMP ND_FRAC
  # True detachment: `nohup ... &` only IGNORES SIGHUP and leaves the job in the
  # launching shell's PROCESS GROUP, which is how V35 was lost on 2026-10-26.
  # macOS ships no setsid(1); python3 is used only to make the syscall.
  if command -v python3 >/dev/null 2>&1; then
    V36_DETACHED=1 nohup python3 -c 'import os,sys; os.setsid(); os.execvp("bash",["bash"]+sys.argv[1:])' "$SELF" >/dev/null 2>&1 &
  else
    echo "WARNING: no python3; falling back to nohup, which does NOT survive a session teardown." >&2
    V36_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  fi
  PID=$!
  echo "$PID" > "logs/${TAG}_${V36_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V36 selection run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V36_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V36_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V36_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export N_OBS="${N_OBS:-800}" M_IMP="${M_IMP:-30}" ND_FRAC="${ND_FRAC:-0.30}"
STAMP="${V36_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"

{
  echo "==========================================================="
  echo "V36 -- outcome-dependent selection   [$(date '+%F %T')]"
  echo "  tag=${TAG}  cells=${CELLS}  n_rep=${N_REP}"
  echo "  n=${N_OBS}  m=${M_IMP}  nd_frac=${ND_FRAC}  cores=${NCORES}"
  echo "  arms: oracle, complete_case, sub_lod2, pipeline_block_fcs, smc_xgrid"
  echo "  design + registered predictions: DESIGN_v36_selection.md"
  if [[ "${PILOT:-0}" == "1" ]]; then
    echo "  *** TIMED PILOT -- measuring per-chunk cost, not scoring gates ***"
  fi
  echo "==========================================================="
} > "$LOG"

# STAGES. One detached job, run SEQUENTIALLY -- concurrent stages would contend
# for the same 22 cores and make every timing in this log meaningless.
#
# STAGES is a ";"-separated list of "n|nd_frac|reps|cells". Each stage writes
# its own result set and its own checkpoint, so a kill costs at most the stage
# in flight and RESUME=1 picks that stage up where it stopped.
#
# MEASURED COST (full-wave pilots, 2026-11-09, 22 cores, m=30, all five arms):
#   n =   800 -> 16.5 s/task -> 0.30 min per wave of 22
#   n =  3200 -> 35.0 s/task -> 0.60 min per wave
#   n = 12800 -> 112.8 s/task -> 1.90 min per wave
# Every number below is measured on a FULL WAVE at the real configuration, not
# extrapolated from a short pilot -- this project's extrapolations have been
# wrong by 3x, 8x and 13x.
if [[ -z "${STAGES:-}" ]]; then
  STAGES="${N_OBS}|${ND_FRAC}|${N_REP}|${CELLS}"
fi

rc=0
STAGE_NO=0
IFS=';' read -r -a STAGE_LIST <<< "$STAGES"
for SPEC in "${STAGE_LIST[@]}"; do
  STAGE_NO=$((STAGE_NO + 1))
  IFS='|' read -r S_N S_ND S_REP S_CELLS <<< "$SPEC"
  STAGE_TAG="${TAG}_n${S_N}_nd${S_ND}"
  { echo ""
    echo "-------------------------------------------------------------------"
    echo "[$(date '+%F %T')] STAGE ${STAGE_NO}/${#STAGE_LIST[@]}"
    echo "  n=${S_N}  nd_frac=${S_ND}  reps=${S_REP}"
    echo "  cells=${S_CELLS}"
    echo "  -> results/${STAGE_TAG}_latest.rds"
    echo "-------------------------------------------------------------------"
  } >> "$LOG"

  OUT_TAG="$STAGE_TAG" N_REP="$S_REP" CELLS="$S_CELLS" \
    N_OBS="$S_N" ND_FRAC="$S_ND" \
    Rscript run_v36_selection.R >> "$LOG" 2>&1 &
  RPID=$!
  # The WORKER pid, not this wrapper's -- run_status.sh tests it with kill -0.
  # Rewritten per stage, so a mid-sweep check reports the stage in flight.
  printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
  wait "$RPID"; rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[$(date '+%F %T')] STAGE ${STAGE_NO} (n=${S_N} nd=${S_ND}) FAILED rc=$rc -- stopping." >> "$LOG"
    break
  fi
  echo "[$(date '+%F %T')] stage ${STAGE_NO} done -> results/${STAGE_TAG}_latest.rds" >> "$LOG"
done

echo "" >> "$LOG"
echo "[$(date '+%F %T')] ALL STAGES finished rc=$rc" >> "$LOG"
