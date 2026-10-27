#!/usr/bin/env bash
# =============================================================================
# V35 -- is the `f^2` exponent a property of the IMPUTATION or of the ANALYSIS
#        MODEL?
# -----------------------------------------------------------------------------
# WHY THIS RUN EXISTS. Two independent routes in this project put the censored
# fraction into the bias as `f^2`: V15's `excess% = -323.4*f^2` (2026-09-02) and
# V34's `f^1.85 +/- 0.20` (2026-10-24). FOUR mechanisms have been proposed and
# refuted -- constant Gamma, Gamma growing in f, direct propagation of the
# per-row imputation error, and (retired 2026-10-26 without a run, see below)
# V15's own evaluation-range conjecture.
#
# P-B (phase1/derive_declared_basis.R, 2026-10-26) localised it. A replica that
# matches V33's imputation in every respect checked against the source -- the
# closed-form conditional prior, the truncated quadratic grid draw, the
# completed-data proper refit, n = 800, m = 10, 3 sweeps, the L-0.5 cold start
# -- and differs ONLY in using OLS on the true basis instead of BKMR agrees with
# the asymptotic fixed point to <= 0.025 in all twelve cells, and still gives
# f-exponents of 0.83-1.33 against the measured 1.46-2.28. Finite n, properness,
# m and sweep count contribute essentially nothing. So the residual belongs to
# the ANALYSIS model's response to degraded curvature information.
#
# WHY V21 IS NOT THIS RUN. V21 (`QLO=0.40`) has been queued as "the theory
# frontier" since 2026-09-02. It is dead on V15's own data and was retired
# WITHOUT being run. The conjecture is `f^2 = f x (fraction of the estimand's
# evaluation range below the LOD)`. At q_lo = 0.25 that range is [-0.674,
# +0.674], so the second factor is IDENTICALLY ZERO for f <= 0.22 -- yet V15
# measured -13.9% +/- 2.0 at f = 0.20 (7.0 sigma) and -29.6% +/- 2.3 at f = 0.30
# (13.0 sigma). At q_lo = 0.40 it is zero for FOUR of the six sweep points.
# Under the only other reading (second factor = f) it does not depend on q_lo at
# all and the test is vacuous. Either way the run buys nothing.
#
# DESIGN. V15's six-point f-sweep, unchanged, with the BKMR arms and the OLS
# arms side by side on the SAME bundle and the SAME seed, so the comparison is
# PAIRED and the only thing that differs is the analysis model. The OLS arm is
# an INSTRUMENT, not a rival: it is handed the generator's own basis, so it is
# correctly specified by construction and its response to information loss is
# known exactly. Whatever BKMR adds shows up as the difference.
#
# REGISTERED PREDICTIONS (2026-10-26, before the run):
#   P1  THE GATE. The log-log slope of |excess%| in f (f >= 0.2) for curv_X1:
#           pipeline_bkmr  ~ 2       (reproducing V15's 1.98 [1.80, 2.17])
#           pipeline_lm    ~ 0.85    (derive_declared_basis.R STAGE 4)
#       and the two CIs must NOT overlap. If they overlap, P-B's localisation is
#       WRONG and the exponent lives in the imputation after all.
#   P2  oracle_lm is FLAT in f and consistent with zero on curv_X1 and int_X1X2
#       -- it is complete data, so any trend in it is a harness fault and
#       invalidates P1. (22-rep pilot: -6.2% to +4.1%, MCSE ~4%. Holds.)
#   P3  pipeline_bkmr reproduces V15 on the CURRENT harness (nd20 -13.9%,
#       nd40 -55.6%, +/- MC error). V15 predates the mice default (v1.7.0) and
#       several other changes; if this fails, V15's curve is not comparable and
#       P1 must be read within this run only.
#   NOT A GATE: the LEVEL difference between the arms. OLS and BKMR need not
#   agree in level -- BKMR's own oracle is +15% biased on curvature -- and only
#   the SHAPE in f is under test. Do not read a level gap as evidence either way.
#
# PRELIMINARY (22-rep full-wave pilot, lm arms only, 2026-10-26):
#   slope 0.64 [0.44, 0.84] against V15's BKMR 1.98 [1.80, 2.17]. Strongly in
#   P1's direction, but 22 reps and UNPAIRED (compared across runs five weeks
#   apart). That is exactly what this run fixes.
#
# COST. MEASURED, not extrapolated. The BKMR half is V15's own identical
# workload on this machine: 2026-09-02 07:58:12 -> 19:08:57 = 11 h 11 m
# (N_REP=100, ITER=3000, NCORES=22, 6 scenarios). The OLS arms were measured by
# a FULL-WAVE pilot (reps = 22 = worker count, all six scenarios): 25 SECONDS,
# so ~2 min at 100 reps. Total ~11 h 15 m. Checkpoints as it goes, so a run
# killed early still yields partial results.
#
# Usage:
#   ./run_v35_analysis_model.sh              # the full paired run (~11 h 15 m)
#   N_REP=10 ./run_v35_analysis_model.sh     # short pilot
#   ARMS=oracle_lm,pipeline_lm N_REP=100 ./run_v35_analysis_model.sh   # ~2 min
#
# Results: results/v35_latest.rds, v35_summary.csv, v35_raw.csv
# Analyse: Rscript analyze_v35.R
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v35}"

if [[ "${V35_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V35_STAMP="$(date +%Y%m%d-%H%M%S)"
  # True detachment. `nohup ... &` only IGNORES SIGHUP -- it leaves the job in
  # the launching shell's PROCESS GROUP, so when a terminal or host app tears
  # that group down the "detached" run dies with it. V35 was lost this way on
  # 2026-10-26: killed 1.6 h into an 11 h run, no R error, no crash report, no
  # OOM, no sleep -- the wrapper and its child went at the same instant.
  # setsid(2) puts the run in its own session so nothing upstream can reach it.
  # macOS ships no setsid(1); python3 is used only to make the syscall.
  if command -v python3 >/dev/null 2>&1; then
    V35_DETACHED=1 nohup python3 -c 'import os,sys; os.setsid(); os.execvp("bash",["bash"]+sys.argv[1:])' "$SELF" >/dev/null 2>&1 &
  else
    echo "WARNING: no python3; falling back to nohup, which does NOT survive a session teardown." >&2
    V35_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  fi
  PID=$!
  echo "$PID" > "logs/${TAG}_${V35_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V35 analysis-model run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V35_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V35_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V35_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
# SEED must match V15's so the BKMR arms are comparable cell for cell.
export SEED="${SEED:-20260828}"
export N_REP="${N_REP:-100}"
export N_OBS="${N_OBS:-800}"
export ITER="${ITER:-3000}"
export M="${M:-10}"
export KNOTS="${KNOTS:-50}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export OUT_TAG="${TAG}"
export SCENARIOS="${SCENARIOS:-mixf10,nd20,mixf30,nd40,mixf50,mixf60}"
export ARMS="${ARMS:-oracle_bkmr,pipeline_bkmr,oracle_lm,pipeline_lm}"

STAMP="${V35_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

{
  echo "=== Track V35: is the f^2 exponent the imputation's or the analysis model's? ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} N_OBS=${N_OBS} M=${M} ITER=${ITER} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
  echo ""
  echo "P1 (GATE): curv_X1 log-log slope in f (f>=0.2), CIs must NOT overlap"
  echo "             pipeline_bkmr ~ 2.0   |   pipeline_lm ~ 0.85"
  echo "           overlap => P-B's localisation is WRONG"
  echo "P2: oracle_lm flat in f and consistent with 0 (else harness fault)"
  echo "P3: pipeline_bkmr reproduces V15 (nd20 -13.9%, nd40 -55.6%)"
  echo "NOT A GATE: the level difference between OLS and BKMR arms."
  echo "MEASURED COST: ~11 h 15 m (BKMR 11 h 11 m from V15; OLS 25 s/full wave)"
  echo "==================================================================="
} > "$LOG"

# Record the worker's PID, not this wrapper's -- see run_v11_ynl.sh.
Rscript run_v3_bkmr.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
echo "[$(date '+%H:%M:%S')] done rc=$rc -> results/${TAG}_latest.rds" >> "$LOG"
