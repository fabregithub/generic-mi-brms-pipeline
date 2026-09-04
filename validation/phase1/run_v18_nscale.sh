#!/usr/bin/env bash
# =============================================================================
# Track V18 -- does the concealment break at larger n?
# -----------------------------------------------------------------------------
# THE QUESTION. V17 found the SHIPPED default (pipeline_bartMI, Y in both
# blocks) biased +4.96% to +10.40% wherever the covariate is a confounder or a
# mediator -- with coverage still 0.93-0.97 and honestly-sized intervals
# (width/SE 0.98-1.13). The offered explanation was that the concealment is a
# SAMPLE-SIZE ACCIDENT: at n = 800 the bias is only 0.35-0.67 empirical SE, and
# since relative bias is ~constant in n while the SE falls as 1/sqrt(n), bias/SE
# must grow and coverage must eventually collapse.
#
# THAT EXPLANATION WAS ARITHMETIC ON ONE CELL, NOT A MEASUREMENT, and it has a
# serious competitor: the bias may be FINITE-SAMPLE -- an artefact of imputing
# from 800 rows -- in which case it shrinks as 1/sqrt(n), bias/SE stays put, and
# coverage holds at every n. The two hypotheses predict opposite things, which is
# what makes this worth a run rather than a paragraph.
#
# IT MATTERS BECAUSE OF WHAT THIS PIPELINE IS FOR. The README's premise is
# large-N, repeated-measures data. A defect that hides at n = 800 and surfaces at
# n = 12800 is the worst possible shape for it, and the one we would ship
# unknowingly.
#
# THE PRIMARY OUTCOME IS RELATIVE BIAS, NOT COVERAGE. Coverage needs a normal
# model and the assumption that width/SE is constant in n; relative bias needs
# neither. Coverage is reported as the consequence.
#
# THE REGISTERED PREDICTION (predict_nscaling.R, closed form; the coverage model
# reproduces all four of V17's measured coverages to within 0.005, which is the
# check that it describes this harness at all):
#
#   relative bias        constant-bias        finite-sample (1/sqrt(n))
#     zr_fork               4.96%  at all n     4.96% -> 1.24% at n=12800
#     zr_pipe               5.39%               5.39% -> 1.35%
#     zrmar_fork            8.52%               8.52% -> 2.13%
#     zrmar_pipe           10.40%              10.40% -> 2.60%
#
#   coverage under constant bias   n=800   n=3200   n=12800
#     zr_fork                      0.939    0.892     0.693
#     zr_pipe                      0.944    0.903     0.731
#     zrmar_fork                   0.927    0.806     0.358
#     zrmar_pipe                   0.936    0.807     0.321
#   under finite-sample bias: coverage stays at its n=800 value at every n.
#
# FIVE n LEVELS, NOT TWO, so the exponent is estimable rather than merely
# testable: regress log|relative bias| on log n. Constant bias predicts slope
# 0, finite-sample predicts -0.5. That is the V15 pattern -- a sweep with a
# slope test beats a two-point comparison, and V15's exponent CI is why.
#
# FALSIFICATION, PRE-DECLARED
#   * PRIMARY: the fitted slope of log|rel bias| on log n. Constant-bias is
#     rejected if the 95% CI excludes 0; finite-sample is rejected if it
#     excludes -0.5. BOTH may be rejected -- an intermediate exponent would say
#     the bias has two components, which no current account predicts.
#   * Relative bias at n = 12800 differing from the constant-bias prediction by
#     more than 2 pp in any of the four biased cells.
#   * CONTROL 1: `oracle` must stay unbiased (|rel bias| < 2%) at every n. It is
#     the check that the estimand and the DGP are not themselves n-dependent; if
#     it drifts, nothing else in the run is readable.
#   * CONTROL 2: `zr_collider` (+0.19% at n=800) and `mcar_z40` (-2.27%) must not
#     develop a fork/pipe-sized bias. mcar_z40 is the informative one -- if its
#     -2.27% shrinks as 1/sqrt(n) while the fork/pipe cells' bias holds, the two
#     have different origins, which is a finding in itself.
#   * width/SE drifting outside 0.9-1.2 at any n voids the coverage projections
#     (not the bias test).
#
# ARM SET IS DELIBERATELY MINIMAL: oracle + pipeline_bartMI. The question is
# about the shipped default's own bias, and the noY arms would triple the cost
# to answer a different question. NOTE (PLAN §10b): an arm's numbers are only
# comparable across runs with an IDENTICAL arm set, because procedures consume
# the task RNG stream in order. So V17's n=800 figures are NOT the baseline here
# -- n = 800 is re-run inside this track with this arm set, and every comparison
# below is within this run.
#
# COST -- MEASURED, not extrapolated. A full-wave pilot (22 reps per cell per
# stage, 660 tasks, 0 errors) at the run's own settings (m = 30, sweeps = 3,
# 22 workers), for 6 cells x 300 reps:
#
#   n =   800    5.9 s/task ->   8 min
#   n =  1600    9.0            12
#   n =  3200   15.2            21
#   n =  6400   28.7            39
#   n = 12800   57.0            78
#   TOTAL ~2.6 h
#
# Cost grows as roughly n^0.75, not n^1 -- a 16x sample costs ~10x, which is why
# five levels are affordable here where two were assumed. n = 12800 alone is half
# the run; STAGES=1,2,3,4 stops at n = 6400 for ~1.3 h and still spans 8x, enough
# for the slope test with a wider CI.
#
# USAGE
#   ./run_v18_nscale.sh                    # all five n levels (~2.6 h)
#   STAGES=1,5 ./run_v18_nscale.sh         # just the ends
#   N_REP=22 ./run_v18_nscale.sh           # full-wave timing pilot
#
# Results: results/v18_n<N>_latest.rds per stage.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v18}"

if [[ "${V18_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V18_STAMP="$(date +%Y%m%d-%H%M%S)"
  V18_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V18_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V18 n-scaling run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V18_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V18_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V18_STAMP}.pid)"
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
export SCENARIOS="${SCENARIOS:-zr_fork,zr_pipe,zrmar_fork,zrmar_pipe,zr_collider,mcar_z40}"
export ARMS="${ARMS:-pipeline_bartMI}"

STAMP="${V18_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 1600 3200 6400 12800)

{
  echo "=== Track V18: does the concealment break at larger n? ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "PRIMARY OUTCOME: relative bias of pipeline_bartMI, not coverage."
  echo "  constant bias  -> 4.96 / 5.39 / 8.52 / 10.40% at EVERY n"
  echo "  finite-sample  -> those / 4 by n = 12800 (1.24 / 1.35 / 2.13 / 2.60%)"
  echo "SLOPE TEST: log|rel bias| on log n. Constant predicts 0, finite-sample"
  echo "  -0.5. Both may be rejected; an intermediate exponent would say the"
  echo "  bias has two components, which no current account predicts."
  echo "COVERAGE (consequence, under constant bias):"
  echo "  zr_fork .939/.892/.693  zr_pipe .944/.903/.731   at n=800/3200/12800"
  echo "  zrmar_fork .927/.806/.358   zrmar_pipe .936/.807/.321"
  echo "CONTROLS: oracle unbiased at every n; zr_collider and mcar_z40 must not"
  echo "  develop a fork/pipe-sized bias."
  echo "NOTE: V17's n=800 numbers are NOT the baseline -- different arm set means"
  echo "  a different RNG stream (PLAN §10b). n=800 is re-run here."
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
