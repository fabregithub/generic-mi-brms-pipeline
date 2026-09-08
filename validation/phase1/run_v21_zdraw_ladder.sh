#!/usr/bin/env bash
# =============================================================================
# Track V21 -- WHAT about a covariate draw has to be right?
# -----------------------------------------------------------------------------
# WHERE THIS SITS. V17 found the shipped default biased +5% to +10.4% under a
# confounder or mediator. V18: it decays as n^-1/3 and never vanishes. V19: only
# BART's decays at all, and "correctly specified in form" is NOT enough --
# micePmm and properZ are correct in form for this Z1 and still carry +4 to +6%
# ASYMPTOTIC bias. V20: the bias is the COVARIATE draw (85-104% of it), not the
# exposure draw. **This track asks the question that stands between that and a
# fix: which PROPERTY of a covariate draw has to be right?**
#
# WHY IT IS TRACTABLE. V20 handed over a known-unbiased anchor -- the exact
# conditional lands within +-0.5% at every n. In this DGP the true conditional of
# Z1 is Gaussian with CONSTANT variance, so V12's heteroscedastic-shape effect
# cannot arise, and the separable candidates are few:
#
#   arm                 what it changes from the anchor
#   z21_exact           nothing -- the true conditional (ANCHOR, ~0)
#   z21_fit_proper      the conditional is ESTIMATED, correct form, proper draw
#   z21_fit_improper    the same fit, MLE plugged in -- isolates PROPERNESS
#   z21_pmm             the same fit, but a DONOR is returned, not a density draw
#   z21_bart            BART -- a nonparametric fit of a linear truth (REFERENCE)
#   z21_bart_inner3     BART with the shipped inner-FCS loop (3 passes)
#
# The X block is held at the shipped `leftcens` draw in every arm: V20 measured
# its share at 0.1-11.7%, so dropping that dimension makes the arms cheap and
# changes nothing being asked.
#
# EACH CANDIDATE HAS A DIFFERENT n-SIGNATURE, which is what makes this decisive
# rather than descriptive:
#
#   estimation error   decays as n^-1/2   (parametric rate, correct form)
#   smoothing          decays as n^-1/3   (V18, V19, V20 all measured this)
#   donor matching     does NOT decay     (a fixed discretisation of the density)
#   properness         moves the INTERVAL, not the point estimate (V4's signature)
#
# THE REGISTERED PREDICTION
#   1. z21_exact stays within 1 pp of zero at every n (re-confirms V20's anchor).
#   2. z21_bart carries the reference bias, ~+4 to +5% at n = 800, decaying at
#      ~-1/3.
#   3. z21_fit_proper is SMALL -- |bias| <= 2 pp at n = 800 -- because a
#      correctly specified linear fit has no smoothing error, only estimation
#      error, and that decays fast.
#   4. z21_fit_improper matches z21_fit_proper on BIAS (within 1.5 pp) but has a
#      SMALLER `b` and lower coverage. That is improper MI's signature from V4,
#      and finding it here would confirm properness is an interval problem, not
#      a bias one.
#   5. z21_pmm's bias does NOT decay with n (slope CI includes 0).
#
# WHAT THIS DOES *NOT* CLAIM, and the distinction matters. `z21_pmm` is a
# hand-rolled pmm on the correct formula for Z1 alone. V19's `micePmm` is the
# pipeline's mice path: it imputes Z1, Z2 AND Y together, with mice's own
# predictor matrix and iteration. They differ in more than the draw, so if the
# ladder does not reproduce micePmm's +4 to +6%, that is INFORMATIVE -- it would
# mean the damage is in the multi-target machinery rather than the draw -- and
# not a contradiction. Registering that in advance is the point.
#
# FALSIFICATION
#   * z21_exact drifting more than 1 pp from zero voids the anchor and with it
#     the whole ladder.
#   * z21_fit_proper exceeding 2 pp at n = 800 would say estimating a CORRECTLY
#     SPECIFIED conditional is itself enough to produce the defect -- which
#     would mean there is probably no shippable fix, only the documented
#     guidance. That is the most consequential single outcome here.
#   * z21_bart failing to reproduce V20's ef_bart_ship (within 1.5 pp) means the
#     instrument changed and nothing is comparable.
#   * If no arm on the ladder reaches the reference bias, the property at fault
#     is not on it, and the next step is the multi-target machinery.
#
# NOTE (PLAN §10b): arms consume the task RNG stream in order, so V20's numbers
# are NOT the baseline -- the ladder is read within this run.
#
# COST -- MEASURED (full-wave pilot, 22 reps per cell per stage, 132 tasks, 0
# errors, at the run's own settings m = 30, sweeps = 3, 22 workers):
#
#   n =   800   18.4 s/task ->  14 min
#   n =  3200   56.5           43 min
#   n = 12800  315.9          239 min      <- 80% of the run
#   TOTAL ~4.9 h
#
# The top level is dominated by the two BART arms, and `z21_bart_inner3` does 3x
# the BART work per sweep by construction. STAGES=1,2 costs ~1 h and still spans
# 4x -- enough for a rough n-signature, with a wide slope CI. Since the
# n-signature is what separates estimation error (-1/2) from smoothing (-1/3)
# from matching (flat), the third level is where the design earns its keep.
#
# USAGE
#   ./run_v21_zdraw_ladder.sh              # all three n levels (~4.9 h)
#   STAGES=1 ./run_v21_zdraw_ladder.sh     # n = 800 only
#   N_REP=22 ./run_v21_zdraw_ladder.sh     # full-wave timing pilot
#
# Results: results/v21_n<N>_latest.rds per stage.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v21}"

if [[ "${V21_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V21_STAMP="$(date +%Y%m%d-%H%M%S)"
  V21_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V21_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V21 covariate-draw ladder detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V21_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V21_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V21_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260908}"
export N_REP="${N_REP:-500}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export SCENARIOS="${SCENARIOS:-zr_fork,zr_pipe}"
export ARMS="${ARMS:-z21_exact,z21_fit_proper,z21_fit_improper,z21_pmm,z21_bart,z21_bart_inner3}"

STAMP="${V21_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 3200 12800)

{
  echo "=== Track V21: what about a covariate draw has to be right? ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "THE LADDER (X block held at the shipped leftcens draw throughout):"
  echo "  z21_exact         the true conditional -- ANCHOR, expect ~0"
  echo "  z21_fit_proper    estimated, correct form, proper draw"
  echo "  z21_fit_improper  same fit, MLE plug-in  -> isolates PROPERNESS"
  echo "  z21_pmm           same fit, donor returned -> isolates MATCHING"
  echo "  z21_bart          BART on a linear truth -> REFERENCE, ~+4 to +5%"
  echo "  z21_bart_inner3   BART + the shipped inner-FCS loop"
  echo ""
  echo "n-SIGNATURES, which is what makes the ladder decisive:"
  echo "  estimation error n^-1/2 | smoothing n^-1/3 | matching flat |"
  echo "  properness moves the INTERVAL, not the point estimate"
  echo ""
  echo "REGISTERED: exact within 1 pp of 0; bart ~+4 to +5% decaying at -1/3;"
  echo "  fit_proper |bias| <= 2 pp at n=800; fit_improper matches fit_proper on"
  echo "  bias but with smaller b and lower coverage; pmm's bias does NOT decay."
  echo ""
  echo "MOST CONSEQUENTIAL OUTCOME: if fit_proper EXCEEDS 2 pp, then estimating"
  echo "  a correctly specified conditional is itself enough to cause the defect"
  echo "  -- and there is probably no shippable fix, only the documented guidance."
  echo ""
  echo "NOT CLAIMED: that this ladder explains V19's micePmm (+4 to +6%)."
  echo "  z21_pmm draws Z1 alone on the correct formula; micePmm imputes Z1, Z2"
  echo "  and Y together with mice's own predictor matrix. A mismatch would mean"
  echo "  the damage is in the multi-target machinery, not the draw."
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
