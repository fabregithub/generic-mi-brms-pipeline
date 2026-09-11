#!/usr/bin/env bash
# =============================================================================
# Track V26 -- does 25 tracks' worth of findings transfer to a NON-GAUSSIAN
#              estimand?
# -----------------------------------------------------------------------------
# THE GAP. Every track V0-V25 estimated a Gaussian linear coefficient by OLS
# (`fit_lm_estimand`). That is what made 500-replication runs affordable, and
# for a Gaussian linear coefficient it loses nothing -- pooled posteriors are
# symmetric and mesokurtic (THEORY.md 6b), so location and scale describe them
# completely. **But the pipeline's real targets are logistic, ordinal `mo()`,
# spline and mixed models**, and the transfer of the whole record to those
# classes has never been tested. It is an assumption, and it is the largest
# open gap in the programme (ROADMAP.md).
#
# WHAT THIS TRACK DOES. Re-runs three cells whose Gaussian answers are already
# known, with the outcome drawn from a LOGISTIC model on the SAME linear
# predictor. The estimand becomes the focal exposure's log-odds-ratio, still
# 0.40. Everything else -- structural coefficients, covariate roles, censoring,
# missingness -- is unchanged, and the same seed gives BYTE-IDENTICAL exposures
# (Y is drawn last), so any difference is the family and nothing else.
#
#   bin_base / mcar_z40   inert covariate      -- the calibration pair
#   bin_fork / zr_fork    confounder           -- V17's headline cell
#   bin_pipe / zr_pipe    mediator             -- V17's other role
#
# Both members of each pair run in the SAME invocation with an identical arm
# set, because PLAN 11 says an arm is only comparable across runs that share
# one -- and V19's "discrepancy" turned out to be partly that.
#
# CURRENCY. A logistic coefficient is on a different scale, so RELATIVE BIAS IS
# NOT COMPARABLE ACROSS FAMILIES. Everything here is read in **bias/SE** and
# coverage, which are. This is the V18 currency argument doing real work rather
# than being a convention.
#
# REGISTERED CRITERIA (PLAN 11: each must be able to fail for the reason it
# names; none may be satisfiable by an uninformative result).
#
#   C1  ORDERING SURVIVES. Within each cell, the arms rank the same by
#       |bias/SE| under both families. Fails if the family changes which
#       imputation strategy is better -- the outcome that would invalidate the
#       most prior guidance. Cannot pass vacuously: the Gaussian cell must
#       itself separate the arms by more than 0.3 in |bias/SE|, or there is no
#       ordering to preserve.
#   C2  SIZE SURVIVES. For each arm, |bias/SE| under binomial is within a
#       factor of 2 of its Gaussian twin. A factor-of-2 bar rather than a tight
#       one because the two families have genuinely different information per
#       observation; what would matter is a defect that VANISHES or DOUBLES.
#   C3  THE ROOT CLAIM SURVIVES. `pipeline_noYx` (Y removed from the exposure
#       draw) is still biased toward the null under binomial. V16 measured
#       -13.1 pp on the Gaussian side; the sign is the claim, not the size.
#   C4  THE SCALE DEFECT'S FAMILY-DEPENDENCE. `pipeline_properZ_ds` minus
#       `pipeline_properZ` is the pre-v1.6.0 covariate-block scale defect,
#       measured on identical data. It was 4.40 pp (Gaussian, confounder).
#       DESCRIPTIVE: this records whether a fixed defect reads the same under a
#       different family, which is the sharpest single test of transfer. No
#       pass/fail, and it must not be reported as one.
#
# MOST CONSEQUENTIAL OUTCOME. If C1 fails, guidance derived from 25 Gaussian
# tracks cannot be quoted for logistic analyses without re-derivation, and the
# claims ledger needs a scope line. If C1 and C2 hold, the record transfers and
# the remaining gap is only Step 6's pooling machinery (see below).
#
# WHAT THIS TRACK DOES *NOT* CLOSE. Step 6's finite-`m` variance correction on
# a support-respecting transform, gated by a bimodality diagnostic, still never
# fires -- this track pools with Rubin's rules in the harness, not through
# Step 6. Exercising that needs real `brms` fits through Steps 3-6 and is a
# separate, smaller piece of work. Recorded so a green result here is not
# mistaken for closing the whole gap.
#
# COST -- measured by full-wave pilots, never extrapolated (CLAUDE.md).
#     n     pilot (24 tasks)   N_REP=500      N_REP=400 (recommended)
#    800          45 s             1.6 h            1.2 h
#   3200         120 s             4.2 h            3.3 h
#  12800         432 s            15.0 h           12.0 h
#   TOTAL                         20.7 h           16.6 h
#
# Measured 2026-09-11 on 22 workers. n = 12800 is 72% of the run.
#
# WHY N_REP = 400 RATHER THAN 500. 500 lands at 20.7 h, which does not fit a
# 20 h unattended window with any margin, and an overrun on the last stage
# wastes the whole stage. At 400 the Monte-Carlo error of a mean bias/SE is
# 1/sqrt(400) = 0.050 against 0.045 at 500 -- immaterial for every gate here
# (C1 needs a 0.3 separation, C2 is a 2x ratio, C3 needs a sign against
# 2*mcse = 0.10). Use 500 only if running over a weekend.
#
#   ~16.6 h :  N_REP=400 ./run_v26_family.sh
#   ~ 4.5 h :  N_REP=400 STAGES=1,2 ./run_v26_family.sh     (settles C1-C3 as levels)
#   ~20.7 h :  ./run_v26_family.sh                          (weekend; N_REP=500)
#
# USAGE
#   ./run_v26_family.sh                  # all stages
#   STAGES=1 ./run_v26_family.sh         # n = 800 only
#   N_REP=4 ./run_v26_family.sh          # full-wave timing pilot
#
# Results: results/v26_n<N>_latest.rds per stage.  Analyse with analyze_v26.R.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v26}"

if [[ "${V26_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V26_STAMP="$(date +%Y%m%d-%H%M%S)"
  V26_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V26_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V26 outcome-family run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V26_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V26_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V26_STAMP}.pid)"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260911}"
export N_REP="${N_REP:-500}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export SCENARIOS="${SCENARIOS:-bin_base,bin_fork,bin_pipe,mcar_z40,zr_fork,zr_pipe}"
export ARMS="${ARMS:-complete_case,leftcens_prestep,pipeline_bartMI,pipeline_noYx,pipeline_noYz,pipeline_properZ,pipeline_properZ_ds}"

# N_REP = 500 comes from the currency. The gates are read in bias/SE, whose
# Monte-Carlo error is about 1/sqrt(R): 0.045 at R = 500. C2's factor-of-2 bar
# needs less, but C1's ordering test needs to separate arms that may sit close
# together, and a logistic coefficient carries less information per observation
# than a Gaussian one -- so this is deliberately not reduced from V23's 500.

STAMP="${V26_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 3200 12800)

{
  echo "=== Track V26: does the record transfer to a non-Gaussian estimand? ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "THE GAP: all 25 prior tracks estimated a GAUSSIAN LINEAR coefficient by"
  echo "  OLS. The pipeline's real targets are logistic, ordinal mo(), spline"
  echo "  and mixed models. Transfer has never been tested -- it is assumed."
  echo ""
  echo "DESIGN: three cells re-run with Y from a LOGISTIC model on the SAME"
  echo "  linear predictor, paired with their Gaussian twins in the same"
  echo "  invocation and the same arm set. Same seed => identical exposures."
  echo ""
  echo "CURRENCY: bias/SE and coverage ONLY. Relative bias is NOT comparable"
  echo "  across families -- a log-OR is on a different scale."
  echo ""
  echo "REGISTERED CRITERIA"
  echo "  C1 arm ORDERING by |bias/SE| survives the family change"
  echo "     (void if the Gaussian cell does not separate arms by > 0.3)"
  echo "  C2 each arm's |bias/SE| is within 2x of its Gaussian twin"
  echo "  C3 pipeline_noYx is still biased TOWARD THE NULL (V16's root claim)"
  echo "  C4 the scale defect's size by family -- DESCRIPTIVE, not a gate"
  echo ""
  echo "MOST CONSEQUENTIAL: if C1 fails, no guidance from the Gaussian record"
  echo "  can be quoted for logistic analyses without re-derivation."
  echo ""
  echo "DOES NOT CLOSE: Step 6's pooling machinery still never fires -- this"
  echo "  track pools with Rubin's rules in the harness. Separate work."
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

{
  echo ""
  echo "[$(date '+%H:%M:%S')] V26 finished (rc=$rc)"
  echo "Analyse: Rscript analyze_v26.R"
} >> "$LOG"
exit $rc
