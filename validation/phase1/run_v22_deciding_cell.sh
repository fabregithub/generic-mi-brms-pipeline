#!/usr/bin/env bash
# =============================================================================
# Track V22 -- THE DECIDING CELL: is a parametric covariate draw a fix, or a
#              different bug?
# -----------------------------------------------------------------------------
# WHAT V21 LEFT UNRESOLVED. V21 attributed the confounder/mediator bias entirely
# to BART's FLEXIBILITY: a correctly specified ESTIMATED linear covariate draw
# came in at +0.26% / +0.58%, against BART's +4.7% / +5.1%. But both V21 cells
# have a LINEAR-GAUSSIAN covariate conditional, so every parametric arm was
# correct by construction. R8 adopted BART precisely because a parametric Z block
# is misspecified when that conditional is non-linear -- V6 measured `mice pmm`
# at -4.70% there. **Trading a 5% bias under linearity for a 5% bias under
# non-linearity is not a fix**, and no cell in the harness could test both sides.
#
# THIS BUILDS THAT CELL. `z_role = "pipe_nl"` is a pipe whose X1 -> Z1 arrow is
# non-linear: Z1 = 1.2*tanh(1.8*logX1) + 0.35*(logX1^2 - 1) + noise. So Z1 is
# genuinely ON an X-Y path AND its conditional has a non-linear mean -- a linear
# imputation model is really misspecified (residual SD 21% worse than the correct
# form) while BART can learn the shape. Verified before adoption: the analysis
# model still recovers b1 (0.4008 at n = 2e5), partial cor(Z1, logX1 | rest) =
# +0.602, and the exact Z1 conditional stays GAUSSIAN (skew 0.001, kurtosis
# 3.001) so the anchor remains closed form.
#
# WHY THE EXPOSURE IS FULLY OBSERVED IN THE TWO PRIMARY CELLS. Not a
# convenience -- required, and found by measurement while building this. With
# logX1 40% censored the ladder's exact-Z anchor sits at **+17%** instead of ~0,
# because `leftcens` draws logX1 from a conditional LINEAR in Z1 while the truth
# has Z1 = g(logX1). The X block is badly misspecified under a non-linear arrow,
# so it swamps the covariate-draw comparison and the ladder loses its zero point.
# With logX1 observed the same anchor is **+1.39% +/- 0.89** -- consistent with
# zero. `zr_pipenl_cens` keeps the censored version as a SECONDARY cell: no valid
# anchor, so it prices the total and flags the X-block interaction as a finding
# in its own right rather than pretending to decompose it.
#
# THE CELLS
#   zr_pipe_nc      linear arrow, exposure observed   <- matched comparator
#   zr_pipenl       NON-LINEAR arrow, exposure observed  <- THE DECIDING CELL
#   zr_pipenl_cens  non-linear arrow, exposure censored  <- secondary, no anchor
#
# THE ARMS (V21's ladder, plus the V19 probe)
#   z21_exact, z21_fit_proper, z21_fit_improper, z21_pmm, z21_bart
#   z21_fit_proper_z2same, z21_pmm_z2same   <- Z2 drawn by the SAME method as Z1
#
# The last two address the other thing V21 left open. V19 measured `micePmm` at
# +2.4 to +3.9% and `properZ` at +4.5 to +6.3% ASYMPTOTIC bias while V21's
# parametric arms sat at +-0.6%; the registered explanation was that V21 draws Z1
# alone with Z2 held exact, while mice imputes Z1, Z2 AND Y together. These arms
# test the first half of that.
#
# THE REGISTERED PREDICTION -- the parametric draw wins on BOTH sides.
#   Derived, not fitted: the analysis conditions on logX1, so a linear draw's
#   error in approximating g(logX1) is absorbed by logX1's own coefficient. What
#   matters for b1 is the part of Z1 carrying information about Y beyond X -- and
#   Y enters the exact conditional LINEARLY, which a linear fit represents
#   exactly. So misspecifying the PRIOR mean should be largely harmless for this
#   estimand, while BART's smoothing of it is not.
#
#     zr_pipenl:  exact ~0 | fit_proper |bias| <= 2 pp | bart >= +3%
#     zr_pipe_nc: the same ordering (V21's result, without censoring)
#
#   A 44-rep pilot corroborates (pipenl: exact -0.02%, fit_proper -0.01%,
#   pmm -0.34%, bart +3.95%) but its MC error is ~1.9 pp, so it cannot settle
#   anything -- it is reported here as corroboration of the mechanism, not as
#   the basis of the prediction.
#
# FALSIFICATION
#   * z21_exact drifting more than 1.5 pp from zero in either primary cell voids
#     the anchor and the ladder with it.
#   * fit_proper exceeding 2 pp in zr_pipenl => the parametric draw IS a
#     different bug, the trade is real, and `docs/covariate-roles.md`'s guidance
#     becomes permanent. **This is the outcome that matters most.**
#   * bart below +3% in zr_pipenl => flexibility stops costing anything once the
#     truth is non-linear, i.e. BART is doing its job and there is nothing to fix.
#   * If the z2same arms move by more than 2 pp, V19's discrepancy is the
#     multi-target machinery and the next step is to add Y as a target.
#
# NOTE (PLAN §10b): arms consume the task RNG stream in order, so V21's numbers
# are not the baseline -- `zr_pipe_nc` is the within-run comparator.
#
# COST -- MEASURED (full-wave pilot, 22 reps per cell per stage, 198 tasks, 0
# errors, at the run's own settings m = 30, sweeps = 3, 22 workers):
#
#   n     pipe_nc  pipenl  pipenl_cens   stage
#     800     5 s     5 s      15 s        9 min
#    3200    23      22        57         39 min
#   12800   290     290       405        373 min
#   TOTAL ~7.0 h
#
# n = 12800 is 89% of the run: BART's cost climbs steeply and there are three
# BART-bearing arms (bart, and the two *_z2same which draw Z2 by BART too when
# z_src is bart -- they do not, but they add a second fit for Z2 in the
# parametric arms). STAGES=1,2 costs ~48 min and still spans 4x, which settles
# the PRIMARY question (does fit_proper stay under 2 pp in zr_pipenl) since that
# is a level, not a slope. The third level is only needed for the n-signatures.
#
# USAGE
#   ./run_v22_deciding_cell.sh              # all three n levels (~7.0 h)
#   STAGES=1 ./run_v22_deciding_cell.sh     # n = 800 only
#   N_REP=22 ./run_v22_deciding_cell.sh     # full-wave timing pilot
#
# Results: results/v22_n<N>_latest.rds per stage.
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

TAG="${OUT_TAG:-v22}"

if [[ "${V22_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V22_STAMP="$(date +%Y%m%d-%H%M%S)"
  V22_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${V22_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V22 deciding-cell run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/${TAG}_${V22_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/${TAG}_${V22_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh ${TAG}"
  echo "Stop:  kill \$(cat validation/phase1/logs/${TAG}_${V22_STAMP}.pid)"
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
export SCENARIOS="${SCENARIOS:-zr_pipe_nc,zr_pipenl,zr_pipenl_cens}"
export ARMS="${ARMS:-z21_exact,z21_fit_proper,z21_fit_improper,z21_pmm,z21_bart,z21_fit_proper_z2same,z21_pmm_z2same}"

STAMP="${V22_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/${TAG}_${STAMP}.log"
mkdir -p results logs

NLEVELS=(0 800 3200 12800)

{
  echo "=== Track V22: the deciding cell ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} SEED=${SEED} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=oracle,${ARMS}"
  echo ""
  echo "THE QUESTION V21 could not answer: is a parametric covariate draw a fix,"
  echo "  or a different bug? Both V21 cells had a LINEAR-GAUSSIAN covariate"
  echo "  conditional, so every parametric arm was correct by construction."
  echo ""
  echo "  zr_pipe_nc      linear arrow, exposure observed  <- comparator"
  echo "  zr_pipenl       NON-LINEAR arrow, observed        <- DECIDING CELL"
  echo "  zr_pipenl_cens  non-linear, exposure censored     <- secondary, NO ANCHOR"
  echo ""
  echo "WHY THE EXPOSURE IS OBSERVED in the primary cells: measured, not assumed."
  echo "  With logX1 40% censored the exact-Z anchor sits at +17% instead of ~0,"
  echo "  because leftcens draws logX1 LINEAR in Z1 while the truth has"
  echo "  Z1 = g(logX1). With logX1 observed the anchor is +1.39% +/- 0.89."
  echo ""
  echo "REGISTERED PREDICTION -- the parametric draw wins on BOTH sides, because"
  echo "  the analysis conditions on logX1, so a linear draw's error in"
  echo "  approximating g(logX1) is absorbed by logX1's own coefficient, while Y"
  echo "  enters the exact conditional linearly and a linear fit gets that exactly."
  echo "    zr_pipenl: exact ~0 | fit_proper |bias| <= 2 pp | bart >= +3%"
  echo ""
  echo "MOST CONSEQUENTIAL OUTCOME: if fit_proper EXCEEDS 2 pp in zr_pipenl, the"
  echo "  parametric draw is a different bug, the trade is real, and the"
  echo "  documented guidance in docs/covariate-roles.md becomes permanent."
  echo ""
  echo "ALSO: z21_*_z2same draw Z2 by the same method as Z1 -- the first half of"
  echo "  V19's discrepancy (micePmm +2.4 to +3.9% asymptotic vs V21's +-0.6%)."
  echo "  A shift over 2 pp points at the multi-target machinery."
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
