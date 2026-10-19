#!/usr/bin/env bash
# =============================================================================
# V33 -- can substantive-model-compatible imputation be used when the surface
#        form is DECLARED rather than known?
# -----------------------------------------------------------------------------
# THE GAP. For a mixture surface the shipped censored-exposure engine must not
# be used: it destroys 57% of the curvature at 40% non-detects, worse than
# LOD/sqrt(2). V3 diagnosed why -- the engine draws the censored exposure from a
# conditional LINEAR in its predictors, so curvature cannot be represented. V9
# then showed the fix works: replacing only that draw closes 89-100% of the
# curvature gap, and the COEFFICIENTS need not be known (plugin mode, which
# estimates the surface each sweep, matches the oracle to 0.2 pp).
#
# But in both V9 modes the FORM came from `dgp_formula(truth)` -- the generator's
# own. Obtaining the form for an unknown surface is what keeps item 07 blocked.
# V33 replaces it with a form the ANALYST DECLARES (`declared_form`, new in
# smc_impute.R; default NULL keeps every pre-V33 number byte-identical, and the
# 16-check harness self-test still passes).
#
# WHY THIS IS THE RIGHT QUESTION FOR A MIXTURE ANALYSIS. The analyst has already
# declared a surface -- it is the mixture model being fitted. So "the form is
# unknown" may be false exactly here, which would make this the one setting
# where SMC is practical.
#
# REGISTERED PREDICTIONS (2026-10-19, before the run). Transferred from the
# COVARIATE-ARROW analogue already measured (SUMMARY.md), where the same
# question was settled for Z1 = g(logX1):
#     correct basis, coefficients fitted ... -0.51%   (as good as knowing them)
#     correct form, wrong constant ......... +3.22%   (degrades gracefully)
#     omit the factor entirely ............. +8.47%
#     correct family, ONE TERM DROPPED ..... +20.71%  (worse than omitting)
#     a straight line ...................... +25.43%
#
#   P1 RECOVERY  `smc_dec_true` recovers the curvature: beats both
#      `pipeline_bkmr` and `sub_lod2_bkmr` on curvature bias, and lands near
#      `oracle_bkmr`. (Also a plumbing check -- it should reproduce V9's plugin
#      arm, since it declares the same form.)
#   P2 INCOMPLETE BASIS IS WORSE THAN NONE  `smc_dec_noint` and
#      `smc_dec_noquad` are WORSE than `smc_dec_lin`. This is the decisive one:
#      if it holds, declared-form SMC is unsafe wherever the basis is uncertain
#      -- which for a PFAS mixture is always -- and the honest finding is that
#      there is no usable method, not that we found one.
#   P3 OVER-DECLARING IS CHEAP  `smc_dec_over` (two spurious terms) is close to
#      `smc_dec_true`. If P2 and P3 both hold the practical rule is "declare
#      generously, never sparsely", which is directly shippable advice.
#
#   FLOOR: the curvature estimand must be measurably non-zero in `oracle_bkmr`,
#   or every ratio below is a ratio of near-zero arms. Seven gates in this
#   programme have been spoiled that way.
#
#   NOT A GATE: which arm "wins". The literature already recommends truncated MI
#   / censored AFT for mixtures (Env Health 2024) and finds qgcomp the most
#   imputation-robust method (2026). V33's contribution is the MECHANISM and its
#   limit, not a leaderboard.
#
# COST -- from V3's own MEASURED header (24 cores, n=800, 50 GPP knots,
# NCORES=22, iter=3000): all 4 arms at M=10 is 13 fits/rep and ~28 min per 22
# reps. V33 runs 3 single-fit arms + 6 MI arms, so at M=5 that is 33 fits/rep:
#     M=5  N_REP=200  1 scenario ... ~10.8 h   <- default
#     M=5  N_REP=100  1 scenario ... ~5.4 h
#     M=10 N_REP=200  1 scenario ... ~20.6 h
# ⚠️ DERIVED, NOT MEASURED. Run the full-wave pilot below before trusting it --
# V3's header calls its own 12.7 h a LOWER BOUND, and V32 was out by 13x from an
# arm-count extrapolation.
#     N_REP=22 SCENARIOS=nd40 bash run_v33_mixture.sh   # one timed full wave
#
# M is the cheap thing to cut (bias is insensitive to m); scale N_REP LAST,
# since bias is the headline and N_REP is what buys precision on it.
#
#   bash run_v33_mixture.sh        # detach; logs/v33_<stamp>.{log,pid}
#   bash run_status.sh v33
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
cd "$(dirname "$SELF")"
TAG=v33
STAMP="$(date '+%Y%m%d-%H%M%S')"
mkdir -p logs results

if [[ -z "${V33_DETACHED:-}" ]]; then
  V33_DETACHED=1 V33_STAMP="$STAMP" nohup bash "$SELF" "$@" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V33 detached (pid $PID)"
  echo "  log:    logs/${TAG}_${STAMP}.log"
  echo "  status: bash run_status.sh ${TAG}"
  exit 0
fi

STAMP="${V33_STAMP:-$STAMP}"
LOG="logs/${TAG}_${STAMP}.log"
export NCORES="${NCORES:-$(( $(sysctl -n hw.ncpu) - 2 ))}"
export N_REP="${N_REP:-200}"
export M="${M:-5}"
export ITER="${ITER:-3000}"
export N_OBS="${N_OBS:-800}"
export SCENARIOS="${SCENARIOS:-nd40}"
export ARMS="${ARMS:-oracle_bkmr,cc_bkmr,sub_lod2_bkmr,pipeline_bkmr,smc_dec_true,smc_dec_noint,smc_dec_noquad,smc_dec_lin,smc_dec_over}"
export OUT_TAG="${OUT_TAG:-${TAG}}"

{ echo "[$(date '+%H:%M:%S')] V33 start -- declared-form SMC on a mixture surface"
  echo "N_REP=${N_REP} M=${M} ITER=${ITER} N_OBS=${N_OBS} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
} > "$LOG"

Rscript run_v3_bkmr.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
{ echo ""; echo "[$(date '+%H:%M:%S')] V33 finished (rc=$rc)"; } >> "$LOG"
exit $rc
