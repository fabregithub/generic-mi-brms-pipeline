#!/usr/bin/env bash
# =============================================================================
# V34 -- does the imputation-induced bias scale with the CENSORED FRACTION?
# -----------------------------------------------------------------------------
# V33 (nd40, 200 reps) measured imputation-induced bias against the ORACLE --
# not against truth, because BKMR is itself +15.1% biased on curvature at these
# settings and every arm inherits that. Worst |bias| over four estimands:
#     smc_dec_over 0.029 | smc_dec_true 0.036 | sub_lod2 0.099 | cc 0.237
#     noint 0.296 | noquad 0.397 | lin 0.460 | pipeline_bkmr 0.501
# Declaring the surface works; the shipped engine destroys half the curvature;
# an omitted basis term is REROUTED onto the terms you kept, not merely lost.
#
# The projection account of that (phase1/DERIVATION_declared_basis.md) says the
# plugin fit of a declared basis is an omitted-variable regression,
#     theta_tilde_ret = theta_ret + Gamma theta_omit,
# so the omitted term is attenuated and the retained ones inflated by
# Gamma*theta_omit -- BOTH to first order LINEAR IN THE CENSORED FRACTION f.
# That is what this run tests.
#
# REGISTERED PREDICTIONS (2026-10-23, before the run):
#   P-C1  bias roughly HALVES from nd40 (40% ND) to nd20 (20% ND), for every
#         arm that has a bias. Specifically, on curv_X1 vs oracle:
#             pipeline_bkmr  -0.501 -> about -0.25
#             smc_dec_lin    -0.460 -> about -0.23
#             smc_dec_noquad -0.397 -> about -0.20
#             smc_dec_noint  +0.296 -> about +0.15
#         A flat or non-monotone result REFUTES the projection account.
#   P-C2  smc_dec_true and smc_dec_over stay near zero at BOTH fractions
#         (|bias| < 0.05). If they drift with f, the declared basis is not in
#         fact congenial and something else is going on.
#   P-C3  nd40_all (ALL THREE exposures censored) is WORSE than nd40 for every
#         biased arm: more of the surface is imputed, so more of it is
#         reconstructed through the declared basis. Direction registered;
#         magnitude not predicted.
#
#   NOT A GATE: sub_lod2_bkmr. Substitution has no declared surface and no
#   Gamma, so the projection argument says nothing about it. Reported for
#   context only -- do not read its f-scaling as evidence either way.
#
#   FLOOR: as V33 -- the oracle must recover a measurably non-zero estimand, or
#   the ratios are undefined.
#
# COST -- nd40 at N_REP=200 measured at 144 min/chunk x 5 = ~12 h (V33).
#   nd20 should be <= that (less censoring, cheaper imputation).
#   nd40_all is UNMEASURED: it censors three exposures instead of one, so its
#   X-block does three draws per sweep. If it runs 1.5-3x, the pair is 28-48 h.
#   The weekend window is ~65 h, and tasks are scenario-major so nd20 finishes
#   FIRST and checkpoints every 44 tasks -- a kill at any point keeps whole
#   chunks. The pilot that would have measured nd40_all was cut short to start
#   this run; that is a deliberate trade of certainty for weekend hours.
#
#   bash run_v34_censoring.sh        # detach; logs/v34_<stamp>.{log,pid}
#   bash run_status.sh v34
#   ESTIMAND=curv_X1 TAG=v34 Rscript analyze_v33.R   (same gates, new data)
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
cd "$(dirname "$SELF")"
TAG=v34
STAMP="$(date '+%Y%m%d-%H%M%S')"
mkdir -p logs results

if [[ -z "${V34_DETACHED:-}" ]]; then
  V34_DETACHED=1 V34_STAMP="$STAMP" nohup bash "$SELF" "$@" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/${TAG}_${STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V34 detached (pid $PID)"
  echo "  log:    logs/${TAG}_${STAMP}.log"
  echo "  status: bash run_status.sh ${TAG}"
  exit 0
fi

STAMP="${V34_STAMP:-$STAMP}"
LOG="logs/${TAG}_${STAMP}.log"
# ⚠️ PIN BLAS TO ONE THREAD PER WORKER. Without this each of the 22 forked
# workers spawns its own BLAS threads (observed 114-160% CPU each on 24 cores)
# and the job thrashes: the 2026-10-19 launch ran >=9 h PER TASK and was killed
# after 18 h having completed zero chunks. Pinned, a 22-rep wave of one SMC arm
# takes 680 s against 411 s uncontended -- a 1.65x contention factor instead of
# ~9x. This is not a tuning preference; it is the difference between feasible
# and not.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export NCORES="${NCORES:-$(( $(sysctl -n hw.ncpu) - 2 ))}"
export N_REP="${N_REP:-200}"
export M="${M:-5}"
export ITER="${ITER:-3000}"
export N_OBS="${N_OBS:-800}"
export SCENARIOS="${SCENARIOS:-nd20,nd40_all}"
export ARMS="${ARMS:-oracle_bkmr,cc_bkmr,sub_lod2_bkmr,pipeline_bkmr,smc_dec_true,smc_dec_noint,smc_dec_noquad,smc_dec_lin,smc_dec_over}"
export OUT_TAG="${OUT_TAG:-${TAG}}"

{ echo "[$(date '+%H:%M:%S')] V34 start -- censoring-fraction sweep"
  echo "N_REP=${N_REP} M=${M} ITER=${ITER} N_OBS=${N_OBS} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
} > "$LOG"

Rscript run_v3_bkmr.R >> "$LOG" 2>&1 &
RPID=$!
printf '%s\n' "$RPID" > "logs/${TAG}_${STAMP}.pid"
wait "$RPID"; rc=$?
{ echo ""; echo "[$(date '+%H:%M:%S')] V34 finished (rc=$rc)"; } >> "$LOG"
exit $rc
