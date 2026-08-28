#!/usr/bin/env bash
# =============================================================================
# Track V3 -- BKMR mixture estimands (roadmap 04, design-plan Phase 2b, R11)
# -----------------------------------------------------------------------------
# THE QUESTION. Phase 1 §7.7 showed linear congenial imputation stays biased on a
# mixture surface -- but measured it on a SCAFFOLD estimand (the `logX1`
# coefficient of a matched linear model), which nobody reports from a BKMR
# analysis. The claims ledger records that verdict as "mechanism demonstrated,
# verdict not manuscript-final". This run replaces the scaffold with the
# estimands people actually report, each with its truth derived ANALYTICALLY
# from the generator.
#
# THE ARMS
#   oracle_bkmr     complete data          -- a GATE on the harness, not the yardstick
#   cc_bkmr         complete case
#   sub_lod2_bkmr   LOD/sqrt(2)            -- standard applied practice
#   pipeline_bkmr   shipped block-FCS, M imputations, BKMR each, Rubin pooled
#
# THE SEVEN ESTIMANDS. overall q75/q50, overall q25/q50, overall q75/q25,
# single-variable X1 at two backgrounds, the X1:X2 interaction, and the curvature
# in X1. The last two are the sharp tests: their analytic truths contain nothing
# but `b_int` and `b_quad`, which a linear imputation conditional cannot
# represent. `overall_q75_q25` is DEGENERATE on this generator (symmetric
# quantiles cancel both non-linear terms) and must not be read as "the mixture
# estimand" -- see ../PLAN_pipeline_validation.md §8.
#
# ACCEPTANCE (full text, and the record of the gate change, in the plan §8)
#   GATE  oracle_bkmr: |rel bias| <= 10% and coverage in [0.90, 0.98] on all seven.
#         A gate failure INVALIDATES the run -- diagnose before reading any arm.
#   MAIN  pipeline_bkmr on int_X1X2 and curv_X1, with Monte-Carlo error, AND the
#         excess over oracle_bkmr on the same estimand.
#         "Materially biased" is |rel bias| > 10% or cov < 0.90.
#   This track is DESCRIPTIVE. It establishes the scope of the mixture
#   restriction; it does not gate a release, and either answer is a result.
#
# THE GATE WAS WIDENED FROM 5% TO 10% BEFORE THE RUN, on measured evidence, not
# after seeing a result we disliked. A 44-rep oracle-only sweep at iter
# 1000/3000/8000 found a STABLE +6-7% bias on curv_X1, overall_q75_q50 and
# singvar_X1_q75 that does not shrink with chain length -- BKMR's own
# finite-sample bias on this surface at n=800, not an MCMC artifact. The original
# 5% bar was set before that floor was known and no procedure could have met it.
# Because the floor exists, the pipeline arm is read as EXCESS OVER THE ORACLE as
# well as raw bias; charging the imputation for error BKMR makes on complete data
# would overstate the failure this track is looking for.
#
# COST (MEASURED on this machine: 24 cores, n=800, 50 GPP knots, NCORES=22)
#   one BKMR fit, unloaded ....... 16-17 s at iter=1000 (full GP is 165 s -- 10x)
#   fits per replication ......... 3 + M   (three single-fit arms, then M for the MI arm)
#
#   MEASURED wall-clock, 22 replications running concurrently on 22 workers:
#     oracle only (1 fit/rep), iter=1000 ....  1.7 min per 44 reps
#     oracle only (1 fit/rep), iter=3000 ....  4.4 min per 44 reps
#     oracle only (1 fit/rep), iter=8000 .... 10.1 min per 44 reps
#     all 4 arms, M=10, iter=1000 ........... 10.7 min per 22 reps
#
#   Note the contention factor: 13 fits that take 17 s each unloaded take ~49 s
#   each with 22 workers competing. Do not cost this run from the unloaded
#   per-fit time -- it is ~3x optimistic.
#
#   DERIVED from those measurements (iter=3000 is 2.6x iter=1000):
#     all 4 arms, M=10, iter=3000 ........... ~28 min per 22 reps
#     registered run, 3 scenarios x N_REP=200 = 600 tasks .......... ~12.7 h
#     two-scenario variant (600 -> 400 tasks) ...................... ~8.5 h
#   Treat 12.7 h as a LOWER BOUND: it assumes nd40_all costs the same per task as
#   nd40, and it does not -- its X block imputes three censored exposures per
#   sweep instead of one. BKMR dominates the per-task cost (13 fits), so the
#   overrun should be modest, but it has not been measured.
#
#   If that is too long, in order of preference:
#     N_OBS=400   roughly halves it (GPP cost is ~linear in n)
#     M=5         0.62x (8 fits/rep instead of 13)
#     SCENARIOS=nd40   halves it; nd40 is the cell with the most censoring
#   Scale N_REP LAST: this track's headline is BIAS, and N_REP is what buys
#   precision on it. M is the cheap thing to cut -- bias is insensitive to m.
#
# BEFORE RUNNING: `Rscript test_v3_harness.R` -- 16 checks, no MCMC, ~1 s. It
# verifies that `h_true()` reproduces the generator's own `eta` exactly, which is
# the assumption every number in this track rests on.
#
# USAGE
#   ./run_v3_bkmr.sh                          # the registered run (see defaults below)
#   N_REP=10 ./run_v3_bkmr.sh                 # a short pilot first -- recommended
#   SCENARIOS=nd40 ARMS=oracle_bkmr ./run_v3_bkmr.sh    # gate only, cheapest check
#   SCENARIOS=nd40,nd40_all ./run_v3_bkmr.sh            # drop nd20, keep the two 40% cells
#   SCENARIOS=nd20,nd40 ./run_v3_bkmr.sh                # the cheaper two-cell set
#   M=30 ./run_v3_bkmr.sh                     # the pipeline's validated m, ~3x the cost
#   KNOTS=0 ./run_v3_bkmr.sh                  # full GP, ~10x slower; only to check the GPP
#   ARMS=oracle_bkmr N_REP=44 ITER=8000 ./run_v3_bkmr.sh   # re-check the oracle floor
#
# Results: results/v3_latest.rds, v3_summary.csv, v3_raw.csv
# Checkpoint after every chunk: results/v3_checkpoint.rds (readable mid-run with
#   summarise_v3(readRDS("results/v3_checkpoint.rds")) -- resume is NOT trusted,
#   see the note in run_v4_variance.R).
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

if [[ "${V3_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export V3_STAMP="$(date +%Y%m%d-%H%M%S)"
  V3_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/v3_${V3_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "V3 BKMR run detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/v3_${V3_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/v3_${V3_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh v3"
  exit 0
fi

# BLAS must stay single-threaded or it fights the R-level fork loop.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260828}"
export N_REP="${N_REP:-200}"
export N_OBS="${N_OBS:-800}"
export ITER="${ITER:-3000}"
export M="${M:-10}"
export KNOTS="${KNOTS:-50}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export OUT_TAG="${OUT_TAG:-v3}"
export SCENARIOS="${SCENARIOS:-nd20,nd40,nd40_all}"
export ARMS="${ARMS:-oracle_bkmr,cc_bkmr,sub_lod2_bkmr,pipeline_bkmr}"

STAMP="${V3_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/v3_${STAMP}.log"
mkdir -p results logs

{
  echo "=== Track V3: BKMR mixture estimands ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} N_OBS=${N_OBS} M=${M} ITER=${ITER} KNOTS=${KNOTS} NCORES=${NCORES}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}"
  echo "GATE: oracle_bkmr must show |rel bias| <= 10% and coverage in [0.90, 0.98]"
  echo "      on all seven estimands, or the run says nothing about the other arms."
  echo "========================================"
} > "$LOG"

Rscript run_v3_bkmr.R >> "$LOG" 2>&1
echo "[$(date '+%H:%M:%S')] done rc=$? -> results/${OUT_TAG}_latest.rds" >> "$LOG"
