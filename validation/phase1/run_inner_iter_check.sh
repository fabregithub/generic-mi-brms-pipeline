#!/usr/bin/env bash
# =============================================================================
# Validation of the block-FCS inner-iteration fix.
# -----------------------------------------------------------------------------
# THE FIX. `run_censored_exposure_block_fcs()` alternates the Z and X blocks for
# `outer_sweeps` iterations. The BART Z-block imputer also ran its own
# `bart_inner_iter` FCS iterations, so the two multiplied: 3 outer x 3 inner = 9
# alternations where 3 were designed, and ~3x the BART fits for nothing.
# `.ce_one_imputation()` now defaults the inner count to 1 when called from
# inside block-FCS (an explicit config value still wins).
#
# THE QUESTION. Fewer alternations *should* leave calibration and bias unchanged,
# because the outer loop still provides three. That is an assumption, and this
# run tests it rather than trusting it.
#
# Scenarios span both shapes the Z block takes inside block-FCS: mcar_z40 /
# nl_mcar_z40 have TWO targets (Z1, Z2, no MID) and combined / nl_combined have
# THREE (Z1, Z2, Y, with MID). The inner-iteration count only bites on the
# multi-target path, so both are worth covering.
#
#   pipeline_bartMI        pipeline's own BART, inherits the new default (inner = 1)
#   pipeline_bartMI_iter3  pipeline's own BART forced to inner = 3 (pre-fix)
#   pipeline_bartHarness   BART via the harness instrument, as V7 actually ran it
#
# The third arm is the PORT CHECK. V7's BART arm overrode
# run_row_level_imputation with the harness instrument, so it measured the
# instrument rather than the code later shipped in 00_common_functions.R. A
# 6-rep check found the two bit-identical; this puts that on the record at 200
# reps. bartHarness and bartMI_iter3 both use inner = 3, so they should agree
# exactly -- any divergence means the port is not faithful.
#
# Everything else is identical, and both arms see the same datasets, so any
# difference is attributable to the inner-iteration count alone.
#
# ACCEPTANCE (registered before the run):
#   * bias:      |bartMI - bartMI_iter3| <= 0.5 pp, paired, in both scenarios
#   * coverage:  within 0.02 of each other
#   * width/SE:  within 0.03 of each other. MEASURED cost is ~48 s per replication
#                for all three arms, so reps are cheap here -- and N_REP is set to
#                1000 specifically so this criterion is meaningful: width/SE has
#                MC error 1/sqrt(2(N-1)), which is 0.050 at 200 reps (useless
#                against a 0.03 band, the V7 mistake) but 0.022 at 1000. See
#                FINDINGS_v7.md for why that band was unattainable there.
#   * runtime:   inner=1 should be materially cheaper; report the ratio
#   * port:      bartMI_iter3 vs bartHarness bit-identical (both inner = 3)
# If those hold, the fix is a free ~3x saving on the censored path and the
# default stands. If bias or coverage moves, the fix is reverted and the outer
# sweeps are doing less work than assumed -- which would itself be worth knowing.
#
#   ./run_inner_iter_check.sh              # 4 scenarios x 1000 reps, ~2 h (measured)
#   N_REP=200 ./run_inner_iter_check.sh    # ~25 min, but width/SE becomes noise
#   SCENARIOS=combined,nl_combined ./run_inner_iter_check.sh   # 3-target cells only
#
# Results: results/inneriter_latest.rds
# =============================================================================
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

if [[ "${II_DETACHED:-0}" != "1" ]]; then
  mkdir -p logs
  export II_STAMP="$(date +%Y%m%d-%H%M%S)"
  II_DETACHED=1 nohup bash "$SELF" >/dev/null 2>&1 &
  PID=$!
  echo "$PID" > "logs/inneriter_${II_STAMP}.pid"
  command -v caffeinate >/dev/null 2>&1 && nohup caffeinate -i -w "$PID" >/dev/null 2>&1 &
  echo "inner-iteration check detached.  PID ${PID}"
  echo "Log:   validation/phase1/logs/inneriter_${II_STAMP}.log"
  echo "Watch: tail -f validation/phase1/logs/inneriter_${II_STAMP}.log"
  echo "Alive? bash validation/phase1/run_status.sh inneriter"
  exit 0
fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd ../.. && pwd)}"
export NCORES="${NCORES:-22}"
export SEED="${SEED:-20260825}"
export N_REP="${N_REP:-1000}"
export M="${M:-30}"
export SWEEPS="${SWEEPS:-3}"
export MARGIN="${MARGIN:-shash}"
export OUT_TAG="inneriter"
export SCENARIOS="${SCENARIOS:-mcar_z40,nl_mcar_z40,combined,nl_combined}"
export ARMS="pipeline_bartMI,pipeline_bartMI_iter3,pipeline_bartHarness"

STAMP="${II_STAMP:-$(date +%Y%m%d-%H%M%S)}"
LOG="logs/inneriter_${STAMP}.log"
mkdir -p results logs

{
  echo "=== inner-iteration fix validation ==="
  echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "N_REP=${N_REP} M=${M} NCORES=${NCORES} SWEEPS=${SWEEPS}"
  echo "SCENARIOS=${SCENARIOS}"
  echo "ARMS=${ARMS}   (bartMI = inner 1 [fixed], bartMI_iter3 = inner 3 [pre-fix])"
  echo "======================================"
} > "$LOG"

Rscript run_v4_variance.R >> "$LOG" 2>&1
echo "[$(date '+%H:%M:%S')] done rc=$? -> results/inneriter_latest.rds" >> "$LOG"
