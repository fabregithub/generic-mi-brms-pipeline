# V8 findings — the inner-iteration fix is half right, and the V7 port is exact

**Run:** 2026-08-28 07:35 → 09:50 JST (**135.3 min**, 4 scenarios × 3 arms × 1000 reps,
4000/4000 tasks, `n_ok` = 1000 in all 12 cells, zero worker failures).
Command: `bash validation/phase1/run_inner_iter_check.sh`. Results:
`results/inneriter_latest.rds`.

Two independent questions, one run:

| # | question | verdict |
|---|---|---|
| **Q1** | Does defaulting the Z block's inner FCS count to 1 inside block-FCS leave bias and calibration unchanged? | ❌ **No — in the 3-target cells.** Fix narrowed, not reverted |
| **Q2** | Does the BART code shipped in `00_common_functions.R` reproduce the harness instrument V7 actually measured? | ✅ **Yes — bit-identical on all 4000 rows** |

Estimand throughout: `b_logX1`, true value 0.4. `m` = 30, `outer_sweeps` = 3,
`margin` = shash.

---

## The arms

| arm | Z-block inner FCS iterations | what it is |
|---|---|---|
| `pipeline_bartMI` | 1 | the pipeline's own BART, inheriting the fix under test |
| `pipeline_bartMI_iter3` | 3 | the same code path, forced to the pre-fix count |
| `pipeline_bartHarness` | 3 | BART via the harness instrument — **as V7 actually ran it** |

All three see the same simulated datasets in every replication, so differences are
attributable to the inner count alone. `bartMI_iter3` and `bartHarness` differ only in
*which implementation* runs, and both use inner = 3 — that pairing is the port check.

---

## Results

| scenario | Z-block targets | arm | rel. bias | width/SE | coverage | `B` |
|---|---|---|---|---|---|---|
| **mcar_z40** | 2 (Z1, Z2) | oracle | −0.327% | 1.0053 | 0.950 | — |
| | | `bartMI` (inner 1) | −1.668% | 1.0539 | 0.957 | 0.00101 |
| | | `bartMI_iter3` | −1.710% | 1.0578 | 0.960 | 0.00101 |
| | | `bartHarness` | −1.710% | 1.0578 | 0.960 | 0.00101 |
| **nl_mcar_z40** | 2 (Z1, Z2) | oracle | −0.415% | 1.0217 | 0.953 | — |
| | | `bartMI` (inner 1) | −1.378% | 1.0630 | 0.958 | 0.00089 |
| | | `bartMI_iter3` | −1.374% | 1.0720 | 0.957 | 0.00090 |
| | | `bartHarness` | −1.374% | 1.0720 | 0.957 | 0.00090 |
| **combined** | 3 (Z1, Z2, **Y** + MID) | oracle | −0.040% | 1.0087 | 0.955 | — |
| | | `bartMI` (inner 1) | **−0.876%** | 1.0562 | 0.959 | 0.00096 |
| | | `bartMI_iter3` | **−0.157%** | 1.0426 | 0.957 | 0.00096 |
| | | `bartHarness` | −0.157% | 1.0426 | 0.957 | 0.00096 |
| **nl_combined** | 3 (Z1, Z2, **Y** + MID) | oracle | +0.611% | 1.0316 | 0.961 | — |
| | | `bartMI` (inner 1) | +0.027% | 1.0610 | 0.960 | 0.00092 |
| | | `bartMI_iter3` | +0.381% | 1.0528 | 0.964 | 0.00091 |
| | | `bartHarness` | +0.381% | 1.0528 | 0.964 | 0.00091 |

---

## Q1 — the inner-iteration fix

### Read the difference **paired**, or you will read noise

The registered criterion is a *paired* bias difference, and that wording carries the whole
result. The arms are run on the same datasets and their estimates correlate at **0.987–0.988**,
so the paired standard error is ~6× smaller than the unpaired one. Reading the arm-level
Monte-Carlo errors instead would have declared every cell inconclusive:

| scenario | unpaired diff | naive MC SE (per arm) | **paired diff** | **paired MC SE** |
|---|---|---|---|---|
| combined | −0.719 pp | ±0.41 | **−0.719 pp** | **±0.064** |
| nl_combined | −0.354 pp | ±0.41 | **−0.354 pp** | **±0.064** |
| mcar_z40 | +0.042 pp | ±0.39 | +0.042 pp | ±0.065 |
| nl_mcar_z40 | −0.004 pp | ±0.38 | −0.004 pp | ±0.060 |

### Against the pre-registered criteria

Percentage points of the true estimand; `bartMI` (inner 1) minus `bartMI_iter3` (inner 3).

| scenario | targets | Δ bias (pp) | 95% CI | *p* | vs 0.5 pp bar |
|---|---|---|---|---|---|
| **combined** | 3 incl. Y | **−0.719** | [−0.845, −0.593] | <1e-6 | ❌ **breach — whole CI outside the bar** |
| nl_combined | 3 incl. Y | −0.354 | [−0.479, −0.229] | <1e-6 | ✅ passes, but a real, resolved effect |
| mcar_z40 | 2, no Y | +0.042 | [−0.085, +0.168] | 0.52 | ✅ indistinguishable from zero |
| nl_mcar_z40 | 2, no Y | −0.004 | [−0.122, +0.114] | 0.95 | ✅ indistinguishable from zero |

| criterion | result | |
|---|---|---|
| bias: \|inner1 − inner3\| ≤ 0.5 pp, paired, in all scenarios | −0.719 pp in `combined` | ❌ |
| coverage: within 0.02 of each other | max \|Δ\| = 0.004 | ✅ |
| width/SE: within 0.03 of each other | max \|Δ\| = 0.0136 | ✅ |
| port: `bartMI_iter3` vs `bartHarness` bit-identical | exact, 4000/4000 rows | ✅ |
| runtime: report the inner1/inner3 ratio | **not measurable from this run** | ⚠️ |

**The bias criterion fails, and the failure is not noise.** In `combined` the entire
confidence interval lies beyond the bar. Coverage and width/SE both pass comfortably — so
this is a **bias-only** failure that the interval diagnostics do not reveal. A run that
checked only calibration would have shipped it.

### Why: the outer loop alternates *between* blocks, not *among* the Z block's targets

The fix rested on an assumption stated in the script header — the outer block-FCS loop
already alternates three times, so the Z block needs no inner iterations of its own. That
assumption is only half true. The outer loop alternates **Z ↔ X**. It does nothing to make
the Z block's own targets condition on *each other*: with `inner = 1`, `run_row_level_imputation_bart()`
does a single pass over its targets, so a target imputed early in the pass is drawn from an
initialisation fill of the ones after it, and never revisited within that call.

The evidence lines up with that mechanism exactly:

- **Two targets → no effect** (+0.042, −0.004 pp; both CIs contain zero). One pass is nearly
  enough when there is almost nothing to chain to.
- **Three targets → a resolved negative effect** in both cells (−0.719, −0.354 pp), and
  `inner = 1` is the more biased arm in each.

The effect scales with the number of targets, which is what an unconverged within-block
chain looks like.

### The fix is narrowed, not reverted

The registered rule for a moved bias was to revert. A plain revert is available and safe,
but it throws away a saving that this run *measured* as free: in the two-target cells,
`inner = 1` costs +0.04 pp of bias with a ±0.13 pp interval, at a third of the BART fits.
So the default is now chosen per sweep from the Z block's actual target set rather than
fixed either way — see **What changed** below. The pre-fix `inner = 3` is restored
everywhere the run found a cost, including every configuration the run did not test.

---

## Q2 — the port check passes exactly

V7's BART arm **overrode** `run_row_level_imputation` with the harness instrument, so V7
measured the instrument rather than the code later shipped in `00_common_functions.R`. That
left V7's conclusions formally unattached to the shipped pipeline. A 6-rep check had found
the two identical; this puts it on the record at 1000 reps per cell.

`pipeline_bartMI_iter3` (shipped code) and `pipeline_bartHarness` (V7's instrument) are
**bit-identical** — `identical() == TRUE`, max \|difference\| = 0 — on all 4000 rows, for
every quantity recorded:

| quantity | identical | max \|diff\| | n |
|---|---|---|---|
| `estimate` | TRUE | 0 | 4000 |
| `se` | TRUE | 0 | 4000 |
| `ci_lo` | TRUE | 0 | 4000 |
| `ci_hi` | TRUE | 0 | 4000 |
| `ubar` | TRUE | 0 | 4000 |
| `b` | TRUE | 0 | 4000 |

**V7's findings therefore transfer to the shipped code without qualification.** This closes
the one methodological gap left open by V7's arm construction.

---

## What changed in the pipeline (v1.5.1)

**`.ce_default_inner_iter(targets, y_var)`** in `00_censored_exposure.R` picks the Z block's
inner FCS count per sweep, from the target set that sweep actually has:

```r
if (length(targets) <= 2L && !(y_var %in% targets)) 1L else 3L
```

`.ce_one_imputation()` calls it inside the sweep loop, right after
`make_row_level_imputation_spec()` — the target set is not known before that. An explicit
`analysis_spec$imputation$bart_inner_iter` still wins, and is read from `analysis_spec`
rather than the mutated `z_spec_as` copy so a per-sweep assignment cannot be mistaken for a
user setting on the next sweep.

**`bart_inner_iter = 3` was removed from the shipped `00_config.R`** (left commented, with
the reason). It had to be: an explicit config value overrides the per-sweep choice, so
leaving it set would have made the new default dead code for anyone using the stock config —
and pinned the two-target path to `inner = 3` for no benefit. Off the censored-exposure path
the fallback in `run_row_level_imputation_bart()` is unchanged at 3.

**Y is excluded from the cheap path deliberately**, beyond what the target count alone would
require. The X block conditions on Y, so an under-conditioned Y feeds straight into the
censored exposure draw — the one path by which a Z-block shortcut can reach the estimand.
See caveat 1: the run cannot separate "3 targets" from "Y is a target", and this predicate
is the reading that stays safe under either explanation.

**Net effect on a stock censored-exposure analysis:** none. That configuration imputes Y
(MID), so it takes the `inner = 3` branch — i.e. the pre-fix behaviour, which is what V7
validated. The saving applies only to blocks with ≤2 non-outcome targets.

### Verification of the change

The new default was checked against the stored 1000-rep results by re-running reps 1-6 of
both scenario shapes with the edited code (`OUT_TAG=v8check`, 0.6 min). Seeding is keyed on
the canonical scenario index (`base_seed + canon_index*100000 + r`), so a subset run
reproduces the full run's replications exactly -- which makes this a bit-identity test, not
a statistical one:

| cell | Z-block targets | new `bartMI` vs stored arm | result |
|---|---|---|---|
| `mcar_z40` | 2, no Y | vs old `pipeline_bartMI` (inner 1) | **identical**, max abs diff = 0 |
| `combined` | 3 incl. Y | vs old `pipeline_bartMI_iter3` (inner 3) | **identical**, max abs diff = 0 |
| `combined` | 3 incl. Y | vs old `pipeline_bartMI` (inner 1) | **differs** (max abs diff = 0.017) -- as required |

Compared on `estimate`, `se`, `ci_lo`, `ci_hi`, `ubar`, `b`. The predicate routes correctly
in both directions: the two-target path is untouched, the three-target path now takes the
validated `inner = 3` behaviour, and it is no longer the arm that failed.
`.ce_default_inner_iter()` was also unit-checked over five target sets, including the two
untested configurations that fall on the safe side.

---

---

## Caveats (scope of this run)

1. **The design confounds "3 targets" with "Y is a target."** Both 3-target cells impute Y
   with MID; neither 2-target cell imputes Y. So the run establishes *that* the cheap path
   is unsafe in the 3-target/with-Y configuration and *that* it is free in the
   2-target/without-Y one, but cannot say which feature is responsible. The shipped
   predicate requires **both** conditions for `inner = 1`, so it is correct under either
   explanation — at the cost of running `inner = 3` in the untested 2-target-with-Y and
   3-target-without-Y configurations. Separating the two needs a cell with three non-outcome
   targets; it is a small follow-up, and only widens the cheap path if run.
2. **The threshold `≤ 2` is a boundary, not a measured optimum.** Two targets were measured
   free and three were measured costly. Four or more were not tested, and are covered only
   because they fall on the safe side.
3. **The runtime saving is asserted, not measured here.** `secs` is logged once per
   replication for all arms together (identical across arms in `inneriter_raw.csv`), so the
   ~3× figure comes from the fit count, not from timing this run. Quantifying it needs
   per-arm instrumentation.
4. **`nl_combined` passes the bar but is not clean.** Its −0.354 pp difference is resolved
   far from zero (*p* < 1e-6). It passes only because the bar is 0.5 pp. Treated as
   corroborating the mechanism, not as a cell where `inner = 1` is safe.
5. **The interval diagnostics were blind to this.** Coverage and width/SE agreed to within
   0.004 and 0.014 while bias moved 0.72 pp. Calibration checks do not substitute for a
   paired bias check on this path.
6. **One estimand, `m` = 30, `outer_sweeps` = 3, shash margin.** The interaction between the
   inner count and `outer_sweeps` was not swept — a longer outer loop might compensate for
   `inner = 1`, and a shorter one might make it worse.
7. **MCAR, additive ERF, linear outcome model** — unchanged limitations inherited from
   V2/V5/V6/V7.
