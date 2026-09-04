# V15 findings — the prediction held: attenuation is quadratic in the censored fraction

**Run:** 2026-09-02 07:58 → 19:09 JST (**670.7 min**, 6 scenarios × 2 arms × 100 reps,
`n_ok` = 100 in all 12 cells of the focal estimand, zero task errors, zero non-finite
estimates). Results: `results/v15_latest.rds`, verdict arithmetic in
[`analyze_v15.R`](analyze_v15.R) → `results/v15_prediction_test.csv`. Criteria registered
in advance: `../PLAN_pipeline_validation.md` §8h.

**This is the first track run the other way round.** V0–V14 registered *acceptance
thresholds* — is the pipeline good enough? V15 registered a *number* for four cells nobody
had measured, derived from a law fitted to two, and could have been simply wrong.

**It was not wrong.** All four unmeasured cells landed inside the pre-declared ±10 pp band,
the worst by 9.3 pp and three of four within 2 pp. Re-fitting the exponent freely gives
**1.90**, and the `f` = 0.60 cell — the sharp one, which required the estimate to flip sign
and overshoot — came in at **−116.5%** against **−116.4%** predicted.

---

## Measured against predicted

`curv_X1`, paired excess over `oracle_bkmr` within replication (so BKMR's own oracle floor
cancels — see FINDINGS_v3.md). Registered law: `excess% = −323.4 · f²`.

| cell | `f` | predicted | measured | miss | MC se | inside ±10 pp |
|---|---|---|---|---|---|---|
| `mixf10` | 0.10 | −3.2% | **−1.3%** | +1.9 | 1.2 | yes |
| `nd20` | 0.20 | −12.9% | −13.9% | −1.0 | 2.0 | *(calibration)* |
| `mixf30` | 0.30 | −29.1% | **−29.6%** | −0.5 | 2.3 | yes |
| `nd40` | 0.40 | −51.7% | −55.6% | −3.9 | 3.3 | *(calibration)* |
| `mixf50` | 0.50 | −80.9% | **−90.2%** | −9.3 | 4.6 | yes |
| `mixf60` | 0.60 | −116.4% | **−116.5%** | −0.1 | 5.2 | yes |

**Falsification test 1 — per-cell tolerance:** 4 unmeasured cells tested, **0 outside**.

**Falsification test 2 — is the exponent 2?** Log–log slope **2.47, 95% CI [1.98, 2.97]**;
a rep-level bootstrap that propagates Monte-Carlo error gives **[1.95, 3.67]**. The CI
contains 2, so the test does not reject — but note it is wide and sits above 2. Fitting `k`
and `p` on the percentage-point scale (the fair comparison; the log fit minimises a
different loss) gives **`k` = 315.5, `p` = 1.90**, RMS residual 3.0 pp against the
registered law's 4.2 pp. **Two is not pinned down to better than roughly ±0.5, but nothing
in this run prefers a different integer.**

**Calibration passed.** `nd20` and `nd40` were re-run from V3's seeds and reproduced within
Monte-Carlo error (−13.9 vs −11.7, −55.6 vs −56.7; both inside 2 × MC se). Appending the new
scenarios rather than inserting them preserved canonical indices 1–3, which is what makes
the six points a single comparable sweep rather than two runs stitched together.

---

## What the sweep shows beyond the exponent

**Curvature is not merely lost — past `f` ≈ 0.55 it inverts.** The pipeline's point estimate
falls monotonically from 0.153 to **−0.003** against a truth of 0.137:

| `f` | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|---|---|---|---|---|---|---|
| mean estimate | 0.153 | 0.141 | 0.118 | 0.078 | 0.034 | **−0.003** |
| reps estimating < 0 | 0% | 0% | 0% | 5% | 19% | **40%** |

At 60% non-detects the fitted exposure–response surface curves **the wrong way** in 40% of
replications. That was pre-declared as the easiest way for the prediction to fail; it
happened instead to be what the prediction required.

**The intervals do not hide it — they widen honestly.** Coverage of `curv_X1` stays at
0.92–0.99 across the whole sweep, including the inverted cell, because the interval grows
with `f` (mean width 0.20 → 0.47, FMI 0.04 → 0.46) rather than because the estimate is
right. **Nominal coverage here is not reassurance:** an interval wide enough to cover a
truth of +0.137 from a centre of −0.003 is an interval that says nothing about the shape.
Anyone reading only coverage from this pipeline at high censoring would conclude the
analysis was sound.

---

## What this revises in the theory

`THEORY.md` §4 listed the acceleration of attenuation with `f` as unexplained, with the
observed 4.85× rise for 2× censoring as its only content. It is now a **measured law over
six points, `excess% ≈ 320 · f²`, with the exponent consistent with 2 and inconsistent with
1** (the CI excludes 1 comfortably). The status changes from *unexplained observation* to
**unexplained law** — which is a sharper target: a derivation must now produce an exponent,
not just a direction.

**The leading candidate is a product of two `f`-linear losses**, and it is a *conjecture*,
not a result of this run:

1. **Rows carry no curvature information.** A fraction `f` of X1 values are draws from a
   conditional linear in X1's predictors, which by clause (A) of the congeniality condition
   cannot represent a quadratic in X1. That fraction is `f`.
2. **The estimand's own support is censored.** `curv_X1` is a contrast at the 25th/50th/75th
   percentiles of X1. As `f` grows the detection limit climbs the X1 distribution, so a
   growing share of the *evaluation range* of the contrast lies below the LOD — where the
   only values are imputed. That share also grows about linearly in `f`.

Their product is `f²`. **This conjecture makes a discriminating prediction**: the second
factor depends on where the contrast is evaluated, so **moving the lower evaluation quantile
should move the curve**. Evaluating at `q_lo` = 0.40 instead of 0.25 should raise the
constant and bend the shape near `f` = `q_lo`, while the first factor alone predicts no
dependence on `q_lo` at all. That is the natural next f-sweep — now queued as V18, behind V16 and V17 — and it reuses this exact harness (the
quantiles are already configurable — `bkmr_estimand_grid()`'s `q_lo`/`q_hi`), and the two
mechanisms predict visibly different sweeps.

---

## Caveats (scope of this run)

- **The constant is still fitted, not derived.** 320 has no derivation; only the exponent
  now has evidence. The law is calibrated to *this* DGP (`b_quad` = 0.15 on the V3 mixture surface),
  `n` = 800, `m` = 10, and BKMR with 50 knots. Nothing here says the constant transfers.
- **`f` is the focal exposure's censored fraction only.** X2 and X3 are fully observed in the
  `mixf*` cells, because `curv_X1` depends on X1 alone — V3 established that censoring all
  three adds little (−56.7% → −58.4%). A sweep in *joint* censoring is a different curve.
- **100 reps, MC error 1.2–5.2 pp per cell.** Adequate against the 17–35 pp gaps between
  candidate laws, and against the ±10 pp gate; not adequate to distinguish `p` = 1.9 from
  `p` = 2.0, which this run does not attempt to do.
- **The exponent CI sits above 2** (point 2.47 on the log scale, 1.90 on the pp scale). If a
  derivation predicts exactly 2 that is compatible with these data; if one predicts 2.5 it is
  also compatible. Distinguishing them needs cells at `f` < 0.10, where the excess is small
  and MC error dominates — a much more expensive run than this one.
- **No claim about the pipeline's fitness changes.** V3 already ruled the mixture path
  unusable under censoring; V15 says only *how fast* it degrades. The remedy remains the
  V14 grid draw, whose mixture-estimand half is untested (roadmap item 07).
