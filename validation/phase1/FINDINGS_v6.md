# V6 findings — forest vs parametric, judged on a DGP that can punish both

**Run:** 2026-08-26 17:16 → 2026-08-27 01:28 JST (**8.2 h**, 2100 tasks, 22 fork workers).
Closes roadmap item **05**. **Zero worker failures, zero task errors, 300/300 reps.**

**Design:** 7 scenarios × 3 arms × 300 reps, `m`=30, seed 20260825. Three diagnostic cells
are run under **both** covariate designs — `mcar_z40`/`nl_mcar_z40`,
`missing_y20`/`nl_missing_y20`, `combined`/`nl_combined` — so the same contrast is
measured with only the covariate structure changed. `base` is the negative control.

**The question.** V5 shipped `properBoot` (bootstrap the forest's training data per
imputation) and it closed only 23% of the calibration gap under heavy missingness. The
obvious alternative — `mice`'s parametric proper draws — calibrated better in V4/V5, but
that evidence came from a DGP whose covariates are entirely linear and Gaussian, so the
parametric model was *correctly specified* and could never be penalised for its one real
risk. This run fixes that: under `z_form = "nonlinear"`, a linear imputation model for `Z1`
leaves **0.387 of explainable R²** unclaimed (against 0.001 under `linear`), while the
oracle still recovers `b_logX1` = 0.392 — the imputation problem got harder, the estimand
did not move.

---

## Results

| scenario | arm | rel. bias | width/SE | `B` | coverage |
|---|---|---|---|---|---|
| **base** (control) | block_fcs | −1.08% | 1.008 | 0.00077 | 0.943 |
| | properBoot | −1.08% | 1.008 | 0.00077 | 0.943 |
| | micePmm | −1.08% | 1.008 | 0.00077 | 0.943 |
| **mcar_z40** | block_fcs | −1.64% | 0.895 | 0.00078 | 0.920 |
| | **properBoot** | −1.57% | **0.966** | 0.00096 | 0.947 |
| | micePmm | **−3.18%** | 1.109 | 0.00111 | 0.947 |
| **nl_mcar_z40** | block_fcs | −1.82% | 0.947 | 0.00072 | 0.940 |
| | **properBoot** | −1.91% | **0.975** | 0.00078 | 0.937 |
| | micePmm | **−4.70%** | 1.063 | 0.00109 | 0.953 |
| **missing_y20** | block_fcs | +2.81% | 0.960 | 0.00070 | 0.917 |
| | properBoot | +2.67% | 0.975 | 0.00077 | 0.937 |
| | **micePmm** | **−0.83%** | 0.983 | 0.00072 | 0.940 |
| **nl_missing_y20** | block_fcs | +2.26% | 0.993 | 0.00071 | 0.943 |
| | properBoot | +2.09% | 1.002 | 0.00081 | 0.947 |
| | **micePmm** | **−1.14%** | 1.034 | 0.00073 | 0.963 |
| **combined** | block_fcs | +2.43% | 0.885 | 0.00068 | 0.910 |
| | properBoot | +2.02% | 0.921 | 0.00081 | 0.927 |
| | **micePmm** | **−0.15%** | 1.011 | 0.00091 | 0.953 |
| **nl_combined** | block_fcs | +3.10% | 0.986 | 0.00067 | 0.943 |
| | properBoot | +2.56% | 1.010 | 0.00079 | 0.957 |
| | **micePmm** | **−0.32%** | 1.084 | 0.00087 | 0.963 |

**The control passes.** In `base` all three arms return *identical* numbers to five
decimals — with no covariate or outcome missingness the Z block is idle, so all three
imputers are no-ops. That is what makes every difference below attributable.

---

## The misspecification penalty is real, and specific to `micePmm`

Changing **only** the covariate structure, same arm, same seeds:

| arm | mcar_z40 → nl_mcar_z40 | penalty |
|---|---|---|
| block_fcs | −1.64% → −1.82% | 0.18 pp |
| properBoot | −1.57% → −1.91% | 0.34 pp |
| **micePmm** | −3.18% → **−4.70%** | **1.52 pp** |

`micePmm` is **4–8× more sensitive** to the non-linearity than either forest arm, and its
bias breaches the pre-registered 3% bar in both `mcar_z40` cells. This is exactly the cost
the old linear DGP could not detect, and it vindicates building item 05 rather than
switching to `mice` on V5's evidence.

## But `micePmm` is better in four of six cells

Paired mean-absolute-error against `properBoot`, same datasets:

| scenario | micePmm vs properBoot | t |
|---|---|---|
| mcar_z40 | better by 0.41 pp | −1.6 (ns) |
| **nl_mcar_z40** | **worse by 0.62 pp** | **+2.5** |
| missing_y20 | better by 0.19 pp | −0.9 (ns) |
| nl_missing_y20 | better by 0.40 pp | −2.0 |
| combined | better by 0.66 pp | −3.3 |
| nl_combined | better by 0.59 pp | −2.8 |

**The only cell `micePmm` loses is the one built to punish it.** There is a clean mechanism:

- In `mcar_z40` the Z block's work is dominated by imputing **`Z1`**, whose true model is
  non-linear → `pmm` is misspecified → bias.
- In `missing_y20` and `combined` the work is dominated by imputing **`Y`**, which *is*
  genuinely linear in `(logX, Z)` → `pmm` is correctly specified → its proper parameter
  draws win outright.

This also explains why `micePmm` survives `combined` (20% covariate missingness) but not
`mcar_z40` (40%): the penalty is dose-dependent on how much `Z1` imputation the Z block has
to do.

One further nuance worth recording: in several cells `micePmm` shows **larger bias but
smaller absolute error** — lower variance partly offsetting a worse centre. Bias was the
registered criterion so the verdict stands, but "worse bias" and "less accurate" are not
the same statement here.

---

## Verdict — criterion 2 met: keep `properBoot`

Judged against the criteria registered in `ROADMAP.md` before the run:

| criterion | outcome |
|---|---|
| 1 — adopt `mice`: calibration edge **and** bias ≤3% | **fails** — bias −4.70%, and it overshoots to 1.06–1.11 rather than holding an edge |
| 2 — keep `properBoot`: `micePmm` bias >3% or materially worse | **met**, in both `mcar_z40` cells |
| 3 — neither suffices: `properBoot` still under-covers | **does not apply** — see below |

**`properBoot` stays the default.** Bounded bias everywhere (max 2.56%), width/SE nearest
nominal in four of six cells, coverage ≥0.927 in every cell, and near-insensitivity to
covariate specification. `micePmm` is more accurate when imputation is outcome-dominated
but can be materially worse when it is covariate-dominated and non-linear — and which
regime real data occupies is not knowable in advance.

### A prediction that did not hold

Before this run the expectation stated in `ROADMAP.md` and in conversation was outcome 3 —
that `properBoot` would still under-cover and neither family would suffice. **It does not
under-cover here.** Coverage is 0.947 / 0.937 / 0.937 / 0.947 / 0.927 / 0.957 across the
six diagnostic cells, passing the 0.92 bar everywhere, while `micePmm` buys its coverage
with over-wide intervals. Recorded because the prediction was on the record.

The *underlying* point does survive in weaker form: neither arm is uniformly best.
`properBoot` carries +2.0 to +2.6% bias in the `combined` cells where `micePmm` sits near
zero; `micePmm` breaks where `properBoot` is steady. A component that is **flexible *and*
properly dispersed** would take the better half of each. Item 05 closes the *choice between
existing options*; it does not close that question.

---

## Two methodological notes

**1. An arm's numbers are only comparable across runs with an identical arm set.**
`properBoot`'s `combined` coverage reads 0.927 here against 0.910 in V5. Same seed, same
datasets — but procedures consume the task's RNG stream *sequentially*, so dropping the
`properZ` / `noMID` arms changed the stream reaching `properBoot`. The datasets are
identical (verified in V4 via the deterministic arms, `max |diff| = 0`); the imputation
draws are not. 0.910 vs 0.927 is ~0.9 SE apart, i.e. noise — but cross-run comparison of a
single arm requires the same arms present.

**2. The non-linear design is not uniformly "harder".** Calibration *improved* for every
arm under it — `block_fcs` width/SE 0.947 / 0.993 / 0.986 against 0.895 / 0.960 / 0.885.
A less predictable `Z1` generates more genuine imputation variability, raising `B`. So
cross-design calibration comparisons are not like-for-like; only the within-design
contrasts carry the argument, which is how the tables above are read.

---

## Caveats (scope of this run)

1. **`pmm` only.** `mice` offers other proper methods (`norm`, `rf`, `cart`). `rf` is
   essentially `properBoot`; `norm` should behave like `pmm`; `cart` is untested and is the
   one plausibly-flexible-and-proper option not tried.
2. **One non-linear form.** `Z1` is a specific tanh + quadratic + interaction construction.
   Other non-linearities may penalise differently.
3. **Additive ERF, MCAR, one estimand** — unchanged limitations from V2/V5.
4. **`Y`'s model is linear by construction**, which is precisely why `pmm` wins the
   outcome-dominated cells. A non-linear outcome model would likely remove that advantage,
   and is untested.
5. **No runtime comparison.** `micePmm` is much cheaper than either forest arm; if it had
   won on accuracy the cost argument would have reinforced it. It did not, so the point is
   moot here — but any future candidate must be costed.
