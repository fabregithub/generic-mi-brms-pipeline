# Which covariates to adjust for, and what the imputation does with them

**Short version.** The pipeline treats every row marked `use_in_model = TRUE` the same way,
but your covariates do not all play the same causal role — and the role decides both whether
adjusting for the variable is correct at all, and how much bias the imputation adds.
Validation measured four things worth knowing before you run:

1. **Adjusting for a collider inverts the exposure effect** — from +0.40 to −0.10 in a
   measured case, on complete data, before imputation enters. No imputation method repairs
   this; it is a model-specification error and the pipeline's `use_in_model = TRUE`
   convention invites it.
2. **Where a covariate is a confounder or a mediator and has missing values, the default
   imputation adds +5% to +10% bias** to the exposure coefficient, and **the coverage of the
   credible interval does not reveal it**. The sign matters as much as the size: the bias is
   **away from the null**, so it *inflates* an exposure–response rather than hiding one.
3. **That bias shrinks slowly with sample size and never disappears** — still 2–4% at
   n = 12,800.
4. **It comes from the covariate imputation, not the exposure imputation** — 85% to 104% of
   it. There is no drop-in fix yet, which is why the guidance below is about reducing your
   exposure to the problem rather than correcting it.

None of this makes the pipeline unusable. It means the covariate list is a modelling
decision, not a bookkeeping one, and that a clean-looking coverage number is not evidence the
adjustment set was right.

---

## The four roles, and what each requires

Write `X` for the exposure, `Y` for the outcome, `Z` for a covariate.

| role | structure | adjust for `Z`? | why |
|---|---|---|---|
| **Confounder** ("fork") | `Z → X`, `Z → Y` | **Yes** | Omitting it leaves the exposure effect confounded |
| **Mediator** ("pipe") | `X → Z → Y` | **Depends on the estimand** | Adjusting gives the *direct* effect; omitting gives the *total* effect. Both are legitimate — say which you mean |
| **Collider** | `X → Z ← Y` | **No** | Adjusting opens a spurious path and biases the estimate, potentially reversing its sign |
| **Precision covariate** | `Z → Y` only | Optional | Improves precision, cannot confound |

`00_variable_dictionary.csv`'s `role` column has one value, `covariate`, for all four. It
controls *handling*, not *identification*: the dictionary cannot tell whether your variable
is a confounder or a collider, and neither can the pipeline. **That judgement is yours, and
it has to be made from subject knowledge or a DAG, not from the data.**

### The collider case is the one that bites hardest

In a measured case where `Z` was caused by both the exposure and the outcome, adjusting for
it moved the exposure coefficient from a true **+0.40** to **−0.10** — a sign reversal — on
*complete data* with no missingness at all. Imputation is irrelevant to it.

The practical trap is that `use_in_model = TRUE` is the natural-looking setting for any
variable you have measured. If a covariate is plausibly a *consequence* of both your exposure
and your outcome — a biomarker measured after exposure that also responds to disease status,
say — set `use_in_model = FALSE` for it.

If you want to keep it for imputation only, that is what `use_as_auxiliary` is for — with one
caveat below.

---

## What the imputation adds when a covariate is causally active

Validation before 2026-09 measured the covariate path only where `Z` was **causally inert**:
independent of the exposures, or generated from them. Under that structure the default is
well behaved (bias ≤1.3%, coverage 0.937–0.963). Under a **confounder or a mediator** with
40% missingness it is not:

| sample size | bias in the coefficient | **bias / standard error** | coverage |
|---|---|---|---|
| n = 800 | +5% to +10% | **0.30 – 0.67** | 0.93–0.95 |
| n = 3,200 | +2.7% to +6.2% | **0.36 – 0.77** | 0.91–0.94 |
| n = 12,800 | +2.0% to +3.6% | **0.50 – 0.67** | 0.87–0.92 |

**Read the middle column, not the left one.** Relative bias shrinks with sample size; the
bias *relative to your standard error* does not — it grows slightly. And the second is what
decides whether the bias can change a conclusion. As a rough guide, a study powered to detect
some effect needs about 2.8 standard errors to do so, so a bias of 0.3 SE consumes roughly a
tenth of that, and a bias of 1.0 SE about a third.

That framing is deliberate: **the effect size you care about is your choice, not a constant.**
A 5% change is a common convention, but a setting with no assumed threshold — ionising
radiation and cancer risk, say — has no such floor, and there the same absolute bias consumes
a much larger share of whatever you are trying to detect.

**And note the direction.** Every measurement of this defect is **positive** — away from the
null. It does not hide an exposure–response, it **manufactures or inflates** one. In a
precautionary or regulatory setting that is the more dangerous of the two directions: a
conservative bias can be disclosed and lived with, an anti-conservative one produces findings
that are not there.

Three things to take from this.

**Coverage does not detect it.** The intervals are honestly sized — they are not too narrow.
They cover because at these sample sizes the bias is still smaller than one standard error.
**A coverage number near 0.95 is not evidence that the point estimate is unbiased**, and in
this failure mode it never will be.

**More data helps only slowly.** The bias falls roughly as `n^−1/3` while the standard error
falls as `n^−1/2`, so the *ratio* grows: coverage degrades as the sample grows, and would
reach 0.80 somewhere above n ≈ 100,000. That is beyond most realistic cohorts, but the
direction is the wrong one.

### Where the bias comes from: the covariate imputation, not the exposure imputation

The `censored_exposure_block_fcs` strategy alternates two blocks — one imputing covariates,
one imputing the censored exposure. Validation swapped each block's draw for the exactly
correct one, in turn, on the same data:

| what was replaced | share of the bias it removes |
|---|---|
| **the covariate draw** | **85% – 104%** |
| the censored-exposure draw | 0.1% – 11.7% |
| both together | the same as the covariate draw alone |

**Essentially all of it is the covariate draw.** Replacing the covariate imputation with the
correct conditional — leaving the exposure imputation exactly as it is — brings the bias to
within ±0.5% of zero at every sample size tested. The two blocks turn out to act
independently here: their combined effect is their sum, to within 0.4 percentage points.

**What that means for you as a user, today:** the imputation of the *covariate* is where this
bias comes from, and there is no drop-in replacement yet.

Later validation narrowed it further: the bias is specifically the **flexibility** of the
default covariate imputer. A draw from a correctly specified *parametric* conditional is
unbiased — and, perhaps surprisingly, stays unbiased even when that conditional is
**misspecified**, because the analysis model already adjusts for the variable carrying the
misspecification and absorbs it. So a parametric covariate draw is a genuine candidate fix,
not a different bug.

**But it is not something you can switch on.** It needs a conditional *form* supplied for
every covariate, which is exactly the burden the flexible default removes. And the two
parametric imputers already in the pipeline were measured at +4% to +6% bias that does not
shrink with sample size — for reasons still not understood, and not shared by the property
itself. **So do not switch `z_imputer` hoping to avoid this.** Keep the default and use the
guidance below.

So the guidance below is what the evidence currently supports. It is deliberately
conservative: reduce the exposure to the problem rather than try to correct it.

### What to do about it

- **Report the exposure effect with this in mind** if a confounder or mediator in your model
  has substantial missingness. A 5% relative bias is often small against the effect size and
  the interval width — but it is systematic, not noise, and it does not average away across
  imputations.
- **Reduce the missingness in causally active covariates** in preference to reducing it
  elsewhere. Missingness in a precision covariate costs almost nothing here (measured at
  ~1 pp); missingness in a confounder costs 20× that.
- **Keep the BART Z-block default** — see below. Since the bias is the covariate draw, the
  covariate imputer is the setting that matters most here, and BART is the only one measured
  whose bias shrinks as the sample grows.
- **A complete confounder is worth more than a large sample.** The bias falls only as
  `n^−1/3`, so collecting more subjects is a weak remedy; filling in the confounder is a
  direct one.
- **Do not read a good coverage number as reassurance** on this point.

---

## The Z-block imputer matters more than it looks, under a causally active covariate

`analysis_spec$imputation$z_imputer` selects how covariates are imputed. Where `Z` is
causally inert, the choice barely matters. Where `Z` is a confounder or a mediator, it
decides whether your bias shrinks with data or not:

| `z_imputer` | bias at n = 800 | bias at n = 12,800 | coverage at n = 12,800 |
|---|---|---|---|
| **`"bart"` (default)** | +4.6 / +4.8% | **+2.1 / +2.0%** | **0.87 / 0.92** |
| `"forest_boot"` | −22.8 / −26.4% | **−22.0 / −23.3%** | **0.00 / 0.00** |
| a parametric alternative | +2.4 / +2.9% | +3.5 / +3.9% | 0.78 / 0.84 |

*(Two covariate structures; more detail in the validation findings linked below.)*

**`"bart"` is the only option measured whose bias shrinks as the sample grows.** Every other
imputer tested carries bias that does not decay — it is simply *there*, at whatever `n` you
have. Keep the default.

### ⚠️ `forest_boot` and the `dbarts` fallback

`"forest_boot"` — the pre-v1.5.0 default — carries **−23% bias that does not decay**, with
coverage falling to **zero** by n ≈ 6,400 under a confounder. Two situations reach it:

1. a configuration written against v1.4.0 that sets `proper_draw = TRUE`, honoured for
   backward compatibility;
2. **`z_imputer = "bart"` when the `dbarts` package is not installed.** The pipeline warns
   and falls back to `forest_boot` so the run completes rather than aborting.

Case 2 is worth guarding against deliberately. **Install `dbarts` before a real analysis**,
and check the run log for the `Z-block imputer:` line — it records which imputer actually
ran:

```bash
grep "Z-block imputer" run_all_stdout.log
```

If it says `forest_boot` and you did not ask for it, install `dbarts` and re-run. The one
consolation is that this failure *is* visible in the diagnostics: coverage collapses, unlike
the subtler bias above.

Reproducing a pre-v1.5.0 analysis is a legitimate reason to set `forest_boot` explicitly —
just do not expect the covariate path to be unbiased while you do.

---

## ⚠️ A non-linear covariate–exposure relationship breaks the censored-exposure draw

Everything above is about imputing the *covariate*. There is a separate and larger problem on
the *exposure* side, and it applies only to the `censored_exposure_block_fcs` strategy.

`leftcens` draws each below-LOD exposure value from a conditional that is **linear in its
predictors**. That is a deliberate design choice — the skew-aware part of it is the *margin*,
not the mean. So if one of your covariates has a genuinely **non-linear** relationship with a
censored exposure, that relationship cannot be represented, and the exposure imputation is
misspecified in a way no amount of data fixes.

Measured: with a covariate related to the censored exposure through `tanh` plus a quadratic,
the exposure coefficient carried **~20% bias, flat across n = 800 to 12,800** — asymptotic,
not a small-sample artefact, and 2–4× the covariate-draw bias described above. In standard-error
units that is **1.2 SE at n = 800 and 2.2 SE at n = 3,200** — comparable to, then larger than,
the effect a study of that size is powered to detect. And it is **away from the null**, so it
inflates the exposure–response rather than hiding it.

**How much censoring it takes.** The mechanism is that the exposure draw fits a *straight line*
to a curved relationship, so what matters is how much of the curve falls below the detection
limit. At **5% non-detects** the censored region sits in the far tail where the relationship is
nearly straight and the effect is negligible (≲2%). It grows steeply with the non-detect rate —
roughly 50× in the unrepresentable-curvature term between 5% and 70% censoring — so **the
practical threshold is the non-detect rate, not the presence of curvature alone**. If your
exposure is 5–10% non-detect, this is unlikely to matter; at 40% or more it is the dominant
error in the analysis.

**What to do:**

- **Check for it.** Plot each censored exposure against your covariates on the scale you
  model them. A visibly curved or saturating relationship is the warning sign.
- **Linearise it if you can.** A transform of the covariate that straightens the relationship
  (log, spline basis expanded into the predictor set) puts it back inside what the draw can
  represent. Passing the expanded terms explicitly via
  `analysis_spec$imputation$censored_exposure$predictors` is the mechanism.
- **If you cannot, treat the exposure–response estimate as biased** and say so. This is not
  a coverage-detectable failure.

This was found in the same run that resolved the covariate-draw question, and it is the
largest open defect in the censored-exposure path. It has not yet been decomposed — see
[`validation/phase1/FINDINGS_v22.md`](../validation/phase1/FINDINGS_v22.md).

---

## `use_as_auxiliary` does not reach the censored-exposure draw

If you mark a variable `impute_target = FALSE`, `use_in_model = FALSE`,
`use_as_auxiliary = TRUE` — "use for imputation only" — that is honoured by the **covariate**
block but **not** by the censored-exposure block. `00_censored_exposure.R` builds the
exposure draw's predictor set from `use_in_model` alone, so an auxiliary variable helps impute
covariates and is silently left out of the exposure imputation.

Measured effect: under 1 percentage point, and its sign differed between missingness
mechanisms — so this is a documentation gap more than a numerical problem. But if you are
relying on an auxiliary variable specifically to inform a censored exposure, name it
explicitly instead:

```r
analysis_spec$imputation$censored_exposure$predictors <- c(
  "outcome", "expo2", "age", "sex", "my_auxiliary_var"
)
```

That overrides the automatic set and is honoured in full.

---

## Where the numbers come from

| finding | source |
|---|---|
| Collider adjustment reverses the exposure effect; the causal role decides the covariate block | [`validation/phase1/FINDINGS_v17.md`](../validation/phase1/FINDINGS_v17.md) |
| The default's bias under a confounder or mediator, and that coverage misses it | [`validation/phase1/FINDINGS_v17.md`](../validation/phase1/FINDINGS_v17.md) |
| How that bias scales with `n` (`n^−1/3`, permanent) | [`validation/phase1/FINDINGS_v18.md`](../validation/phase1/FINDINGS_v18.md) |
| The imputer comparison, and `forest_boot`'s −23% | [`validation/phase1/FINDINGS_v19.md`](../validation/phase1/FINDINGS_v19.md) |
| Why `Y` belongs in every imputation block | [`validation/phase1/FINDINGS_v16.md`](../validation/phase1/FINDINGS_v16.md) |
| That the bias is the **covariate** draw, not the exposure draw | [`validation/phase1/FINDINGS_v20.md`](../validation/phase1/FINDINGS_v20.md) |
| That it is the imputer's **flexibility**, and a parametric draw is a real candidate fix | [`validation/phase1/FINDINGS_v21.md`](../validation/phase1/FINDINGS_v21.md), [`FINDINGS_v22.md`](../validation/phase1/FINDINGS_v22.md) |
| The ~20% asymptotic bias when a covariate relates non-linearly to a censored exposure | [`validation/phase1/FINDINGS_v22.md`](../validation/phase1/FINDINGS_v22.md) |

**What is not yet known:** what a shippable replacement for the covariate draw would be. The
bias is located (the covariate imputation) and its cause identified (the imputer's
flexibility), and a parametric draw is now known to work on both sides of the trade that
stood against it. What is missing is a way to get a conditional form for every covariate
without asking the analyst for it — and an explanation of why the two parametric imputers
already present carry +4% to +6% bias when the *property* does not. Until both are settled the
guidance above is the whole of the remedy.
