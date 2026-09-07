# Which covariates to adjust for, and what the imputation does with them

**Short version.** The pipeline treats every row marked `use_in_model = TRUE` the same way,
but your covariates do not all play the same causal role — and the role decides both whether
adjusting for the variable is correct at all, and how much bias the imputation adds.
Validation measured three things worth knowing before you run:

1. **Adjusting for a collider inverts the exposure effect** — from +0.40 to −0.10 in a
   measured case, on complete data, before imputation enters. No imputation method repairs
   this; it is a model-specification error and the pipeline's `use_in_model = TRUE`
   convention invites it.
2. **Where a covariate is a confounder or a mediator and has missing values, the default
   imputation adds +5% to +10% bias** to the exposure coefficient, and **the coverage of the
   credible interval does not reveal it**.
3. **That bias shrinks slowly with sample size and never disappears** — still 2–4% at
   n = 12,800.

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

| sample size | bias in the exposure coefficient | coverage |
|---|---|---|
| n = 800 | **+5% to +10%** | 0.93–0.95 |
| n = 3,200 | +2.7% to +6.2% | 0.91–0.94 |
| n = 12,800 | **+2.0% to +3.6%** | 0.87–0.92 |

Two things to take from this.

**Coverage does not detect it.** The intervals are honestly sized — they are not too narrow.
They cover because at these sample sizes the bias is still smaller than one standard error.
**A coverage number near 0.95 is not evidence that the point estimate is unbiased**, and in
this failure mode it never will be.

**More data helps only slowly.** The bias falls roughly as `n^−1/3` while the standard error
falls as `n^−1/2`, so the *ratio* grows: coverage degrades as the sample grows, and would
reach 0.80 somewhere above n ≈ 100,000. That is beyond most realistic cohorts, but the
direction is the wrong one.

### What to do about it

- **Report the exposure effect with this in mind** if a confounder or mediator in your model
  has substantial missingness. A 5% relative bias is often small against the effect size and
  the interval width — but it is systematic, not noise, and it does not average away across
  imputations.
- **Reduce the missingness in causally active covariates** in preference to reducing it
  elsewhere. Missingness in a precision covariate costs almost nothing here (measured at
  ~1 pp); missingness in a confounder costs 20× that.
- **Keep the BART Z-block default** — see below.
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

**What is not yet known:** which part of the imputation produces the confounder/mediator
bias — the exposure draw, the covariate draw, or their alternation. Until that is settled
there is no fix to offer, only the guidance above.
