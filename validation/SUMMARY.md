# What we learned, and what the pipeline now does about it

**2026-10-13.** A synthesis of 32 validation tracks, the theory they produced, and the changes
that theory justified in the pipeline. Written to be readable without the other documents;
every number is a measurement with a named source.

Companions: [`THEORY.md`](THEORY.md) (the mechanism and its algebra),
[`README.md`](README.md) (the status index), [`ROADMAP.md`](ROADMAP.md) (what is still open),
[`phase1/FINDINGS_v24.md`](phase1/FINDINGS_v24.md) (the decisive run).

---

## 1. The problem, in one paragraph

An environmental exposure is often reported only as "below the limit of detection" for a large
share of samples — 30–60% is routine. Those rows are not missing at random in the ordinary
sense; they are *censored*, and we know something exact about them (the value lies below a known
bound). The pipeline imputes them, and the question this whole programme answers is: **when is
that imputation good enough to leave the exposure–response coefficient alone, and when does it
quietly bend it?**

The answer turned out to depend on something the imputation code cannot see: **the causal role
of the covariates**. Two analyses with identical data, identical missingness and identical
models can differ by 22 percentage points in their exposure coefficient depending on whether a
covariate is a *cause* of the exposure or an *effect* of it.

---

## 2. The theory, in the order it was derived

### 2.1 The starting point is Bayes, not a heuristic

For a censored value `v`, the distribution to impute from is

> `p(v | Y, rest) ∝ p(Y | v, rest; θ) · p(v | rest)`

— the outcome's likelihood times the exposure's own conditional. This is an identity, not a
modelling choice. Its consequence is **congeniality**: an imputation that ignores `Y` is drawing
from the wrong distribution, and the exposure–response coefficient attenuates. Measured
directly, with `Y` removed from one imputation block at a time: **−13.1 percentage points**
(Track V16).

### 2.2 The bias has exactly two parts

Writing `Δ` for the imputation's error in the censored cells,

> `bias = T₁ + T₂`,  `T₁ = −β₀[(W′W)⁻¹W′Δ]₁`,  `T₂ = [(W′W)⁻¹W′ε]₁`

verified to ~1e-17 (`THEORY.md` §0c). `T₁` is **attenuation** — the imputed exposure is a
shrunken version of the truth. `T₂` is **outcome borrowing** — the imputation has copied some of
the outcome's own noise into the exposure. They have opposite signs, and *congeniality is the
`Y`-loading at which they cancel*. This is why "use `Y` in the imputation" is not a slogan: too
little `Y` leaves `T₁`, too much creates `T₂`.

### 2.3 Bayes plus the DAG gives the actual target

The identity above hides which variables belong in `p(v | rest)`. Factorising along the causal
graph makes it explicit:

> `p(x | rest, x ≤ L)` ∝ `p(x | Pa(X))` · `p(Y | x, Pa(Y))` · `∏_{C ∈ Ch(X)} p(C | x, ·)` · `1{x ≤ L}`
>
> parents · the outcome · **every covariate the exposure causes** · the censoring bound

The third factor is the one nobody writes down. It is empty when every covariate is a *cause* of
the exposure — which is why the confounder case has always worked. It is not empty when a
covariate is an *effect* of the exposure.

### 2.4 Direction is DAG-aware, and that is the practically important part

A bias toward the null is conservative; one away from the null invents an association. The
measured direction is not a single fact about the method — it depends on the graph:

| what is wrong | direction | why |
|---|---|---|
| `Y` omitted from the imputation | **toward null** | pure `T₁`; `T₂ = 0` because nothing borrowed `Y` |
| an informative covariate omitted | **toward null** | same reason; 6/6 cells, three sample sizes (V24) |
| a covariate's curved relationship forced through a straight line | **away from null** | `T₂` dominates; +22% measured |
| a collider adjusted for | **sign reversal** | +0.40 became −0.10 on *complete* data (V17) |

---

## 3. What was measured, condensed

Only the facts that change a decision. Sources in brackets.

**The exposure draw is exact more often than we assumed.** At 40% non-detects, `n` = 3200:
confounder **+0.20%**, linear mediator **−0.02%** (direct) and **+0.19%** (total), collider
**+0.16%**. The single regression `leftcens` performs *is* the exact conditional wherever the
joint is Gaussian. [V24 §2b]

**It is badly biased in exactly one configuration: a mediator whose relationship with the
exposure is curved.** **+22.0%** (direct) and **+17.9%** (total), `bias/SE` 2.4 and 6.1. This
does **not** shrink with sample size. [V22, V24 §2b]

**Dropping an informative covariate from the exposure draw costs 6–7%, permanently.** Flat in
`n`, so coverage *falls as the study grows*: 0.945 → 0.830 → **0.535** for a direct effect and
0.830 → 0.480 → **0.040** for a total effect across `n` = 800 → 12,800. A large study is more
exposed, not less. [V24 §2]

**The covariate block, not the exposure block, carried the confounder-case bias** — 85–104% of
it — and within that, BART's *flexibility* was the whole story (100%); properness moved the
interval, not the point. [V20, V21]

**A parametric covariate draw is not a trade.** It wins on both sides, including where it is
misspecified, because the analysis model absorbs the error. [V22]

**So the default became a parametric draw.** The objection that had blocked it for months — a
measured −4.70% where the covariate's conditional is non-linear [V6] — was re-run on exactly
that ground and **did not reproduce**: mice came within 0.09 SE of bart and on the right side of
it, at every `n`. Relative bias fell 10.21% → **2.95%** in the worst cell, and mice was better or
equal in every cell at every `n`. `z_imputer` now defaults to `"mice"` (v1.7.0). [V27]

**A derived-vs-measured ratio taken at one `n` can be measuring the sample size.** The chain
below is an *asymptotic* statement; the measured effect climbs with `n` toward it. The
"fork residual" — recorded for a month as the chain closing on one causal role (1.011) and
missing on the other (1.147) — is an artefact of both having been measured at `n = 800`. At
`n = 12,800` the fork reads **0.992** and the pipe **0.953**: the fork agrees *better*. ⚠️ Pre-V28, the **end-to-end** comparisons — a closed form against a *pipeline-measured* effect — were all taken at `n = 800`. The component-level links are unaffected, being verified at `N = 300000` on both sides. [V28]

**The outcome-family attenuation is exactly the outcome's Fisher information — but that is only
half of a transfer.** The cumulative-logit generalisation predicted `SE_gauss/SE_ord` in 24 of 24
cell×`n` comparisons to ≤0.001, including a 3-category outcome split 85/10/05 that carries *less*
information than a balanced binary one. ⚠️ A `bias/SE` ratio is `(bias ratio)×(SE ratio)` though,
and V26's assumption that a defect's *absolute* bias is the same under both families — never
tested until now — is false: it runs **0.26 to 2.2**. **So ~0.474× is not a transferable
constant**, and whether a Gaussian-measured defect is conservative elsewhere is not yet knowable.
[V30]

**But there is no transfer law, and the apparent localisation was an identification artefact.**
With five arms spanning surviving fraction 0.069–1.000, a single pair of channel scales per
family is **rejected** for binomial and ordinal outcomes (χ²/df 17.8 and 4.4) and the borrowing
scale is arm-dependent (p = 0.004, 0.008). Under a Gaussian outcome the instrument is essentially
exact (χ²/df = 0.24, both scales within 5% of 1). ⚠️ V31's claim that attenuation transfers
untouched while borrowing fails is **withdrawn**: its arms were too ill-conditioned to separate
the two scales, which are collinear when every arm cancels heavily. [V32]

**What V31 established, before that correction.** Recording the imputation error `d`
directly shows the bias is `borrow + atten`; under a Gaussian outcome the prediction is
essentially exact (χ²/df = 1.1 across four arms and three missingness fractions), the
**attenuation** channel transfers across outcome families untouched, and only the **borrowing**
channel fails — over-predicted 1.6× under ordinal and 2.4× under binomial. It is not a question
of expansion order: that was registered as a prediction and refuted, each arm's error being flat
across a 2.6× range of `|d|`. [V31]

**And the residual that did survive is governed by `s²`, not by the causal role.** Moving
`Var(logX1)` with the arrow held fixed, the chain's `derived/measured` ratio rises in *both*
roles, and one line fits all six cells — the role adds nothing (offset p = 0.79, slope p = 0.21),
crossing 1 at `s²` ≈ 1.42. The fork/pipe contrast was a parameterisation artefact: only the
fork's arrow moves `s² = sd_x² + delta_zx²`. The chain is accurate near `s²` ≈ 1.4 and degrades
in both directions, ~0.16 in ratio per unit `s²`. [V29]

**Bias/SE, not relative bias, is the right currency.** Bias decays like `n^−1/3` while the SE
decays like `n^−1/2`, so bias/SE *grows* like `n^(1/6)`: an arm that just passes at `n` = 800
fails at `n` = 12,800. Every gate is stated in bias/SE. [V18, V23]

**The quantity that governs the defect is not estimable from data.** It lives below the LOD. A
natural spline returns exactly zero there by construction, a quadratic is 4.6× low, a cubic
false-positives on a null. So this can never be a diagnostic computed from a dataset — only a
sensitivity analysis. [`THEORY.md` §0b]

---

## 4. The algorithm

Add the missing factor, and draw from the product on a grid (or, equivalently, by reweighting
draws that already cover the first two factors):

```
per outer sweep:
  draw θ from the analysis model's posterior
  fit p(x | Pa(X)) interval-censored, using ALL rows      # not an lm on observed rows
  for each censored cell i:
    grid x over (L_i − 8 s_i, L_i]                        # anchored on the combined Gaussian
    log w = log p(x | Pa_i)
          + log p(Y_i | x, Pa(Y)_i; θ)
          + Σ_{C ∈ Ch(X)} log p(C_i | x, ·)               # THE MISSING TERM
    normalise, inverse-CDF sample
```

Validated across all four cases that need it, simultaneously, at three sample sizes ×
200 replicates: **`|bias/SE|` ≤ 0.13 everywhere**, coverage nominal, indistinguishable from
having no missing data — against 1.2 to 8.2 for the draw without the term. [V24 §1]

Three details that are easy to get wrong and each produce believable wrong numbers:

- The outcome factor conditions on `Pa(Y)`, **not** "all covariates". Under a collider the
  covariate is a *child* of `Y`, so including it means conditioning on a collider inside the
  imputation.
- The parent factor must be fitted **interval-censored** on all rows. An `lm` on the above-LOD
  rows is a response-truncated regression — measured σ 24% low.
- The grid must be anchored per cell on the parent × outcome Gaussian. A fixed `L − 6σ` window
  clipped the target for extreme-`Y` cells.

---

## 5. What the pipeline now does

### 5.1 Auxiliary and out-of-model covariates now reach the exposure draw *(fix, no assumption)*

The automatic predictor set was `outcome + use_in_model covariates + other exposures`. A
mediator in a **total-effect** analysis is marked `use_in_model = FALSE` — correctly — and so
**dropped out of the exposure imputation too**. Same for `use_as_auxiliary = TRUE`, which was
documented as a sub-1-percentage-point curiosity. It is not: that is the 6–7%-flat-in-`n` cell
above, with coverage reaching 0.040.

Now `use_as_auxiliary = TRUE` and `role = "auxiliary"` both enter the set. **This costs nothing
and assumes nothing** — conditioning *linearly* on the covariate is sufficient whenever the
relationship is linear. `00_censored_exposure.R`.

### 5.2 A preflight that says whether you are in the unsafe cell *(diagnostic, no assumption)*

`01_validate_config.R` now reports, per (exposure, covariate) pair, the strength of their
association and how much of it a straight line misses, alongside the non-detect rate. Curvature
*above* the LOD is estimable, so unlike the defect itself this needs no assumption. It stays
silent on a linear DGP and fires on a curved one.

It is a screen, not a guarantee: it cannot see below the LOD, and it cannot tell a confounder
from a mediator — only the analyst can.

### 5.3 The covariate block conditioned on the wrong scale *(fix, v1.6.0 — behaviour change)*

The block-FCS loop returned the censored exposure to the **data** scale after each exposure
draw, so every subsequent covariate block regressed on `X` while the analysis model — and the
exposure draw itself — worked in `log X`. Congeniality clause (A) is a statement about the
*functional* of a variable, not only about which variables, so this was a genuine violation:
a model linear in `X` cannot represent a conditional linear in `log X`.

**Who was affected: only LINEAR covariate imputers.** Regression trees split on order
statistics, so a tree ensemble is invariant to any strictly monotone transform of a predictor
— `exp()` is monotone, and BART never saw the defect. Measured pre/post fix:

| covariate imputer | change |
|---|---|
| `bart` (the default), `forest` | +0.02 / −0.06 pp — none |
| `mice` pmm, parametric proper | **−3.16 to −5.36 pp** |

So the shipped default and V17/V20/V23 are unaffected; what needed restating was V19's two
parametric arms. The magnitude is derived in `THEORY.md` §1a(ii) — retained fraction
`σ̃²/(e^{s²} − 1 − τ²)` verified to ≤0.002, and the full chain predicts one cell to 1.1%.

### 5.4 The child factor does **not** ship, and the reason is measured

The validated algorithm is a **grid sampler**: build the product of the three factors on a grid
over the censored interval and inverse-CDF sample it. That lands at −0.09% to −0.63% in every
case, at three sample sizes.

What was implemented instead — and then removed — asked `leftcens` for K candidate draws and
reweighted them by the child factor. Attractive because it leaves `leftcens`'s skew-aware margin
untouched, and reimplementing that margin caused four separate defects earlier in this programme.
Measured against the grid sampler on identical data:

| | linear pipe | `pipe_nl` (the case the feature exists for) |
|---|---|---|
| weight `p(C \| x)` | +7.61% | +27.37% |
| weight `p(C \| x, Y)` ← the correct weight | **+0.46%** | +32.53% |
| grid sampler | −0.60% | −0.41% |

Two findings, in order.

**The weight is `p(C | x, Y)`, not `p(C | x)`.** The proposal already conditions on `Y`, and
since `p(Y|x,C)·p(C|x) = p(Y|x)·p(C|Y,x)`, the residual factor carries `Y`. That correction is
real — it fixes the linear pipe outright.

**But it does not help where the feature is needed.** With the child removed from the predictor
set (required, or its information is counted twice — measured at +24.48%), `leftcens` regresses
`x` on `Y` *linearly*, while marginalising the child out makes `Y = b₁x + γ₁·g(x) + …`, genuinely
non-linear in `x`. So the **proposal itself** is misspecified, and importance reweighting can only
*add* a missing factor, never repair a wrong proposal — that needs the full target/proposal ratio,
and `leftcens` does not expose its density. The clean contrast is the table: same code, same
weight, linear proposal works and non-linear proposal does not.

**What a correct implementation needs.** The grid sampler, whose outcome factor conditions on the
child — where `Y` *is* linear in `x`. That means evaluating the **analysis model's** linear
predictor as a function of a candidate exposure value, per row. Feasible, but not a small
addition: the analysis model here is an arbitrary `brms` formula, and the margin would have to be
preserved by working on the shash-`z` scale through `leftcens`'s exported `fit_shash_margin()` /
`x_to_z()` / `z_to_x()`. Tracked in `ROADMAP.md`.

**Shipping the reweighting version behind a warning was the wrong instinct**, and worth recording
as such: a feature that is off by default but wrong when enabled is worse than no feature. That
is this track's own finding — a wrong child model is worse than none — applied to the
implementation rather than to the user.

---

## 6. The honest limit of the fix

The child factor must be evaluated **below** the detection limit, where there is no exposure data
to fit it from. V24 tested only fitted *polynomials* against the true arrow, and every polynomial
was **worse than omitting the term**:

| assumed shape | bias | coverage |
|---|---|---|
| none (omit the factor) | +8.8% | 0.375 |
| linear | +25.8% | 0.000 |
| quadratic | +31.1% | 0.000 |
| cubic | +37.4% | 0.000 |
| the true arrow, coefficients known | **−0.2%** | **0.945** |

That left a question a shipped feature lives or dies on: the true arrow's *coefficients* can
never be known, so is the feature useful at all? **Tested separately (n = 3200, 200 reps), and
the answer is yes — what has to be right is the *basis*, not the numbers:**

| declared shape | bias | `bias/SE` |
|---|---|---|
| correct basis, coefficients **fitted from the observed rows** | **−0.51%** | −0.07 |
| correct form, wrong internal constant (`tanh(x)` for `tanh(1.8x)`) | **+3.22%** | +0.44 |
| none (omit the factor) | +8.47% | +1.09 |
| correct family, one term dropped | +20.71% | +2.33 |
| a straight line | +25.43% | +2.67 |

Three consequences, and they set the interface:

1. **Declaring the shape is enough.** Coefficients estimated above the LOD extrapolate correctly
   when the basis is right — as good as knowing the answer.
2. **A roughly-right shape degrades gracefully.** Getting an internal constant wrong cost 3.7
   percentage points, not 25. So this is usable by someone who knows the mechanism only
   approximately.
3. **An *incomplete* basis is worse than nothing.** Dropping a term from the right family gave
   +20.7% against +8.5% for omitting the factor. The failure mode is missing structure, not
   insufficient flexibility — which is why every fitted polynomial fails, and why the declared
   form must never be chosen by goodness of fit. A good in-range fit is not evidence about an
   extrapolation: all five rows above fit the observed region acceptably and differ by 26
   percentage points below the LOD.

So when the grid sampler is built, its interface takes a **declared formula**
(`~ tanh(1.8 * .x) + I(.x^2)`), not a flexibility setting, and the defensible use is still to run
more than one form and report the spread.

**And a zero-assumption mitigation exists that we would have missed.** Under a curved mediator,
conditioning on the covariate *linearly* is worse than **dropping it**: +22.0% against +9.07%.
Removing such a covariate from `censored_exposure$predictors` roughly halves the bias at no
assumption cost. Not a fix — 9% is still 9% — but strictly better than the default, and free.

---

## 7. What is still open

Grouped by what kind of work would close it.

### Theory not yet closed

| | state |
|---|---|
| The `u²` law's exponent and constant | still a **fit**, not derived (`THEORY.md` §4). V23's exponent CI was [0.99, 2.28] — no power. V25 narrowed the role: `u` **modifies** the defect, it does not generate it |
| The link from `‖r‖` to `T₂` | why a misspecified conditional *over-borrows* by an amount `u²` happens to track. This is the step that would fix the exponent |
| The scale defect's residual, **after V29** | **resolved as to shape, open as to mechanism.** The error is a monotone function of `s²` crossing zero near `s²` ≈ 1.4 (slope +0.161 ± 0.033 in ratio per unit `s²`); the causal role is irrelevant (p = 0.79 / 0.21), so there is no "fork residual". But it is **not** the `exp`-linearisation: that must vanish as `s² → 0`, and the ratio instead falls to 0.709 at `s²` = 0.489. At least two `s²`-dependent effects of opposite sign, or one `s²`-dependent mis-scaling — the truncation factor `Var(L\|L≤c)/σ_x²` is the suspect, since the LOD is a quantile of `L`. Separately, the chain under-predicts the pipe by a flat ~5% at `n` = 12,800 |
| ~~Why logistic defects attenuate by ~0.474~~ | ✅ **derived 2026-10-11** — it is Fisher information, not imputation: `SE_g/SE_l = σ·√E[p(1−p)]`, predicting each cell to ≤0.016 with no free parameters. **Generalises to any GLM** as `σ·√E[w]`, so V26's conservative-transfer result need not be re-measured per family. `THEORY.md` §6c |

### Build work

| | state |
|---|---|
| The DAG-factorised exposure draw in the pipeline | algorithm **validated** (`\|bias/SE\|` ≤ 0.13 across four cases, three `n`); a cheap SIR implementation was built, measured and **removed** — it cannot work in the curved case it exists for. A faithful version needs the analysis model's linear predictor as a function of a candidate exposure value |
| Both blocks missing at once | depends on the above; `dag_impute_datasets()` refuses it rather than returning `NA` |
| ~~A shippable covariate-block fix~~ | ✅ **delivered 2026-10-13 (v1.7.0)** — it is `z_imputer = "mice"`, now the default. V22 removed the objection in principle; V27 tested it on the ground that had blocked it (a non-linear covariate conditional, where V6 measured a parametric draw at −4.70%) and it did not reproduce. Better or equal in every cell at every `n`; worst-cell relative bias 10.21% → **2.95%** |
| Item 07 | blocked — gates **A07-1/2/3** all currently fail; its grid draw inverts to +142% under a non-linear covariate arrow |

### Validation coverage

| | state |
|---|---|
| ~~Posterior predictions (Step 7) pool as a plain mixture~~ | ✅ **settled 2026-10-11 — the correction applies.** The mixture is right for the predictive *shape*, but the interval is centred on an estimate of the infinite-`m` centre, so Rubin's `B/m` term still belongs. Measured coverage: **0.9446** uncorrected at `m` = 3 against nominal 0.95, 0.9508 by `m` = 20. Effect is ~0.5 pp at `m` = 3 and nil by `m` = 10, so it ships for consistency rather than impact |
| ~~Step 10's monotonic OR tables bypassed the MI machinery~~ | ✅ **closed 2026-10-10.** Read the same draws as Step 6 but summarised them unweighted and uncorrected, so the OR table was narrower than the parameter table beside it (0.8–8.5% depending on FMI and `m`; <1% at `m` ≥ 50). Both now share `mi_pool_draws()`. The `mo()` arithmetic itself was correct and is now pinned by `test/test_step10_mo.R` |
| ~~Step 6's pooling machinery~~ | ✅ **closed 2026-10-10 — and it was wrong in two independent ways.** The finite-`m` correction was **nullified** by a denominator mismatch (missing Rubin's total by up to 22% at `m` = 3, an interval ~12% too narrow), and the **logit branch was unreachable** so probability-scale parameters were corrected on the log scale. Both fixed and now exact; first test at `test/test_step6_pooling.R`. **Behaviour change — pooled intervals get wider** |
| Ordinal `mo()`, spline, mixed families | untested. V26 covered binomial only; the mechanism suggests the same direction, which is a conjecture |
| The shipped default under a confounder or mediator | documented and attributed (+5% to +10.4%), **no fix** |
| The MAR conjecture | open, small — does a flexible covariate draw recover `Y` from proxies? |
