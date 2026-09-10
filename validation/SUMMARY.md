# What we learned, and what the pipeline now does about it

**2026-09-09.** A synthesis of 24 validation tracks, the theory they produced, and the changes
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

### 5.3 The child factor does **not** ship, and the reason is measured

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

- **Both blocks missing at once.** Every V24 cell has fully observed covariates, to isolate the
  exposure draw from the covariate-draw ladder V21 settled. The realistic case needs the child
  factor inside the block alternation.
- **Where the declared shape comes from.** We can say a wrong shape is harmful and the right one
  is exact. We cannot tell a user how to choose, beyond §3's finding that the data cannot.
- **No arrow-absent control in V24.** The claim that "unrepresentable curvature" modifies rather
  than generates the defect currently rests on a comparison across two instruments. One extra
  scenario would close it.
- **Every result is on a Gaussian linear estimand.** All 24 tracks use OLS for the estimand, and
  the pipeline's real targets are logistic, ordinal `mo()`, spline and mixed models. Step 6's
  pooling machinery — the transform and its bimodality gate — **has never fired**, because no
  validation cell produces a non-Gaussian pooled posterior. This is arguably the largest
  remaining gap in the programme.
