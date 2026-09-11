# Why the pipeline fails where it fails

**One condition explains thirteen validation tracks.** This document states it, proves the
boundary exactly, and maps each empirical finding onto a specific clause of it.

It is deliberate about the limits: the boundary — *when* linear-Gaussian imputation is
correct — is derived and numerically exact. The *magnitudes* when it is violated are not
derived, and two attempts to derive them are recorded below as failures so nobody repeats
them. Every number quoted as theory has been checked against the Monte-Carlo results, and
where the two disagree the disagreement is stated rather than smoothed.

Written 2026-09-01, after V13. Referenced from the findings files, which each explain
themselves locally but none explain the set.

---

## 0. The whole account — and what in it is an identity, a claim, or a guess

**Added 2026-09-09, corrected the same day.** Twenty-three tracks and the framework behind
them had been implicit. Writing it down exposed that its three parts have *very* different
epistemic status, and an earlier draft of this section blurred them. They are separated here
because the difference decides what is worth testing.

### The identity — not testable, and not tested

For a missing value `v`,

> `p(v | Y, rest) ∝ p(Y | v, rest; θ) · p(v | rest)`,  so
> `ℓ(v) ≡ log p(v | Y, rest) = log p(Y | v, rest; θ) + log p(v | rest) + const`

**This is Bayes' theorem.** It is true by construction, it has never been "validated", and
validating it would be a category error. Every use of it in this document is a use of an
identity, not of a finding. What Meng's congeniality argument adds is the *observation* that an
imputation drawing from anything else is uncongenial — also not an empirical claim.

### The empirical claim — what family the pipeline actually draws from

This *is* testable, and it is where an earlier draft of this section was wrong.

| block | what it draws from | representable class |
|---|---|---|
| **covariate** (BART, or mean + homoscedastic Gaussian) | a Gaussian on the **raw** scale | log-density quadratic in `v` |
| **exposure** (`leftcens::impute_censored_conditional`) | shash margin → **transform to a latent normal scale** (`x_to_z`) → interval-censored Gaussian AFT (`survreg`) → truncated normal draw → `z_to_x` | log-density quadratic in **`z(v)`**, the shash transform — **not** in `v` |

**The earlier draft stated one projection, onto `span{1, v, v²}`, for both blocks. That is
wrong for the exposure block**, whose representable class is quadratic in a *monotone
transform* of `v`. The transform is not incidental: Phase 1 §7.5 measured coverage falling to
0.47 when it is replaced by a plain Gaussian on skewed exposures, so it is doing real work and
enlarges the representable class in a way the raw-scale statement misses.

Verified by reading the installed function, not from memory: `fit_shash_margin` → `x_to_z` →
`survreg` → `rnorm_trunc` → `z_to_x`.

### The conjecture — why the bias should be second order

Let `q*` be the member of the block's representable family that the fit converges to. **It is
the KL projection**, because these are maximum-likelihood fits under misspecification — *not*
the L² projection of the log-density, which is what the earlier draft asserted. The
first-order condition at `q*` is that the true distribution is orthogonal to the fitted
family's **score**,

> `E_true[ ∂/∂β log q_β ]|_{β*} = 0`

i.e. moment matching, not log-density projection. The bias in `θ̂` is a smooth functional of
the imputation distribution, and its derivative along directions *inside* the family vanishes
at `q*` — so the leading term is **quadratic in the distance from the family**.

**That argument may well be right, but three things about it are unsettled and one of them is
a substitution I have not justified:**

1. **The norm is unspecified.** "Distance from the family" is a distributional distance (KL,
   Hellinger, …). Nothing here says which, and the exponent could depend on the choice.
2. **`u` is a proxy, not the object.** V23's `u` is the *L² residual of the covariate arrow
   `g` after a linear fit below the LOD* — neither a KL distance nor an L² log-density
   residual. It predicted four of six cell levels to within 1.5 pp from one calibration point,
   which is striking and **unexplained**. Why a heuristic proxy tracks the theory's quantity so
   well is an open question, not a result.
3. **The exponent is untested and the constant is unfixed.** V23 measured the exponent at
   **1.638, CI [0.994, 2.282]** — containing 1 as readily as 2. Refitting the constant with the
   exponent held at 2 gave **249**, against the registered 300, and the two largest-`u` cells
   were over-predicted, so the one-parameter form saturates.

**Status: a motivated form with an unspecified norm, a proxy standing in for its argument, an
exponent consistent with anything between linear and quadratic, and a constant that moved 17%
on refitting.** That is weaker than "derived", and the earlier draft of this section read as
though it were stronger.

### The proposition that *is* well tested — propagation

A non-zero misspecification for `v` moves estimand `θ` only if `v` is conditioned on by the
analysis (§1b), **and** lies on an open exposure–outcome path (§1b), **and** the
misspecification is not absorbable by a coefficient the analysis already fits (V22). All three
required. V17 measured a factor of ~22 between the first two holding and not; V22 found a case
where the third fails and a misspecified draw is nonetheless unbiased for the focal
coefficient.

**This is the part of the framework with real evidence behind it**, and it is structural rather
than numerical — which is why it survived the currency change in §4b intact while every
magnitude had to be re-read.

### What the account buys, and what it does not

It buys an enumeration: every defect measured anywhere in this project is a way for the fitted
family to miss the target, and the ways are structured — a non-quadratic outcome term, a
non-representable `p(v | rest)`, or a scale on which neither is quadratic. It says what a fix
must do: enlarge the family, shrink the distance, or break propagation.

It does not buy magnitudes. Every number in this document is conditional on a chosen distance
from the family, and §0b shows **that distance is not identifiable from data**. So the project
has measured a response function it cannot locate a real study on — and the response function's
own exponent and constant are, after V23, still unfixed.

---

## 0c. What a wrong draw actually does: an exact two-term decomposition

**Derived and verified 2026-09-09.** §0 gave the identity and a conjecture about magnitudes.
This section is neither: it is an **exact algebraic consequence** of the identity for the
estimator this pipeline actually computes, verified numerically to machine precision. It is
the first equation in this document that is derived, checked, and not fitted.

### The decomposition

Let the analysis model be correct,

> `Y = β₀X + γ₀′Z + ε`

with `X` left-censored at `L`. The imputation supplies `X̃`, equal to `X` on observed rows, and
write

> `Δ = X̃ − X`   (zero on observed rows, non-zero only where `X ≤ L`)

Substituting `X = X̃ − Δ` into the normal equations for `W = (X̃, Z, 1)` gives, with **no
asymptotics and no approximation**,

> `b̂ = (β₀, γ₀′)′ − β₀ (W′W)⁻¹W′Δ + (W′W)⁻¹W′ε`

so the bias in the focal coefficient is **exactly the sum of two terms**:

> **T₁ = −β₀ · [(W′W)⁻¹W′Δ]₁**  — the **imputation-error** term
> **T₂ = +[(W′W)⁻¹W′ε]₁**  — the **outcome-borrowing** term

`T₂` exists because a `Y`-aware imputation makes `X̃` a function of `Y`, and `Y` contains `ε`.
It is zero if and only if the imputation ignores `Y`.

### What each term does — measured, not asserted

`b₀ = 0.40`, 40% left-censored, shipped `leftcens` draw, 120 reps at `n` = 2,000. The residual
`bias − T₁ − T₂` is ~10⁻¹⁷ in every row, as it must be:

| cell | bias | `T₁` (attenuates) | `T₂` (borrows) |
|---|---|---|---|
| precision | +0.0010 | −0.0578 | +0.0589 |
| **precision, `Y` omitted** | **−0.0602** | −0.0610 | **+0.0008** |
| fork | −0.0017 | −0.0631 | +0.0613 |
| pipe | +0.0049 | −0.0710 | +0.0759 |
| **`pipe_nl`** | **+0.0933** | −0.0088 | **+0.1021** |
| `pipe_nl`, `Y` omitted | −0.0092 | −0.0151 | +0.0059 |

**`T₁` is negative in every cell** — attenuation toward the null, behaving exactly like
classical measurement error on the censored rows.

**`T₂` is positive and vanishes when `Y` is dropped** (+0.0008, +0.0059 against +0.06 to +0.10
when `Y` is used).

### `T₂` is proportional to the imputation's `Y`-loading

Scaling the fitted `Y` coefficient of a *correctly specified* draw by `k` (so `k` = 1 is the
congenial draw), `n` = 3,000, 100 reps:

| `k` | bias | `T₁` | `T₂` | `T₂/k` |
|---|---|---|---|---|
| 0.0 | −0.0430 | −0.0419 | −0.0011 | — |
| 0.5 | −0.0208 | −0.0471 | +0.0263 | 0.053 |
| **1.0** | **+0.0005** | −0.0550 | +0.0555 | 0.056 |
| 1.5 | +0.0196 | −0.0659 | +0.0855 | 0.057 |
| 2.0 | +0.0347 | −0.0796 | +0.1143 | 0.057 |

`T₂/k` is constant to two figures. So:

> **Congeniality is the `Y`-loading at which `T₁ + T₂ = 0`.**

The net bias crosses zero at `k` = 1 (+0.0005). Congeniality is not "the draw happens to be
unbiased" — it is the exact point at which outcome-borrowing offsets censoring attenuation.
That is a sharper statement than §1's clause list, and it is derived from the identity rather
than measured against it.

### The direction rule, and the correction it forces

| what the imputation does | `T₂` | net direction | measured |
|---|---|---|---|
| ignores `Y` | 0 | **toward the null** — pure attenuation | −0.043 (−10.8%) at `k` = 0 |
| under-borrows (`k` < 1) | small | toward the null | −0.021 at `k` = 0.5 |
| **congenial** (`k` = 1) | offsets `T₁` | **unbiased** | +0.0005 |
| over-borrows (`k` > 1) | dominates | **away from the null** | +0.035 at `k` = 2 |
| wrong conditional *while using* `Y` | dominates | **away from the null** | +0.093 in `pipe_nl` |

**So "imputing left-censored `X` incorrectly gives bias toward the null" is true only of the
`T₂` = 0 case** — the classical result, and the one Phase 1's H1 and V16 measured (−13.1 pp).
Get the conditional wrong *while still conditioning on `Y`* and the sign reverses, because the
draw over-borrows from the outcome. Every away-from-null entry in the §11b direction ledger is
this: a `T₂` that `T₁` no longer cancels.

In `pipe_nl` the mechanism is visible in the split — `T₁` collapses to −0.0088 while `T₂` rises
to +0.1021. The non-linear covariate arrow does not make the draw noisier; it makes it lean on
`Y` for information the arrow should have supplied.

### A falsifiable consequence, already observed

Item 07's grid draw proposes from the exposure prior and **reweights by `p(Y | x, rest)`** —
which is maximal outcome borrowing with no covariate-arrow information at all. The
decomposition therefore predicts a runaway `T₂` and a large away-from-null bias in exactly the
cells where the arrow matters. V23 measured **+43.8%, +94.1%, +142.1%** as `u` rises, with
coverage 0.000 (`phase1/FINDINGS_v23.md`). That was observed before this decomposition was
written; it is now predicted by it.

### What this changes about the fix

The target is **not** "make `Δ` small". `T₁` and `T₂` are both large in every well-behaved
cell (≈ ±0.06) and cancel; a fix that shrank `Δ` alone would break the cancellation. The target
is `T₁ + T₂ = 0`, i.e. **the correct `Y`-loading**, which is what congeniality delivers and
what a misspecified conditional gets wrong in either direction.

**This also bounds what any diagnostic could do.** `T₂` depends on `ε`, which is unobservable,
and `T₁` on `Δ`, which requires the unobserved `X`. Neither term is estimable from data — a
sharper version of §0b's non-identifiability, now with a reason rather than an experiment.

---

## 1. The condition

For a variable `v` with missing values, the imputation must draw from the distribution
implied jointly by the analysis model and `v`'s own model:

> **p(v | Y, rest) ∝ p(Y | v, rest; θ) · p(v | rest)**

This is Meng's congeniality requirement, specialised. The first factor is the **substantive
model's likelihood**; the second is the **imputation model for `v`**. An imputation that
draws from anything else is uncongenial, and the resulting estimator is biased even when
every component is individually well fitted.

The pipeline instead draws `v` from a **linear-Gaussian regression of `v` on (Y, rest)** —
`leftcens`'s interval-censored Gaussian AFT for the exposure block, and
`mean + homoscedastic Gaussian noise` for the covariate block.

**Those two coincide if and only if both of the following hold:**

| clause | requirement |
|---|---|
| **(A)** | the outcome model is **linear in `v`** |
| **(B)** | `p(v | rest)` is **Gaussian** |

> **Scale correction, 2026-09-09.** For the **exposure** block both clauses are on the
> **shash-transformed** scale, not the raw one: `leftcens` fits its Gaussian AFT after
> `x_to_z()`, so what must hold is that the outcome is linear in `z(v)` and `p(z(v) | rest)`
> is Gaussian. That is a *weaker* requirement than the raw-scale version — it is why the shash
> margin recovers skewed exposures where a plain Gaussian gives coverage 0.47 (Phase 1 §7.5)
> — and every clause-(B) statement about the exposure block in §§2–4 should be read on that
> scale. For the **covariate** block the raw-scale statement is correct as written.

Under (A) and (B) the product of two Gaussians is Gaussian, its mean is linear in `Y`, and
its variance is free of `Y` — exactly what a linear-Gaussian regression represents. Violate
either and the target acquires a non-linear mean, a `Y`-dependent variance, skew, or all
three, and no amount of flexibility in modelling the *mean* recovers it.

---

## 1b. When a violation reaches the estimand

**Added 2026-09-09.** Clauses (A) and (B) say when the *draw* is right. They say nothing about
when a wrong draw *matters*, and for seventeen tracks that gap was invisible because every
scenario happened to sit on the same side of it. V17 measured it, and it belongs beside the
condition rather than in the findings catalogue.

> **Propagation.** A violation of the condition for variable `v` biases estimand `θ` only if
> `v` is (i) **conditioned on by the analysis model**, and (ii) **on an open path between the
> exposure and the outcome** in the data-generating process.

Measured (V17, the Z-block `Y`-omission penalty on the focal exposure coefficient):

| `v` | adjusted for? | on an open X–Y path? | penalty |
|---|---|---|---|
| precision covariate | yes | **no** | **+1.2 pp** |
| confounder (fork) | yes | yes | **+25.5 pp** |
| mediator (pipe) | yes | yes | **+28.2 pp** |
| collider | **no** *(correctly omitted)* | yes | **−0.2 pp** |

A factor of ~22 in the same arm, one arrow apart. **This is why fifteen tracks measured the
covariate block and found small numbers** — they were all in the first row.

**Note what the criterion is not.** An earlier version of this section keyed on the
covariate's *association with the exposure*. The mediator refutes that: a pipe has no
confounding path, yet produces the largest penalty of any role. Association is not the
criterion; **adjusted-for-and-on-a-path** is.

**And (i)–(ii) are necessary, not sufficient.** V22 found a covariate that satisfies both —
adjusted for, on a pipe — where a *misspecified* imputation draw is nonetheless unbiased for
the focal coefficient (−0.14%), because the misspecification lies in a direction the analysis
absorbs: the error in approximating `g(logX1)` is taken up by `logX1`'s own coefficient. So a
third question follows the first two: **is the violation orthogonal to what the analysis
already adjusts for?** If it is, it does not reach that estimand — though it would reach an
estimand depending on `v`'s own coefficient, for which no such absorption exists.

**Consequence for reading this document.** Every magnitude in §§2–4 is an estimand-specific
number. The condition is a property of the imputation; whether violating it costs anything is
a property of the imputation *and* the DAG *and* the estimand. All three have to be named.

### 1a(ii). Congeniality is a statement about the SCALE too, not only the variables

**Added 2026-09-10.** Clause (A) requires the imputation model to be able to represent the
analysis model's conditional. That is usually read as a statement about which *variables* the
imputation conditions on. It is equally a statement about which **functional** of them.

If the analysis model is linear in `log X` and a block of the imputation regresses on `X`
itself, clause (A) fails — a linear predictor in `X` cannot represent a conditional that is
linear in `log X`, at any sample size. No new mechanism is involved; it is the same clause, and
the consequence is the same attenuation-plus-borrowing of §0c.

**This was a live defect in the shipped path, not a hypothetical.** `00_censored_exposure.R`
returned the censored exposure to the data scale *inside* the block-FCS sweep loop, so every
subsequent covariate block regressed on `exp(log X)` while the analysis model — and the DGP —
were linear in `log X`. Measured: both covariate-block arms sat **3–5 pp** above their harness
equivalents, with two entirely different imputation engines; moving the conversion out of the
loop closed all four gaps to **≤ 0.20 pp** and moved the harness arms and the oracle by
**0.00 pp**. (`ROADMAP.md`; `phase1/derive_zblock_scale.R`.)

**How much is lost, derived.** Write `L = logX`, `W` for the block's other predictors, and split
`L = μ_W + L̃` with `τ² = Var(μ_W)`, `σ̃² = Var(L̃)`, `s² = τ² + σ̃²`. Then `R = e^L = A·B` with
`A = e^{μ_W}` independent of `B = e^{L̃}`. Residualising `R` on `W` and applying Stein's identity
to `Cov(L̃, e^{L̃}) = σ̃²e^{σ̃²/2}`:

> `Cov(L̃, R̃) = E[A]E[B]·σ̃²`,  `Var(R̃) = E[A]²E[B]²(e^{s²} − 1 − τ²)`
>
> ⟹ **`Corr(L, R | W)² = σ̃² / (e^{s²} − 1 − τ²)`**

the fraction of the explainable *partial* variance a linear predictor in `R` can retain. It
nests the marginal case: `τ² = 0` gives `s²/(e^{s²} − 1)`.

| cell | `s²` | `τ²` | derived | measured | marginal form |
|---|---|---|---|---|---|
| `zr_fork` | 1.360 | 0.566 | **0.3408** | 0.3386 | 0.4696 |
| `zr_pipe` | 1.000 | 0.420 | **0.4466** | 0.4473 | 0.5818 |
| `zr_collider` | 1.000 | 0.314 | **0.4882** | 0.4869 | 0.5818 |

No free parameters; agreement ≤0.002 across three different `(s², τ²)` pairs, with the collider
out-of-sample. **Conditioning on `W` is what the marginal form was missing** — `Y` carries `log X`
directly, so the governing quantity was always a partial correlation.

**From the retained fraction to the block's error.** The fitted coefficient on `R̃` is
`c = βσ̃² / [E[A]E[B](e^{s²} − 1 − τ²)]`, where `β` is the true partial coefficient of `L`. So
`Cov(D, L̃) = βσ̃²f − βσ̃² = −βσ̃²(1 − f)`, and the imputation error's `L`-loading is

> **`−β(1 − f)`**

with `f` entering exactly once. Verified: −0.2208 / −0.2239 / −0.3061 derived against −0.2168 /
−0.2211 / −0.3016 measured (~2%, with a consistent sign — the linearisation's second-order
term). Carrying the exposure block too, an error `D` moves the censored draw by `λD`, giving the
fixed point `D = a₀L̃/(1 − a₀λ)`; measured `a₀λ = −0.088`, an **8%** shrinkage.

**The exposure block's response, and why truncation matters.** An error `D` in the covariate
moves the censored exposure draw. The draw is `L | rest ~ N(m, σ_x²)` **truncated above at the
LOD** `c`, so with `α = (c − m)/σ_x` and `h = φ/Φ`, `E[L | L ≤ c] = m − σ_x h(α)` and
`h′ = −αh − h²`, giving

> `∂E[L | L ≤ c] / ∂m = 1 − αh − h² = ` **`Var(L | L ≤ c) / σ_x²`**

— the truncated-to-untruncated *variance ratio*. Truncation therefore attenuates the block's
sensitivity, and heavily: measured **0.460** and **0.465** at 40% non-detects, so `λ_eff` is less
than half `λ`. Ignoring this over-predicts the bias by ~30%.

**The completed chain.** `f = σ̃²/(e^{s²} − 1 − τ²)` → error `L`-loading `−β(1 − f)` →
`λ_eff = λ·Var(L|L≤c)/σ_x²` → fixed point `D = a₀L̃/(1 − a₀λ_eff)` → the imputation **replaces**
the covariate (`Ẑ = m_wrong + e_wrong`, not truth-plus-noise) → the `L`-coefficient shift in the
contaminated analysis. Against the measured wrapper effect:

| cell | derived | measured | ratio |
|---|---|---|---|
| `zr_pipe` | **+5.42% ± 0.03** | +5.36 pp | **1.011** |
| `zr_fork` | **+5.05% ± 0.03** | +4.40 pp | 1.147 |

**One cell closed to 1%, the other 15% high** — and the 15% is outside the derivation's own
Monte-Carlo error, so it is a real residual rather than noise.

**The fork residual is NOT a linearisation error (tested 2026-09-11, refuted).** The obvious
suspect was the chain's linearisation of `exp`: the fork has `Var(logX)` = 1.36 against 1.00, so
`e^{s²}` is 3.90 against 2.72, and link 1's miss was already +0.004 in a consistent direction.
That gives a falsifiable prediction — the agreement must *degrade* as `s²` grows. Swept over
`s²` = 0.52 to 2.05 against an independent two-block simulation, the ratio moves the other way:
**0.841 → 0.879 → 0.900 → 0.908**, i.e. steadily *toward* 1. The residual is therefore not the
`exp` linearisation, and the fork's larger `s²` is not what distinguishes it.

> **A caution on that sweep, which also applies to anyone extending it.** The two-block
> simulation used for the comparison does **not** reproduce the pipeline's measured +4.40%
> (it reads +5.76% at the matching `s²`, 31% high), so it can refute a hypothesis about the
> closed form's *shape* but cannot adjudicate the pipeline's *level*. Its first version was worse
> still (+17%), because the exposure block was fitted with an `lm` on the above-LOD rows — a
> response-truncated regression, `σ` ~24% low, which is the same defect §6b documents for the
> DAG draw's parent factor. It was caught only because one point of the sweep had an
> independently measured value to check against; without that anchor, "roughly flat in `s²`"
> would have read as a clean refutation for the wrong reason. **The fork residual remains
> unexplained**, with the `exp`-linearisation hypothesis eliminated.
> `phase1/derive_fork_residual.R`.

> **Three wrong turns, recorded so they are not retried.** (i) That the exposure-block feedback
> explained an earlier 2× overshoot — it is an 8% effect. (ii) That the draw's inflated residual
> should be added as noise on top of the truth — that is the *classical measurement-error*
> structure, and it double-counts the residual; an imputation **replaces** the value. (iii) That
> Rubin pooling averages the noise contribution away — it does not, because attenuation is a bias
> in each imputation's estimate, not variance (measured identical at `m` = 1, 5 and 30). What
> actually closed most of the original gap was finding two errors in the projection step:
> regressing the error on a design holding the *true* covariate rather than the imputed one, and
> multiplying by the missing fraction when the error is already zero on observed rows.
> `ROADMAP.md`.

### 1b(i). The adjustment set is an estimand choice, and it changes what propagates

**Added 2026-09-09.** Condition (i) says "conditioned on by the analysis model" as though the
adjustment set were given. It is not — it follows from the estimand, and the DAG decides which:

| DAG | estimand | `Z` in the analysis model? |
|---|---|---|
| **fork** (confounder) | the causal effect | **must be** — otherwise the backdoor stays open |
| **pipe** (mediator) | **total** effect (`X→Y` plus `X→Z→Y`) | **must not be** |
| **pipe** | **direct** effect | **must be** |
| collider | any | **must not be** |

Every pipe cell in V17–V23 used the *direct* effect, so the total-effect case — a legitimate
and common target — had never been tested. Measured (40% censored, `m` = 20, 200 reps,
`Z1` always **in** the imputation):

| cell | estimand | truth | rel. bias | `bias/SE` | coverage |
|---|---|---|---|---|---|
| pipe | direct | 0.400 | +0.54% | 0.06 | 0.975 |
| pipe | **total** | 0.700 | **+0.05%** | **0.01** | 0.955 |
| `pipe_nl` | direct | 0.400 | +22.75% | 1.86 | 0.600 |
| `pipe_nl` | **total** | 0.830 | **+17.85%** | **3.29** | **0.100** |

**Two results, and they point opposite ways.**

**Where the covariate draw is the problem, the estimand choice protects you.** In the linear
pipe the total effect is essentially exact (+0.05%, `bias/SE` 0.01), and dropping `Z1` from the
*imputation* as well costs the total effect nothing (+0.03%) while costing the direct effect
**+5.07%**. That is condition (i) doing exactly what it says: an analysis that does not
condition on `Z` cannot be reached through `Z`.

**Where the exposure draw is the problem, it does not.** In `pipe_nl` the total effect is
+17.85% — smaller in *relative* terms than the direct effect's +22.75%, but **`bias/SE` 3.29
against 1.86, and coverage 0.100 against 0.600**. Worse in the currency that decides (§4b),
because the total effect has the smaller standard error.

**So condition (i) governs covariate-draw violations only.** A violation in the **exposure**
draw reaches *every* estimand involving `X`, whatever the adjustment set, because `X` is in
every model. Choosing the total effect is a real mitigation for one class of defect and no
mitigation at all for the other.

**A prediction of mine was refuted here, and the reasoning error is worth recording.** I
predicted the total effect would be biased *toward the direct effect* whenever the imputation
conditions on `Z`, on the grounds that the imputer then "knows" the part of `Y` that `Z`
explains. Measured: +0.05%. The error was conflating *an imputer richer than the analyst's
model* — which classically affects the **variance**, and does so here, mildly and
conservatively (coverage 0.955–0.975, width/SE ≈ 1.01) — with *an imputer that is wrong*, which
affects the point estimate. If the draw is from the correct posterior given observed data, the
completed data has the correct joint distribution and **every** functional of it is
consistently estimable, including one the imputer's own model does not name.

---

## 2. The boundary, derived and verified

### Covariate block

With `μ(z) = a + γ₁z`, `z ~ N(0,1)`, and `Y | z ~ N(μ(z), σ²)`:

```
log p(z | Y)  =  −(Y − a − γ₁z)² / 2σ²  −  z²/2  +  const
```

which is quadratic in `z`, hence Gaussian, with

> **precision = 1 + γ₁²/σ²**,  so  **SD = (1 + γ₁²/σ²)^(−1/2)**

For the harness (`γ₁` = 0.5, `σ` = 1): **SD = (1.25)^(−1/2) = 0.8944**, at every `Y`.

| | value |
|---|---|
| theory | **0.8944** |
| numerical integration, `Y` ∈ {−2,−1,0,1,2} | **0.8944** at all five |
| V12 measured | **0.894** |

### Exposure block

With `Y | x ~ N(a + b₁x, σ²)` and `x | other exposures ~ N(m, s²)`:

> **SD = (1/s² + b₁²/σ²)^(−1/2)**

For the harness, exchangeable `ρ` = 0.4 across three exposures gives
`s² = 1 − s₁₂S₂₂⁻¹s₂₁ = 0.7714`, so `s` = 0.8783, and with `b₁` = 0.4, `σ` = 1:

| | value |
|---|---|
| theory | **0.8287** |
| numerical integration, five values of `Y` | **0.8287** at all five |
| V9's probe, exposure prior SD `s` measured empirically | 0.8774 vs analytic 0.8783 |

**Both blocks' safe regime is therefore fully explained**, not merely observed. This is why
V1 and V2 found the engine sound: their DGP satisfies (A) and (B), so the shipped linear
conditional *is* the correct conditional, and no bias should be expected — nor was any found
beyond 2.7%.

---

## 3. Every finding is a named violation

| track | finding | clause violated | in which block |
|---|---|---|---|
| **V1, V2** | bias ≤0.8%, coverage 0.957; holds across ten scenarios | **none** — (A) and (B) both hold | — |
| **V4, V5** | intervals too narrow; `B` understated 8–59% | neither: `θ` not drawn from its posterior (a *properness* failure, orthogonal to congeniality) | Z |
| **V7** | BART fixes properness, best bias in all cells | as V4 — flexibility of the *mean* was never the issue | Z |
| **V3, V9** | mixture curvature destroyed (57%); SMC restores 89–100% | **(A)** — outcome quadratic in the exposure | **X** |
| **V10** | MAR no worse than MCAR (slightly better) | **none** — MAR changes *which* rows are missing, not the conditional's form | — |
| **V11** | ≈3 pp for **every** imputer, spread 0.44 pp | **(A)** in the covariate: all arms model the mean differently and share the wrong shape | Z |
| **V12** | shape accounts for 32–49%; SD 0.63→1.22, skew →−1.34 | **(A)** ⇒ conditional non-Gaussian ⇒ homoscedastic draw wrong | Z |
| **V13** | exposure draw dominates (3.1–3.4 of 3.9 pp) | **(A)** violated harder in X than in Z | **X** ≫ Z |
| **V14** | importance sampling fails; `μ(x) = Y` has a far root | **(A)** ⇒ the target is not where a linear conditional puts it | X |

Three consequences fall straight out of the table and were each discovered the hard way:

**Why V11's arms degraded equally.** BART, `mice pmm`, forest and bootstrapped forest differ
only in how they model the conditional **mean**. Clause (A)'s failure corrupts the
conditional's **shape**. A defect they share cannot be distinguished by a property they
differ in — hence 0.44 pp of spread across arms that differ enormously in flexibility.

**Why V10 found nothing.** MAR alters the missingness *mechanism*, not `p(v | Y, rest)`.
Neither clause mentions the mechanism, so no bias was predicted and none was found. This is
the one track whose null result the theory explains rather than merely records.

**Why fixing one block backfired (V13).** Bias contributions from the two blocks add
approximately, and here they carry **opposite signs**. Correcting Z alone moved the total
from −2.19% to −3.62% because it removed an error that was partly cancelling X's. Nothing
about congeniality guarantees that a partial fix is an improvement.

---

## 3b. What is *not yet validated* — the condition itself

**Added 2026-09-03.** Sections 1–3 state the congeniality condition and prove its safe
boundary; §5 catalogues thirteen tracks as instances of it. Every one of those tracks tests
a **refinement** — whether the conditional has the right functional form (V3, V9, V11), the
right shape (V12), the right margin (Phase 1 §7.5), the right block (V13, V14). **None of
them varies the condition's own antecedent**: every arm in every track had `Y` in the
imputation model. The foundation is the one thing that has never been moved.

The project's cited evidence is Phase 1's **H1** — `leftcens_prestep` attenuating the focal
coefficient −5.6% at 20% non-detects and −14.4% at 40%, coverage 0.92 → 0.82. That is real
and its direction and monotonicity support the claim, but the arm it compares against
differs in **three** ways at once: no `Y`, **no `Z` covariates either**, and a **different
estimator** (`leftcens::gsimp_mi`, a copula MI of the exposure matrix, not an
interval-censored conditional regression). Attributing the effect to `Y` is an inference
from that design, not a measurement.

Two cells in the same document cut the other way. On the **mixture** surface the no-`Y`
pre-step (+3.8% / +4.0%) **beat** the `Y`-aware linear draw (+7.7% / +19.6%, coverage 0.62);
under **skew 0.75** the no-`Y` copula (−1.2% / −5.7%) beat the `Y`-aware Gaussian (−8.6% /
−14.2%, coverage 0.47). Neither refutes the condition as stated here — it is a
**conjunction**, and in both cells a different clause is the one being broken. But both
refute the claim in its loose form, *"include `Y` or you attenuate"*, and that loose form is
how the condition tends to get quoted. **Including `Y` with the wrong functional form or the
wrong margin is worse than omitting it.**

> **TESTED, 2026-09-03.** V16 measured it: removing `Y` from both blocks of the shipped
> engine costs **−13.1 pp** (derived −14.6), from the exposure block alone **−14.4**
> (derived −15.8). **The antecedent holds**, isolated inside one estimator for the first
> time, so Phase 1's H1 number was indeed dominated by the `Y` omission — which could only
> be known by separating it. Two things did not survive: the **positive leak** (predicted
> +1.6 pp, measured +0.09 ± 0.13 — the blocks are additive here), and the loose form
> *"include `Y` or you attenuate"*, which now stands only for the **exposure** block.
> `phase1/FINDINGS_v16.md`.

**V16 tests the antecedent directly** (`PLAN_pipeline_validation.md` §8i): one estimator, `Y`
removed from one block at a time. Its prediction is *derived* rather than fitted — **−14.9
pp** when `Y` is omitted everywhere (the configuration a real analysis produces), and
**zero** for the covariate block alone in cells where `Z ⊥ X`, because with no confounding
path there is nothing to under-adjust. If that second figure lands, it establishes the
**weaker half** of a scope qualifier the condition does not currently carry. The stronger
half needed a DGP in which a covariate causes an exposure, which this harness did not have —
drawing the DAG is what surfaced that, at the cost of a registered prediction withdrawn
before the run. **V17 built the missing structures, and the derivation (2026-09-03,
`phase1/predict_roles.R`) states the qualifier sharply:**

**MEASURED, 2026-09-03** (500 reps; the derived predictions are in brackets). The
propagation rule these numbers established has since been promoted to **§1b**, beside the
condition itself, because it governs how every other magnitude in this document is read:

| covariate `Z` | adjusted for? | on an X–Y path? | penalty, MCAR | penalty, MAR-on-`Y` |
|---|---|---|---|---|
| precision | yes | no | **+1.2 pp** *(≈0)* | **+0.9 pp** *(+0.1)* |
| fork (confounder) | yes | yes | **+25.5 pp** *(+26.7)* | **+26.2 pp** *(+32.4)* |
| pipe (mediator) | yes | yes | **+28.2 pp** *(+30.4)* | **+28.9 pp** *(+36.1)* |
| collider | **no** | yes | **−0.2 pp** *(−0.7)* | **+2.1 pp** *(+3.1)* |
| mixed | yes | yes | **+22.6 pp** *(+19.5)* | **+22.0 pp** *(+23.8)* |

**The condition binds a covariate's draw exactly when that covariate is adjusted for *and*
lies on an open X–Y path** — whether the path is confounding or mediating. Both other cells
are near zero: a covariate off any path contributes nothing to `β₁` however badly it is
imputed, and a covariate on a path that the analysis correctly omits never reaches `β₁` at
all.

**The missingness mechanism turns out to be almost nothing.** The derivation predicted
MAR-on-`Y` would add **+4 to +6 pp** on-path — omitting `Y` there drops the variable the
mechanism depends on, so the covariate imputation is invalid and not merely inefficient.
**Measured, MAR − MCAR is +0.70, +0.73, −0.58 and +2.31 pp** (fork, pipe, mixed, collider):
three of four within ±1 pp of zero, one with the wrong sign, and only the collider clearly
positive. That prediction is **refuted**.

The conjecture for why — untested, and discriminable with a `mice pmm` Z block — is that the
covariate draw recovers most of the `Y`-relevant information from the exposures and the other
covariate, which a flexible imputer can exploit and the linear stand-in could not.

**So the qualitative claim survives in its stronger form: the path structure dominates the
mechanism.** Not merely by a margin — the mechanism term is within noise of zero. In the
precision row the `Z` imputation is genuinely invalid under MAR and `β₁` still does not move.
A biased covariate distribution cannot reach the focal estimand along a path that is not
there, and whether the covariate is MCAR or MAR barely registers next to whether the path
exists.

This is **not** the qualifier §8i proposed before the derivation. That version keyed on the
covariate's *association with the exposure*, and it would have got the pipe wrong: a mediator
has no confounding path, yet produces the largest penalty of any role. Association is not the
criterion; adjusted-for-and-on-a-path is.

**What it costs the *pipeline*, not just the record.** The paired design cancels the
reference arm's own bias — which is what makes those contrasts clean, and what hid this until
a gate was added for it. `pipeline_bartMI` is the **shipped configuration**, `Y` in both
blocks, and under a causally active covariate it is biased **+5.0% (fork), +5.4% (pipe),
+8.5% and +10.4% under MAR**, against ≤2.3% in every precision cell. The oracle is unbiased in
all four, so this is the imputation's bias and not the design's. Coverage sits at 0.93–0.97
with honestly-sized intervals — the **fourth** independent instance of coverage failing to
detect a large bias (V3, V8, V15, V17) — and it covers only because at `n` = 800 the bias is
0.35–0.67 empirical SE. Relative bias is roughly constant in `n` while the SE falls as
`1/√n`, so that concealment looked like a sample-size accident.

> **MEASURED, 2026-09-03 (V18) — and the projection was wrong.** Over five `n` levels to
> 12,800, the bias decays as **n^−0.335** (CI [−0.378, −0.291]), rejecting both registered
> accounts: constant bias (slope 0) and a finite-sample artefact (slope −0.5). Coverage at
> n = 12,800 is **0.867–0.923**, not the 0.32–0.73 projected here. The coverage *model* was
> sound — fed measured inputs it reproduces all 20 points to within 0.024 — the
> constant-relative-bias premise was not.
>
> **The defect is still permanent.** Bias falls as n^−0.335 while the SE falls as n^−0.50, so
> bias/SE grows as **n^0.17** and coverage is unbounded below — just 3× more slowly in the
> exponent than claimed. At every sample size tested the shipped default carries 2–10% bias
> that coverage does not reveal. `phase1/FINDINGS_v18.md`.

**What it costs the record.** V12, V13 and V16 all reason about the Z block's contribution,
and all of them measured it in the precision structure — the one row of that table where the
Z block barely matters (+1.2 pp measured). V13's headline, that a Z-only fix moves bias *away* from
truth, is a statement about an inert covariate and may not survive a confounder or a
mediator. None of it is wrong; all of it is narrower than it reads.

**The blocks are not separable, and the derivation says by how much.** Because block-FCS
alternates, `Y` reaches the exposure block through a `Y`-informed covariate draw even when it
is removed from the exposure block's own predictor set. The predicted leak is **+1.4 pp**,
and its sign is the interesting part: removing `Y` from the exposure block *while leaving it
in the covariate block* is predicted to be **worse** (−15.9 pp) than removing it from both
(−14.9 pp) — a `Y`-informed `Z̃` partly absorbs `Y`, and adjusting for it alongside a
`Y`-blind `X̃` over-adjusts. If that holds, **partial compliance with the condition is worse
than none**, which is a stronger and more useful statement than the condition alone makes,
and belongs in §5's taxonomy as a general property of block-wise imputation rather than a
quirk of this engine.

---

## 3c. A derived mechanism, and a reason for the quadratic form

**Added 2026-09-08.** Until now this document recorded two measured exponents it could not
produce (§4). V22's ~20% censored-exposure defect looked like a third. It is not — it is
derivable, and the derivation supplies something the earlier laws never had: a **reason for
the square**.

**Why the defect must exist.** `leftcens` draws a censored exposure from a conditional
**linear in its predictors**. Where a covariate satisfies `Z = g(X)` with `g` non-linear, the
true conditional for a censored `X` contains

> `p(Z | X) = N(Z; g(X), s²)`, whose log contributes `−(Z − g(X))² / 2s²`

— non-linear in `X` wherever `g` is curved, and therefore outside what a linear-Gaussian draw
can represent. This follows from reading the conditional `leftcens` fits; it needs no
simulation.

**Why it is asymptotic.** The fitted coefficients converge to the best *linear* approximation
of a non-linear target. That approximation stays wrong at every `n`, so the bias does not
decay — which is what V22 measured (+20.8 / +20.3 / +20.3% across a 16× range of `n`).

**What governs its size.** Not the amplitude of `g` but **how much of `g` a straight line
cannot express where the censored mass sits**: the density-weighted residual SD of `g` after
its best linear fit below the LOD. Write it `u` (`ef_unrep_curvature()` in
`phase1/R/dgp.R`). A 5-setting probe separates the two: a `tanh`-only arrow with `sd(g)` =
1.20 — as large as V22's — gives **+0.09%** bias, because `tanh` saturates and is nearly
linear-fittable over the censored region, while V22's arrow at `u` = 0.262 gives **+22.3%**.

**Why `u` enters SQUARED.** `u` is the residual of a density-weighted L² projection of `g`
onto the linear span, so it is orthogonal to that span, and a first-order expansion of the
bias functional would pair the misspecification with the score and annihilate it — leaving a
second-order leading term:

> ⚠️ **This argument is looser than it reads; see §0.** The fit converges to the **KL**
> projection, not the L² projection of a log-density, so the orthogonality that actually holds
> is to the fitted family's *score* (moment matching). And `u` is the L² residual of `g`
> alone — neither a KL distance nor a log-density residual — so it is a **proxy standing in
> for the theory's quantity**, not the quantity itself. That it predicts four of six cell
> levels to 1.5 pp is empirically striking and theoretically unexplained.

> **bias ≈ k · u²**

This is a heuristic argument rather than a proof.

> **TESTED, 2026-09-09 (V23) — and the test had no power.** The law's *levels* hold well:
> four of six predictions within 1.5 pp from a single calibration point, the null cell at
> −0.80%, and the bias flat to ≤0.27 pp across a 4× range of `n`, confirming it is asymptotic.
> But the exponent came back at **1.638 with a 95% CI of [0.994, 2.282]** — containing 2 and 1
> alike, on six points. **So the orthogonality argument is neither confirmed nor refuted.**
> The point estimate leans below 2, and the two largest-`u` cells are over-predicted (the
> biggest by 8.8 pp, with a freely fitted constant of 249 rather than 300), so **the
> one-parameter `u²` form saturates** and something is missing at large `u`.
>
> Distinguishing 1.64 from 2 needs more `u` levels at the low end, where the bias is small and
> Monte-Carlo error dominates — materially more expensive, and worth doing only if this
> argument is load-bearing for something else. Until then it stands as a *motivated* form
> rather than a demonstrated one. `phase1/FINDINGS_v23.md`.

**What it may explain beyond itself.** §4's `f²` law for curvature attenuation has been frank
curve-fitting since V15. If the orthogonality argument is right, any misspecification
expressible as an L²-projection residual should enter the bias quadratically — and the censored
fraction's effect may be one. That would turn `f²` from a fitted shape into a consequence.
**Not claimed here**; recorded as the first candidate mechanism this document has had for a
rate rather than a magnitude — and after V23, still a candidate, because the exponent test
could not tell 1 from 2.

---

## 4. What is *not* derived

**No closed form for the magnitudes when (A) fails.** Under a non-linear outcome the
log-conditional is quartic; its moments require numerical integration. Two attempts to
approximate them were tested against the Monte-Carlo values and **both failed**:

| attempt | prediction | truth | verdict |
|---|---|---|---|
| **Laplace** (mode + curvature): `precision = 1 + [μ'(z*)² − (Y−μ(z*))·2b_zq]/σ²` | SD *falls* with `Y`: 0.669 → 0.577 | SD *rises*: 0.626 → 1.216 | **wrong in sign** — a mode-based quadratic cannot represent a strongly skewed target (skew −1.34) |
| **Second root**: `μ(z) = Y` is quadratic, so inflation should track the far root's prior weight | far root's relative weight *falls* 2.3e−01 → 2.4e−02 as `Y` rises | SD *rises* over the same range | **wrong direction** — the far root's weight moves opposite to the spread |

A third tempting prediction also fails: that curvature attenuation should be roughly
**proportional to the censored fraction**, since imputed rows carry no curvature by
construction. V3 measured 20% non-detects → 11.7% excess and 40% → 56.7% — markedly
superlinear.

**V15 (2026-09-02) replaced that gap with a law.** A six-point sweep in `f` confirmed the
registered prediction `excess% = −323.4·f²` in all four previously unmeasured cells (worst
miss 9.3 pp against a ±10 pp gate), and a freely fitted exponent gives **1.90**, with the
log–log slope CI **[1.98, 2.97]** — excluding 1, containing 2. So the exponent is now
evidence, not a guess:

| `f` | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|---|---|---|---|---|---|---|
| excess | −1.3% | −13.9% | −29.6% | −55.6% | −90.2% | −116.5% |

*(§3c now supplies a candidate mechanism for the quadratic form, which may reach `f²`
as well. Read it before treating the exponents below as wholly unexplained.)*

**A second measured rate joined it on 2026-09-03.** V18 found the shipped default's bias
under a causally active covariate decaying as **n^−1/3** — also measured, also underived.
The theory now carries **two exponents it cannot produce**: `f²` for curvature attenuation in
the censored fraction, and `n^−1/3` for imputation bias in the sample size. The V18 one carried a
conjecture — that `n^−1/3` is the Z-block imputer's own nonparametric convergence rate
showing up in the estimand — and **V19 refuted it** (2026-09-07). A second nonparametric
imputer (`properBoot`) shows **no decay at all** (slope −0.029 [−0.087, +0.028]), and two
*correctly specified* parametric imputers sit at **≈0** rather than the predicted −1/2. Family
membership predicts nothing; **the split is BART versus everything else**. So both exponents
remain measured and underived, and this document still has no mechanism for a rate.

What V19 established instead is not theoretical but it is the more useful half: **BART is the
only Z-block imputer whose bias decays with `n`**, and the only one whose coverage survives to
n = 12,800 under a causally active covariate (0.87–0.92, against 0.78–0.84 for `micePmm`,
0.61–0.62 for `properZ`, and **0.00** for `properBoot`). `phase1/FINDINGS_v19.md`.

**What is still not derived is the exponent itself, and the constant.** The status improved
from *unexplained observation* to **unexplained law**, which is a sharper target: a
derivation must now produce a power, not merely a direction. The leading conjecture is a
product of two `f`-linear losses — (i) the fraction of rows whose exposure is a draw from a
conditional linear in its own predictors, hence carrying no curvature by clause (A), and
(ii) the fraction of the *estimand's evaluation range* that lies below the detection limit,
since `curv_X1` is a contrast at fixed quantiles of X1 and the LOD climbs that distribution
as `f` grows. Their product is `f²`. **This is a conjecture with a discriminating test**:
factor (ii) depends on where the contrast is evaluated and factor (i) does not, so moving
`q_lo` from 0.25 to 0.40 should move the curve under the conjecture and leave it alone
otherwise. See FINDINGS_v15.md; that comparison is the natural next prediction.

V15 also found that the damage is not confined to attenuation: at `f` = 0.60 the fitted
surface **inverts**, curving the wrong way in 40% of replications, while nominal coverage
stays at 0.92 because the interval widens to 2.7× its `f` = 0.10 width. **Coverage is not a
guard against this failure mode** — a point worth carrying into any guidance about
censored mixtures.

### Revised by V14 (2026-09-01)

V9 found that supplying the *form* of the outcome model and **estimating** its coefficients
matched an oracle that knew them, to 0.2 pp — the basis for saying "the functional form is
what matters, not the parameters". **V14 qualifies that.** In the scalar non-linear-outcome
setting the same substitution costs **0.53 pp on average and 1.36 pp at worst** — an order of
magnitude more.

The two measurements are not in conflict (V9: BKMR mixture estimands; V14: a scalar
coefficient), but the general claim is **setting-dependent** and should not be carried
forward unqualified. What survives is the weaker, still-useful statement: getting the form
right recovers most of the penalty (83% here, 89–100% in V9), and getting the parameters
right recovers the rest.

**Not addressed at all:** why the shape effect has a fixed *sign* (V13's −1.3 to −1.9 pp
regardless of the exposure draw); the +15–17% BKMR oracle floor on `curv_X1` (V3), shown not
to be chain length or the GPP approximation; and why MAR is *easier* than MCAR (V10).

---

## 4b. Two rates, not one — and why the currency of a claim matters

**Added 2026-09-09.** Everything this document measures is a *bias*. Every decision made from
it compares that bias to a *standard error*. Those two move at different rates in `n`, so a
claim stated in one currency can be silently wrong in the other — and for seventeen tracks
this project stated all of its bars in the wrong one.

> bias decays as **n^−1/3** (V18, replicated in V19, V20, V21)
> the standard error decays as **n^−1/2**
> therefore **bias / SE grows as n^(1/6)** — measured at n^0.10 to n^0.17

**A relative-bias bar therefore gets easier to pass exactly as the bias becomes more
consequential.** What the old 10% bar permitted, using each run's own measured `emp_se`:

| `n` | emp. SE | 10% of `b₁` | in SE units | share of a 2.8-SE detectable effect |
|---|---|---|---|---|
| 800 | 0.054 | 0.040 | 0.74 | 26% |
| 3,200 | 0.028 | 0.040 | 1.45 | 52% |
| 12,800 | 0.015 | 0.040 | 2.67 | **95%** |

At `n` = 12,800 a result passing the old bar could carry a bias nearly as large as the effect
the study is powered to find. **The convention changed on 2026-09-08**
(`PLAN_pipeline_validation.md` §11): report `bias/SE`, gate at ≤ 0.3, and leave the
conversion to a share of the detectable effect to the reader — **because the target effect is
a choice, not a constant.** A 5% change is a common epidemiological convention; a
no-threshold setting, such as ionising radiation and cancer risk, has no such floor and the
same absolute bias consumes a much larger share of whatever is being detected.

**This is a theoretical point, not only an administrative one.** It says that a defect's
*consequence* has a different exponent from its *magnitude*, so any account that predicts one
without the other is incomplete. §4's `f²` and `n^−1/3` are magnitude laws; the decision-
relevant quantity is their ratio against `n^−1/2`, and nothing here derives that ratio from
first principles.

### Direction, and why it is not a detail

A magnitude cannot be acted on without a sign. The three cases differ in kind:

| direction | meaning | how it can be handled |
|---|---|---|
| **toward the null** | conservative — risks missing a real effect | disclosable; a stated attenuation can be lived with |
| **away from the null** | anti-conservative — inflates or invents an exposure–response | **dangerous**, especially where any positive finding drives policy |
| **sign flips** with a design parameter | neither | **no single correction and no single caveat covers it** |

The third is not hypothetical: V15 found the curvature estimate **inverting** past `f` ≈ 0.55,
so the same pipeline both attenuates and reverses depending on the non-detect rate.

**And the record has an asymmetry worth stating plainly.** Across 118 shipped-arm cells in
V16–V23, 75 are away from the null and 43 toward it — and the split is not random. **The
defects this project has closed were the conservative ones** (an omitted `Y`, −13%; the
improper Z block; `forest_boot`, −23%). **Both defects still open are anti-conservative**:
the covariate-block defect at +4.4% to +10.4%, and the censored-exposure defect at ~+20%
asymptotic. The full ledger is `PLAN_pipeline_validation.md` §11b.

That ordering should drive what gets fixed next. A conservative defect can be documented; an
anti-conservative one manufactures findings.

---

## 5. What the condition is good for

Even without magnitudes, (A) and (B) are **decision-useful**, because they are checkable
before running anything:

1. **Is the analysis model linear in every variable that has missing values?** If a variable
   appears in a spline, an interaction, `mo()`, a quadratic, or a mixture surface *and* has
   missing values, clause (A) is violated for that variable and the imputation is
   uncongenial. This is a property of the formula, readable at configuration time.
2. **Is the variable's own conditional plausibly Gaussian?** For skewed exposures it is not,
   which is why `leftcens`'s shash margin exists and why V1 found a Gaussian conditional
   gives coverage 0.47.
3. **Are the imputation model's parameters drawn from a posterior?** Orthogonal to both
   clauses, and the subject of V4/V5/V7.

The pipeline satisfies all three for **additive exposure–response functions with a linear
outcome model**, which is precisely the scope the README claims. Every documented failure
lies outside it, and the two open build tracks map onto the two clauses: roadmap **07** is
clause (A) in the exposure block, and the closed **08** was clause (A) in the covariate
block.

---

## 6. How this document is meant to be used from now on

**Adopted 2026-09-01** (`PLAN_pipeline_validation.md` §7b): tracks are designed forwards
from here — state the mechanism, derive a number for a cell nobody has measured, run, and
revise the theory against the result. The first such track is
[V15](PLAN_pipeline_validation.md#8h-track-v15--the-attenuation-shape-first-prediction-led-track),
which registered `excess% = 323.4·f²` and the conditions that would refute it.

That law was **frank curve-fitting to the two points in §4** — the mechanistically motivated
candidates all fit worse. It was registered anyway, labelled as such, because a
weak-but-falsifiable prediction moves the theory and a post-hoc description does not.

**Outcome: it survived** (FINDINGS_v15.md, 2026-09-02) — four for four inside the gate, one
cell to within 0.1 pp. The cycle's first pass therefore worked in the less useful of the two
directions: the theory was not corrected, it was promoted, and the open question moved from
*what shape?* to *why that exponent?*. The value of registering the prediction was not that
it was right; it was that a wrong number would have been visible as such.

**Predictions that fail stay in this document** with the superseding result beside them, on
the same principle by which findings are corrected rather than overwritten.

**Three things must be named for any magnitude quoted from here**, and §1b and §4b are why:
the **estimand** (a violation reaches some and not others), the **currency** (`bias/SE`, not
relative bias — they move at different rates in `n`), and the **direction** (conservative,
anti-conservative, or sign-flipping). A number without all three is not actionable, and every
figure in this document registered before 2026-09-08 was quoted without the second.

---

## 6b. What a correct draw would be: the DAG-factorised target, and an algorithm

**Derived and design-tested 2026-09-09.** Everything before this section says where the
pipeline goes wrong. This one says what right would look like. It is Bayes plus the DAG's
factorisation, and nothing else.

### The target

For a censored exposure `X` with parents `Pa(X)` and children `Ch(X)` among the modelled
variables (`Y` is always a child of `X`):

> **p(x | rest, x ≤ L) ∝ p(x | Pa(X)) · ∏_{C ∈ Ch(X)} p(C | x, Pa(C)\{X}) · 1{x ≤ L}**

That is just Bayes applied to the joint the DAG implies. Written out per role:

| DAG | `Pa(X)` | `Ch(X)` | target ∝ |
|---|---|---|---|
| **fork** (`Z→X`, `Z→Y`) | `{Z}` | `{Y}` | `p(x \| Z) · p(Y \| x, Z)` |
| **pipe** (`X→Z→Y`) | `∅` | `{Z, Y}` | `p(x) · ` **`p(Z \| x)`** ` · p(Y \| x, Z)` |
| **collider** (`X→Z←Y`) | `∅` | `{Z, Y}` | `p(x) · ` **`p(Z \| x, Y)`** ` · p(Y \| x)` |

**What `leftcens` implements is the fork form**, and only that: a regression of `x` on
`(Y, all covariates)` treats every covariate as a *parent*. It is exact when the joint is
Gaussian-linear — which is why the `fork` cell measures **+0.28%** — and it has **no child
factor at all**, which is why a pipe with a non-linear `X→Z` arrow fails.

### The estimand does not change the target

This is worth stating because it is not obvious and it simplifies the design. The target above
is the correct posterior **given observed data**, so the completed data has the correct joint
distribution and *every* functional of it is consistently estimable — including the total
effect, whose analysis model does not name `Z`. §1b(i) measured this: +0.05% for the total
effect from an imputation that conditions on `Z`.

So the estimand enters **only** in the analysis model applied afterwards. What it costs is
variance, not bias: a target richer than the analyst's model gives mildly conservative
inference (coverage 0.955–0.975, width/SE ≈ 1.01), which is Meng's classical result.

| | imputation target | analysis model |
|---|---|---|
| fork | the fork form | `Y ~ X + Z` |
| pipe, **direct** | the pipe form | `Y ~ X + Z` |
| pipe, **total** | the pipe form — *unchanged* | `Y ~ X` |
| collider | the collider form | `Y ~ X` |

### The algorithm

`X` is scalar and bounded above by `L`, so a grid inverse-CDF draw is exact to grid resolution
and needs no proposal, no weights and no acceptance step:

```
per outer sweep t:
  draw θ⁽ᵗ⁾ from the analysis model's posterior          # properness (V4, V21)
  draw φ⁽ᵗ⁾ from each child model's posterior
  fit the interval-censored model of X on Pa(X)          # what leftcens already does
  for each censored cell i:
    1. m_i, s_i  ← parent fit for row i
    2. grid  x_g  over (L_i − 6 s_i, L_i],  G points
    3. log w_g = log p(x_g | Pa_i)                            # parents
               + log p(Y_i | x_g, Z_i; θ⁽ᵗ⁾)                  # the outcome child
               + Σ_{C ∈ Ch(X)} log p(C_i | x_g, ·; φ⁽ᵗ⁾)      # THE MISSING TERM
    4. normalise w, inverse-CDF sample x̃_i
```

Only **step 3's third line** is new. It is empty for a fork, so the fork case needs no change.

**This is not item 07.** Item 07 proposes from the exposure prior and *reweights by
`p(Y | x, rest)`* — the same missing child factor, but with the outcome term amplified rather
than balanced, which is why V23 measured it at **+142%** with coverage 0.000. In §0c's terms it
maximises `T₂` while leaving `T₁` alone. The grid form above balances the factors instead of
reweighting by one of them.

### Design test: with the child factor known, the defect goes away

`pipe_nl`, 40% censored, `m` = 20, 200 reps, `g` supplied (so this tests the *design*, not the
estimation of `g`):

| child factor | estimand | rel. bias | `bias/SE` | coverage |
|---|---|---|---|---|
| off | direct | +6.50% | 0.96 | 0.965 |
| **on** | direct | **−0.26%** | **0.04** | 0.995 |
| off | total | +4.73% | 1.62 | 0.815 |
| **on** | total | **−0.26%** | **0.08** | 0.980 |

*(Reference: shipped `leftcens` gives +19.5% and +16.9% in these cells. The "off" rows are
better than that because this grid uses the exact parent factor and the true `θ`, so they
isolate the child factor's own contribution.)*

**The design is right.** All the remaining risk is in obtaining `p(Z | x)` below the LOD, which
§0b shows is **not identifiable** — so it has to be a declared assumption with a sensitivity
axis, not an estimate. That is the honest shape of the fix: the algorithm is exact given a
child model, and the child model is an assumption.

### `Z` is needed in the imputation even when it is absent from the model

**Added 2026-09-09.** A natural simplification is: if `Z` is not in the analysis model, the
imputation does not need it either — so only the fork and the pipe-direct cases need the new
child factor. **That is wrong, and it fails in the cases where the current code is furthest
off.**

`Z` belongs in the target iff `Z` is **observed and informative about `X`** — which has
nothing to do with the analysis model. In a pipe `Z` is a *child* of `X`, so `p(Z | x)` is
observed information about `x`, and dropping it means drawing from the wrong posterior no
matter what the analyst then fits.

**And the non-linearity does not disappear when `Z` is marginalised out — it moves.** With
`Z₁ = g(X₁) + ν` and `Y = b₁X₁ + γ₁Z₁ + …`, marginalising `Z₁` gives

> `Y = b₁X₁ + γ₁·g(X₁) + … + (γ₁ν + ε)`

so once `Z₁` is no longer conditioned on, **`Y` is non-linear in `X₁`** through `γ₁·g(X₁)`. A
linear draw is misspecified either way: **clause (B)** if you condition on `Z`, **clause (A)**
if you do not. Measured in `pipe_nl` with the total effect as the estimand:

| `Z₁` in the exposure draw | rel. bias | `bias/SE` | coverage |
|---|---|---|---|
| in | +16.87% | 2.91 | 0.145 |
| out | +14.31% | 2.45 | 0.320 |

Both badly wrong. Dropping `Z` buys 2.6 pp and costs the structure that would have let the
child factor fix it — the DAG-factorised draw with the child factor **on** gives **−0.26%** for
this same total-effect estimand.

**So the algorithm is needed in four cases, not two:**

| DAG / estimand | `Z` in the analysis? | child factor needed? |
|---|---|---|
| fork | yes | no — `Z` is a *parent*; current code is already correct (+0.28%) |
| pipe, **direct** | yes | **yes** — `p(Z \| x)` |
| pipe, **total** | **no** | **yes** — `p(Z \| x)`, unchanged |
| collider | **no** | ~~**yes**~~ — **corrected by V24: available, not necessary.** The `Y`-only draw is already unbiased there (+0.06 to +0.16% across three `n`), because `Z₁` is off every `X`–`Y` path and absent from the correct analysis model, so omitting it is *congenial*. The child factor buys 1.4 ± 0.1% of the SE — not worth an unidentifiable assumption |

Of the two rows where `Z` is absent from the analysis model, **one** is where the child factor
is indispensable (pipe-total) and one is where it turned out to be optional (collider) — see
the corrected row above. The distinction is whether `Z` is on an `X`–`Y` path, not whether it
is in the model.

### All four cases at once, and the one place the algorithm turns harmful

> **Superseded in scale by Track V24 (2026-09-09), which ran the same design through the real
> runner at `n` = 800 / 3200 / 12800 × 200 reps and confirmed every result below, with two
> additions the design test could not see — the current draw class is **asymptotically** biased
> wherever the covariate is a mediator (6–7% flat in `n`, coverage → 0.040), and a *more
> flexible* child model is *monotonically worse*. Read
> [`phase1/FINDINGS_v24.md`](phase1/FINDINGS_v24.md) first; this section is the derivation and
> the single-`n` design test it was built from.**

**Added 2026-09-09.** The four-cases table above was assembled from cells measured
separately, which cannot rule out a sign being a property of one cell rather than of the
role. `phase1/derive_dag_four.R` runs all four simultaneously — the two axes are the
covariate's role and, for a pipe, the estimand — with the two **linear** pipe cells included
as controls that separate the *factorisation* from the *non-linearity*. `n` = 2000, LCR 40%,
`m` = 20, **200 reps**, arms paired on data and imputation seed within a rep.

`Δ SE` is the paired change in the pooled SE against the `Y`-only draw, with its Monte-Carlo
error — `bias/SE` cannot show an efficiency change, since an arm that sharpens the draw moves
numerator and denominator together.

| case | arm | rel. bias | `bias/SE` | coverage | `Δ SE` (paired) |
|---|---|---|---|---|---|
| fork | no `Y` | −17.34% | 2.42 | 0.270 | +1.7 ± 0.2% |
| fork | **`Y` only** (= DAG-correct) | **+0.08%** | **0.01** | 0.965 | — |
| pipe, direct | no `Y` | −23.90% | 2.94 | 0.095 | +0.5 ± 0.2% |
| pipe, direct | `Y` only | −6.00% | 0.74 | 0.905 | — |
| pipe, direct | **+ child, linear** | **−0.07%** | **0.01** | 0.960 | +4.3 ± 0.1% |
| pipe, total | `Y` only | −6.86% | 1.54 | 0.670 | — |
| pipe, total | **+ child, linear** | **−0.24%** | **0.06** | 0.950 | −4.6 ± 0.1% |
| `pipe_nl`, direct | no `Y` | −19.23% | 1.94 | 0.530 | +0.4 ± 0.2% |
| `pipe_nl`, direct | `Y` only | +9.04% | 0.92 | 0.870 | — |
| `pipe_nl`, direct | + child, **linear** | **+26.11%** | **2.17** | 0.415 | +22.3 ± 0.2% |
| `pipe_nl`, direct | **+ child, true** | **−0.11%** | **0.01** | 0.985 | −9.3 ± 0.2% |
| `pipe_nl`, total | `Y` only | +5.65% | 1.48 | 0.695 | — |
| `pipe_nl`, total | + child, **linear** | **+21.94%** | **5.74** | **0.000** | −0.2 ± 0.2% |
| `pipe_nl`, total | **+ child, true** | **+0.00%** | **0.00** | 0.950 | +1.0 ± 0.2% |
| collider | no `Y` | −14.94% | 2.12 | 0.435 | +1.6 ± 0.1% |
| collider | `Y` only | −0.14% | 0.02 | 0.945 | — |
| collider | + child, true | −0.13% | 0.02 | 0.950 | −1.4 ± 0.1% |

Four things follow, and only the first was expected.

**1. The factorisation is right in all four cases simultaneously.** With the child model
supplied, `|bias/SE|` is at most **0.02** across fork, `pipe_nl`-direct, `pipe_nl`-total and
collider, with coverage 0.95–0.985 — against 1.9–2.9 for the incongenial draw in the same
cells. This is not a vacuous pass: every cell has real bias to remove.

**2. A *misspecified* child factor is worse than no child factor at all.** In `pipe_nl` the
linear child form takes `bias/SE` from 0.92 → **2.17** (direct) and 1.48 → **5.74** (total,
coverage **0.000**). The linear pipe controls show this is misspecification and not the factor
itself: there the same linear form is correct by construction and it *improves* both cells to
`|bias/SE|` ≤ 0.06.

> **This is the constraint on the fix, and it is sharper than "the child model is an
> assumption".** An assumption that degrades gracefully could ship with a default. This one
> does not: getting `p(Z | x)` wrong below the LOD is *more* damaging than omitting the term,
> so the algorithm cannot ship with a fitted-by-default child model. It has to be a
> **sensitivity procedure over a declared form** — and a run that declares the form wrongly is
> worse off than one that never used the algorithm. §0b's non-identifiability result is
> therefore not a caveat on this fix, it is the fix's binding constraint.

**3. Under a collider the child factor is an *efficiency* term, not a bias fix.** The `Y`-only
draw is already at −0.14% there, because with `Z₁` absent from the analysis model and off every
`X`–`Y` path, ignoring `Z₁` entirely is *congenial*. Adding `p(Z₁ | x, Y)` buys **1.4 ± 0.1%**
of the SE — real, but small against an unidentifiable assumption. **The collider is the one
case of the four where the honest recommendation is not to use the algorithm.** The four-cases
table above overstates the collider row: the child factor is *available* there, not
*indispensable*.

**4. `Δ SE` moves in both directions, and the sign tracks the estimand.** For the *direct*
effect the child factor widens intervals (+4.3% linear pipe); for the *total* effect it narrows
them (−4.6%). Consistent with §1b(i): the total-effect analysis discards `Z₁`, so the
information the child factor recovers has nowhere else to enter.

### What V24 added, at scale

**Two results needed three `n` levels, so the single-`n` design test could not reach them.**

**The draw the pipeline uses today is *asymptotically* biased wherever the covariate is a
mediator.** `dag_Yonly` — parents plus `p(Y | x, Pa(Y))`, no child factor, the class every
track from V0 to V23 used — has **flat** relative bias and therefore `bias/SE` growing as √n:

| cell | rel. bias (800 / 3200 / 12800) | `bias/SE` | coverage | exponent |
|---|---|---|---|---|
| pipe, direct | −7.48 / −6.18 / −6.18% | 0.58 → 0.97 → 1.94 | 0.945 → 0.830 → 0.535 | 0.43 |
| pipe, total | −7.37 / −6.90 / −6.68% | 1.04 → 1.97 → 3.80 | 0.830 → 0.480 → **0.040** | 0.47 |
| `pipe_nl`, direct | +7.52 / +9.25 / +8.82% | 0.48 → 1.19 → 2.28 | 0.940 → 0.780 → 0.375 | 0.56 |
| `pipe_nl`, total | +5.10 / +5.59 / +5.76% | 0.84 → 1.86 → 3.83 | 0.895 → 0.565 → **0.050** | 0.55 |
| fork | ±0.4% | ≤ 0.07 | 0.925–0.980 | — |
| collider | ≤ +0.71% | ≤ 0.08 | 0.925–0.975 | −0.03 |

**The two linear `pipe` rows are the important ones.** §3c and V22 attributed the ~20%
censored-exposure defect to the arrow's *non-linearity*, through `u²`. These cells have a
**linear-Gaussian** arrow, so `u` = 0 by construction — and they still carry 6–7%, flat over a
16× range in `n`.

`u` therefore cannot be what *generates* this defect, because **`u` = 0 in two situations that
differ by 6 pp**: when the arrow is *absent* (V23's `cv000`, `g ≡ 0`, measured −0.80%) and when
the arrow is *present but linear* (V24's `dag_pipe_*`, −6.18 / −6.68%). Both have zero
unrepresentable curvature; only the second has a `p(Z₁ | x)` that depends on `x` at all. So the
generator is **whether the child factor carries information**, and `u²` governs how much worse
things get once it also cannot be represented — a modifier, not the mechanism.

**Confirmed within one instrument (2026-09-10).** The comparison above crossed instruments —
`cv000` on the pipeline arms, `dag_pipe_*` on V24's harness draw — so a paired control was
built: `dag_pipe_null`, which is the pipe cell with `delta_xz = 0`. `Z₁` is still present,
observed and in the analysis model; it simply carries no information about the exposure.
Setting the arrow to zero leaves the RNG consumption order untouched, so this cell sees
**byte-identical exposures** to `dag_pipe_dir` (verified: correlation 1.000). Three cells, one
arm set, 200 reps, `n` = 800:

| arrow | `u` | `Z₁` informative? | `dag_Yonly` `bias/SE` (n = 800 → 3200) | rel. bias |
|---|---|---|---|---|
| **absent** (`delta_xz` = 0) | 0 | **no** | **−0.13 → +0.05** ± 0.08 | −1.39% → +0.28% |
| **linear** | 0 | yes | **−0.60 → −0.86** | −7.72% → −5.45% |
| **curved** (`pipe_nl`) | > 0 | yes | +0.45 → +1.26 | +7.14% → +9.77% |

`u` = 0 in the first two rows and they differ by 6 pp, so **`u` cannot be the generator**. What
separates them is whether `p(Z₁ | x)` depends on `x` at all. Under an absent arrow all four
arms agree to within 0.05 pp — as they must, since the child factor then contributes nothing —
and they stay at zero at **both** sample sizes, while the other two cells grow in `bias/SE` the
way an asymptotic bias against a shrinking SE must.

*(One artefact worth naming, because it misleads at the smaller size: every oracle in these
cells sits near −1.2% at `n` = 800 and near +0.4% at `n` = 3200. That is the estimator's own
finite-sample bias, not the imputation's — which the two-level run makes visible and a
single-level run would have hidden inside the arm numbers.)*

Cross-checking in the other direction: `dag_noY` is biased in **all three** cells at both sizes
(−1.42/−2.61 absent, −1.96/−3.68 linear, −1.30/−2.40 curved) — confirming that the `Y`-omission
bias of §1 is about `Y` and is independent of the covariate structure.

**A more flexible child model is monotonically worse.** At `n` = 12800, `pipe_nl` direct:
no child +8.82%, linear +25.75%, quadratic +31.07%, cubic +37.38%, true **−0.21%**; total:
+5.76 / +22.27 / +29.04 / +41.66 / **+0.04%**, coverage 0.000 for every fitted form.

> This is §0b's non-identifiability turned into a *sign*. The usual defence — fit something
> flexible and let the data decide — is **exactly backwards** below an LOD, because there is no
> data there to decide with and a richer basis extrapolates further wrong. It is why the fix
> cannot carry a fitted default, and why "declared form plus sensitivity axis" is not
> conservatism but the only correct interface.

**Where this is not yet tested.** The covariates are fully observed in every cell above
(`mcar_frac` = 0), deliberately — all three factors condition on `Z`, so a missing `Z₁` is a
different question, and leaving it missing would reintroduce V21's covariate-draw ladder as a
confound. `dag_impute_datasets()` refuses missing covariates rather than returning `NA`. The
combined case is open (`ROADMAP.md`).


### Why point estimates, and where that choice runs out

**Added 2026-09-09.** Every verdict in this project compares a point estimate, an empirical SE
and a coverage rate — not whole posteriors. For a Gaussian linear estimand that loses nothing,
and it can be checked exactly without MCMC, because the conjugate posterior for a coefficient
is closed form. Pooling `m` = 20 exact posteriors and comparing to the complete-data posterior:

| cell | location shift | scale ratio | pooled skew / kurtosis | complete-data skew / kurtosis |
|---|---|---|---|---|
| fork | −0.0015 | 1.24 | +0.034 / 2.94 | +0.000 / 3.05 |
| pipe | +0.0028 | 1.20 | +0.035 / 2.95 | +0.000 / 3.05 |
| `pipe_nl` | **+0.0945** | **1.58** | +0.030 / 2.95 | +0.000 / 3.05 |

**The pooled posterior is symmetric and mesokurtic in every cell** — indistinguishable in shape
from the complete-data posterior. Location and scale describe it completely, and those are
exactly what bias, `emp_se` and coverage measure. So for this estimand class, comparing whole
posteriors would add nothing, and MCMC would add nothing but cost: `fit_lm_estimand()` is what
makes 500 reps × 7 cells × 4 arms a five-hour run rather than a months-long one.

**But that argument is specific to a Gaussian linear model, and the pipeline is not one.** Its
real targets are logistic, ordinal via `mo()`, splines and mixed models, where the posterior
for a coefficient is skewed and location-plus-scale is *not* a complete summary. Two
consequences, both unaddressed:

1. **Every validation result in this project is on a Gaussian linear estimand.** Transfer to
   the pipeline's actual model classes is an assumption, not a measurement.
2. **Step 6's posterior-pooling machinery has never been exercised.** It applies the finite-`m`
   variance correction on a support-respecting transform, *gated by a bimodality diagnostic* —
   machinery that only acts when the pooled posterior is non-Gaussian. In 23 tracks that
   condition has never been met, so the gate has never fired and the transform has never been
   tested against a known answer.

**That is a larger gap than any remaining cell in the covariate programme**, and it is the one
place where "compare the whole posterior" would be the right instrument rather than a more
expensive version of the same answer.

### A shortcut that does not work

I hypothesised that if `p(Z | x)` cannot be modelled, *dropping* `Z` from the exposure draw
might beat conditioning on it linearly. **Refuted**, measured by changing only `leftcens`'s
predictor set:

| cell | estimand | `Z1` in the draw | `Z1` out |
|---|---|---|---|
| `pipe_nl` | direct | +19.48% (`bias/SE` 1.53) | **+32.32% (2.95)** |
| `pipe_nl` | total | +16.87% (2.91) | +14.31% (2.45) |
| `fork` | direct | +0.28% (0.03) | **+4.11% (0.56)** |

Dropping it is clearly worse for the direct effect and for the fork, and only marginally better
for the total effect. There is no shortcut: the child factor has to be modelled, or its absence
disclosed.

---

## 7. What would make this stronger

In rough order of value:

- **A magnitude result for clause (A)'s failure**, even asymptotic. The superlinear
  attenuation in the censored fraction is the sharpest open question, because it says the
  penalty for uncongeniality grows faster than the amount of missing data.
- **A test of the additive-decomposition assumption** used in §3's third consequence. V13's
  three steps summed to the total *exactly*, which is suggestive but is one data point at
  one set of parameter values.
- **A derivation of the sign** of the shape effect, which would say whether a covariate-only
  fix is always harmful or only harmful here.
- **Extending the condition to `mo()` and splines**, where "linear in `v`" needs care: a
  spline basis is linear in its coefficients but not in `v`, and clause (A) concerns `v`.
