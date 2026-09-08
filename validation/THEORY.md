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

Under (A) and (B) the product of two Gaussians is Gaussian, its mean is linear in `Y`, and
its variance is free of `Y` — exactly what a linear-Gaussian regression represents. Violate
either and the target acquires a non-linear mean, a `Y`-dependent variance, skew, or all
three, and no amount of flexibility in modelling the *mean* recovers it.

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

**MEASURED, 2026-09-03** (500 reps; the derived predictions are in brackets):

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

**Why `u` enters SQUARED.** `u` is by construction orthogonal, in the density-weighted L²
sense, to the span of the linear predictors — it is the residual of that projection. A
first-order expansion of the bias functional pairs the misspecification with the score, and
orthogonality annihilates that term. The leading contribution is therefore **second order**:

> **bias ≈ k · u²**

This is a heuristic argument rather than a proof, and V23
(`PLAN_pipeline_validation.md` §8p) is its test — with the log–log slope's CI excluding 2 as
the condition that would kill it.

**What it may explain beyond itself.** §4's `f²` law for curvature attenuation has been frank
curve-fitting since V15. If the orthogonality argument is right, any misspecification
expressible as an L²-projection residual should enter the bias quadratically — and the censored
fraction's effect may be one. That would turn `f²` from a fitted shape into a consequence.
**Not claimed here**; recorded as the first candidate mechanism this document has had for a
rate rather than a magnitude.

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
