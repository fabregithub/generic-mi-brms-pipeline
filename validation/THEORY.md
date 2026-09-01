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
superlinear. **The mechanism by which attenuation grows faster than the imputed fraction is
not understood.**

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
which registers `excess% = 323.4·f²` and the conditions that would refute it.

That law is **frank curve-fitting to the two points in §4** — the mechanistically motivated
candidates all fit worse. It is registered anyway, labelled as such, because a
weak-but-falsifiable prediction moves the theory and a post-hoc description does not.

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
