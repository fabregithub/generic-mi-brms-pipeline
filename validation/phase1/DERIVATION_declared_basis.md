# Why a declared basis loses exactly what it omits — a projection argument

**Proposed 2026-10-23 to explain V33. Derivation only; the test is specified in §5 and has
not been run.**

## 1. What needs explaining

Imputation-induced bias (vs oracle), V33, 200 reps, nd40:

| arm | curv_X1 | int_X1X2 | overall | singvar |
|---|---|---|---|---|
| `smc_dec_true` | +0.036 | +0.004 | 0.000 | −0.001 |
| `smc_dec_over` | +0.029 | +0.004 | 0.000 | 0.000 |
| `smc_dec_noint` (drops `x₁x₂`) | **+0.296** | **−0.132** | +0.010 | −0.052 |
| `smc_dec_noquad` (drops `x₁²`) | **−0.397** | **+0.101** | −0.051 | +0.025 |
| `smc_dec_lin` (drops both) | −0.460 | +0.027 | −0.136 | −0.057 |
| `pipeline_bkmr` (linear conditional) | −0.501 | +0.025 | −0.126 | −0.046 |

Three facts to account for:

1. **Omit a basis function and its own coefficient is attenuated.** Drop `x₁x₂` → interaction
   −13.2%. Drop `x₁²` → curvature −39.7%.
2. **The retained terms are INFLATED, not merely left alone.** Drop `x₁x₂` → curvature
   **+29.6%**. Drop `x₁²` → interaction **+10.1%**. The bias has the *opposite sign* on the
   terms you kept.
3. **The shipped engine and the explicitly-linear SMC arm agree** — −0.501 vs −0.460 on
   curvature, a 1.0σ difference. They are doing the same thing.

## 2. Setup

True surface `h(x) = θ'φ(x)` on a basis `φ = (φ_ret, φ_omit)`. Here
`φ = {x₁, x₂, x₃, z₁, z₂, x₁x₂, x₁²}` with `θ_int = 0.25`, `θ_quad = 0.15`.

The SMC draw for a censored `x₁` is

> `p(x₁ | rest, Y) ∝ p(x₁ | rest) · p(Y | x₁, rest; h_dec) · 1{x₁ ≤ L}`

and the **declared surface `h_dec` enters only through the Y factor**. In `plugin` mode its
coefficients are fitted each sweep from the rows where `x₁` is observed.

## 3. The key step — plugin fitting of a declared basis is an omitted-variable regression

Regressing `Y` on `φ_ret` alone over the observed region gives not `θ_ret` but

> **`θ̃_ret = θ_ret + Γ θ_omit`,  `Γ = E[φ_ret φ_ret']⁻¹ E[φ_ret φ_omit']`**

`Γ` is the regression of the omitted basis functions on the retained ones. So the declared
surface is not "the truth minus a term" — it is **the truth with the omitted term's effect
projected onto the retained basis**.

## 4. What the imputation then does

The censored draw is made consistent with `h_dec`, so on censored rows the imputed `x₁`:

* carries **no** `Y`-information through `φ_omit` — the `(φ_omit, Y)` association there is
  generated only by `p(x₁|rest)`, not by the outcome. Hence **attenuation of `θ̂_omit`**;
* carries **inflated** `Y`-information through `φ_ret`, by exactly `Γθ_omit` — the omitted
  term's signal has been rerouted. Hence **inflation of `θ̂_ret`**.

Writing `f` for the censored fraction and `a, b ∈ (0,1]` for the share of each coefficient's
identification coming from censored rows:

> `θ̂_omit ≈ θ_omit (1 − a f)`  — attenuated
> `θ̂_ret  ≈ θ_ret + b f · Γ θ_omit` — inflated, **in the direction `Γ` points**

This is ordinary omitted-variable bias, displaced one level: not in the analysis model's
covariates, but **in the imputation model's basis**. It is also congeniality in scale
(`THEORY.md` §1a(ii)) made quantitative — clause (A) fails for the specific functional that was
dropped, and the failure has a computable size.

## 5. What this predicts, and how to test it

**P-A (sign and target).** Omitting `φ_j` attenuates `θ̂_j` and inflates the retained terms with
which `φ_j` is correlated over the censored region. ✅ *Consistent with all four V33 rows, and
the signs are right in both directions.* Not yet a quantitative test.

**P-B (magnitude).** The inflation equals `b f Γ θ_omit`. `Γ` is computable from the DGP with
no free parameters — regress `x₁x₂` on `{x₁, x₂, x₃, z₁, z₂, x₁²}` over `x₁ ≤ L`, and
symmetrically for `x₁²`. **This is the derivation to do next**, in the style of
`derive_scale_retention.R`: compute `Γ`, predict +29.6% and +10.1%, compare.

**P-C (scaling in `f`).** Both effects are **linear in the censored fraction** to first order.
The censoring sweep now being set up (`nd20` at 20% vs `nd40` at 40%) tests this directly:
*`smc_dec_noint`'s curvature inflation should roughly halve, +29.6% → ≈ +15%.* A flat or
non-monotone result refutes the projection account.

**P-D (already confirmed, and it is the strongest support).** The shipped engine draws from a
conditional **linear in its predictors** — which is exactly `smc_dec_lin`'s declared basis. The
projection argument therefore says they must agree. Measured: **−0.501 vs −0.460 on curvature,
1.0σ.** Two independent implementations of the same declared basis land on the same number.

## 6. What this does NOT explain

**LOD/√2.** Substitution imputes a constant, so there is no declared surface and no `Γ`. Its
profile is different in kind (curvature −1.6% but interaction +9.9%), and the projection
argument says nothing about it. A constant on censored rows removes all within-censored
variation in `x₁`, so `x₁x₂` and `x₁²` become degenerate there; that is classical measurement
error, needing a separate treatment. **Do not stretch this argument to cover it.**

**Why `smc_dec_over` is marginally the best arm.** Two spurious terms cost nothing measurable
(+0.029 worst case vs +0.036 for the true basis, 0.2σ apart). Under the projection account
extra terms have `θ_omit = 0`, so `Γθ_omit = 0` and no rerouting occurs — they add variance,
not bias. That matches, but with 0.2σ separating the two arms it is consistent rather than
demonstrated.

## 7. If it holds

The practical rule stops being "use SMC" and becomes **"declare generously; an omitted term is
not neutral, it is rerouted onto whatever it correlates with."** That is a sharper and more
useful statement than the covariate-arrow rule it replaces — which V33 refuted for the outcome
surface (P2 NOT MET).
