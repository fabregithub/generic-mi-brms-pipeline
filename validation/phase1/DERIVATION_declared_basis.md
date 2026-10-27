# Why a declared basis loses exactly what it omits — a projection argument

**Proposed 2026-10-23 to explain V33. Status after testing (2026-10-26).**

* **P-A (signs and targets): stands.** Omitted terms are rerouted onto what they correlate
  with, at every censoring level.
* **P-B (magnitude): MET on both registered targets, and 7 of 12 cells overall.** Computed in
  `derive_declared_basis.R` — but only once the account is corrected. See §5 P-B and §3a.
* **P-C (scaling in `f`): refuted, and three repairs have now failed** (`Γ` constant, `Γ`
  growing, per-row error propagation). See P-C, P-C′, P-C″.
* **Where the residual lives is now localised.** A harness-faithful replica that differs from
  V33 **only in using OLS on the true basis instead of BKMR** reproduces the asymptotic
  prediction to ≤0.025 in all 12 cells. So the unexplained part — the curvature over-prediction
  and the whole `f`-exponent gap — is a property of **the analysis model's response to
  information loss**, not of the imputation. See §5 P-B, "What this localises".

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

## 3a. Correction (2026-10-26): `Γ` is a FIXED POINT, and the redraw does most of the work

Two things in §3 and §4 are wrong as written, both found by reading the implementation rather
than by arguing:

**(i) The plugin fit is not over the observed region.** `R/smc_impute.R:238` is
`fit <- stats::lm(fm, data = w)` — over the **current completed data**, imputed rows included.
So `h_dec` is not a one-shot omitted-variable regression on the rows where `x₁` is observed. It
is a **fixed point** of

> `h_dec → draw censored x₁ given h_dec → refit h_dec on all rows → …`

That feedback is absent at `f` = 0 and strengthens as `f` grows, because the imputed rows
increasingly determine the very surface used to impute them. (Measured: the iteration converges
by **sweep 2**, so the harness's three sweeps from a cold start at `L − 0.5` reach it. Not a
confound.)

**(ii) §4 ignores that `x₁` itself is redrawn, and that is most of the answer.** §4 reasons as
if the censored rows keep their `x₁` and only their `Y`-information is rerouted. Computing that
literally (`derive_declared_basis.R` STAGE 2, one-shot projection, `x₁` at its true values)
gives, for `smc_dec_noint` at `f` = 0.40, an interaction bias of **−0.756** against a measured
**−0.132** and a full prediction of **−0.169**. The written account overstates by a factor of
about **4.5**. The redraw is not a detail that refines the projection — it **substantially
self-corrects** it, because the drawn `x₁` is chosen to be consistent with `h_dec`, which
partially restores the association the projection destroyed.

The account is only quantitatively right when both corrections are made. §4 below should be
read as the qualitative mechanism (which is what P-A tests and what survives), not as the
calculation.

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
the signs are right in both directions.* Also holds at nd20 and nd40_all (V34): every sign is
preserved as the fraction changes. Not yet a quantitative test.

**P-B (magnitude) — TESTED 2026-10-26. BOTH REGISTERED TARGETS MET; 7 of 12 cells MET.**
`phase1/derive_declared_basis.R`, an independent implementation (own DGP, own censored draw,
own plugin fit, own analysis model — nothing but `MASS` shared with `R/`). Four stages:
`Γ` over the observed region; the one-shot projection §4 literally predicts; the **fixed point
with the `x₁` redraw** (`n` = 200,000, 4 reps); and a **harness-faithful** replica
(`n` = 800, `m` = 10, 3 sweeps, proper `β` draw, cold start, the same closed-form conditional
prior and truncated-quadratic grid draw the harness uses) whose **only** difference from V33 is
an OLS analysis model on the true basis instead of BKMR.

There are no free parameters: §4's `a` and `b` are gone, because the fixed point determines the
completed-data design outright.

**The two registered numbers:**

| target | measured | predicted | σ | |
|---|---|---|---|---|
| curvature inflation, omit `x₁x₂`, `f` = 0.40 | **+0.296** ± 0.036 | **+0.273** ± 0.010 | 0.6 | ✅ MET |
| interaction inflation, omit `x₁²`, `f` = 0.40 | **+0.101** ± 0.020 | **+0.124** ± 0.004 | 1.1 | ✅ MET |

**All twelve cells** (predicted = asymptotic fixed point; MET is |σ| < 2):

| estimand | arm | `f` | measured | predicted | σ | |
|---|---|---|---|---|---|---|
| `int_X1X2` | `lin` | 0.20 | +0.056 | +0.058 | 0.1 | ✅ |
| `int_X1X2` | `lin` | 0.40 | +0.027 | +0.026 | 0.1 | ✅ |
| `int_X1X2` | `noquad` | 0.20 | +0.055 | +0.069 | 0.7 | ✅ |
| `int_X1X2` | `noquad` | 0.40 | +0.101 | +0.124 | 1.1 | ✅ |
| `int_X1X2` | `noint` | 0.20 | −0.049 | −0.080 | 1.7 | ✅ |
| `int_X1X2` | `noint` | 0.40 | −0.132 | −0.169 | 2.1 | ~ marginal |
| `curv_X1` | `noint` | 0.20 | +0.061 | +0.109 | 1.5 | ✅ |
| `curv_X1` | `noint` | 0.40 | +0.296 | +0.273 | 0.6 | ✅ |
| `curv_X1` | `noquad` | 0.40 | −0.397 | −0.457 | 2.7 | ~ marginal |
| `curv_X1` | `noquad` | 0.20 | −0.144 | −0.252 | 3.9 | ❌ |
| `curv_X1` | `lin` | 0.40 | −0.460 | −0.580 | 4.8 | ❌ |
| `curv_X1` | `lin` | 0.20 | −0.102 | −0.328 | 7.2 | ❌ |

**The pattern is not scattered.** Every interaction cell is MET or marginal. Every curvature
cell is MET where the declared basis **retains** `x₁²` (`noint`), and over-predicts where it
**omits** it (`noquad`, `lin`) — always in the same direction: **the measurement shows less
curvature damage than the linear algebra says survives.**

**`f`-scaling, as exponents over `f` = 0.20 → 0.40:**

| arm | curvature measured | curvature predicted | interaction measured | interaction predicted |
|---|---|---|---|---|
| `noint` | 2.28 | 1.33 | 1.44 | 1.07 |
| `noquad` | 1.46 | 0.86 | 0.87 | **0.84** |
| `lin` | 2.17 | 0.83 | **−1.04** | **−1.19** |

The account reproduces the **interaction's** scaling, including the counterintuitive case:
`smc_dec_lin`'s interaction bias **shrinks** as censoring doubles (+0.056 → +0.027), a negative
exponent, and the prediction gives −1.19 against a measured −1.04 without being told to. A
mechanism that reproduces a sign-reversing non-monotonicity it was not fitted to is doing real
work. It does **not** reproduce curvature's scaling in the quadratic-omitting arms.

**What this localises — the first narrowing of the `f²` question.** STAGE 4 matches the
harness's imputation in every respect that was checked against the source (the true closed-form
conditional prior from `.smc_exposure_prior()`, the truncated quadratic grid draw, the
completed-data proper refit, `n` = 800, `m` = 10, 3 sweeps, the `L − 0.5` cold start) and
differs only in the analysis model. It agrees with the asymptotic STAGE 3 to **≤ 0.025 in all
twelve cells**. Finite `n`, properness, `m` and sweep count therefore contribute essentially
nothing, and **the entire residual — the curvature over-prediction and the whole `f`-exponent
gap — is attributable to BKMR's response to degraded curvature information, not to the
imputation.** Both of the project's `f²` findings (V15 and V34) were measured with BKMR as the
analysis model. That is now a specific, cheap, testable common cause rather than a coincidence.

⚠️ **This does not say BKMR is wrong.** The estimands are contrasts of the fitted surface, and
BKMR's own oracle is +15.1% biased on curvature at these settings. The claim is narrower: where
the imputation has destroyed curvature information, a shrinking GP and an OLS fit on the true
basis respond *differently*, and the difference carries the extra power of `f`.

**Internal check.** STAGE 1 recomputes `Γ` independently and reproduces P-C′: `Γ` ≈ 0.40 → 0.41
(flat) for `noint`, 0.51 → 0.29 (decreasing) for `noquad` across `f` = 0.05 … 0.50.

**P-C (scaling in `f`) — TESTED BY V34, REFUTED. The repair offered for it was then also refuted (P-C′).**
Registered: both effects **linear** in the censored fraction, so halving `f` halves the bias
(`smc_dec_noint` +29.6% → ≈ +15%).

Measured over nd20 (f = 0.20) vs nd40 (f = 0.40), 200 reps:
**pooled exponent `p` = 1.85 ± 0.20** — against linear (`p`=1) **4.2σ, refuted**; against
quadratic (`p`=2) 0.7σ, consistent. `smc_dec_noint` came in at **+6.1%**, not +15%.

**The repair I proposed:** that §4 treated `Γ` as constant, whereas `Γ` depends on the region
and hence on `f`, so `f · Γ(f)` would give `f²`. **P-C′ below tested that and it is false** —
`Γ` is flat or decreasing in `f`. Both the original prediction and its repair are refuted; see
P-C′ for what survives.

This is also the **second independent route to `f²`** in the project: V15 found curvature
attenuation quadratic in the censored fraction (`excess% = −323.4·f²`) by a different
mechanism. See `THEORY.md` §1a(vi).

**P-C′ — TESTED 2026-10-26, REFUTED.** Is `Γ(f)` linear in `f`? Computed from the design at
`f` = 0.05 … 0.50 (n = 400,000):

| | `Γ(0.05)` | `Γ(0.20)` | `Γ(0.40)` | `Γ(0.50)` |
|---|---|---|---|---|
| `noint` (omits `x₁x₂`, read on `x₁²`) | 0.400 | 0.403 | 0.405 | 0.407 |
| `noquad` (omits `x₁²`, read on `x₁x₂`) | 0.507 | 0.389 | 0.315 | 0.287 |

**`Γ` is flat in one case and decreasing in the other.** The censored-region extrapolation
residual is likewise flat to slightly decreasing (0.39 → 0.29; 0.56 → 0.46). The product
`f · residual(f)` therefore scales with exponent **0.88** and **0.98** — essentially linear —
against V34's measured **1.85 ± 0.20**. The account under-predicts the exponent by ~4.8σ.

Also corrected here: §3 is right that `Γ` is an expectation over the **observed** region (it is
the plugin fit's omitted-variable bias); §5's earlier wording said "censored". Neither reading
produces `f²`.

**Status of the account after this test.** §4's qualitative content stands — the signs and
targets are right at every censoring level (P-A). Its quantitative content does not: it cannot
produce the measured `f`-scaling. **`f²` is an empirical regularity with no mechanism**, which
is where V15's `excess% = −323.4·f²` has also sat since 2026-09-02. Two unexplained `f²`
findings by different routes is now a question in its own right, not a coincidence to note in
passing.

**P-C″ — the direct measurement, done 2026-10-26. Also negative.** Rather than propose a third
closed form, the chain was measured: imputed `x₁` error per censored row, and its propagation
into the curvature regressor, at `f` = 0.05 … 0.50 (n = 4,000, 6 reps, no BKMR).

Declared `lin` (V34 bias exponent 2.17 ± 0.45):

| `f` | mean \|d\| | `curv_err` = mean(2x₁d + d²) | `f × curv_err` |
|---|---|---|---|
| 0.05 | 0.353 | −0.098 | −0.005 |
| 0.20 | 0.467 | −0.115 | −0.023 |
| 0.40 | 0.567 | −0.073 | −0.029 |
| 0.50 | 0.618 | −0.042 | −0.021 |

**Exponents over `f` 0.20 → 0.40: mean\|d\| 0.28, `curv_err` −0.66, `f × curv_err` 0.34.**
The per-row error grows only weakly with `f`; the induced regressor error is flat to shrinking.
Nothing here produces 1.85. (Control: declared `true` gives `curv_err` ≈ 0 at every `f`, so the
measurement is sound.)

**What the three failures have in common.** Constant `Γ`, growing `Γ`, and direct error
propagation are all arguments about **magnitudes**. A mean shift in `x₁²` does not determine a
fitted coefficient — coefficient bias is set by `Cov(error, Y)` and by **leverage**. All three
were therefore answering the wrong question, which is why each came out flat.

**The remaining candidate — and P-B has since sharpened it.** The hypothesis stated here was
that the extra power lives in the **fitting** step rather than the imputation: censored rows lie
in the low-`x₁` tail, which is where a quadratic is identified, so as `f` grows the corrupted
rows take over more of the curvature information and tail leverage grows faster than linearly.
P-B (above, same day) supports the *location* and changes the *object*. STAGE 4 reproduces the
harness's imputation exactly and still under-predicts the exponent, so the residual is in the
analysis model — but the analysis model is **BKMR's shrinking GP**, not an OLS hat matrix.
Computing OLS leverage is therefore no longer the test. **The test is to re-run V33/V34's two
scenarios with an `lm` analysis arm on the true basis and see whether the exponent falls to
≈ 0.85.** `derive_declared_basis.R` STAGE 4 already predicts that it will.

**RESOLVED 2026-10-27 (V35): the exponent is the ANALYSIS MODEL's.** V15's six-point sweep with
an OLS analysis arm, paired on identical imputed datasets at 100 reps/cell, gives slope
**0.75 [0.69, 0.81]** against BKMR's **2.01 [1.80, 2.21]**, and over all six points
**0.87 [0.75, 0.99]** — containing STAGE 4's registered **0.85**. The prediction made here was
right. `phase1/FINDINGS_v35.md`.

**Standing caution.** Three mechanisms had been proposed and refuted for this one
exponent before it was localised. The qualitative account (§4) has survived every test; the quantitative one has failed
every test. `f²` should be cited as an **empirical regularity** — twice measured, by V15 and
V34, through different mechanisms — until something predicts it.

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
