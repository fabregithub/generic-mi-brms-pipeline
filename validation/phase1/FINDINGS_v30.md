# V30 — the ordinal attenuation factor is exact; the thing it was bolted to is not

**Question.** V26 measured defects reading ~0.474× as large in `bias/SE` under a binary outcome,
and `THEORY.md` derived that factor as Fisher information. Does the no-free-parameter
generalisation to a cumulative-logit outcome predict the measured attenuation?

**Outcome, in two parts.**

1. **The information derivation is exact** — `SE_gauss/SE_ord` matched in **24 of 24** cell×`n`
   comparisons, every one within **0.001**.
2. **It does not predict the `bias/SE` ratio**, because that ratio is
   `(bias ratio) × (SE ratio)` and the first factor is **not 1**. V26 assumed it was; V30 is the
   first test of that assumption, and it ranges **0.26 to 2.2**.

So **V26's 0.474 is not a transferable constant.** It is the median of a product with one exactly
known factor and one that is not constant at all.

## 1. Design

14 cells × 8 arms × 300 reps × `n` ∈ {3200, 12800}, 20.3 h. Each family cell paired to its
Gaussian twin by `seed_as` on byte-identical exposures, in one invocation with an identical arm
set (V26's arms exactly, for comparability). A ratio is formed only where the **Gaussian twin has
a real defect**, `|bias/SE| > 0.5` — registered up front this time rather than corrected
mid-analysis as in V26.

Cells sweep `K` ∈ {2,3,5,7} balanced, plus — the point of the design — **unbalanced splits at
fixed `K = 3`** (70/20/10 and 85/10/05), which move the prediction *without* moving `K`. That span
(0.173) is wider than the whole balanced `K` sweep (0.110), and it makes a counterintuitive claim:
a 3-category outcome split 85/10/05 carries **less** information than a balanced binary one.

## 2. The derivation: exact, including the correction found before launch

| | `n` = 3200 | `n` = 12800 |
|---|---|---|
| worst `SE_gauss/SE_ord` miss over 12 cells | **0.0009** | **0.0008** |

Every cell, both `n`, same sign convention, no free parameters — including `ord3v` at 0.3350
measured against 0.3358 derived, i.e. **below** its own binary cell (0.4456), as predicted.

Two errors in the naive `√E[w]` form were caught by a complete-data pre-flight *before* the run,
and both had to be fixed to get here:

* the information for `β` is `E[w·x̃²]`, not `E[w]·E[x̃²]` — `w` peaks mid-range while `η` rises in
  `X₁`, so the covariance is negative;
* `polr` **estimates** the `K−1` cutpoints and, unlike the Gaussian twin's nuisance parameters,
  they are not orthogonal to `β`, so the variance is the (1,1) element of the inverted full
  `(β, τ)` information. Flat in `n` (+0.0073/+0.0061/+0.0071 at `n` = 2k/8k/32k), hence asymptotic.

**P4 MET at both `n`.** The unbalanced fan falls and crosses below binary in both roles:

| role | balanced | 70/20/10 | 85/10/05 | its own `K=2` |
|---|---|---|---|---|
| fork (`n` = 12800) | 0.670 | 0.519 | **0.400** | 0.623 |
| pipe (`n` = 12800) | 0.610 | 0.463 | **0.356** | 0.568 |

## 3. What actually fails: bias is not family-invariant

`THEORY.md` states the assumption explicitly — *"If a given imputation defect produces the same
**absolute** bias in the coefficient under both families — plausible, since the structural
coefficient is the same number on both scales"*. It was never tested. Measured
`|bias_fam| / |bias_gauss|`:

| arm | range over cells and `n` |
|---|---|
| `pipeline_bartMI` | 0.24 – 0.33 |
| `pipeline_noYx` | 0.33 – 1.06 |
| `leftcens_prestep` | 0.50 – 1.07 |
| `pipeline_noYz` | 0.80 – 0.91 |
| `pipeline_properZ_ds` | **1.13 – 2.16** |

It varies by arm **and** by cutpoint configuration, and drifts upward with `n` — which is exactly
the signature in the gate misses, every one positive and larger at `n` = 12800. A defect's
absolute size is a property of the defect's interaction with the outcome's likelihood, not a
constant of the coefficient's scale.

**Consequence for the record.** Anywhere "the Gaussian record transfers at ~0.474×" is asserted,
it should instead read: the **standard-error** part of the transfer is exactly computable from
the outcome's Fisher information, and the **bias** part is not known and is not constant. The
conservative-transfer reading survives only where the bias factor happens to be ≤ 1, which is most
arms but emphatically not `pipeline_properZ_ds`, where it reaches 2.16.

## 4. ⚠️ Three of the four gates were mis-specified, and P2 could not have passed

**P1 and P3 tested a composite against a prediction for one of its two factors.** They were
written as though the bias factor were 1 — the very assumption the run disproved. Their "NOT MET"
verdicts say nothing about the derivation; they are a restatement of §3.

**P2 required strict monotonicity in `K`, where the last step is unmeasurable.** The predicted
`K=5 → 7` step is **0.011** against a between-arm SE of **0.08–0.16**. It could not have been met
except by luck, and it failed on a 0.002 dip. The trend it was meant to test is plainly present
(fork 0.623 → 0.670 → 0.677 → 0.675; pipe 0.568 → 0.610 → 0.624 → 0.626). **The verdict is
uninformative; the bar is not moved after the fact, and this is recorded instead.**

This is the **seventh** gate in this programme spoiled by a bar that cannot fail for the reason it
names. The six before it were all denominators that could approach zero, and the standing fix —
a floor on the quantity being compared — worked here: the 0.5 floor was registered up front and
caused no trouble. **This one is a different failure mode and needs its own rule: a gate on a
TREND must have steps larger than the noise on a step.** P4, designed with a 0.173 span against
the same noise, is what a usable gate looks like.

## 4a. An attempt at the bias half, and why it failed (2026-10-15)

**The first-order structure.** The analysis solves an M-estimating equation and the imputation
delivers `x̂ = x + d`. Expanding to first order in `d`:

> `b = I⁻¹ ( E[d·ε] − E[w·x·d'] θ₀ )`,  `I = E[w x x']`, `ε = y − μ(η)`, `w = μ'(η)`

The two numerator terms are exactly `THEORY.md`'s `T₂` (outcome borrowing) and `T₁`
(attenuation), and they carry the outcome family differently: if `I_f ≈ E[w]·I_g` the weights
**cancel** in `T₁` but not in `T₂`. Writing `A` and `B` for those parts on the Gaussian scale,
and using the validated `SE_gauss/SE_fam = r = √E[w]`:

> **M1**  `bias_fam = A/r² + B`  — if `E[d·ε]` is family-invariant
> **M2**  `bias_fam = A/r + B`   — if `E[d·ε]` scales with `sd(ε) = √w`

Both are testable on data already in hand: V30 gives **seven** `r` values per arm per role,
spanning `1/r²` from 1 to 8.9, with the Gaussian anchored at `r = 1`.

**Both are rejected.** Weighted by each point's Monte-Carlo error, with 5 df:

| | arms failing goodness-of-fit | worst `χ²/df` |
|---|---|---|
| M1 (`1/r²`) | 12 of 14 arm×role | **73.6** |
| M2 (`1/r`) | 12 of 14 | **101.1** |

The only passes are `complete_case`, which has no imputation defect to transfer — a null, not a
success. Neither model consistently beats the other.

**Why, and what it costs.** The step that fails is the same one that failed in the SE derivation:
factoring `E[w]` out of `E[w x x']` and `E[w x d']`. There it was fixable, because the full
information matrix is computable from the design. Here it is not, because `E[w x d']` and
`E[d·ε]` depend on the **imputation error distribution `d`**, which is arm-specific and is not
recoverable from summary statistics. A second caveat: the biases are ~25% of the estimand, so a
first-order expansion in `d` may simply be inadequate.

**So the decomposition stands as structure and fails as a predictor.** It says correctly that the
transfer depends on an arm's mix of borrowing and attenuation — which is why the ratio ranges
0.26–2.2 rather than sitting at a constant — but it cannot be collapsed onto any function of the
information ratio alone. **Closing this needs the harness instrumented to record `d` per
imputation**, which is a new run, not a re-analysis. Nothing further should be fitted to these
seven points; with two free parameters and a rejected model, anything that fit would be curve
fitting.

## 5. What is now open

* **The bias transfer has no theory.** The information half is closed exactly; the other half is
  wide open, and it is the half that decides whether a defect measured under a Gaussian outcome
  understates or overstates the same defect under a binary or ordinal one.
* **`pipeline_properZ_ds` roughly doubles its absolute bias** under a non-Gaussian outcome
  (1.13–2.16). That arm is the pre-v1.6.0 scale defect, already fixed in shipped code, so this is
  not a live exposure — but it shows the bias factor can be well above 1, which is the
  anti-conservative direction.
* **V26's logistic table should be re-read** against the exact form: it used the known-threshold
  approximation, and the exact values are 0.4459/0.4488 for `bin_fork`/`bin_pipe`.

## 6. Reproduce

```bash
bash run_v30_ordinal.sh        # 14 cells x 8 arms x 300 reps x n in {3200, 12800}; 20.3 h measured
Rscript analyze_v30.R
Rscript derive_ordinal_information.R   # the derivation and its two pre-launch checks
```

Cells `ord{3,5,7}_{fork,pipe}`, `ord3{u,v}_{fork,pipe}` in `R/robustness.R` (axis `"y_family"`);
results in `results/v30_n{3200,12800}_*`.
