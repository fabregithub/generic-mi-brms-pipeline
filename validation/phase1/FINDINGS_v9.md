# V9 findings — the V3 diagnosis is confirmed causally, and the fix works

**Run:** 2026-08-29 08:00 → 18:42 JST (**641.5 min**, 2 scenarios × 4 arms × 7 estimands ×
100 reps = 200 tasks, `n_ok` = 100 in all 56 cells, zero task errors, zero non-finite
estimates). Results: `results/v9_latest.rds`. Registered criteria:
`../PLAN_pipeline_validation.md` §8b.

**Headline.** V3 inferred that the shipped X block destroys curvature *because* its
conditional is linear in the predictors. Replacing only that component — drawing the
censored exposure from the correct conditional under the same surface — **closes 89–100%
of the gap**. Across all fourteen cells the mean absolute excess over the oracle falls
from **17.2% to 0.9%**, a factor of ~18.

**And the realistic version matches the oracle version.** Estimating the surface's
coefficients from the data each sweep, rather than knowing them, costs **0.2 percentage
points**. What matters is getting the functional *form* right, not the parameters.

---

## The control: bit-identical to V3

V9 shares V3's seeding rule and canonical scenario indices, so `oracle_bkmr` and
`pipeline_bkmr` should reproduce V3 exactly on the overlapping replications — not merely
agree within Monte-Carlo error. They do, on **all 2,800 shared rows**:

| scenario | arm | identical | max abs diff |
|---|---|---|---|
| `nd40` | `oracle_bkmr` | TRUE | 0 |
| `nd40` | `pipeline_bkmr` | TRUE | 0 |
| `nd40_all` | `oracle_bkmr` | TRUE | 0 |
| `nd40_all` | `pipeline_bkmr` | TRUE | 0 |

The registered control asked only for agreement within MC error. Exact reproduction is
stronger, and it means the SMC arms are being compared against the *same* pipeline
behaviour V3 measured, on the *same* data.

---

## Primary result — `curv_X1`, paired excess over the oracle

| scenario | `pipeline_bkmr` | `smc_oracle_bkmr` | `smc_plugin_bkmr` |
|---|---|---|---|
| `nd40` | **−55.62%** (±3.33) | **+5.90%** (±2.78) | **+6.15%** (±2.98) |
| `nd40_all` | **−59.18%** (±2.83) | **−0.26%** (±2.00) | **−0.10%** (±2.50) |

*(paired on identical datasets; ± is Monte-Carlo error)*

**Read as a recovery:**

| scenario | gap closed by `smc_oracle` | by `smc_plugin` | cost of estimating vs knowing |
|---|---|---|---|
| `nd40` | **89%** | 89% | 0.2 pp |
| `nd40_all` | **100%** | 100% | −0.2 pp |

**Against the registered criterion** (`smc_oracle` within ±10 pp of zero): `nd40_all`
passes outright, with its whole confidence interval inside the bar. `nd40` passes on the
point estimate (+5.90%) but its interval reaches 11.4%, so it clears the bar without
clearing it by a margin. See the caveats — that +5.90% is 2.1σ from zero among fourteen
comparisons, which is not compelling evidence of a real residual.

---

## The effect is not confined to curvature

Excess over the oracle, all seven estimands:

| scenario | estimand | `pipeline` | `smc_oracle` | `smc_plugin` |
|---|---|---|---|---|
| `nd40` | `curv_X1` | −55.62 | +5.90 | +6.15 |
| | `int_X1X2` | +1.52 | −0.96 | −0.26 |
| | `overall_q25_q50` | −14.98 | −1.21 | −1.14 |
| | `overall_q75_q25` | −13.02 | −0.12 | +0.14 |
| | `overall_q75_q50` | −12.14 | +0.36 | +0.71 |
| | `singvar_X1_q50` | −9.83 | −0.41 | −0.03 |
| | `singvar_X1_q75` | −5.22 | −0.75 | −0.19 |
| `nd40_all` | `curv_X1` | −59.18 | −0.26 | −0.10 |
| | `int_X1X2` | −11.41 | +0.97 | +1.17 |
| | `overall_q25_q50` | −0.30 | −0.35 | −0.32 |
| | `overall_q75_q25` | −12.52 | −0.54 | −0.58 |
| | `overall_q75_q50` | −17.95 | −0.62 | −0.70 |
| | `singvar_X1_q50` | −15.85 | −0.63 | −0.69 |
| | `singvar_X1_q75` | −11.55 | −0.42 | −0.28 |

| arm | mean \|excess\| | median | worst cell |
|---|---|---|---|
| `pipeline_bkmr` | **17.22%** | 12.33% | −59.18% |
| `smc_oracle_bkmr` | **0.97%** | 0.58% | +5.90% |
| `smc_plugin_bkmr` | **0.89%** | 0.45% | +6.15% |

V3 found the pipeline biased on nearly every estimand, not only the sharp ones — the
overall mixture effect was off by −12% to −18%. **SMC removes essentially all of it.**
Twelve of the fourteen SMC cells sit within ±1.2 pp of the oracle.

Coverage is unproblematic throughout (mean 0.951 oracle, 0.974 pipeline, 0.981
`smc_oracle`, 0.974 `smc_plugin`). The SMC arms run slightly conservative rather than
anti-conservative.

---

## What this settles, and what it does not

**Settled: the V3 diagnosis was right, and it was right about the mechanism.** The failure
is caused by the *functional form of the imputation conditional*, not by censoring as such,
not by the block-FCS alternation, and not by anything downstream. Changing one component
and nothing else recovers the estimands.

**Settled: knowing the parameters is not the hard part.** `smc_plugin` estimates the
surface's coefficients from the current completed data each sweep and draws them from
their posterior, and it matches `smc_oracle` to 0.2 pp. Roadmap item 07 does not need
oracle knowledge — it needs the right *form*.

**Not settled: how to get the form when the surface is unknown.** Both SMC arms are given
the generator's functional form (`dgp_formula()`). That is the whole point of a mechanism
test, and it is also its limit. A real BKMR analysis does not know the surface — that is
what BKMR is for. This run shows the target is reachable *if the form is available*; it
does not show how to obtain the form. That gap is item 07's actual content, and it is now
much better defined: an SMC-FCS-style scheme would need to condition on a *fitted* surface
rather than a known parametric one.

---

## Caveats (scope of this run)

1. **`smc_oracle` at `nd40` is +5.90% (p = 0.036), not zero.** It is 2.1σ among fourteen
   comparisons, so it is weak evidence at best, and the harder cell (`nd40_all`, all three
   exposures censored) shows −0.26%. If the residual were real, the harder cell should be
   worse rather than better; the two differ by only ~1.8σ, so the most likely reading is
   noise. It is recorded rather than rounded away.
2. **The SMC arms are given the correct functional form.** See above — this is a mechanism
   test, not a proposed implementation.
3. **`smc_*` are harness instruments, not pipeline code.** They do not touch
   `00_censored_exposure.R`. Nothing here changes what the shipped pipeline does, and the
   mixture restriction in the README stands unchanged.
4. **Two scenarios, 100 replications.** `nd20` was dropped as the least informative cell.
   MC error on `curv_X1` excess is ~2–3 pp, which resolves a 57 pp effect easily but leaves
   a few pp of residual indistinguishable from zero.
5. **The oracle floor is still present and still unexplained** (V3 caveat 2). It cancels in
   the paired excess, which is why every number here is reported that way; raw relative
   bias for the SMC arms on `curv_X1` is +18.7% / +14.8%, almost all of which is the floor.
6. **One mixture surface**, `n` = 800, correctly specified BKMR analysis model, `m` = 10.
   The *direction* of the result should generalise; the completeness of the recovery is
   specific to a surface whose form the imputer was handed.
