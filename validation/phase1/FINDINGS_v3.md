# V3 findings — the mixture verdict holds, and is worse than the scaffold showed

**Run:** 2026-08-28 14:16 → 2026-08-29 02:50 JST (**753.8 min**, 3 scenarios × 4 arms ×
7 estimands × 200 reps = 600 tasks, `n_ok` = 200 in all 84 cells, zero task errors, zero
non-finite estimates). Command: `bash validation/phase1/run_v3_bkmr.sh`. Results:
`results/v3_latest.rds`. Closes roadmap item **04** = design-plan **Phase 2b** = **R11**.

**Headline.** Linear congenial imputation does not merely fail to help on a mixture
surface — on the estimands a BKMR analysis actually reports, it is **the worst of the four
arms tested**, worse on average than the LOD/√2 substitution it was built to replace. The
damage concentrates on **curvature**, where it destroys 57% of the effect.

| arm | mean \|paired excess over oracle\| across 21 cells | median | worst cell |
|---|---|---|---|
| `sub_lod2_bkmr` | **6.4%** | 4.6% | −20.5% |
| `cc_bkmr` | 10.5% | 7.5% | −30.6% |
| **`pipeline_bkmr`** | **14.6%** | 11.7% | **−58.4%** |

---

## Read the gate first — it failed, and here is the diagnosis

The pre-registered rule was that an `oracle_bkmr` failure invalidates the run and must be
diagnosed before any other arm is read. It failed: **6 of 21 oracle cells breach**, and the
substantive breach is `curv_X1` at **+15.1% to +16.9%** against a 10% bar, with Monte-Carlo
error of only ~2.3% — not noise.

**The 10% bar was mis-set, by me, on an underpowered calibration.** It was chosen from a
44-replication oracle sweep that put `curv_X1` at 6.61%. Reps 1–44 of this run reproduce
that number **exactly** (6.61%), so the two runs agree and nothing drifted — the 44-rep
estimate simply carried 5.27% Monte-Carlo error. Rolling 44-rep windows across these 200
replications range from **5.2% to 26.1%**; the calibration landed near the bottom of that
range by luck.

| reps used | `curv_X1` oracle rel. bias | its MC error |
|---|---|---|
| 1–44 (the calibration) | 6.61% | ±5.27% |
| 1–100 | 12.83% | ±4.00% |
| 1–200 (this run) | **15.12%** | ±2.43% |

This is the V7 lesson recurring: **a criterion set from a measurement noisier than the
criterion itself is not a criterion.** It is recorded here rather than quietly widened.

### Why the run is still interpretable

The gate exists to detect a *broken harness* — wrong estimand extraction, wrong quantile
convention, a bad GPP approximation, too short a chain. Three pieces of evidence say the
harness is sound and the breach is confined:

1. **`int_X1X2` oracle bias is 0.54–1.12%** across all three scenarios. If the contrast
   construction, the point ordering or the quantile convention were wrong, the pure
   interaction contrast could not come out essentially unbiased.
2. **`curv_X1` is the hardest contrast by construction** — a second difference, so it
   differences away the most signal, and its truth (0.1365) is the smallest in the grid.
   A modest absolute error (~0.021) is therefore a large *relative* one. Its MC error is
   also the largest of the seven, for the same reason.
3. **The MAIN analysis is immune to the floor by design.** Arms see the same datasets in
   each replication, so the headline is the *paired excess over the oracle*, in which
   BKMR's own finite-sample bias cancels. This was decided before the run precisely because
   a floor was known to exist — it just turned out to be larger than measured.

**The oracle floor is a property of BKMR at n = 800 on this surface, not of the pipeline.**
Every pipeline number below is reported as excess over it.

---

## The sharp estimands

`int_X1X2` and `curv_X1` are the tests that matter: their analytic truths contain nothing
but `b_int` and `b_quad`, the two parameters a linear imputation conditional cannot
represent. Paired excess over the oracle, 200 reps, all p < 1e-4 unless noted:

| scenario | `curv_X1` excess | `int_X1X2` excess |
|---|---|---|
| `nd20` (X1 censored 20%) | **−11.7%** | +5.7% |
| `nd40` (X1 censored 40%) | **−56.7%** | +2.6% |
| `nd40_all` (all three censored 40%) | **−58.4%** | **−12.5%** |

**Curvature is destroyed, and the damage scales with censoring.** At 20% non-detects the
pipeline loses 12% of the curvature; at 40% it loses **57%**. Censoring all three exposures
rather than one adds almost nothing beyond that (−56.7% → −58.4%), which locates the damage
in the *focal* exposure's own imputation rather than in the mixture as a whole — `curv_X1`
depends only on X1.

**The interaction is far more robust**: +5.7% and +2.6% where only X1 is censored. It breaks
down (−12.5%) only when X2 is censored too, which is what the estimand's dependence on X1
*and* X2 predicts. Interaction survives because the X-block conditions on the other
exposures linearly and that is enough to carry a product term; curvature in X1 cannot
survive a conditional that is linear in X1's own predictors.

This is exactly the §7.7 mechanism, now measured on real estimands instead of a scaffold
coefficient — and it is **localised**: the failure is curvature, not mixtures in general.

---

## The uncomfortable comparison

**The pipeline is worse than LOD/√2 substitution in 12 of 21 cells**, and much worse on
curvature:

| scenario | `curv_X1` paired excess: `sub_lod2` | `pipeline` |
|---|---|---|
| `nd20` | +4.4% | −11.7% |
| `nd40` | −1.9% *(p = 0.47, indistinguishable from zero)* | −56.7% |
| `nd40_all` | −2.0% *(p = 0.42)* | −58.4% |

On the curvature estimand, **substitution is unbiased relative to the oracle and the
congenial engine is not.** The mechanism is coherent rather than paradoxical: substitution
does not impute at all — it inserts a constant, which distorts the marginal distribution
but leaves the surface's shape to be estimated from the observed rows. The block-FCS X block
*does* impute, from a conditional that is **linear in the predictors**, and so actively
imposes a linear relationship on precisely the rows where curvature would have to show up.
Imputing with the wrong functional form is worse than not imputing.

`cc_bkmr` (complete case) is the worst on `int_X1X2` (−12.2% to −30.6%), as expected: it
discards the low-exposure region where the interaction is identified.

---

## Coverage does not reveal any of this

| scenario | `curv_X1` pipeline rel. bias | coverage | width/SE | FMI |
|---|---|---|---|---|
| `nd20` | +5.2% | 0.980 | 1.12 | 0.084 |
| `nd40` | **−41.5%** | **0.945** | 1.50 | 0.239 |
| `nd40_all` | **−42.2%** | **0.960** | 1.63 | 0.250 |

**A −42% point bias sits behind near-nominal coverage.** The pooled intervals are 1.5–1.6×
wider than the estimates' actual sampling spread, and that inflation is enough to keep the
truth inside them. Anyone validating this path on coverage alone would conclude it was
working.

This is the second time in two tracks that calibration diagnostics have been blind to a real
bias — V8 found the same thing on the inner-iteration shortcut (0.72 pp of bias with
coverage moving < 0.004). **Coverage is not a substitute for a bias check against known
truth.**

---

## What this settles

The claims ledger recorded the mixture verdict as *"mechanism demonstrated, verdict not
manuscript-final"*, resting on a scaffold estimand (`b_logX1`, +19.6%, coverage 0.62). It is
now final, on the reported estimands, with analytic truth:

- **The mixture restriction in the README is correct and should be stated more firmly.**
  The censored-exposure engine must not be used for BKMR-style mixture analyses.
- **It is a curvature failure specifically**, not a blanket mixture failure. Interaction
  estimands survive single-exposure censoring at +2.6% to +5.7%.
- **It is not fixed by preferring the "principled" option.** On curvature the naive
  substitution beats the congenial engine, decisively and reproducibly.
- **The fix remains substantive-model-compatible imputation** — putting the exposure–response
  surface into the imputation model. V3 was scoped as confirmatory and deliberately did not
  attempt it. It is now the sole remaining item for the mixture path.

---

## Caveats (scope of this run)

1. **The gate failed as written.** The harness is argued sound on the three grounds above,
   not on the criterion having been met. A reader who rejects that argument should treat the
   `curv_X1` results as provisional — though the *paired* comparison would still stand,
   since the oracle floor cancels in it.
2. **The oracle floor on `curv_X1` (+15–17%) is unexplained.** It was shown not to be chain
   length (flat from 3000 to 8000 iterations) and not the GPP approximation (50 knots
   matched the full GP). Narrowing the evaluation quantiles reduced *absolute* error on all
   seven estimands but not relative error. It is a real property of BKMR on this surface at
   `n` = 800 and remains uncharacterised beyond that.
3. **`m` = 10, not the pipeline's validated 30.** Chosen because bias is insensitive to `m`
   and cost is linear in it. The `fmi` (0.084–0.250) and interval-width figures are
   therefore on thinner evidence than the bias figures, and the width/SE inflation reported
   above should not be quoted as a precise calibration result.
4. **`overall_q75_q25` is degenerate on this generator** — symmetric quantiles cancel both
   non-linear terms, so it is numerically identical on a mixture and an additive surface. It
   is reported for completeness and tests nothing about the mixture.
5. **One mixture surface**, one `rho`, `n` = 800, MCAR covariates (none missing here), a
   correctly specified BKMR analysis model, and a single `b_int` / `b_quad` pair. The
   *direction* of the curvature failure should generalise; its magnitude is specific to this
   generator.
6. **`iter` = 3000 with the second half kept.** Adequate by the chain-length sweep, but no
   formal convergence diagnostic (R-hat across parallel chains) was run per fit; `bkmrhat`
   would be the tool if that is wanted.
