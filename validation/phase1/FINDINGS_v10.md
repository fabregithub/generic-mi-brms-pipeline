# V10 findings — MAR does not degrade the engine; it is slightly easier than MCAR

**Run:** 2026-08-29 (**87.0 min**, 12 scenarios × 1000 reps = 12,000 tasks, `n_ok` = 1000
in every cell, zero task errors, zero non-finite estimates). Results:
`results/v10_latest.rds`. Registered criteria: `../PLAN_pipeline_validation.md` §8c.

**Why this ran.** Every claim the root README makes about the additive path — V1's bias
≤0.8% and coverage 0.957, V2's ten-scenario sweep, the V4–V8 variance work — rests on
**MCAR** covariates. MCAR is the easy case and the unrealistic one. This is the first test
of the mechanism multiple imputation is actually built for.

**Headline.** The engine is not degraded by MAR. At matched missing fractions it is
**less** biased under MAR than under MCAR — −0.17% vs −1.67% at 40%. Coverage is
0.948–0.967 throughout. **The registered ±1 pp criterion is nonetheless breached**, because
the difference (+1.50 pp) exceeds it — in the *favourable* direction. That is reported as a
breach rather than reinterpreted.

---

## Control: the injector changed only the mechanism

`mar_z40_mcarctl` sets the logit coefficients to zero, which reduces the MAR injector to
MCAR. It must reproduce `mcar_z40`:

| scenario | rel. bias | coverage | width/SE | `b` |
|---|---|---|---|---|
| `mcar_z40` | −1.668% | 0.957 | 1.054 | 0.00101 |
| `mar_z40_mcarctl` | −1.781% | 0.949 | 1.047 | 0.00100 |
| **difference** | **−0.113 pp** (MC error 0.558) | −0.008 | | |

**Control passes.** The new code path changes the mechanism and nothing else.

---

## Main result

| cell | rel. bias | MC error | coverage | width/SE | FMI |
|---|---|---|---|---|---|
| `mcar_z20` | −1.116% | 0.392 | 0.948 | 1.023 | 0.340 |
| `mar_z20` | **−0.437%** | 0.387 | 0.963 | 1.049 | 0.351 |
| `mcar_z40` | −1.668% | 0.394 | 0.957 | 1.054 | 0.376 |
| `mar_z40` | **−0.167%** | 0.402 | 0.967 | 1.055 | 0.391 |

| comparison | Δ bias | bar | Δ coverage | bar | verdict |
|---|---|---|---|---|---|
| 20%: MAR vs MCAR | +0.678 pp | ±1.0 | +0.015 | ±0.03 | **passes** |
| 40%: MAR vs MCAR | **+1.502 pp** | ±1.0 | +0.010 | ±0.03 | **breaches the bias bar** |

The 40% breach is 2.7σ against its own Monte-Carlo error, so it is real rather than noise.
**Its direction is favourable**: MAR sits at −0.17% and MCAR at −1.67%, so the MAR cell is
*closer to unbiased*. The registered bar was two-sided and is recorded as breached; what
follows explains the direction rather than excusing it.

### It is not an artifact of realising less missingness

An obvious confound: if the MAR cells simply lost fewer values, lower bias would follow
trivially. Measured over 200 replications:

| cell | realised Z1 missingness |
|---|---|
| `mcar_z40` | 0.4000 (sd 0.0000) |
| `mar_z40` | 0.4012 (sd 0.0137) |
| `mar_z40_strong` | 0.4017 (sd 0.0108) |
| `mar_z40_mcarctl` | 0.4008 (sd 0.0158) |

The intercept solve works: fractions match to ~0.1 pp. And the MCAR cells themselves give
the sensitivity — bias moves from −1.116% to −1.668% across a 20 pp change in missingness,
i.e. **−0.028 pp of bias per pp of missingness**. A 0.12 pp difference in realised
missingness therefore accounts for **0.003 pp** of the 1.50 pp gap. The confound is ruled
out.

---

## Dose: no degradation as the mechanism strengthens

| cell | rel. bias | coverage | FMI |
|---|---|---|---|
| `mar_z40` | −0.167% | 0.967 | 0.391 |
| `mar_z40_strong` (doubled coefficients) | +0.252% | 0.966 | 0.412 |
| difference | +0.418 pp (MC error 0.574) | −0.001 | |

Doubling the logit coefficients moves nothing. A null at one strength could be a weak
mechanism; a null across a doubling is much harder to dismiss.

---

## Stress cells

| cell | rel. bias | coverage | width/SE | FMI |
|---|---|---|---|---|
| `base` (no covariate missingness) | −0.478% | 0.940 | 0.977 | 0.313 |
| `nl_mcar_z40` | −1.378% | 0.958 | 1.063 | 0.348 |
| `nl_mar_z40` | **+0.097%** | 0.959 | 1.073 | 0.371 |
| `combined` | −0.157% | 0.957 | 1.043 | 0.324 |
| `mar_combined` | −0.410% | 0.959 | 1.040 | 0.365 |
| `nl_mar_combined` | −0.945% | 0.958 | 1.049 | 0.341 |

The non-linear pair repeats the pattern — `nl_mar_z40` at +0.10% against `nl_mcar_z40` at
−1.38%, a +1.48 pp difference, again favourable and again outside a ±1 pp bar. The
combined cells (MAR covariates *plus* 20% missing outcomes, so MID fires as well) sit at
−0.41% and −0.95% with coverage 0.958–0.959: **the hardest cells in this run behave.**

---

## What this settles

- **The additive-path claims survive their first non-MCAR test.** Bias stays within ±1% of
  truth in every MAR cell, coverage within 0.948–0.967.
- **The MCAR-based claims are, if anything, conservative.** Every headline number in the
  root README was measured under the mechanism that turns out to be *harder* for this
  engine. Nothing needs weakening; the scope note can honestly widen from MCAR to MAR.
- **Doubling the mechanism strength changes nothing**, so this is not a knife-edge result.

## Caveats (scope of this run)

1. **The registered ±1 pp criterion is breached at 40%** (+1.50 pp) and in the non-linear
   pair (+1.48 pp). Both breaches are favourable in direction, but the bar was two-sided
   and is recorded as failed rather than rewritten after the fact. The substantive
   conclusion — MAR does not degrade the engine — does not depend on the bar.
2. **Why MAR is *easier* is not established.** The plausible story is that missingness
   driven by `Y`, which the Z block conditions on, makes the missing values more
   predictable than a random subset does — but that was not tested, and a mechanism this
   run cannot explain should not be leaned on. What *is* established is the direction and
   that it is not a missing-fraction artifact.
3. **MAR, not MNAR.** Missingness depends on observed `Y` and `logX2` only. Verified at
   n = 40,000: conditioning on the drivers, the true `Z1`'s coefficient is +0.0015
   (p = 0.91), while its marginal association is p < 1e-150. MNAR is **out of scope by
   decision** — it is a study-design question, not an analysis one.
4. **One MAR form** — logistic in two standardised drivers, both of them covariates the
   imputation model already conditions on. A mechanism driven by something *outside* the
   imputation model would be the adversarial case and was not run.
5. **Additive ERF, scalar estimand `b_logX1`, `pipeline_bartMI` only.** Nothing here speaks
   to the mixture path, which V3/V9 govern.
