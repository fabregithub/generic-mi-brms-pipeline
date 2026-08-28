# V7 findings — BART closes R8, and the pilot-`m` question closes with it

**Run:** 2026-08-27 08:18 → 23:14 JST (**~15 h**, five chained phases, zero worker
failures, zero task errors). Closes roadmap item **06** and requirement **R12**.

Three phases, one command (`run_r8_overnight.sh`):

| phase | design | question |
|---|---|---|
| **P1** | 7 scenarios × 4 arms × 300 reps, `m`=30, `ntree`=50 | Does BART close R8? |
| **P2** | 3 non-linear cells × 2 arms × 300 reps, `ntree`=**200** | Is the P1 verdict tuning-fragile? |
| **P3** | `combined`/`nl_combined`, `properBoot`, `m` ∈ {10, 20, 50} × 200 reps | Where does interval width stabilise? |

---

## P1 — the decision grid

| scenario | arm | rel. bias | width/SE | `B` | coverage |
|---|---|---|---|---|---|
| **base** (control) | all four arms | −1.076% | 1.008 | 0.00077 | 0.943 |
| **mcar_z40** | block_fcs | −1.64% | 0.895 | 0.00078 | 0.920 |
| | properBoot | −1.57% | 0.966 | 0.00096 | 0.947 |
| | micePmm | −3.18% | 1.109 | 0.00111 | 0.947 |
| | **bartMI** | **−1.29%** | 1.083 | 0.00101 | **0.957** |
| **nl_mcar_z40** | block_fcs | −1.82% | 0.947 | 0.00072 | 0.940 |
| | properBoot | −1.91% | 0.975 | 0.00078 | 0.937 |
| | micePmm | **−4.70%** | 1.063 | 0.00109 | 0.953 |
| | **bartMI** | **−0.84%** | **1.009** | 0.00090 | 0.947 |
| **missing_y20** | block_fcs | +2.81% | 0.960 | 0.00070 | 0.917 |
| | properBoot | +2.67% | 0.975 | 0.00077 | 0.937 |
| | micePmm | −0.83% | 0.983 | 0.00072 | 0.940 |
| | **bartMI** | **−0.24%** | **0.998** | 0.00084 | 0.937 |
| **nl_missing_y20** | block_fcs | +2.26% | 0.993 | 0.00071 | 0.943 |
| | properBoot | +2.09% | 1.002 | 0.00081 | 0.947 |
| | micePmm | −1.14% | 1.034 | 0.00073 | 0.963 |
| | **bartMI** | **−0.36%** | 1.033 | 0.00084 | 0.950 |
| **combined** | block_fcs | +2.43% | 0.885 | 0.00068 | 0.910 |
| | properBoot | +2.02% | 0.921 | 0.00081 | 0.927 |
| | micePmm | −0.15% | 1.011 | 0.00091 | 0.953 |
| | **bartMI** | **+0.56%** | **1.000** | 0.00096 | 0.950 |
| **nl_combined** | block_fcs | +3.10% | 0.986 | 0.00067 | 0.943 |
| | properBoot | +2.56% | 1.010 | 0.00079 | 0.957 |
| | micePmm | −0.32% | 1.084 | 0.00087 | 0.963 |
| | **bartMI** | **+0.29%** | 1.064 | 0.00091 | 0.963 |

**The control passes.** In `base` all four pipeline arms are **bit-identical** (−1.076%,
1.00818, `B` 0.00077, coverage 0.943) — with nothing for the Z block to impute, every
imputer is a no-op. `n_ok` = 300 in all 35 cells.

### Against the pre-registered criteria

| criterion | result | |
|---|---|---|
| bias ≤1% in `combined` / `nl_combined` | +0.56% / +0.29% | ✅ |
| bias ≤2% in `nl_mcar_z40` (misspecification cell) | −0.84% | ✅ |
| width/SE 0.97–1.03 in all six diagnostic cells | 1.083, 1.033, 1.064 outside | ❌ |
| `base` control identical across arms | bit-identical | ✅ |
| `n_ok` = 300 everywhere | 35/35 cells | ✅ |

**BART fails criterion 3 as written** — three cells 3–8% too wide. But see the next
section: the criterion itself was unattainable, and the "failure" is not statistically
distinguishable from perfect calibration.

### No arm met the band — and only one arm's misses were real

| | band misses | direction | significant? | max \|bias\| | min coverage |
|---|---|---|---|---|---|
| `properBoot` | 2 — 0.921, 0.966 | **too narrow (anti-conservative)** | **yes** (z = −2.10 in `combined`) | 2.67% | 0.927 |
| `bartMI` | 3 — 1.083, 1.033, 1.064 | too wide (conservative) | no (max z = +1.87) | **1.29%** | **0.937** |

`bartMI` has better bias in **all six cells** and a higher coverage minimum. And once the
metric's own noise is accounted for (next section), *its* band misses vanish while
`properBoot`'s does not.

### Criterion 3 was unattainable: the band was tighter than the metric's own noise

`width/SE` divides the mean CI width by the **empirical SD of the estimates**. The standard
error of an SD is `sd/sqrt(2(n-1))`, so at 300 replications the metric carries ~4.1%
Monte-Carlo error. The registered band was ±3%. **A ±3% band on a quantity with ±4% noise
cannot be reliably met by any procedure**, including a perfectly calibrated one.

Testing each deviation against that noise:

| cell | arm | width/SE | MC SE | z vs 1.00 |
|---|---|---|---|---|
| mcar_z40 | bartMI | 1.083 | 0.044 | +1.87 (ns) |
| nl_combined | bartMI | 1.064 | 0.043 | +1.46 (ns) |
| nl_missing_y20 | bartMI | 1.033 | 0.042 | +0.79 (ns) |
| nl_mcar_z40 | bartMI | 1.009 | 0.041 | +0.22 (ns) |
| combined | bartMI | 1.000 | 0.041 | −0.00 (ns) |
| missing_y20 | bartMI | 0.998 | 0.041 | −0.05 (ns) |
| **combined** | **properBoot** | **0.921** | 0.038 | **−2.10** |

**Not one of BART's deviations is significant.** The only statistically real deviation from
1.00 in the whole table is `properBoot` being **too narrow** in `combined` — i.e. the
anti-conservatism R8 exists to eliminate.

So the honest reading of criterion 3 is not "BART overshoots" but "**the criterion could not
discriminate, and where it could, it favoured BART**". What can be claimed is that BART's
interval width is *consistent with calibration, bounded to within roughly ±8%*; not that it
is calibrated. Establishing ±1% would need ~4,800 replications (16×), which is not
obviously worth it.

**Two candidate artefacts were tested and refuted before reaching this conclusion:**
a `t`-vs-`z` bias in the 3.92 divisor (`t` = 1.965–1.974 throughout, moves the ratio ≤0.008),
and light-tailed sampling distributions (`emp95/SD` = 0.95–1.02, excess kurtosis −0.43 to
+0.51 — near-normal). The metric is unbiased; it is simply noisy.

**The design lesson.** `width/SE` is an excellent *mechanism* diagnostic — it is how V4
localised the shortfall to `B` rather than `Ubar` — but a poor *acceptance criterion*,
because its noise floor exceeds the effect size of interest. **Coverage** is the acceptance
criterion: MC SE 0.013 at 300 reps, and it is the quantity that actually matters. Criterion
3 should have been stated on coverage, and the symmetric band was a second, smaller error
on top of that.

### Two refuted explanations, kept for the record

Before landing on the noise-floor explanation above, two candidate artefacts were tested
and both failed:

1. **A `t`-vs-`z` bias in the divisor.** `width/SE` divides by 3.92 (= 2 × 1.96) while
   Rubin intervals use a `t` critical value that grows with `B`, which would penalise
   high-`B` arms like BART. **False:** with `m` = 30 the Rubin df are large, `t` ranges
   1.965–1.974 across every arm and cell, and the correction moves the ratio by ≤0.008
   (`bartMI` `mcar_z40` 1.083 → 1.076).
2. **Light-tailed sampling distributions.** If the estimates were platykurtic, a calibrated
   95% interval would genuinely exceed 3.92·SD. **False:** the empirical central-95% range
   divided by 3.92·SD is 0.95–1.02, and excess kurtosis runs −0.43 to +0.51 — near-normal
   throughout.

So the metric is neither biased nor mis-scaled. It is simply **too noisy to support the
band it was used with**, which is the conclusion above.

### Cost

P1 ran four arms over the same 7 scenarios × 300 reps that V6 ran with three, taking
**8.5 h against 8.2 h** — so the BART arm cost roughly **4% of total runtime**, far less
than either forest arm. The single-target fast path (one fit serving all `m` draws) covers
the `missing_y20` cells; the multi-target cells run `m` chains with no shortcut.

---

## P2 — the verdict is not tuning-fragile

`ntree` 50 → 200 (the `dbarts` and literature default), non-linear cells only:

| cell | bias @50 | bias @200 | width/SE @50 | width/SE @200 |
|---|---|---|---|---|
| nl_mcar_z40 | −0.84% | −0.18% | 1.009 | 1.012 |
| nl_missing_y20 | −0.36% | −0.44% | 1.033 | 1.049 |
| nl_combined | +0.29% | +0.79% | 1.064 | 1.067 |

Bias stays under 0.8% at both settings; width/SE moves by ≤0.016. **`properBoot` is
identical at both** — the control confirming the tuning knob touches only the BART arm.

Two conclusions. The P1 verdict does not depend on an arbitrary tree count — which is what
this phase was for. And the width readings are essentially unchanged at 200 trees
(≤0.016), consistent with the noise-floor reading above: there is no tuning-driven
calibration gain to chase, so `ntree` = 50 stands, at 4× less cost.

---

## P3 — R12 closed: `m` = 30 confirmed

Interval width for the shipped default (`properBoot`), 200 reps per point:

| `m` | `combined` CI width | `nl_combined` CI width | FMI |
|---|---|---|---|
| 10 | 0.2104 | 0.2116 | 0.306 / 0.309 |
| 20 | 0.2072 | 0.2068 | 0.306 / 0.299 |
| **30** | **0.2048** | **0.2053** | 0.304 / 0.298 |
| 50 | 0.2052 | 0.2052 | 0.302 / 0.294 |

**Width stabilises by `m` ≈ 30.** Going from 10 to 30 narrows intervals by ~2.7% (the
`(1 + 1/m)·B` term shrinking); going from 30 to 50 changes them by <0.2%. FMI is flat at
≈0.30 across `m`, as a population quantity should be, and the classical `m ≈ 100 × FMI`
rule independently gives 30.

**Recommended default: `m` = 30** — which is what the pipeline already ships. Design-plan
Phase 5 / requirement R12 is closed, with the empirical curve and the FMI rule agreeing.

---

## What changed in the pipeline (v1.5.0)

`analysis_spec$imputation$z_imputer` selects the Z-block imputer, defaulting to **"bart"**:

| value | behaviour |
|---|---|
| `"bart"` | BART posterior draws. **Default since v1.5.0.** Needs `dbarts` |
| `"forest_boot"` | Bootstrapped forest — the v1.4.0 default |
| `"forest"` | Plain miceRanger — improper; pre-v1.4.0 behaviour |

**Back-compatibility is explicit.** `proper_draw` is still honoured: `TRUE` maps to
`"forest_boot"`, `FALSE` to `"forest"`. So a config written against v1.4.0 keeps the
imputer it was validated with, and only configs specifying *neither* move to the new
default. `z_imputer` wins when both are set. Verified against all seven selector
combinations.

**Graceful degradation:** if `"bart"` is selected and `dbarts` is absent, the pipeline
warns loudly and falls back to `"forest_boot"` rather than failing an upgraded install.
The imputer actually used is logged every run.

This is a **behaviour change** — analyses with missing covariates or outcomes will get
slightly wider intervals and slightly different point estimates than under v1.4.0.
Analyses must cite the pipeline version used.

---

## Caveats (scope of this run)

1. **Interval width is consistent with calibration, not proven calibrated.** The apparent
   3–8% overshoot is within the diagnostic's own ±4% Monte-Carlo error (max z = 1.87) and
   should not be treated as an established defect. It is bounded to roughly ±8%; tightening
   that bound needs ~16× the replications. No tuning is warranted on this evidence.
2. **`micePmm` remains better in some cells** — bias −0.15% in `combined` against BART's
   +0.56%. BART wins on the *worst case*, not on every case.
3. **Linear outcome model.** `Y` is linear in `(logX, Z)` by construction, which is why
   parametric arms do well in the outcome-dominated cells. Untested under a non-linear
   outcome.
4. **MCAR, additive ERF, one estimand** — unchanged limitations from V2/V5/V6.
5. **One non-linear covariate form** (tanh + quadratic + interaction). Other
   non-linearities may penalise differently.
6. **P3 used `properBoot`, not BART.** The `m` curve was measured on the v1.4.0 imputer,
   which was the shipped default when the run was designed. BART's larger `B` could shift
   the required `m` upward slightly; re-measuring on BART is a small open follow-up.
