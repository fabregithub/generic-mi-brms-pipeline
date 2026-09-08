# Roadmap — one defect, three tracks

**As of 2026-08-31**, after V11 tested the imputer default against a non-linear outcome. Items 01 and 02 are complete
(partially — see below); item 05 is new and now carries the residual gap.

Five validation tracks are closed. The censored-exposure engine is validated for additive
exposure–response functions. One real defect surfaced — belonging to the pipeline as a
whole, not to the feature that exposed it — and has now been **partially** fixed and
shipped in v1.4.0. The residual gap, plus the two deferred design-plan phases, are the
three tracks that remain.

Companion documents: [`PLAN_pipeline_validation.md`](PLAN_pipeline_validation.md) (the
validation plan, tracks V0–V4) and
[`PLAN_leftcensored_exposure_integration.md`](PLAN_leftcensored_exposure_integration.md)
(the original design plan, Phases 1–5).

---

## Where this sits in the two numbering schemes

The project carries **two** independent numbering schemes, which is easy to trip over.
This roadmap straddles both:

| Roadmap item | Design plan (Phases 1–5) | Validation plan (V0–V4) |
|---|---|---|
| ~~01 — Proper Z-block draw~~ ✅ | **not in either plan** — emerged from V4 | shipped v1.4.0 |
| ~~02 — Re-validate the fix~~ ✅ | not in either plan | done as the V5 `properBoot` arm |
| 05 — Non-linear covariate DGP | not in either plan — emerged from V5 | **new, blocks the imputer choice** |
| 03 — Pilot `m` / FMI | **Phase 5** | V4's original scope, undone |
| ~~04 — BKMR estimands~~ ✅ | **Phase 2b** | Track V3, done 2026-08-29 |
| 07 — SMC imputation | not in either plan — emerged from V3 | **the remaining fix** |

**Design-plan Phase 2 was the decision gate** and was resolved on 2026-08-13: build the
linear congenial component for additive/per-analyte ERFs; mixtures need
substantive-model-compatible imputation. Phases 1, 2, 3 and 4 are complete. Phases **2b**
and **5** were deferred and are now roadmap items 04 and 03.

> **A note on how the defect got through.** The design plan states the governing principle
> in §4: *draw the regression coefficients (and margin parameters) from their posterior per
> imputation, so between-imputation variance is honest.* That standard was applied
> rigorously to the **X block** the plan was building — and never audited on the
> pre-existing **Z block** it inherited. V4 found exactly the failure the principle exists
> to prevent. Worth remembering when the next component is bolted onto an existing path.

---

## What closed

| Track | What it settled | Headline |
|---|---|---|
| **V0** | Steps 5–11 and the Quarto report on the censored path, never previously run | 78 s, 5/5 criteria |
| **V1** | The *shipped* engine against known truth, replacing prototype-only evidence | bias ≤0.8%, coverage 0.957 |
| **V2** | Robustness: all exposures censored, n = 80,000, skew, correlation, MCAR, MID | 9/10 pass, 63.7 CPU·h |
| **V4** | Attribution of V2's under-coverage, one component at a time | 1,500 tasks, 5.25 h |
| **V5** | The candidate fix, validated on shipped code | partial: 66% / 23% of gap closed, 8.05 h |

Findings: [`phase1/FINDINGS_v1.md`](phase1/FINDINGS_v1.md),
[`phase1/FINDINGS_v2.md`](phase1/FINDINGS_v2.md),
[`phase1/FINDINGS_v4.md`](phase1/FINDINGS_v4.md),
[`phase1/FINDINGS_v5.md`](phase1/FINDINGS_v5.md).

---

## The defect, and what remains of it

**`miceRanger` row-level imputation was improper multiple imputation** in Rubin's sense:
it did not draw imputation-model parameters from a posterior, so between-imputation
variance `B` was understated and intervals came out too narrow, worsening with imputation
load. Point estimates were never affected (|bias| ≤ 2.71% across all ten V2 scenarios).

**v1.4.0 fixes most of it.** `proper_draw` now defaults to `TRUE`, bootstrapping the
imputation model's training data once per imputation:

| scenario | before | after | nominal |
|---|---|---|---|
| 40% MCAR covariates | 0.933 | **0.947** | 0.95 |
| 20% missing outcome | 0.927 | 0.937 | 0.95 |
| heavy combined missingness | 0.893 | **0.910** | 0.95 |

**What remains.** Under heavy combined missingness coverage reaches ~0.91, not 0.95 — only
23% of that shortfall closed. The reason is structural: a random forest is *already* bagged,
so an outer bootstrap shifts it only slightly. Measured on identical data, the
bootstrap-on-forest inflated `B` by +31–33% where a parametric posterior draw inflated it
by +41–59%.

Ruled out along the way: raising `m` (no effect, proved on byte-identical data at `m` = 30
vs 100) and MID (disabling it makes calibration *worse* and adds −6 to −7% bias).

**Scope: this was never a censored-exposure defect.** The Z block is the pipeline's general
row-level imputation path, so every analysis with substantial missingness was affected —
and every one benefits from the fix.

---

## What comes next

Ordered by dependency, then value. **Item 05 is now the blocker** — it gates the choice
between the two families of imputation fix. Items 01 and 02 are retained below, marked
done, so the reasoning stays readable.

### 01 — Make the Z block a proper draw · ✅ **PARTIALLY DONE (v1.4.0)**

> **Shipped and validated 2026-08-26.** `proper_draw` now defaults to `TRUE`:
> the Z block bootstraps its training data once per imputation. Closes **66%** of the
> shortfall at 40% covariate missingness (coverage 0.933 → 0.947), improves every case,
> never over-corrects, best bias of any variant. **But only 23% under heavy combined
> missingness** — coverage 0.893 → 0.910, still short of nominal. A forest is already
> bagged, so an outer bootstrap injects less parameter uncertainty than a parametric
> posterior draw. See [`phase1/FINDINGS_v5.md`](phase1/FINDINGS_v5.md). The residual gap
> is now item **05** below.

<details>
<summary>Original scope (for the record)</summary>


The only known defect in the pipeline, sitting in the general imputation path — so the fix
improves every analysis with missing covariates or outcomes. Candidate approaches,
cheapest first:

- `mice`-style proper draws (`norm` / `logreg`) for the imputation models
- a bootstrapped random forest, preserving miceRanger's flexibility while restoring
  parameter uncertainty
- miceRanger's own settings, if any configuration delivers a genuinely proper draw

**Cost is part of the choice.** V2 measured the Z block at ~90 imputation fits per
replication, dominating pipeline runtime. A candidate that calibrates correctly but
triples the cost is not automatically better than a cheaper, slightly less clean one — that
is a judgement call, not a purely statistical one.

</details>

### 02 — Re-validate the fix · ✅ **DONE (V5)**

> Done as part of 01: the `pipeline_properBoot` arm exercised the shipped fix across
> 5 scenarios × 300 reps on the same seeds as V4, giving a paired before/after.
> The harness worked exactly as designed — no new code was needed.

<details>
<summary>Original scope (for the record)</summary>


No new harness needed. The V4 `pipeline_properZ` arm already **is** the "swap the imputer
and re-measure" scaffold: point a candidate at it and width/SE, coverage and the
`Ubar`/`B` split come out directly, on the same seeds as every prior run.

**Target:** width/SE ≈ 1.00 across `mcar_z40`, `missing_y20` and `combined` without the
overshoot the validation instrument showed, and bias inside 3%.

</details>

### 05 — Extend the DGP with non-linear covariates · ✅ **DONE (V6)**

> **Resolved 2026-08-27.** 7 scenarios × 3 arms × 300 reps, 8.2 h, zero failures.
> `mice pmm`'s bias breaches the 3% bar exactly where the non-linear design punishes it
> (−3.18% → **−4.70%**, a penalty **4–8× larger** than either forest arm's), so
> **`properBoot` stays the default**: bounded bias everywhere (max 2.56%), calibration
> nearest nominal, coverage ≥0.927 in every cell. But `pmm` is *better* in four of six
> cells — it wins wherever imputation is outcome-dominated (`Y` is genuinely linear) and
> loses only where covariate imputation is heavy and non-linear. Neither is uniformly
> best. The prediction that `properBoot` would still under-cover **did not hold**.
> Full reading: [`phase1/FINDINGS_v6.md`](phase1/FINDINGS_v6.md).
>
> **What remains open:** a component that is flexible *and* properly dispersed. See the
> BART hypothesis in [`INTEGRATION_SUMMARY.md`](INTEGRATION_SUMMARY.md) §2.

<details>
<summary>Original scope and pre-registered criteria (for the record)</summary>


The residual gap from 01 raises a question the current simulation **cannot answer**: every
covariate in it is linear and Gaussian, so a parametric imputation model is correctly
specified and is never penalised for misspecification. That is exactly why we did not
switch to `mice`'s `pmm`/`norm` despite its better calibration (0.973 vs 0.878 in the worst
cell) — and why the choice cannot be made on current evidence.

Add non-linear and interacting covariates to `R/dgp.R`, then re-run with a `pmm` arm
alongside `properBoot`. Only then does the comparison mean something: the forest keeps its
flexibility advantage, and the parametric draw has to earn its calibration advantage
against a model it can get wrong.

**Until this is done, `properBoot` is the defensible default** — it improves calibration
without giving up flexibility, and carries the best bias on record.

**Built 2026-08-26.** `simulate_complete(z_form = "nonlinear")` makes `Z1` a saturating,
curved, interacting function of the *non-focal* exposures and `Z2`. Verified: a linear
imputation model leaves **0.387 of explainable R²** unclaimed under `nonlinear` against
**0.001** under `linear`, while the oracle recovers `b_logX1` = 0.392 in *both* — so the
imputation problem got harder and the estimand did not move. Scenarios `nl_mcar_z40`,
`nl_missing_y20`, `nl_combined`; `mice` arm in
[`phase1/R/mice_impute.R`](phase1/R/mice_impute.R).

**Decision criteria, registered before the run.** Judged on the *non-linear* cells, where
both families can lose something:

| outcome | reading | action |
|---|---|---|
| `micePmm` keeps its calibration edge **and** bias stays within 3% | parametric proper draws win; the forest's flexibility was not worth its under-dispersion | adopt `mice` for the Z block |
| `micePmm` bias exceeds 3%, or materially exceeds `properBoot`'s | misspecification cost is real — the calibration edge was an artefact of the old linear DGP | keep `properBoot`; question closed |
| `properBoot` still under-covers **and** `micePmm` is biased | neither family is adequate alone | a flexible model *with* proper draws is required — a genuinely new component, not a settings change |

The third outcome is a real possibility and the most consequential, so it is named in
advance rather than discovered as a disappointment.

</details>

### 06 — A flexible *and* properly dispersed imputer · ✅ **DONE (V7, v1.5.0)**

> **Resolved 2026-08-28.** BART adopted as the `z_imputer` default. Best bias of any
> variant in **all six** diagnostic cells (max 1.29% against `forest_boot`'s 2.67%),
> coverage 0.937–0.963, robust at both 50 and 200 trees, and ~4% of pipeline runtime.
> R8 closes. Residual: intervals 3–8% **wide** in three cells — conservative, and now
> tracked as its own (low-priority) item in
> [`INTEGRATION_SUMMARY.md`](INTEGRATION_SUMMARY.md) §2.
> Full reading: [`phase1/FINDINGS_v7.md`](phase1/FINDINGS_v7.md).

<details>
<summary>Original scope and criteria (for the record)</summary>


V6 closed the choice between the two existing families but not the underlying problem:
the forest under-disperses (bagging damps the bootstrap), and the parametric draw is
misspecified under non-linear covariates. The better half of each would be an imputer that
is non-parametric **and** genuinely Bayesian, so posterior draws are proper by construction
rather than approximated by resampling.

**Leading candidate: BART** (Bayesian additive regression trees) — `dbarts` / `bartMachine`;
`sequential BART` (Xu, Daniels & Winterstein) was designed for multiple imputation
specifically. **This is a hypothesis, not a result.** The V4/V6 harness measures any
candidate directly, so the test is cheap once the arm is written.

**Built 2026-08-27** (`dbarts` 0.9.34). `phase1/R/bart_impute.R`, arm `pipeline_bartMI`.
`dbarts::bart()` returns posterior draws of the conditional mean (`yhat.test`) and matching
residual-SD draws (`sigma`), so imputation *j* is
`yhat.test[j, ] + rnorm(n, 0, sigma[j])` — a posterior *predictive* draw, proper by
construction. Binaries fit as factors (probit form). One fit serves all `m` draws when
there is a single target with complete predictors (`missing_y20`); otherwise `m` FCS chains
run, with no shared-fit approximation. A failed fit is fatal rather than falling back to
something improper.

**Decision criteria, registered before the run.** BART must beat **both** incumbents, not
just the shipped reference — otherwise it is not the resolution:

| test | target | incumbent to beat |
|---|---|---|
| bias, `combined` / `nl_combined` | ≤1% | `properBoot` +2.02% / +2.56% |
| bias, `nl_mcar_z40` (misspecification cell) | ≤2% | `micePmm` −4.70% |
| width/SE, all six diagnostic cells | 0.97–1.03 | `properBoot` 0.92–1.01 |
| `base` negative control | identical to other arms | — |
| `n_ok` | 300 in every cell | — |

**If any test fails, R8 stays ⚠️ partial and the gap is documented, not papered over.**
Early signal from the 8-rep smoke: `B` is largest for `bartMI` in all three cells tested
(0.00113 / 0.00110 / 0.00094 against `properBoot`'s 0.00085 / 0.00081 / 0.00076) and bias
in `nl_mcar_z40` is −0.73% against `micePmm`'s −4.73% — the right direction on both axes,
but 8 reps is noise. The cell to watch is `nl_combined`, where three targets force the
`m`-chain path and the smoke read −3.20%.

</details>

### 03 — Pilot `m` and document FMI · ✅ **DONE (V7 P3)**

> **Closed 2026-08-28.** Width stabilises by **`m` ≈ 30**: going 10 → 30 narrows intervals
> ~2.7%, 30 → 50 changes them <0.2%. FMI flat at ≈0.30, and `m ≈ 100 × FMI` independently
> gives 30. The shipped default was already right. Design-plan **Phase 5** closes.

<details>
<summary>Original scope (for the record)</summary>


Most of the way there already: FMI is now reported per arm and lands between **0.26 and
0.44** across scenarios. What remains is turning that into a recommended default `m` for
the censored path and cross-checking it against the auto-increment loop's stopping rule.

**Do this after 01** — a proper draw changes `B`, which changes FMI.

</details>

### 04 — BKMR estimands · ✅ **DONE (V3, 2026-08-29)** — design-plan Phase 2b

Resolved, negatively and more strongly than expected. On the estimands a BKMR analysis
reports, measured against analytic truth at 200 replications per cell, the shipped
censored-exposure engine is the **worst of four arms** — mean 14.6% absolute paired excess
over the oracle, against 6.4% for LOD/√2 substitution.

It **destroys 57% of the curvature** at 40% non-detects. Interaction estimands survive
censoring of a single exposure and break only when a second is censored too. The failure is
therefore *curvature*, not mixtures in general.

Evidence: [`phase1/FINDINGS_v3.md`](phase1/FINDINGS_v3.md). The mixture restriction in the
root README has been strengthened accordingly.

---

### 07 — Substantive-model-compatible imputation · *the mixture-path fix*

V3 established what is broken and why: the X block draws the censored exposure from a
conditional **linear in (Y, other X, Z)**, which cannot represent a surface with curvature
in the focal exposure, and so imposes linearity on exactly the rows where curvature would
appear. Imputing with the wrong functional form is worse than not imputing — LOD/√2
substitution beats the congenial engine on curvature.

> **V14 tested a shippable method and it mostly works.** `.ce_smc_x_grid()` — a grid draw
> given only the analysis *formula* — removes **83%** of the non-linear-outcome penalty and
> is the best-calibrated arm measured (width/SE 1.004, coverage 0.950). Two conditions
> before building: it needs a **missing-`Y` path** (the control failed where multiple
> censored exposures meet missing outcomes), and 83% is not 100% — estimating the outcome
> model costs 0.5–1.4 pp. The **mixture** half of the claim is still untested and needs the
> BKMR harness (~10 h). See [`phase1/FINDINGS_v14.md`](phase1/FINDINGS_v14.md).

> **V13 widened this item.** The `leftcens` linear conditional is not only the mixture
> defect — it accounts for **3.1–3.4 pp of the ≈3.9 pp non-linear-outcome penalty** too
> (V13). Fixing it therefore helps *additive* analyses with missing covariates, not just
> mixture ones, which materially raises the value of this track.

**V9 has now validated that target empirically**, which removes most of the risk from this
track. Drawing the censored exposure from the correct conditional under the same surface
closes **89–100%** of the curvature gap and cuts the mean absolute excess over the oracle
from 17.2% to 0.9%. Two results narrow the work sharply:

- **Oracle parameters are not needed.** The plug-in arm estimates the surface's
  coefficients from the current completed data each sweep, draws them from their posterior,
  and matches the oracle arm to **0.2 pp**. What matters is the functional *form*.
- **Obtaining the form is the whole remaining problem.** Both V9 arms were *handed* the
  generator's formula. A real BKMR analysis does not know the surface — that is what BKMR
  is for.

So the work is no longer "does SMC help?" but "how is the surface represented in the
imputation conditional when it is being estimated non-parametrically?" A plausible route is
an SMC-FCS-style step conditioning on the *fitted* BKMR surface at the current iteration,
which would make the imputation and the analysis model genuinely congenial. V9's sampler
(`phase1/R/smc_impute.R`) already handles a one-dimensional truncated draw against an
arbitrary mean function, so the machinery for the draw itself exists.

For a **known parametric** surface the fix is essentially done and demonstrated.

**Until it exists, the mixture path is documented as unusable rather than approximate.**
That is now stated in the root README, the claims ledger, and `FINDINGS_v3.md`.

This is a build track, not a validation track. When it exists, V3's harness re-runs against
it unchanged — the estimands, analytic truth, arms and criteria are all in place, and a new
arm is one function.

---

### 08 — A shape-aware Z-block draw · ❌ **CLOSED AS SCOPED (V13, 2026-09-01)**

**The defect.** Every row-level imputer in the pipeline draws a covariate as
*(fitted conditional mean) + homoscedastic Gaussian noise*. Under an outcome linear in the
covariates that is exactly right: `p(Z1 | Y, X)` is Gaussian with constant SD (0.894) and
zero skew at every `Y`. Under a **non-linear** outcome the same conditional becomes
heteroscedastic and skewed — SD running 0.626 → 1.216 across `Y`, skew reaching −1.338 —
and the Gaussian draw is then wrong in a way no amount of flexibility in the *mean* model
can repair.

**The evidence that it is real.** V11 measured a ≈3 pp bias penalty under a non-linear
outcome for **all four** imputers — BART, `mice pmm`, forest, bootstrapped forest — with a
spread of only 0.44 pp, while the oracle arm (which imputes nothing) moved 0.01 pp. Arms
that differ enormously in how they model the conditional mean degrade almost identically,
which is what a shared defect in the *noise* looks like and not what a mean-modelling
problem looks like.

**Why this is not item 07.** Item 07 is the **X** block drawing the exposure from a
conditional with the wrong functional *form* — a mixture-path problem, and research-scale.
This is the **Z** block drawing a covariate with the right mean and the wrong *shape*. They
are different defects; neither fix closes the other.

| | 07 — X-block form | **08 — Z-block shape** |
|---|---|---|
| affects | mixture / BKMR analyses | **any analysis with missing covariates and a non-linear outcome** |
| measured cost | 57% of the curvature destroyed | ≈3 pp of bias, every imputer |
| the fix | research: obtain the surface when it is estimated non-parametrically | bounded: sample the conditional in `run_row_level_imputation_bart()` instead of adding Gaussian noise |

> **CLOSED. Do not build this.** [V13](phase1/FINDINGS_v13.md) separated the Z draw from
> the exposure draw and found that a Z-only fix **moves bias away from truth** in the
> configuration the pipeline actually ships. The shape correction that V12 identified — and
> that this item was written around — is *actively harmful* on its own.

**What V13 measured.** Decomposing the ≈3.9 pp penalty, paired, with the three steps summing
exactly to the total:

| step | `ynl_mcar_z40` | `ynl_mar_z40` |
|---|---|---|
| Z conditional **mean** | +1.85 pp | +2.14 pp |
| Z draw **shape** | **−1.43 pp** | **−1.32 pp** |
| **Exposure draw** | **+3.41 pp** | **+3.14 pp** |
| total | +3.83 pp | +3.95 pp |

**Why the shape fix hurts.** Its effect is a consistent negative shift in signed bias
(−1.3 to −1.9 pp) whichever exposure draw is used. With the *exact* exposure draw the bias
sits positive (+1.69%), so that shift helps — which is all V12 could see. With the *shipped*
exposure draw the bias is already negative (−2.19%), so the same shift carries it further
away, to −3.62%. **The Gaussian Z draw and the linear exposure draw were partly cancelling**;
removing one error exposes the other.

**A Z-only fix recovers 11–21%** of the gap, and only because the mean correction outweighs
the shape correction working against it.

**Where this leaves the work.** The exposure draw is the priority at 3.1–3.4 pp of ~3.9 —
and that is the *same* `leftcens` linear conditional that V3/V9 found destroys mixture
curvature. **The two defects share a cause**, which is the opposite of what V12's framing
suggested. Item 07 therefore absorbs this work rather than competing with it.

**Superseded by item 07.** The cheap-and-wide argument for doing this first does not
survive V13: the cheap part (the Z draw) is the part that does not work alone.

**The unattributed half is the obvious next question.** An arm with the exact Z draw but
the *shipped* exposure draw would separate the mean and X-block contributions — a small
addition to the existing harness, and worth doing before building anything, since it says
whether a Z-only fix is enough.

**Sketch of the fix.** `run_row_level_imputation_bart()` currently does
`yhat + rnorm(n_mis, 0, sigma)`. A shape-aware version would sample from the conditional
implied by the fitted outcome model rather than assume it Gaussian — the machinery exists
in `phase1/R/smc_impute.R` (`.smc_draw_z`) as a validated one-dimensional grid sampler. The
open design question is how to obtain that conditional inside the pipeline, where the
outcome model is a `brms` fit rather than a known formula.


---

## Dependencies

```
DONE (v1.4.0)
  [01 Proper Z-block draw] ──▶ [02 Re-validate] ──▶ partial: 66% / 23% of gap closed

DONE
  [05 Non-linear covariate DGP] ──▶ [06 BART imputer] ──▶ R8 CLOSED (v1.5.0)
  [03 Pilot m / FMI] ──▶ R12 CLOSED: m = 30 confirmed

DONE
  [04 BKMR estimands] ──▶ mixture path CONFIRMED UNUSABLE (V3, 2026-08-29)
  [V9 mechanism test]  ──▶ cause CONFIRMED: the linear conditional (2026-08-29)
                           SMC closes 89-100% of the gap; form matters, not params

REMAINING
  [08 Shape-aware Z-block draw]                 CLOSED as scoped (V13): harmful alone
                                                the exposure draw dominates -> see 07
  [07 Substantive-model-compatible imputation]  research-scale; mixture path only
                                                V3's harness re-runs against it unchanged

FOLLOWS
  [03 Pilot m / FMI]                                    (B has changed, so FMI has too)

INDEPENDENT
  [04 BKMR estimands]                                   (install bkmr first)
```

**One build track remains: item 07 (substantive-model-compatible imputation).** V13 closed item 08 as scoped and showed the exposure draw — item 07's territory — accounts for most of the non-linear-outcome penalty as well as the whole mixture failure. The two defects turned out to share a cause. Everything else is loose ends plus the low-priority interval-overshoot refinement.

**V15 does not change that.** It measured how fast the item-07 defect grows with censoring
(quadratically) rather than testing any fix, so the build queue is untouched — but it does
raise the value of finishing 07: at 40% non-detects the defect costs 56% of the curvature,
and at 60% it reverses its sign.

---

## Loose ends

| Item | State | Why it matters |
|---|---|---|
| `resume` kills the V4 runner | Unexplained | Reproducibly killed the parent after one chunk; a fresh checkpoint runs fine. Now off by default. Partial results remain recoverable straight from the checkpoint, so it costs nothing operationally — but an unexplained crash is worth understanding. |
| ~~Translated READMEs~~ | ✅ **closed 2026-09-01** | Deleted in `1c01b2e` rather than re-translated: they predated the censored-exposure feature and mentioned none of it. Readers are pointed at the English vignettes from the root README. |
| ~~MAR covariate missingness~~ | ✅ **done (V10, 2026-08-29)** | MAR does not degrade the engine — it is slightly *easier* than MCAR (−0.17% vs −1.67% at 40%), across a doubling of mechanism strength and in the non-linear and missing-outcome cells. MNAR is **out of scope by decision** — a study-design question, not an analysis one. |
| ~~The root claim was never isolated~~ | ✅ **closed 2026-09-03 (V16)** | The condition every track since V1 refines — that an imputation of `v` must draw from `p(v \| Y, rest)` — has never been tested on its own: every arm in fifteen tracks had `Y` in the imputation model. Phase 1's H1 is the cited evidence, but its no-`Y` arm also drops the `Z` covariates *and* is a different estimator, and two cells in the same document have the no-`Y` arm **beating** the `Y`-aware one. **Measured: −13.1 pp** when `Y` is omitted from both blocks of the shipped engine, isolated inside one estimator for the first time (`phase1/FINDINGS_v16.md`). The registered *leak* prediction failed — the blocks are additive in that structure. Original design: **−14.9 pp** for the realistic case (omitted everywhere), **zero** for the Z block where `Z ⊥ X`, and a **+1.4 pp leak** — because the blocks alternate, `Y` reaches one through the other, and the mixed configuration is predicted to be *worse* than omitting it from both. `phase1/run_v16_noy.sh`, ~1.3 h. |
| **`use_as_auxiliary` is half-honoured in the censored-exposure path** | ⚠️ **measured 2026-09-03 (V17b); fix still unmade** | `docs/variable-dictionary.md` documents `FALSE / FALSE / TRUE` as "use as imputation predictor only". The Z block honours it (`make_row_level_imputation_spec()` selects on `use_in_model \| use_as_auxiliary \| impute_target`); the X block does not (`00_censored_exposure.R`'s `auto_preds` reads `use_in_model` only). So an auxiliary helps impute covariates and is silently dropped from the **censored exposure draw** — the one draw the strategy exists for. Verified directly. Compounding it, **no validation track has ever exercised `auto_preds`**: every arm since V1 passes an explicit `ce$predictors`. Measured: **−0.30 ± 0.10 pp** (MCAR) and **+0.60 ± 0.10 pp** (MAR) — real, under 1 pp, and opposite in sign between mechanisms. `auto_preds` also ran for the first time in any track, without error. |
| ~~No DGP where a covariate causes an exposure~~ | ✅ **closed 2026-09-03 (V17)** | Four roles built and run (fork, pipe, collider, mixed), MCAR and MAR. **The covariate's role decides the covariate block**: the `Y`-omission penalty is +1.2 pp when `Z` is inert and +25 to +28 pp when `Z` is adjusted for *and* on an X–Y path. See `phase1/FINDINGS_v17.md`. Two follow-ups opened below. |
| **Shipped default is biased under a confounder or mediator** | ⚠️ **documented and attributed 2026-09-07; no fix yet** | `pipeline_bartMI` with `Y` in both blocks — the recommended configuration — carries +4.96% (fork), +5.39% (pipe), +8.52% and +10.40% (both under MAR), against ≤2.3% in every precision cell. The oracle is unbiased in all four, so it is the imputation's bias, not the design's. **Coverage does not detect it** (0.93–0.97, width/SE 0.98–1.13) because at `n` = 800 the bias is only 0.35–0.67 empirical SE. V18 measured how it scales: **n^−1/3**, so it shrinks with data but never washes out (2.0–3.6% still present at n = 12,800). Now documented for users in `docs/covariate-roles.md`. **V20 attributed it: the covariate draw carries 85–104% of it**, the exposure draw 0.1–11.7%. No fix yet — see the two rows below. |
| ~~Does the concealment break at larger `n`?~~ | ✅ **closed 2026-09-03 (V18)** | **Measured over five `n` levels to 12,800: the bias decays as n^−0.335 (CI [−0.378, −0.291])**, rejecting both the constant-bias account (slope 0) and the finite-sample one (−0.5). **My projection was wrong** — coverage at n = 12,800 is 0.867–0.923, not 0.32–0.73. But the defect is permanent: bias/SE still grows, as n^0.17, and 2–10% bias is present at every `n` tested. `phase1/FINDINGS_v18.md`. |
| ~~Does the bias exponent track the Z-block imputer?~~ | ✅ **closed 2026-09-07 (V19) — refuted** | V18's n^−1/3 is the classical nonparametric convergence rate for a Lipschitz function in 1-D, and the Z block is BART. In the `fork`/`pipe` cells the true `Z₁` conditional is **linear-Gaussian**, so a parametric imputer is correctly specified there and should converge at −1/2 — while a misspecified one should give 0. **Refuted**: a second nonparametric imputer (`properBoot`) shows no decay at all (−0.029 [−0.087, +0.028]) and both *correctly specified* parametric arms sit at ≈0, not −1/2. Family predicts nothing — **the split is BART versus everything else**. `THEORY.md` still has no mechanism for a rate. What the run established instead: **BART is the only Z-block imputer whose bias decays with `n`**. `phase1/FINDINGS_v19.md`. |
| **`pipeline_properBoot` is 17× the cost of every other Z-block arm** | Open, small — measured 2026-09-04 | At n = 800, 22 reps, one cell, 22 workers: `bartMI` 9 s, `micePmm` 9 s, `properZ` 5 s, **`properBoot` 156 s**. It is ~80% of a four-arm run's cost and is what put V19's first design at ~64 h. Not a correctness issue, but it silently prices `properBoot` out of any multi-arm sweep, and the cause is unknown — the loader caches per variant (measured: ≤0.01 s after the first load), so it is real computation in `run_row_level_imputation_proper()`. Worth an hour with a profiler before the next track that wants it. |
| **MAR conjecture: does a flexible Z draw recover `Y` from proxies?** | Open, small | V17 predicted MAR-on-`Y` would add +4 to +6 pp on-path and measured +0.7 / +0.7 / −0.6 / +2.3 — refuted. The conjecture is that BART recovers most `Y`-relevant information from the exposures and the other covariate. Discriminating run: the same cells with a `mice pmm` Z block, which should show a larger MAR term if the conjecture holds. |
| **`forest_boot` carries −23% asymptotic bias under a confounder** | ✅ **documented 2026-09-07** | Measured over five `n` levels: bias −22% to −29% that **never decays** (slope −0.029), coverage falling to **0.000** by n = 6,400. Reachable two ways today: a v1.4.0-era config setting `proper_draw = TRUE`, and — more importantly — `z_imputer = "bart"` when `dbarts` is not installed, where `run_row_level_imputation()` falls back to `forest_boot`. That fallback warns with `immediate. = TRUE` and logs the imputer, and coverage catches this failure unlike every other one on record, so **do not change the fallback** — completing with a warning beats aborting. Now said in `docs/covariate-roles.md` (with the `grep "Z-block imputer"` check), in `docs/censored-exposures.md`'s corrected v1.5.0 note, and in the root README's cautions. |
| ~~Attribute the bias to a block~~ | ✅ **closed 2026-09-07 (V20) — it is the covariate draw** | Fixing the Z block removes **85–104%** of the bias; fixing the exposure draw **0.1–11.7%**; interaction ≤0.40 pp. Every registered gate passed. `phase1/FINDINGS_v20.md`. Two consequences below. |
| **Item 07 is NOT the fix for the confounder/mediator bias** | Clarified 2026-09-07 (V20) | Item 07 has been the only open build track and the presumed remedy for everything. It is a substantive-model-compatible **exposure** draw, and it removes 0.1–11.7% of a bias that lives in the **covariate** block. It stays the right fix for the mixture failure (V3/V9) and the non-linear-outcome penalty (V13/V14) — both X-block defects. Two different defects in two different blocks had been folded into one queue. |
| **What a shippable Z-block fix would be** | ⚠️ **the objection is gone (V22); not yet a proposal** | V21's ladder attributes **100%** of the bias to BART's *flexibility*: a correctly specified **estimated** linear draw is +0.26% / +0.58%, and so are improper and donor-matched versions. Properness moves the interval (`b` −12–16%), not the bias; the inner-FCS loop contributes nothing. **But both cells have a linear-Gaussian covariate conditional, so every parametric arm is correct by construction** — and R8 adopted BART precisely because a parametric Z block is misspecified when that conditional is non-linear (V6: `mice pmm` at −4.70%). **V22 removed that objection**: the parametric draw wins on *both* sides (−0.14% vs BART's +2.84% under non-linearity), because a misspecified prior mean is absorbed by the coefficient the analysis already adjusts for. So a parametric covariate draw is now supported by evidence on both sides of the one trade that stood against it. **Still not a proposal**: it needs the analyst to supply a conditional form, and V19's implementations (`micePmm`, `properZ`) remain unexplained. `phase1/FINDINGS_v21.md`, `phase1/FINDINGS_v22.md`. |
| ~~Build a non-linear covariate conditional ON an X–Y path~~ | ✅ **built and run 2026-09-08 (V22)** | `z_role = "pipe_nl"` plus three scenarios. **Verdict: there is no trade** — the parametric draw is −0.14% where the covariate conditional is non-linear, against BART's +2.84%, and wins in both cells. It is unbiased even while *misspecified*, because the analysis absorbs the error into `logX1`'s coefficient — estimand-specific. `phase1/FINDINGS_v22.md`. |
| **⚠️ Item 07's grid exposure draw may be HARMFUL under a non-linear covariate arrow** | Open — **22-rep probe 2026-09-08; must be confirmed before item 07 ships** | A quick probe of `zr_pipenl_cens` with a real shipped baseline: `pipeline_bartMI` **+24.9%**, `pipeline_micePmm` +14.4%, exact-Z anchor +17.5% — and **`smc_xgrid` +93.7% with coverage 0.000**. The mechanism is plausible: the grid draw proposes from the exposure prior and reweights by `p(Y \| x, rest)`, so it never uses the `Z1 = g(logX1)` information at all, discarding the (linear, partially-right) covariate conditioning `leftcens` does use. **Item 07 is the only open build track and this is the one cell where it inverts.** 22 reps is not a result — but it is enough to say item 07 must not ship without this cell in its acceptance set. |
| **~20% ASYMPTOTIC bias in the censored-exposure draw under a non-linear covariate arrow** | ⚠️ **derived and registered as V23 2026-09-08; run pending (~4.8 h)** | With `logX1` censored and `Z1 = g(logX1)` non-linear, the exact-`Z` anchor sits at **+20.8 / +20.3 / +20.3%** across n = 800/3,200/12,800 — flat, so asymptotic. The X block alone carries it: `leftcens` draws `logX1` from a conditional **linear in its predictors** by construction (Phase 1 §7.5 chose the shash *margin* for skew, not a non-linear mean), so a covariate with a non-linear relationship to the censored exposure is outside what it can represent. **No cell before V22 had one.** 2–4× the covariate-draw defect this whole line of work has been chasing, and it lands on the strategy's entire purpose. Caveat: that cell has **no valid anchor and no shipped baseline** (`pipeline_bartMI` was not in the arm set), so it is a signal to investigate rather than a decomposed result. **V23 is that track, and it is derivation-led**: the defect follows from `leftcens` fitting a conditional linear in its predictors, it is asymptotic because the fit converges to the best linear approximation of a non-linear target, and its size is governed by `u` = the density-weighted residual of `g` after that linear fit — entering **squared** because `u` is an L²-projection residual and the first-order term vanishes (`THEORY.md` §3c). Registered law **bias% = 300·u²** over seven curvature settings. `phase1/run_v23_curvature.sh`. |
| **Why V19's parametric imputers are biased when the property is fine** | Open — narrowed 2026-09-07 (V21) | V19 measured `micePmm` at +2.4 to +3.9% and `properZ` at +4.5 to +6.3%, both *asymptotic*; V21's parametric arms sit at ±0.6%. The boundary was registered before V21 ran, so this is a **located discrepancy, not a contradiction**: `z21_pmm` draws `Z1` alone on the correct formula with `Z2` exact, while `micePmm` imputes `Z1`, `Z2` **and `Y`** together with mice's own predictor matrix and iteration. **V22 eliminated the `Z2` draw** — drawing it by the same method as `Z1` moves the estimate ≤0.17 pp. Two candidates remain: **`Y` as an imputation target**, and mice's own predictor matrix. Extending the ladder to impute `Y` is the next step and is cheap. |
| **`use_as_auxiliary` fix in `auto_preds`** | ⚠️ **documented 2026-09-07; code fix still unmade** | V17b measured the defect at **−0.30 ± 0.10 pp** (MCAR) and **+0.60 ± 0.10 pp** (MAR): detectable, under 1 pp, and opposite in sign between mechanisms. Adding `use_as_auxiliary` to `00_censored_exposure.R`'s `auto_preds` is two lines and makes the code match `docs/variable-dictionary.md`. Do it for consistency; do **not** present it as an accuracy fix. **The gap is now documented** in `docs/variable-dictionary.md` and `docs/covariate-roles.md`, with the `censored_exposure$predictors` override as the workaround. |
| Why the exponent is 2 | Open — the theory frontier | V15 confirmed `excess% ≈ 320·f²` over six points (exponent CI [1.98, 2.97], free fit 1.90), so the shape is measured but not derived. The conjecture in [`THEORY.md`](THEORY.md) §4 is a product of two `f`-linear losses, the second being the censored share of the *estimand's own evaluation range*. It is discriminable: re-running the V15 sweep with `q_lo` = 0.40 instead of 0.25 should move the curve under the conjecture and leave it unchanged otherwise. Same harness, same cost (~11 h), no new code: `QLO=0.40 ./run_v15_fsweep.sh` (the runner already reads `QLO`/`QHI`). **Queued as V21 — the theory frontier, but behind the remediation work.** |
| Curvature sign inversion at high censoring | Open, **user-facing** | V15 found the fitted curvature crosses zero near `f` ≈ 0.55 (40% of reps negative at `f` = 0.60) while coverage stays 0.92 on a 2.7× wider interval. Whatever guidance ships about censored mixtures should say that nominal coverage does not detect this, and should probably refuse the analysis above some `f`. |
| V8 confound: target count vs `Y`'s presence | Open, small | V8's 3-target cells all impute `Y` and its 2-target cells do not, so it cannot say which feature makes the inner-iteration shortcut unsafe. `.ce_default_inner_iter()` requires both conditions, which is safe but conservative. A cell with three non-outcome targets would separate them and could widen the cheap path. |
| V8 runtime saving unmeasured | Open, small | The ~3× figure comes from the fit count, not from timing: `secs` is logged per replication for all arms together. Needs per-arm instrumentation in the harness. |
| External testers | Standing goal | Have other users run the pipeline on their own data before broadening release claims. Each analysis pattern that trips someone up should become a bundled example, entering permanent regression coverage. |

Longer-term ideas (sensitivity-analysis runner, multi-outcome support, DAG-based
confounder selection, export hardening) remain in `NOTES.md` under *Ideas for future
development*.

---

## One thing not to do

**Do not adopt the V4 validation instrument as the fix.**
`.v4_proper_row_imputation()` is a deliberately simple linear/logistic draw, built to
isolate a mechanism rather than to impute well. It has no non-linearity and no
interactions — adequate for the simulation's linear-Gaussian covariates, inadequate in
general.

It also **overshoots**: width/SE reaches 1.085 where the shipped engine sat at 0.899,
i.e. intervals roughly 9% too *wide*, and its bias is marginally worse (−3.18% against
−2.11%). It proves where the variance goes missing. It is not the replacement.
