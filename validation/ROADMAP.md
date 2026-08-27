# Roadmap — one defect, three tracks

**As of 2026-08-26**, after the v1.4.0 imputation fix. Items 01 and 02 are complete
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
| 04 — BKMR estimands | **Phase 2b** | Track V3, not started |

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

### 06 — A flexible *and* properly dispersed imputer · *new; carries the residual gap*

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

### 03 — Pilot `m` and document FMI · *small; design-plan Phase 5*

Most of the way there already: FMI is now reported per arm and lands between **0.26 and
0.44** across scenarios. What remains is turning that into a recommended default `m` for
the censored path and cross-checking it against the auto-increment loop's stopping rule.

**Do this after 01** — a proper draw changes `B`, which changes FMI.

### 04 — BKMR estimands · *independent; design-plan Phase 2b*

The last piece before the mixture verdict is manuscript-final. Today the claim that linear
congenial imputation fails for mixtures rests on a single local coefficient rather than the
estimands that matter: the overall mixture effect q25→q75, pairwise interactions, and the
response surface at representative exposure profiles.

Each estimand needs its truth computed **analytically from the generator**, not estimated
from an oracle fit — that is the real work in this track. Requires `bkmr`, **not currently
installed**.

---

## Dependencies

```
DONE (v1.4.0)
  [01 Proper Z-block draw] ──▶ [02 Re-validate] ──▶ partial: 66% / 23% of gap closed

DONE (V6)
  [05 Non-linear covariate DGP] ──▶ choice closed: properBoot stays; neither uniformly best

CRITICAL PATH (residual gap)
  [06 Flexible + properly dispersed imputer]   BART hypothesis; harness already exists

FOLLOWS
  [03 Pilot m / FMI]                                    (B has changed, so FMI has too)

INDEPENDENT
  [04 BKMR estimands]                                   (install bkmr first)
```

Item 06 now carries the residual gap; 03 is unblocked and cheap; 04 remains independent.

---

## Loose ends

| Item | State | Why it matters |
|---|---|---|
| `resume` kills the V4 runner | Unexplained | Reproducibly killed the parent after one chunk; a fresh checkpoint runs fine. Now off by default. Partial results remain recoverable straight from the checkpoint, so it costs nothing operationally — but an unexplained crash is worth understanding. |
| Translated READMEs | Stale | `docs/README.{de,es,fr,ja}.md` contain **zero** mentions of the censored-exposure strategy or its `leftcens` dependency — they predate the feature. Needs a translation pass, not a patch. |
| MAR covariate missingness | Unswept | Every robustness scenario used MCAR. MAR-on-covariates is the more realistic mechanism and is untested. |
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
