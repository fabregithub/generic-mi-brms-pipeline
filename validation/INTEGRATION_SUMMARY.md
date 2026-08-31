# Left-censored exposure integration — what was needed, what exists, what remains

**As of 2026-08-29**, after V9 confirmed the mixture failure's cause, V3 closed the mixture question, and v1.5.1 narrowed the Z-block inner-iteration default (V8), v1.5.0 adopted BART for the Z block (R8 closed) and the
pilot-`m` question closed (R12). One-page summary of the whole effort, from the
congeniality argument through to the current open questions.

---

## A note on the table design

The obvious shape — *item · implemented · not yet · blocks · solution* — has three problems
worth avoiding:

1. **"Implemented" and "not yet" are one field, not two.** A single `Status` column says the
   same thing without inviting half-filled rows.
2. **The most important column is missing: evidence.** For this project the interesting
   question is never just "is it built?" but "how do we know it works?". Making `Evidence` a
   column turns the table into a gap detector: **a row with implementation but no evidence
   is an untested claim; a row with a claim but no implementation is a stale one.** That has
   been the failure mode this whole effort kept catching.
3. **Status tracking and blocker analysis are different jobs.** Squeezing "theoretical
   resolution" into a status table makes both cramped — a blocker needs room to explain
   *why* it is blocked and what would unblock it.

So: **three tables.** A traceability matrix (§1) answering *how do we know*; an open-items
table (§2) answering *what is blocked and what would resolve it*; and a claims ledger (§3)
answering *what may we state publicly, and with what scope limit*. The third exists because
the documentation repeatedly risked outrunning the evidence, and a ledger makes that
visible.

---

## 1. Traceability matrix — requirement → design → implementation → evidence

| # | Requirement | Design decision | Implementation | Evidence | Status |
|---|---|---|---|---|---|
| R1 | The censored exposure must be imputed **conditional on the outcome** (congeniality), or the exposure–response coefficient is attenuated | Two-engine block-FCS; the X-block draw takes `Y` as a predictor | `leftcens::impute_censored_conditional()`; `00_censored_exposure.R` | `FINDINGS.md` H1: no-`Y` pre-step attenuates −5.6% → −14.4%, coverage 0.92 → 0.82. `FINDINGS_v1.md`: shipped engine tracks the oracle (bias ≤0.8%) | ✅ done |
| R2 | The **estimand must be preserved** — the coefficient on the modelling scale, not a proxy | Matched analysis model in the harness (`dgp_formula`); log-scale handling decoupled from the dictionary | `.ce_exposure_bounds()`, `log_scale` config flag | Oracle recovers `b_logX1` = 0.392–0.400 across every design tested | ✅ done |
| R3 | **Skew-robust** margins — real exposure distributions are right-skewed and the shape is undiagnosable under censoring | sinh-arcsinh (shash) margin, interval-censored MLE, latent-normal transform | `leftcens::fit_shash_margin()`, `margin = "shash"` | `FINDINGS.md` skew sweep: Gaussian/tobit breaks (bias −14%, coverage **0.47**); shash holds (≤2.7%, ~0.95) | ✅ done |
| R4 | **Interval censoring**: three-tier ND / detected-not-quantified / quantified | Per-row `X_lo`/`X_hi` bound columns (leftcens interval form) | `.ce_exposure_bounds()`; `examples/censored_exposure/` | Two-exposure example with a three-tier `expo2` recovers both coefficients; V2 `censor_all` | ✅ done |
| R5 | **Multiple** censored exposures | X-block loops over `exposure_vars` | `.ce_one_imputation()` | V2: 600 tasks at `cens=3` (every exposure censored), zero failures, bias ≤2.71% | ✅ done |
| R6 | **Missing outcomes** must not bias the fit | Multiple-imputation-then-deletion (von Hippel 2007): impute `Y` for the X-block predictor step, delete those rows before fitting | `mid_delete_imputed_y`, runtime assertion | V2: MID fired on 600 tasks with exact row counts. **V4: disabling MID makes calibration worse and adds −6 to −7% bias** — MID is load-bearing | ✅ done |
| R7 | **MAR covariates** imputed jointly with the censored exposures | Alternating block-FCS: miceRanger Z-block + leftcens X-block, several outer sweeps | `run_row_level_imputation()` reused inside the block loop | V2 `mcar_z20`/`mcar_z40`/`combined` exercise the alternation under Monte Carlo | ✅ done |
| R8 | **Proper MI** — imputation-model parameters drawn from a posterior so between-imputation variance `B` is honest | X-block: proper by construction. Z-block: **BART posterior draws** — non-parametric *and* fully Bayesian | `run_row_level_imputation_bart()`, `z_imputer = "bart"` (default v1.5.0) | V4 attributed the defect (`B` understated 8–59%). V5/V6: bootstrap-on-forest closes only part of it; `mice pmm` misspecified. **V7: BART has the best bias in all six cells (max 1.29% vs 2.67%), coverage 0.937–0.963, tuning-robust, ~4% of runtime** | ✅ done *(intervals run 3–8% wide — conservative; see caveat)* |
| R9 | **Production scale** (n ≈ 80k, wide covariate sets) | Reduced X-block predictor set (cost ~`k²`, not `n`) | `censored_exposure$predictors` | Scale test: X-block ~1.6–2.4 min at `k` ≤ 50. V2 `large_n` at **n = 80,000**: bias −0.09%, coverage 0.960 | ✅ done |
| R10 | **Must not break existing analyses** — config compatibility | Opt-in strategy dispatch; new behaviour behind a flag; numbers versioned, not frozen | One `else if` in `03_impute.R`; `proper_draw` flag; two over-strict Step-1 checks downgraded to warnings | Regression test gave **bit-identical** coefficients before the v1.4.0 default flip; previously-valid configs verified to still pass | ✅ done |
| R11 | **Mixture / BKMR** exposure–response surfaces | *Undecided* — linear congenial draw shown insufficient | none (`bkmr` / `bkmrhat` / `aBKMR` now installed) | `FINDINGS.md` §7.7: linear congenial imputation is itself biased on a mixture surface (+7.7% → +19.6%, coverage → 0.62) | ❌ **not started** |
| R12 | A **documented `m`** for this path | FMI-based, per PLAN §6, cross-checked against a measured width-vs-`m` curve | `m` = 30 shipped default | **V7 P3:** width stabilises by `m` ≈ 30 (10 → 30 narrows ~2.7%; 30 → 50 changes <0.2%). FMI flat at ≈0.30, and `m ≈ 100 × FMI` independently gives 30 | ✅ done |

---

## 2. Open items — what blocks them, and the theoretical resolution

| Item | Nature of the block | Theoretical resolution | Effort |
|---|---|---|---|
| ~~**Residual variance gap** (R8)~~ — **RESOLVED v1.5.0** | Was: forest under-disperses (bagging damps the bootstrap); parametric draw misspecified under non-linear covariates | **BART** — non-parametric *and* fully Bayesian, so posterior draws are proper by construction. Adopted as the `z_imputer` default. Bias ≤1.29% in every cell, coverage 0.937–0.963, robust at 50 and 200 trees | Done. `phase1/FINDINGS_v7.md` |
| ~~**interval overshoot**~~ — **NOT AN ESTABLISHED DEFECT** | The apparent 3–8% overshoot sits inside the diagnostic's own Monte-Carlo error. `width/SE` divides by the empirical SD of the estimates, whose SE is `sd/sqrt(2(n-1))` ≈ **4.1%** at 300 reps — so the ±3% acceptance band was tighter than the metric's noise floor and could not be met by anything. **None of BART's deviations is significant (max z = +1.87);** the only real deviation in the table is `properBoot` being *too narrow* (z = −2.10). Two alternative explanations (a `t`-vs-`z` divisor bias; light-tailed estimates) were tested and refuted | Nothing to fix on this evidence. BART's width is *consistent with calibration, bounded to ~±8%* — not proven calibrated. Tightening to ±1% needs ~4,800 reps (16×), which is likely not worth it. **Lesson: `width/SE` is a mechanism diagnostic, not an acceptance criterion — coverage is** (MC SE 0.013) | None warranted |
| ~~**Redundant BART sweeps on the censored path**~~ — **RESOLVED v1.5.1, mostly negative** | Was: block-FCS calls the Z block with `m`=1 per outer sweep, and BART then ran its own `bart_inner_iter`=3 iterations, so the path did ~3× the fits it needed — *if* the outer loop's alternation substitutes for the inner one | **It does not, in general.** The outer loop alternates Z ↔ X, not among the Z block's own targets. V8 (4×1000 reps, paired): with 3 targets including an imputed Y, `inner`=1 cost **−0.72 pp** of bias (bar 0.5 pp, whole CI outside); with 2 non-outcome targets it was free (+0.04 pp, CI ±0.13). Coverage and width/SE moved <0.02 in every cell — **a bias-only defect the calibration diagnostics could not see**. Shipped as `.ce_default_inner_iter()`: `inner`=1 only on the measured-free path, pre-fix 3 elsewhere. **No change for a stock censored-exposure analysis** (it imputes Y, so it takes the 3 branch) | Done. `phase1/FINDINGS_v8.md` |
| ~~**Mixture estimands** (R11)~~ — **RESOLVED 2026-08-29 (V3), negatively** | Was: the mixture verdict rested on a scaffold coefficient, not the reported estimands | **Confirmed and worse.** On real BKMR estimands with analytic truth, the shipped engine is the **worst of four arms** (mean 14.6% paired excess over the oracle vs 6.4% for LOD/√2). It destroys **57% of the curvature** at 40% non-detects; interaction survives single-exposure censoring. Imputing from a conditional that is *linear in the predictors* is worse than not imputing. **The mixture path must not be used.** Fix = substantive-model-compatible imputation, unbuilt | Done. `phase1/FINDINGS_v3.md` |
| **Substantive-model-compatible imputation** *(target validated by V9)* | The only remaining fix for the mixture path. The X block draws the censored exposure from a conditional linear in (Y, other X, Z); a mixture surface with curvature in X1 cannot be represented by it, and V3 measured the cost | Put the exposure–response *surface* into the imputation model, or use a joint model with the surface. Tractable in principle for a known parametric surface; the general case is research, not engineering | Large |
| ~~superseded row below, kept for its still-valid package notes~~ | **One block now, not two.** `bkmr` 0.2.2, `bkmrhat` 1.1.7 and `aBKMR` 0.1.0 are installed, so the fitting side is solved. What remains: each estimand's truth must be derived **analytically** from the generator rather than estimated from an oracle fit — otherwise "bias" is contaminated by the reference's own error. (Note `aBKMR` is a speed/memory variant of `kmbayes`, *not* a censored-data extension — none of the three handles censoring) | Substantive-model-compatible imputation: put the exposure–response *surface* into the imputation model, or use a joint model with the surface. Until then the mixture path should not be used. Tractable: `h(X) = b'logX + b_int·x₁x₂ + b_quad·x₁²`, so the q25→q75 overall effect is a closed-form difference given the MVN quantiles | Medium — the analytic truth derivation is the work; `bkmrhat::kmbayes_parallel` covers the MCMC cost |
| **Pilot `m` / FMI** (R12) | Nothing blocking; needs writing up | FMI → required `m` by the standard rule, cross-checked against the auto-increment loop's stopping rule | Small — FMI is already measured |
| ~~**MAR covariate missingness**~~ — **RESOLVED 2026-08-29 (V10)** | Was: every scenario to date was MCAR, so every additive-path claim rested on the easy mechanism | **MAR does not degrade the engine — it is slightly easier.** −0.17% vs MCAR's −1.67% at 40% missing; coverage 0.948–0.967; unchanged across a doubling of mechanism strength and in the non-linear and missing-outcome cells. The MCAR-based claims are therefore conservative. **MNAR untested** | Done. `phase1/FINDINGS_v10.md` |
| **Non-linear outcome model** | `Y` is linear in `(logX, Z)` by construction, which is *why* `pmm` wins the outcome-dominated cells | Adding curvature to the outcome model would test whether that advantage survives. Currently an untested confound in the V6 verdict | Small — one DGP option |
| **`resume` kills the V4 runner** | Cause never established. Reproducibly killed the parent after one chunk in both detached and foreground runs; a fresh checkpoint runs fine; the checkpoint's structure looks correct | Unknown. Off by default with a warning; partial results are recoverable directly from the checkpoint, so it costs nothing operationally | Unknown |
| **Translated READMEs** | `docs/README.{de,es,fr,ja}.md` contain **zero** mentions of the censored-exposure strategy or `leftcens` — they predate the feature | Needs a translation pass, not a patch | Medium |

---

## 3. Claims ledger — what may be stated, and its scope limit

| Claim | Scope limit | Evidence |
|---|---|---|
| Imputing a censored exposure **without** the outcome attenuates the exposure–response coefficient by up to **−14%** with coverage falling to 0.82 | Additive ERF; replicated on two independent seeds | `FINDINGS.md`, `FINDINGS_v1.md` |
| The shipped censored-exposure engine recovers the coefficient with **bias <1% and coverage 0.957** at up to 40% non-detects | **Additive / per-analyte ERFs only**; complete covariates | `FINDINGS_v1.md` |
| That holds with **every** exposure censored, at **n = 80,000**, under right-skew, and at exposure correlation 0.8 | Max bias 2.71% across ten scenarios | `FINDINGS_v2.md` |
| A **tobit/Gaussian** conditional is unsafe for skewed exposures — coverage collapses to **0.47** | Skew 0.75; motivates the shash default | `FINDINGS.md` |
| For **mixture surfaces**, linear congenial imputation remains biased (**+19.6%**, coverage 0.62) | Scaffold estimand, not true BKMR estimands — mechanism demonstrated, verdict not manuscript-final | `FINDINGS.md` §7.7, `FINDINGS_v1.md` |
| Pipeline intervals were **too narrow** wherever much was imputed, because the imputation step was improper (`B` understated 8–59%) | Applies to the **general** imputation path, not just censored exposures | `FINDINGS_v4.md` |
| v1.4.0 largely corrected this; **v1.5.0 closes it**: BART posterior draws give bias ≤1.29% in every cell and coverage **0.937–0.963** | Intervals run 3–8% **wide** in some cells — conservative, so precision is understated rather than confidence overstated | `FINDINGS_v5.md`, `FINDINGS_v6.md`, `FINDINGS_v7.md` |
| The imputer choice is **tuning-robust**: BART's bias stays ≤0.8% and width/SE moves ≤0.016 between 50 and 200 trees | One non-linear covariate form tested | `FINDINGS_v7.md` |
| **`m` = 30 is the right default** for this path — width stabilises there, and the FMI rule agrees | Measured on the v1.4.0 imputer; BART's larger `B` could shift it slightly upward | `FINDINGS_v7.md` |
| Raising `m` does **not** fix interval width | Proved on byte-identical data at `m` = 30 vs 100 | `FINDINGS_v2.md` |
| MID is **necessary** — disabling it worsens calibration and adds −6 to −7% bias | | `FINDINGS_v4.md` |
| A **bootstrap** makes a forest imputer proper only partially, because bagging damps it | Measured against a parametric posterior draw on identical data | `FINDINGS_v5.md`, `FINDINGS_v6.md` |
| The mixture failure is **caused by the imputation conditional's functional form**: replacing only that draw closes **89–100%** of the curvature gap and cuts mean excess over the oracle from 17.2% to 0.9% | Mechanism test with the surface's form supplied; harness instrument, not shipped code | `FINDINGS_v9.md` |
| **Oracle parameters are not required** — estimating the surface each sweep and drawing from its posterior matches the oracle version to 0.2 pp | The functional *form* was still given; obtaining it for an unknown surface is unsolved | `FINDINGS_v9.md` |
| The additive-path claims **hold under MAR**, not only MCAR: bias −0.17% vs MCAR's −1.67% at 40% covariate missingness, coverage 0.948–0.967 | One MAR form, driven by variables the imputation model already conditions on. MNAR untested | `FINDINGS_v10.md` |
| For **mixture surfaces the censored-exposure engine must not be used**: on real BKMR estimands it is the worst of four arms, destroying **57% of the curvature** at 40% non-detects — worse than LOD/√2 substitution | Curvature specifically; interaction estimands survive single-exposure censoring (+2.6% to +5.7%). One mixture surface, n = 800 | `FINDINGS_v3.md` |
| **Imputing with the wrong functional form is worse than not imputing** — on curvature, LOD/√2 is indistinguishable from the oracle (p = 0.42) while the congenial engine loses 57% | Measured against analytic truth, paired on identical datasets | `FINDINGS_v3.md` |
| **Coverage does not detect these failures.** A −42% point bias sat behind 0.945–0.960 coverage, because pooled intervals ran 1.5–1.6× wider than the true sampling spread | Observed twice independently (V3, V8) | `FINDINGS_v3.md`, `FINDINGS_v8.md` |
| The **V7 evidence applies to the shipped code**, not only to the validation harness: the arm V7 measured and `00_common_functions.R`'s BART path are **bit-identical** | 4000 paired replications, all six recorded quantities, max diff exactly 0 | `FINDINGS_v8.md` |
| Inside block-FCS, one Z-block FCS pass is **not** a substitute for the outer sweeps when the block has several targets: bias moved **−0.72 pp** with 3 targets, while coverage and width/SE moved <0.02 | Establishes the 3-target-with-Y and 2-target-without-Y configurations only; the design confounds target count with Y's presence | `FINDINGS_v8.md` |
| `mice pmm` is **not** a safe drop-in: bias −4.70% when covariate imputation is heavy and non-linear | It *is* better where imputation is outcome-dominated (bias ≈0 vs +2.0–2.6%) — the choice is regime-dependent | `FINDINGS_v6.md` |

**Out of scope by decision (2026-08-29):** **MNAR** missingness. It is a study-design
question rather than an analysis one — no imputation engine can recover a mechanism that
depends on the unobserved values, and the answer lies in how the data were collected. Not
listed as a gap.

**Not claimable yet:** behaviour under a non-linear *outcome* model; that BART's
intervals are correctly *calibrated* (they are conservative, not calibrated); *why* the
inner-iteration shortcut fails — target count and the outcome's presence are confounded in
V8's design.

---

## Reading guide

- **Status index and navigation:** [`README.md`](README.md)
- **What to do next, in dependency order:** [`ROADMAP.md`](ROADMAP.md)
- **Why the design is what it is:** [`PLAN_leftcensored_exposure_integration.md`](PLAN_leftcensored_exposure_integration.md)
- **Pre-registered acceptance criteria:** [`PLAN_pipeline_validation.md`](PLAN_pipeline_validation.md)
- **Evidence:** `phase1/FINDINGS.md`, then `FINDINGS_v1`, `v2`, `v4`, `v5`, `v6`
