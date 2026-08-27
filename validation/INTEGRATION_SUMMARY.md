# Left-censored exposure integration — what was needed, what exists, what remains

**As of 2026-08-27**, after the BART arm was built and the BKMR stack installed. One-page
summary of the whole effort, from the congeniality argument through to the current open
questions.

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
| R8 | **Proper MI** — imputation-model parameters drawn from a posterior so between-imputation variance `B` is honest | X-block: proper by construction. Z-block: bootstrap the training data per imputation | `run_row_level_imputation_proper()`, `proper_draw` (default TRUE, v1.4.0) | V4 attributed the defect (`B` understated 8–59%). V5: fix closes 66% at 40% MCAR, 23% under heavy load. V6: `properBoot` kept over `mice pmm`, judged on a DGP that can punish both. **BART arm built, run pending** | ⚠️ **partial** |
| R9 | **Production scale** (n ≈ 80k, wide covariate sets) | Reduced X-block predictor set (cost ~`k²`, not `n`) | `censored_exposure$predictors` | Scale test: X-block ~1.6–2.4 min at `k` ≤ 50. V2 `large_n` at **n = 80,000**: bias −0.09%, coverage 0.960 | ✅ done |
| R10 | **Must not break existing analyses** — config compatibility | Opt-in strategy dispatch; new behaviour behind a flag; numbers versioned, not frozen | One `else if` in `03_impute.R`; `proper_draw` flag; two over-strict Step-1 checks downgraded to warnings | Regression test gave **bit-identical** coefficients before the v1.4.0 default flip; previously-valid configs verified to still pass | ✅ done |
| R11 | **Mixture / BKMR** exposure–response surfaces | *Undecided* — linear congenial draw shown insufficient | none (`bkmr` / `bkmrhat` / `aBKMR` now installed) | `FINDINGS.md` §7.7: linear congenial imputation is itself biased on a mixture surface (+7.7% → +19.6%, coverage → 0.62) | ❌ **not started** |
| R12 | A **documented `m`** for this path | FMI-based, per PLAN §6 | none | FMI measured at 0.26–0.44 (implying `m` ≈ 29–38); current default `m` = 30 looks about right | ❌ **not written up** |

---

## 2. Open items — what blocks them, and the theoretical resolution

| Item | Nature of the block | Theoretical resolution | Effort |
|---|---|---|---|
| **Residual variance gap** (R8). Under heavy covariate imputation, coverage reaches ~0.93–0.95 rather than nominal, and `properBoot` carries +2.0–2.6% bias where `mice pmm` sits near zero | Not a settings problem. A random forest is **already bagged**, so an outer bootstrap shifts it only slightly — bagging damps the very uncertainty the bootstrap injects (`B` +31–33% vs +41–59% for a parametric posterior draw). But the parametric draw is misspecified under non-linear covariates (bias −4.70%) | **BART** — non-parametric like a forest, fully Bayesian, so posterior draws are proper *by construction* rather than approximated by resampling. **Built 2026-08-27** (`dbarts` 0.9.34, `phase1/R/bart_impute.R`, arm `pipeline_bartMI`): imputation *j* is `yhat.test[j, ] + rnorm(n, 0, sigma[j])`, a posterior predictive draw. Criteria registered in `ROADMAP.md` item 06 — it must beat **both** incumbents | **Awaiting its run** (~7 h, 152 CPU-h measured). 8-rep smoke: `B` largest of all arms, and bias −0.73% in the misspecification cell against `micePmm`'s −4.73% — right direction, but noise. Watch `nl_combined`, where 3 targets force the `m`-chain path |
| **Mixture estimands** (R11) | **One block now, not two.** `bkmr` 0.2.2, `bkmrhat` 1.1.7 and `aBKMR` 0.1.0 are installed, so the fitting side is solved. What remains: each estimand's truth must be derived **analytically** from the generator rather than estimated from an oracle fit — otherwise "bias" is contaminated by the reference's own error. (Note `aBKMR` is a speed/memory variant of `kmbayes`, *not* a censored-data extension — none of the three handles censoring) | Substantive-model-compatible imputation: put the exposure–response *surface* into the imputation model, or use a joint model with the surface. Until then the mixture path should not be used. Tractable: `h(X) = b'logX + b_int·x₁x₂ + b_quad·x₁²`, so the q25→q75 overall effect is a closed-form difference given the MVN quantiles | Medium — the analytic truth derivation is the work; `bkmrhat::kmbayes_parallel` covers the MCMC cost |
| **Pilot `m` / FMI** (R12) | Nothing blocking; needs writing up | FMI → required `m` by the standard rule, cross-checked against the auto-increment loop's stopping rule | Small — FMI is already measured |
| **MAR covariate missingness** | Every scenario to date is MCAR | MAR-on-`Z` is the realistic mechanism; the injector exists (`inject_mcar_covariates`) and needs a MAR sibling | Small |
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
| v1.4.0 **largely corrects** this: 66% of the shortfall closed at 40% covariate missingness (coverage 0.933 → 0.947) | Partial. Under heavy combined missingness coverage reaches ~0.93, not 0.95 | `FINDINGS_v5.md`, `FINDINGS_v6.md` |
| Raising `m` does **not** fix interval width | Proved on byte-identical data at `m` = 30 vs 100 | `FINDINGS_v2.md` |
| MID is **necessary** — disabling it worsens calibration and adds −6 to −7% bias | | `FINDINGS_v4.md` |
| A **bootstrap** makes a forest imputer proper only partially, because bagging damps it | Measured against a parametric posterior draw on identical data | `FINDINGS_v5.md`, `FINDINGS_v6.md` |
| `mice pmm` is **not** a safe drop-in: bias −4.70% when covariate imputation is heavy and non-linear | It *is* better where imputation is outcome-dominated (bias ≈0 vs +2.0–2.6%) — the choice is regime-dependent | `FINDINGS_v6.md` |

**Not claimable yet:** anything about mixture/BKMR estimands; a recommended `m` for this
path; behaviour under MAR (as opposed to MCAR) missingness; behaviour under a non-linear
outcome model.

---

## Reading guide

- **Status index and navigation:** [`README.md`](README.md)
- **What to do next, in dependency order:** [`ROADMAP.md`](ROADMAP.md)
- **Why the design is what it is:** [`PLAN_leftcensored_exposure_integration.md`](PLAN_leftcensored_exposure_integration.md)
- **Pre-registered acceptance criteria:** [`PLAN_pipeline_validation.md`](PLAN_pipeline_validation.md)
- **Evidence:** `phase1/FINDINGS.md`, then `FINDINGS_v1`, `v2`, `v4`, `v5`, `v6`
