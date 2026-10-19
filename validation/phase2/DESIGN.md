# Phase 2 validation harness — design sketch

**Goal.** Extend the Phase-1 harness so the scored quantity is the **recovered
mixture HBEB benchmark** (not a scalar coefficient), on the birth-weight DGP with a
known non-linear surface (`manuscript_demo/PROVENANCE.md`). This is what licenses the
main-text real-data mixture HBEB. Maps to `../../NOTES.md` §2 and the integration
plan's **Phase 2b** (estimand extension) + **Phase 3** (surface-congenial engine).

Status: design only. Nothing here is built yet.

---

## 0. The one architectural principle

**One inversion function, applied to `h` for the truth and to each procedure's `ĥ`
for the estimate.** The benchmark estimand is *defined by the inversion pipeline*, so
truth and every estimate must go through the *same* code:

```
            ┌── true surface h ───────────────► invert ──► BENCHMARK_true
 DGP truth ─┤                                                    │  score
            └── population exposure dist f(x) ──► invert ──►      ▼  (bias/RMSE/coverage)
 procedure ── impute → fit surface ĥ → predict ĥ(x) ► invert ──► BENCHMARK_hat
```

Consequence: the inversion module is **shared** and should be the *same code we ship
for the real data* — extend the HBEB paper's `…/hbeb-framework-nakayamasf/R/burden_pipeline.R`
to the mixture case rather than writing a parallel inversion. Validating a re-implementation
would not validate the pipeline.

---

## 1. What changes vs Phase 1 (module map)

Phase-1 files: `phase1/R/{dgp,censoring,procedures,metrics}.R` + `run_phase1.R`.
Phase-2 mirrors the layout under `phase2/R/`. Reuse verbatim where marked.

| Module | Phase-1 | Phase-2 change |
|---|---|---|
| `dgp.R` | scalar `b_logX1`; `Y` Gaussian; p=3 | **birth-weight grams**; 7 PFAS; the PROVENANCE surface `h` (main + `u²` curvature + 2 interactions); **GFR confounder `R`**; three-tier censoring inputs; `make_truth()` also emits the **true benchmark** |
| `censoring.R` | single left-censored analyte | **multi-analyte** + **DNQ interval tier**; `limit_scale`; emits `leftcens` `_lo/_hi` columns |
| `procedures.R` | return a pooled scalar | return **completed datasets** (`imputed_list`); analysis+inversion moved to a shared step |
| **`analysis.R`** *(new)* | — | fit the **mixture surface** (qgcomp / BKMR) on a completed dataset → an object exposing `predict_h(x)` + draws |
| **`benchmark.R`** *(new)* | — | the shared inversion (extend `burden_pipeline.R` to mixtures); truth and estimate both call it |
| `metrics.R` | scalar bias/rmse/coverage | score the **benchmark** (headline) + surface-recovery + intermediate estimands; pool across imputations |
| `run_phase2.R` | grid × reps | grid gains the **surface/correlation** axes; per-engine compute budget |

Reused as-is: `rubin_pool()`, `one_row()`, the shash draw
`impute_censored_conditional()`, the driver skeleton, `inject_mcar_covariates()`.

---

## 2. Estimands (what gets scored)

Scored top-to-bottom; **the benchmark is the headline**, the rest are diagnostics that
explain *why* a procedure's benchmark is off.

1. **Mixture HBEB benchmark (headline).** For each burden target `T` and each
   **reduction direction** `d`, a scalar benchmark = how far the mixture must move along
   `d` to meet `T`. Directions keep a (p−1)-dim contour scorable as a few scalars:
   - `d_prop` — proportional reduction of all analytes (the overall-mixture direction);
   - `d_top`  — reduce only the highest-weight / highest-PIP analyte;
   - `d_pc1`  — reduce along the principal exposure axis.
   Targets `T`: relative reductions {10, 25, 50}% (counterfactual-robust; primary) and a
   couple of absolute LBW rates (sensitivity). Report the benchmark as (a) the mixture
   shift and (b) the implied component concentrations.
2. **Overall mixture effect** — joint q25→q75 change in mean birth weight (grams).
3. **Single-exposure effects** — `h` varying `X_j`, others at median (recovers `β_j`).
4. **Interactions / curvature** — the `u_PFOS·u_PFOA`, `u_PFOS·u_PFHxS`, `u_PFOS²` terms.
5. **Surface recovery** — integrated squared error of `ĥ − h` over the exposure support,
   reported **separately for the low-dose region** (below the analyte LODs) — this is where
   the benchmark lives and where partial identification bites.

`make_truth()` computes 1–4 exactly from the known `h` + DGP exposure distribution, via
the shared `benchmark.R`, and stores them in the `*_truth.rds`.

---

## 3. Procedures = imputation strategies (analysis+inversion is shared)

Cleaner separation than Phase-1: a **procedure produces `m` completed datasets**; a shared
`analyze_and_invert()` maps each to the estimands of §2 and pools. Oracle / complete-case /
substitution are degenerate "imputations" (`m = 1`). This mirrors the production pipeline
(`imputed_list → analysis`), so the harness exercises the real data path.

| # | Procedure | What it is | Expectation |
|---|---|---|---|
| 1 | `oracle` | uncensored complete data | upper bound |
| 2 | `substitution` | half-limit (MDL/2) fill | biased (reference-paper method) |
| 3 | `impute_blind` | `leftcens` pre-step, **no Y** | biased (uncongenial) |
| 4 | `impute_linear` | Y-aware **linear** shash draw (`impute_censored_conditional`) | OK additive; **biased on the surface (§7.7)** — the trap, now on the *benchmark* |
| 5 | `impute_surface` | **surface-congenial** (SMC) — §4 below | the proposed fix: recovers the benchmark |
| 6 | `complete_case` | drop censored rows | ~unbiased slope, loses low-dose region → benchmark inefficient/■ |
| (R) | `joint_brms` | joint model with the surface, subsample | gold-standard reference (§4) |

Procedures 1–4, 6 already exist in miniature in `phase1/R/procedures.R` (retarget from the
scalar to producing completed datasets). **5 and R are the new build.**

---

## 4. The surface-congenial engine (the research contribution + main risk)

Phase-1 proved procedure 4 fails on the surface. Two principled fixes; build one as
production, keep the other as the reference.

**(A) SMC-FCS censored draw — production candidate. → BESPOKE BUILD (resolved 2026-08-14).**
Draw the censored `X` from `p(X | Y, Z, X_obs)` **under the substantive (surface) model**
(Bartlett 2015), so the imputation is congenial with the non-linear `h`.

> **`smcfcs` capability check (v2.0.2, verified against the installed package).** smcfcs has
> **no native left-/interval-censored (LOD) covariate support**: its covariate `method`
> options are `norm`/`logreg`/`poisson`/`latnorm`/…, where **`latnorm` is classical
> *measurement error*** (requires an `errorProneMatrix` of ≥2 replicate error-prone measures)
> — a different mechanism from censoring below a detection limit. No `lod`/`detection`
> argument exists, and the API exposes no hook to substitute a truncated/censored draw for a
> covariate. **However**, its rejection sampler *does* handle interaction + quadratic
> substantive models (`ex_lininter`, `ex_linquad`, `ex_logisticquad`, `ex_coxquad`) — i.e.
> the congeniality side we need. Conclusion: smcfcs solves the *surface-congeniality* half
> but not the *censoring* half, and cannot be a black-box engine here.

**Decision: build procedure 5 bespoke, fusing the two published mechanisms —**
- **censoring + skew (proposal):** `leftcens`'s interval-censored sinh-arcsinh draw
  (`impute_censored_conditional`) supplies a bound-respecting, skew-correct *proposal* for the
  censored `X` on the latent scale;
- **surface congeniality (rejection):** the SMC-FCS accept step (Bartlett 2015) — accept the
  proposed `X` with probability ∝ `p(Y | X, Z)` under the current substantive surface. This is
  the exact synthesis smcfcs cannot do off-the-shelf; smcfcs's source is a **reference
  implementation of the rejection sampler to mirror, not a runtime dependency.**
- Home: a `leftcens` variant (Phase-3 package item) or a harness-local procedure first.

> **Congeniality nuance to record.** SMC-FCS makes the draw congenial with a *parametric*
> substantive model. Our analysis engine is BKMR (non-parametric). So procedure 5 is congenial
> with a **flexible parametric surrogate** of the surface (interactions + splines) used as the
> imputation's substantive model, *not* with the BKMR kernel itself. The harness must therefore
> check that BKMR-on-procedure-5-completed-data recovers the benchmark — that residual gap
> (parametric-surrogate vs kernel) is itself a §7.7-style question, and the joint-model
> reference (B) is the fully-congenial yardstick.

**(B) Joint `brms` model with the surface — gold-standard reference (subsample).**
`bf(Y ~ h(mi(X_1..X_p)) + Z) + bf(X_j | mi() + cens(lod) ~ Z)`, congenial by construction.
> **Known risk (from Phase-1).** The *naive* `mi() + cens()` encoding collapses to
> **LOD-substitution** in the installed brms (censored rows carry the LOD, no latent
> per-obs parameter is created) and is biased. The reference therefore needs a **correct
> latent-variable encoding** (NA at censored cells + a truncated latent, or custom Stan),
> validated on the additive DGP first (where the answer is known) before it is trusted on
> the surface. Treat (B) as reference-on-a-subsample, not production.

**De-risk order:** settle `smcfcs` censoring/skew support → validate (B)'s latent encoding
on the additive DGP → run 5 vs R on the mixture surface at moderate `n`. Whichever recovers
the **benchmark** against (B) becomes the production mixture engine.

---

## 5. Analysis engines (shared downstream)

`analysis.R` fits the mixture surface on a completed dataset and returns an object with
`predict_h(x)` + uncertainty draws, consumed by `benchmark.R`.

**BKMR is the benchmark engine, not qgcomp.** The benchmark is a continuous concentration;
qgcomp **quantizes** exposures and cannot return one, so it is at most a quantized
robustness comparator — never the inversion engine (see `../../NOTES.md` §3).

- **A-BKMR (primary; accelerated, continuous surface).** Guo et al. 2026, *The Innovation*,
  doi:10.1016/j.xinn.2025.101248 (`aBKMR`): GPP + matrix decomposition; continuous ERF curves
  + 95% CI that invert to the benchmark; scales past JECS n. Benchmark inverts **per posterior
  draw** → a benchmark posterior; pool across imputations by stacking within-imputation draws
  + the finite-`m` correction in `06_posterior_summary.R`. Integration: Python backend
  (`reticulate` / subprocess).
- **Full `bkmr` (fidelity reference / fallback).** Bobb 2015, native GPP (`knots=`), on a
  tractable n — the yardstick the accelerated surface is checked against (§7), and the
  published fallback.
- **F-BKMR** — pure-R `bkmr`-API drop-in alternative (arXiv:2411.10858); adopt if it clears
  peer review or if the A-BKMR Python dependency is costly.
- **qgcomp / WQS** — optional quantized comparators only.

Matched-model principle preserved: every procedure feeds the **same** engine + inversion,
so differences are attributable to censored-exposure handling alone.

---

## 6. Metrics

Extend `summarise_phase2()`:
- **Benchmark:** bias, rel_bias, RMSE, **interval coverage** per (procedure × T × direction
  × scenario). Coverage of the benchmark is the headline number.
- **Surface recovery:** ISE(`ĥ−h`) overall and **low-dose-restricted**.
- **Intermediates:** overall mixture effect, single-exposure, interaction terms (as Phase-1
  but vectorised).
- MC error reported per cell (as Phase-1's `n_ok`).

---

## 7. Grid, engines, compute

- **Design factors** (PROVENANCE §7): ND fraction {real, increased via `limit_scale`} ×
  inter-analyte correlation {low, real, high} × surface {additive, interaction+curvature}.
  Target regime where procedure 5 ≫ 4: **weak correlation + heavy censoring + interactions**.
- **`n`:** moderate (2000–5000) for the full grid; confirm at one large `n`
  (imputation FMI is ~n-invariant — E8 / Phase-1 scale test).
- **Reps:** ~200 for the fast configurations (small MC error); fewer for the heaviest
  A-BKMR / subsample cells, enough to separate procedures.
- **BKMR engine — decided (2026-08-14; see `../../NOTES.md` §3).** Primary = **A-BKMR** (Guo
  et al. 2026, *The Innovation*, doi:10.1016/j.xinn.2025.101248; `aBKMR`), a published,
  continuous-surface accelerator that scales past JECS n. Fallback / fidelity yardstick =
  **full `bkmr`** native GPP (`knots`; Bobb 2015). Watch = **F-BKMR** (arXiv:2411.10858;
  pure-R `bkmr`-API drop-in). Not variational-HS as primary (VB underestimates variance vs
  the honest-benchmark goal).
- **Low-dose fidelity check (first-class, not optional).** Every accelerator approximates full
  BKMR and is least accurate in sparse low-density exposure regions — exactly the low-dose,
  censored region the benchmark is read from. Before trusting A-BKMR in the grid, check its
  surface **vs full `bkmr`** on a tractable n, scored on the **low-dose surface (ISE) and the
  recovered benchmark**, not global R². This reuses §6's surface-recovery metric and doubles
  as a paper result (the engine is certified where it matters).
- **Parallelism:** over reps (as Phase-1, `NCORES`); BKMR cells checkpoint per cell.

---

## 8. Build order (each gates the next)

1. **`benchmark.R`** — extend `burden_pipeline.R` to a mixture joint counterfactual +
   inversion; unit-test that truth-in = truth-out on the known `h`. *(Foundation; also used
   by the real-data main text.)* **DONE (2026-08-14)** — `phase2/R/benchmark.R`: surface
   contract (`make_surface`, predict-capable posterior draws), directional counterfactuals
   (`dir_proportional`/`dir_analyte`/`dir_pc1`), the mixture burden `Δμ = E[h(X_obs)] −
   E[h(X_cf)]`, `residual_cases`, `benchmark_mixture` (uniroot), and `mc_benchmark_mixture`.
   `benchmark_selftest()` passes: (A) single-pollutant parity vs the shipped
   `burden_pipeline.R` (cases 2.6e-9, benchmark 3e-5), (B) truth-in=truth-out round trip on a
   known non-linear surface (≤3e-6). Run: `BENCHMARK_SELFTEST=1 Rscript R/benchmark.R`.
2. **`dgp.R` + `censoring.R`** — the birth-weight surface DGP + multi-analyte three-tier
   censoring; emit `*_truth.rds` incl. the true benchmark (via step 1). **DONE (2026-08-14)** —
   `phase2/R/dgp.R` (`make_truth`/`simulate_complete`/`h_predict`/`true_surface`/
   `true_benchmark`) + `phase2/R/censoring.R` (`inject_censoring` three-tier ND/DNQ/quant →
   `leftcens` `_lo`/`_hi`; `inject_mar_covariates`). Smoke tests (`DGP_SMOKETEST=1 Rscript
   R/dgp.R`) confirm vs PROVENANCE: medians exact, BW mean/SD 3027/378 & LBW 8.05%,
   grounded correlations matched, PFHxS weakest, GFR confounding ratio 1.86≈2×, detection
   targets hit (PFOS .990…PFTrDA .785), `limit_scale` increases censoring. *(`true_benchmark`
   computes the truth in-memory; persisting `*_truth.rds` is the driver's job, step 7.)*
   Tunable knobs left in `make_truth`: realised `sdlog` ~0.1 above target and PFHxS corr ~0.4
   (shared determinants add spread/correlation) — documented, not blocking.
3. **`analysis.R`** — the **A-BKMR** surface engine + `predict_h`/draws (full `bkmr` on
   small n as the correctness check); wire `analyze_and_invert()`.
4. **Retarget procedures 1–4, 6** to produce completed datasets; reproduce the §7.7 trap
   **on the benchmark** (procedure 4 biased on the surface) — the first headline result.
5. **Procedure 5 (surface-congenial)** + reference R; the de-risk order of §4.
6. **Scale the A-BKMR grid** (after the §7 low-dose fidelity check) — the reduced BKMR grid
   at the larger n / subsample.
7. **`metrics.R` + `run_phase2.R`** — full grid, tables/figures for the Supplementary
   (Fig 3 benchmark recovery; Fig 4 robustness surface).

Steps 1–4 already produce a publishable result (substitution/naive/linear-congenial bias the
benchmark; the trap is real for HBEB). Step 5 is the fix and the contribution. Step 6 is the
BKMR showcase.

---

## 9. Risks / open questions

- **`smcfcs` censoring support — RESOLVED (2026-08-14): none.** smcfcs 2.0.2 has no
  censored/LOD covariate method (`latnorm` = measurement error). Procedure 5 is therefore a
  **bespoke** `leftcens`-proposal + SMC-rejection draw (§4), not an smcfcs call. New residual
  risk: the parametric-surrogate-vs-BKMR congeniality gap (§4), measured in the harness.
- **Joint-brms latent encoding** — the naive form is a known trap (§4B). Validate on the
  additive DGP before trusting it as the reference.
- **Benchmark pooling under BKMR** — combine within-imputation posterior draws with
  between-imputation variance; reuse the pipeline's finite-`m` correction, but check its
  symmetry assumption holds for the (possibly skewed) benchmark posterior.
- **Partial identification** — the low-dose surface is not point-identified; the benchmark
  for aggressive targets inherits it. Report the modelled-tail sensitivity as a first-class
  result, not a footnote (connects to the HBEB "no extrapolation beyond observed range" point).
- **Keep the harness auditable** in the `leftcens`/`validation` style (seed-fixed,
  reproducible tables) as Phase-1 did.
