# Validation — status index

> ## 🔑 Magic word: `ORIENT`
>
> Type **`ORIENT`** as the entire first message in a new Claude Code session. It means:
>
> *"Read this file, check for a validation run in flight or new results, check
> `git status`, then tell me in a few lines where we are, what is running or just
> finished, and what the next step is. Do not re-explain the project history, and do not
> commit anything."*
>
> Nothing else needs typing. Everything required is in this file, `INTEGRATION_SUMMARY.md`,
> `ROADMAP.md`, and `phase1/FINDINGS_*.md`.

**Single source of truth for what has been validated, what is open, and what is next.**

> **New: [`THEORY.md`](THEORY.md)** — the thirteen tracks turn out to be instances of a
> single condition. The imputation must draw from `p(v | Y, rest) ∝ p(Y | v, rest; θ)·p(v | rest)`,
> and the pipeline's linear-Gaussian draw equals that **iff** the outcome is linear in `v`
> **and** `p(v | rest)` is Gaussian. That boundary is derived and numerically exact (theory
> 0.8944 against 0.894 measured); the magnitudes when it is violated are **not** derived, and
> three failed attempts are recorded there so they are not retried.
Keep this file current: the validation study is the evidence base for the claims the
README makes about bias, coverage and scope, so a stale entry here silently becomes a
wrong claim in the documentation.

Last updated: **2026-09-07** · after V21 — the bias is **BART's smoothing, 100% of it**: a correctly specified *estimated* covariate draw is unbiased (+0.3 / +0.6%). So a fix exists in principle — but it is a **trade** (flexibility for correct specification) and this run tests only the side where parametric wins. User guidance in [`docs/covariate-roles.md`](../docs/covariate-roles.md)

---

## Current status

### Design plan — Phases 1–5

Tracks the *build* of the censored-exposure feature.
Document: [`PLAN_leftcensored_exposure_integration.md`](PLAN_leftcensored_exposure_integration.md)

| Phase | What | Status |
|---|---|---|
| 1 | Confirm the ERF bias; test the congenial fix | ✅ done 2026-08-13 |
| 2 | Decision gate — build or not | ✅ resolved: build linear for additive; mixtures need SMC |
| 2b | Estimand extension (BKMR estimands) | ✅ **done 2026-08-29** — V3; mixture verdict now manuscript-final |
| 3 | Build the §4 component in `leftcens` | ✅ core shipped in `leftcens` 0.9.0 |
| 4 | Wire block-FCS into the pipeline | ✅ built + hardened |
| 5 | Pilot `m`, document FMI | ✅ **done 2026-08-28** — `m` = 30 confirmed (V7 P3) |

### Validation plan — Tracks V0–V21

Tracks the *validation of the shipped code*.
Document: [`PLAN_pipeline_validation.md`](PLAN_pipeline_validation.md)

| Track | Question | Status | Headline |
|---|---|---|---|
| V0 | Do Steps 5–11 + Quarto work on the censored path? | ✅ passed | 78 s, 5/5 criteria |
| V1 | Is the *shipped* engine unbiased with nominal coverage? | ✅ passed | bias ≤0.8%, coverage 0.957 |
| V2 | Does that hold outside the tested corner? | ✅ 9/10 + 1 finding | max bias 2.71%; one under-coverage cell |
| V3 | Does the mixture verdict hold on real BKMR estimands? | ✅ **resolved** | **Yes, and worse.** Pipeline destroys 57% of curvature and is the **worst of four arms**; LOD/√2 beats it. Interaction survives single-exposure censoring |
| V4 | Where does the V2 under-coverage come from? | ✅ resolved | improper MI in the Z block understates `B` by 8–59% |
| V5 | Does the candidate fix work? | ✅ partial | `proper_draw` default TRUE in v1.4.0; 66% of gap closed at 40% MCAR, 23% under heavy combined missingness |
| V6 | Forest+bootstrap or parametric+posterior? | ✅ resolved | `properBoot` stays: `mice pmm` bias −4.70% under non-linear covariates. Neither uniformly best |
| V7 | Does BART close R8? And what `m`? | ✅ resolved | **BART adopted (v1.5.0)**: best bias in all 6 cells, coverage 0.937–0.963, tuning-robust. `m` = 30 confirmed |
| V8 | Is the Z-block inner-iteration shortcut safe? Is V7's BART arm the shipped code? | ✅ resolved | Shortcut **unsafe with ≥3 targets** (−0.72 pp bias, bar 0.5); kept only where measured free (v1.5.1). V7 port **bit-identical** |
| V9 | Is the linear conditional really what destroys curvature? | ✅ resolved | **Yes.** Replacing only that component closes **89–100%** of the gap; mean excess 17.2% → 0.9%. Estimating the surface rather than knowing it costs 0.2 pp |
| V10 | Do the additive-path claims hold under **MAR**, not just MCAR? | ✅ resolved | **Yes — MAR is *easier*.** Bias −0.17% vs MCAR's −1.67% at 40%; coverage 0.948–0.967. Registered ±1 pp bar breached, favourably |
| V11 | Does the imputer choice survive a **non-linear outcome**? | ✅ resolved | **Yes** — `bartMI` never significantly beaten. But **every** arm loses ≈3 pp of bias: a pipeline-wide scope limit, not an imputer question |
| V12 | Is that penalty the **shape** of the Z-block draw? | ✅ resolved | **Partly** — shape costs 1.3–1.9 pp at heavy covariate missingness (*p* < 1e-4), but only **32–49%** of the penalty. The rest is the conditional mean and/or the exposure draw |
| V13 | Which block causes the rest of that penalty? | ✅ resolved | **The exposure draw** (3.1–3.4 pp of ~3.9). A **Z-only fix moves bias away from truth** — item 08 closed as scoped, item 07 widened |
| V14 | Does a **shippable** exposure draw (formula only) match the exact one? | ✅ resolved | **Mostly** — removes **83%** of the penalty and is the **best-calibrated arm** (width/SE 1.004). Fails ±1 pp in 1 of 3 cells; needs a missing-`Y` path |
| V15 | **Prediction:** curvature attenuation = `−323.4·f²`. Right? | ✅ **held** | **Yes** — 4 of 4 unmeasured cells inside ±10 pp, one to 0.1 pp; free exponent **1.90**, slope CI [1.98, 2.97]. Past `f` ≈ 0.55 the surface **inverts** while coverage stays 0.92 |
| V16 | **The root claim:** does omitting `Y` from the imputation attenuate — and in which block? | ✅ resolved | **Yes, −13.1 pp**, isolated inside one estimator for the first time (derived −14.6). The **exposure** block carries all of it; the covariate block **+1.2 pp**. The registered *leak* prediction **failed** — the blocks are additive here |
| V17 | Does the covariate's **causal role** change any of this? fork · pipe · collider · mixed | ✅ resolved | **It decides it.** The covariate block's penalty is **+1.2 pp** when `Z` is inert and **+25 to +28** when `Z` is adjusted for *and* on an X–Y path. MAR-vs-MCAR adds ≈0. **New: the shipped default is +5% to +10.4% biased there, coverage 0.93–0.97** |
| V18 | Does that bias wash out at larger `n`? | ✅ resolved | **No — and neither registered account survives.** Bias decays as **n^−0.335** (CI [−0.378, −0.291]), not 0 (constant) and not −0.5 (finite-sample). Coverage 0.87–0.92 at n = 12,800 against my projected 0.32–0.73 — **the V17 alarm was overstated**, the defect is permanent |
| V19 | Does that decay rate track the Z-block **imputer**? | ✅ resolved | **No — family predicts nothing**, but only `bartMI` decays at all (−0.311). `properBoot` **−23% asymptotic, coverage → 0**; `properZ` and `micePmm` flat at +4 to +6%. **The v1.5.0 BART default is vindicated on a question it was never tested against** |
| V20 | **Which block** carries that bias? | ✅ resolved | **The covariate draw.** Fixing the Z block removes **85–104%**; fixing the exposure draw **0.1–11.7%**; interaction ≤0.40 pp. **All registered gates passed** — the first since V15. So **item 07 is not the fix for this**, and a shippable Z-block fix is the open question |
| V21 | **Which property** of the covariate draw has to be right? | ✅ resolved | **Flexibility, and nothing else.** BART's smoothing carries **100%**; a correctly specified estimated draw is **+0.3 / +0.6%**, and so are improper and donor-matched versions. Properness moves the *interval* (`b` −12–16%), not the bias. **A fix exists in principle — but it is a trade, and the deciding cell (non-linear covariate conditional, on-path) does not exist yet** |

> **Two numbering schemes coexist** — design-plan Phases and validation-plan Tracks. They
> are not the same sequence and do not map one-to-one. The mapping table is in
> [`ROADMAP.md`](ROADMAP.md); consult it before answering "which phase is this?".

---

## Open items

| # | Item | Where |
|---|---|---|
| ~~01~~ | ~~Proper Z-block draw~~ — **shipped v1.4.0**, partial (see V5) | [`phase1/FINDINGS_v5.md`](phase1/FINDINGS_v5.md) |
| ~~02~~ | ~~Re-validate the fix~~ — **done** as the V5 `properBoot` arm | [`phase1/FINDINGS_v5.md`](phase1/FINDINGS_v5.md) |
| ~~05~~ | ~~Non-linear covariate DGP~~ — **done (V6)**, choice closed | [`phase1/FINDINGS_v6.md`](phase1/FINDINGS_v6.md) |
| ~~06~~ | ~~Flexible + properly dispersed imputer~~ — **BART adopted v1.5.0**, R8 closed | [`phase1/FINDINGS_v7.md`](phase1/FINDINGS_v7.md) |
| ~~03~~ | ~~Pilot `m` / FMI~~ — **closed**, `m` = 30 confirmed | [`phase1/FINDINGS_v7.md`](phase1/FINDINGS_v7.md) |
| ~~04~~ | ~~BKMR estimands~~ *(= Phase 2b, = Track V3)* — **closed 2026-08-29**; superseded by open item **07** below | [`phase1/FINDINGS_v3.md`](phase1/FINDINGS_v3.md) |
| ~~08~~ | ~~Shape-aware Z-block draw~~ — **closed as scoped (V13)**: harmful on its own; the exposure draw dominates and belongs to 07 | [`phase1/FINDINGS_v13.md`](phase1/FINDINGS_v13.md) |
| **07** | **Substantive-model-compatible imputation** — the mixture-path fix. **V9 validated the target**: the fix works and needs the right functional *form*, not oracle parameters. What remains is obtaining that form when the surface is unknown | [`phase1/FINDINGS_v9.md`](phase1/FINDINGS_v9.md) |
| — | ~~MAR covariate missingness~~ — **done (V10)**; the MCAR-only scope gap is closed | [`phase1/FINDINGS_v10.md`](phase1/FINDINGS_v10.md) |
| — | Interval overshoot (BART runs 3–8% wide) — low priority, conservative direction | [`INTEGRATION_SUMMARY.md`](INTEGRATION_SUMMARY.md) §2 |

Plus loose ends — the unexplained `resume` failure, stale translated READMEs, unswept MAR
covariate missingness, external testers — all listed in [`ROADMAP.md`](ROADMAP.md).

### Earlier run, kept for its still-open follow-ups — V8, 2026-08-28 (135 min)

`bash validation/phase1/run_inner_iter_check.sh` — 4 scenarios × 1000 reps × 3 arms,
`results/inneriter_latest.rds`. Written up in [`phase1/FINDINGS_v8.md`](phase1/FINDINGS_v8.md).

1. **The inner-iteration fix — failed as shipped, now narrowed.** Defaulting the Z block's
   inner FCS count to 1 inside block-FCS assumed the outer loop's alternation substitutes
   for it. It does not: the outer loop alternates Z ↔ X, not among the Z block's own
   targets. With **3 targets** (Z1, Z2, imputed Y) the paired bias difference was
   **−0.72 pp** against a 0.5 pp bar, with the whole CI outside it; with **2 non-outcome
   targets** it was free (+0.04 pp, CI ±0.13). Coverage and width/SE moved <0.02 in every
   cell — **a bias-only defect the calibration diagnostics could not see.** Shipped in
   v1.5.1 as `.ce_default_inner_iter()`, which picks the count per sweep from the target
   set; `bart_inner_iter` is now left unset in `00_config.R` so that choice is live.
   **A stock censored-exposure analysis is unaffected** — it imputes Y, so it takes the
   pre-fix `inner = 3` branch.
2. **The port check — passed exactly.** V7's BART arm overrode `run_row_level_imputation`
   with the harness instrument. That instrument and the shipped `00_common_functions.R`
   path are **bit-identical** across all 4000 replications on `estimate`, `se`, `ci_lo`,
   `ci_hi`, `ubar` and `b` (max difference exactly 0). **V7's conclusions transfer to the
   shipped code without qualification.**

Open follow-ups from V8, both small: the design confounds "3 targets" with "Y is a target"
(a 3-non-outcome-target cell would separate them and could widen the cheap path), and the
~3× runtime saving is inferred from the fit count, not measured — `secs` is logged per
replication, not per arm.

### ✅ Last completed run — V21, 2026-09-07 (316 min)

3 `n` levels × 2 cells × 6 arms × 500 reps = 3,000 tasks, zero errors.
[`phase1/FINDINGS_v21.md`](phase1/FINDINGS_v21.md) · verdict arithmetic:
[`phase1/analyze_v21.R`](phase1/analyze_v21.R)

**V20 located the bias in the covariate draw; V21 asked which *property* of that draw has to
be right.** A ladder from V20's known-unbiased exact anchor to BART, changing one property at
a time.

| arm | `zr_fork` @ 800 → 12,800 | `zr_pipe` | share of the span |
|---|---|---|---|
| `exact` *(anchor)* | −0.05% → −0.01% | +0.31% → +0.21% | 0% |
| `fit_proper` | **+0.26%** → −0.03% | **+0.58%** → +0.24% | 6.6% / 5.8% |
| `fit_improper` | +0.01% → −0.02% | +0.41% → +0.26% | 1.3% / 2.1% |
| `pmm` | +0.07% → −0.04% | +0.64% → +0.23% | 2.5% / 7.0% |
| **`bart`** | **+4.73%** → +1.31% | **+5.11%** → +1.64% | **100%** |

**It is the flexibility, and nothing else.** A correctly specified *estimated* linear draw is
unbiased; so is an improper one, and so is a donor-matched one. Properness behaves exactly as
V4 predicted — it moves the **interval** (`b` down 12–16%, coverage 0.95 → 0.94), not the
point estimate. The inner-FCS loop contributes nothing. Every registered gate passed, and
V21's BART arm matches V20's to ≤0.10 pp.

**A fix exists in principle** — that was the consequential registered outcome, and it went
the useful way: had a correct estimated draw exceeded 2 pp, there would be no fix at all.

**But it is a trade, and this run tests one side of it.** Both cells have a linear-Gaussian
covariate conditional, so every parametric arm is correct *by construction*. R8 adopted BART
precisely because a parametric Z block is misspecified when that conditional is non-linear —
V6 measured `mice pmm` at −4.70% there. Trading a 5% bias under linearity for a 5% bias under
non-linearity is not a fix. **The deciding cell — a non-linear covariate conditional on an
X–Y path — does not exist in the harness** and has to be built before anything is shipped.

**One discrepancy, located rather than explained away.** V19 measured `micePmm` and `properZ`
at +2.4 to +6.3% asymptotic, both parametric; here every parametric arm is ±0.6%. The driver
registered that boundary in advance: `z21_pmm` draws `Z1` alone on the correct formula, while
`micePmm` imputes `Z1`, `Z2` **and `Y`** together with mice's own predictor matrix. So the
*property* suffices and the *implementations* fail for some other reason — now the narrowed
open question.

*(Cost: 316 min against 4.9 h predicted.)*

## Document map

| Document | What it is | Read it when |
|---|---|---|
| [`INTEGRATION_SUMMARY.md`](INTEGRATION_SUMMARY.md) | **Start here.** Requirement → design → implementation → evidence traceability; open blocks with theoretical resolutions; the claims ledger | Orienting, writing methods text, or checking what may be claimed |
| [`THEORY.md`](THEORY.md) | **Why the pipeline fails where it fails.** One congeniality condition, derived and numerically exact, that explains all thirteen tracks — plus an explicit list of what is *not* derived | Before designing a fix, or when a new result looks surprising |
| [`ROADMAP.md`](ROADMAP.md) | What is left to do, in dependency order | Deciding what to work on next |
| [`PLAN_leftcensored_exposure_integration.md`](PLAN_leftcensored_exposure_integration.md) | The design/theory document: congeniality argument, two-engine block-FCS, estimand gates | Understanding *why* the feature is built the way it is |
| [`PLAN_pipeline_validation.md`](PLAN_pipeline_validation.md) | The validation plan: tracks, pre-registered acceptance criteria, gates | Before running a validation track, to see its criteria |
| [`phase1/README.md`](phase1/README.md) | How to run the Monte-Carlo harness | Running or extending the simulation study |
| [`phase1/FINDINGS.md`](phase1/FINDINGS.md) | Phase 1: the ERF bias confirmed, skew sweep, scale test | The original evidence for the whole approach |
| [`phase1/FINDINGS_v1.md`](phase1/FINDINGS_v1.md) | V1: the shipped engine validated against known truth | Citing bias/coverage for the censored path |
| [`phase1/FINDINGS_v2.md`](phase1/FINDINGS_v2.md) | V2: robustness sweep + the under-coverage finding | Checking whether a claim holds under some condition |
| [`phase1/FINDINGS_v4.md`](phase1/FINDINGS_v4.md) | V4: the under-coverage attributed to improper MI | Working on the Z-block fix |
| [`phase1/FINDINGS_v5.md`](phase1/FINDINGS_v5.md) | V5: the fix validated — partial; why the bootstrap under-delivers on a forest | Citing current interval behaviour, or picking an imputer |
| [`phase1/FINDINGS_v6.md`](phase1/FINDINGS_v6.md) | V6: forest+bootstrap vs `mice pmm` on a DGP that can punish both | Choosing an imputer, or reading the regime-dependence |
| [`phase1/FINDINGS_v7.md`](phase1/FINDINGS_v7.md) | V7: BART closes R8; the pilot-`m` curve closes R12 | Citing current interval behaviour, or the `m` recommendation |
| [`phase1/FINDINGS_v8.md`](phase1/FINDINGS_v8.md) | V8: the Z-block inner-iteration shortcut is unsafe with several targets; V7's port verified bit-identical | Touching block-FCS sweep counts, or citing V7 as evidence about the *shipped* code |
| [`phase1/FINDINGS_v3.md`](phase1/FINDINGS_v3.md) | V3: the mixture verdict on real BKMR estimands — curvature destroyed, pipeline worst of four arms | Before anyone points the censored path at a mixture analysis; citing the scope limit |
| [`phase1/FINDINGS_v9.md`](phase1/FINDINGS_v9.md) | V9: the mechanism test — replacing the linear conditional restores the estimands | Designing the SMC fix (item 07), or citing *why* the mixture path fails |
| [`phase1/FINDINGS_v10.md`](phase1/FINDINGS_v10.md) | V10: MAR covariate missingness — the additive claims survive their first non-MCAR test | Citing scope beyond MCAR, or before claiming anything about MNAR |
| [`phase1/FINDINGS_v11.md`](phase1/FINDINGS_v11.md) | V11: the imputer default tested against a non-linear outcome; the ≈3 pp scope limit | Before quoting additive-path bias figures, or revisiting `z_imputer` |
| [`phase1/FINDINGS_v12.md`](phase1/FINDINGS_v12.md) | V12: the Z-draw shape isolated — confirmed, and quantified at ~40% of the penalty | Before building the shape fix, or citing what it would buy |
| [`phase1/FINDINGS_v13.md`](phase1/FINDINGS_v13.md) | V13: the three-way decomposition — the exposure draw dominates and a Z-only fix backfires | Before building any single-component fix |
| [`phase1/FINDINGS_v14.md`](phase1/FINDINGS_v14.md) | V14: the shippable grid draw — 83% recovery, best calibration, and the missing-`Y` gap | Before building item 07 |
| [`phase1/FINDINGS_v15.md`](phase1/FINDINGS_v15.md) | V15: the attenuation law — quadratic in `f`, confirmed by prediction; sign inversion at high censoring | Before quoting how bad censoring is, or advising on a censored mixture analysis |
| [`phase1/FINDINGS_v16.md`](phase1/FINDINGS_v16.md) | V16: the root claim isolated — `−13.1 pp`, and the leak prediction refuted | Before citing Phase 1's H1, or claiming the blocks are separable |
| [`phase1/FINDINGS_v17.md`](phase1/FINDINGS_v17.md) | V17: the covariate's causal role decides the covariate block, and the shipped default is +5% to +10.4% biased under a confounder or mediator | **Before quoting any covariate-block result (V4, V5–V7, V10, V12, V13)** — all were measured with an inert covariate |
| [`phase1/FINDINGS_v18.md`](phase1/FINDINGS_v18.md) | V18: that bias decays as `n^−1/3` — permanent, but slower than V17 warned | Before quoting V17's coverage projection, or advising on large-`n` use |
| [`phase1/FINDINGS_v19.md`](phase1/FINDINGS_v19.md) | V19: only BART's bias decays; `forest_boot` is −23% asymptotic with coverage → 0 | **Before changing `z_imputer`, or when `dbarts` is missing** |
| [`phase1/FINDINGS_v20.md`](phase1/FINDINGS_v20.md) | V20: the bias is the **covariate draw** — 85–104% of it — not the exposure draw | **Before treating item 07 as the fix for the confounder bias** |
| [`phase1/FINDINGS_v21.md`](phase1/FINDINGS_v21.md) | V21: it is BART's **smoothing**; a correctly specified estimated draw is unbiased — but the fix is a trade, untested on the side that matters | **Before proposing a parametric Z block as the fix** |

**Runners** live in `phase1/`: `run_phase1.sh` (Phase 1), `run_v1_pipeline.sh` (V1),
`run_v2_robustness.sh` (V2), `run_v4_variance.sh` (V4, V5, V6, V7 and V8 — the arm/scenario
set is chosen by the `ARMS`/`SCENARIOS` environment variables), `run_inner_iter_check.sh` (V8), `run_v3_bkmr.sh` (V3). `run_status.sh` reports which
runs are alive. `results/` and `logs/` are
git-ignored — regenerate rather than commit them.

---

## Keeping this current

When a validation run completes or a finding changes:

1. **Write or update the findings file** for that track (`phase1/FINDINGS_*.md`) — table,
   reading against the pre-registered criteria, and an explicit caveats section naming
   what the run does *not* establish.
2. **Update the status tables above**, and the open-items list.
3. **Update [`ROADMAP.md`](ROADMAP.md)** if the remaining work or its ordering changed.
4. **Update the root `README.md`** if the *scope* of a validated claim moved — the README
   is where users read those claims, and it must not outrun the evidence.
5. **Correct, don't overwrite.** Where a later run overturns an earlier conclusion, mark
   the superseded passage as withdrawn and say why, leaving the original reasoning visible
   (see the withdrawn two-mechanism section in `FINDINGS_v2.md`). The record of how a
   conclusion was reached — and where it went wrong — is part of the evidence.

Acceptance criteria are written **before** a run, in the plan. If a criterion turns out to
be wrong, change it explicitly and say why; do not quietly re-read a result against a
softer bar.
