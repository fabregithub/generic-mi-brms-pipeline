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
Keep this file current: the validation study is the evidence base for the claims the
README makes about bias, coverage and scope, so a stale entry here silently becomes a
wrong claim in the documentation.

Last updated: **2026-08-31** · after V11 — the `bart` default survives a non-linear outcome; a new pipeline-wide scope limit

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

### Validation plan — Tracks V0–V11

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
| **08** | **Shape-aware Z-block draw** — every imputer draws covariates as *mean + homoscedastic Gaussian*, which is wrong under a non-linear outcome (≈3 pp, all arms). **Pending V12**; if confirmed, cheaper and wider-reaching than 07 | [`ROADMAP.md`](ROADMAP.md) |
| **07** | **Substantive-model-compatible imputation** — the mixture-path fix. **V9 validated the target**: the fix works and needs the right functional *form*, not oracle parameters. What remains is obtaining that form when the surface is unknown | [`phase1/FINDINGS_v9.md`](phase1/FINDINGS_v9.md) |
| — | ~~MAR covariate missingness~~ — **done (V10)**; the MCAR-only scope gap is closed | [`phase1/FINDINGS_v10.md`](phase1/FINDINGS_v10.md) |
| — | Interval overshoot (BART runs 3–8% wide) — low priority, conservative direction | [`INTEGRATION_SUMMARY.md`](INTEGRATION_SUMMARY.md) §2 |

Plus loose ends — the unexplained `resume` failure, stale translated READMEs, unswept MAR
covariate missingness, external testers — all listed in [`ROADMAP.md`](ROADMAP.md).

### ✅ Last completed run — V8, 2026-08-28 07:35 → 09:50 (135 min)

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

### ✅ Last completed run — V11, 2026-08-31 (10.5 h)

2,400 tasks, 300 reps in all 40 cells, zero errors.
[`phase1/FINDINGS_v11.md`](phase1/FINDINGS_v11.md)

`z_imputer = "bart"` was chosen in V6/V7 on a design where the outcome is **linear in
`(logX, Z)` by construction** — so a parametric imputer was correctly specified against the
outcome and could never be penalised for it. V7 named that confound and left it untested.

**The default survives.** `bartMI` is never significantly beaten in any non-linear-outcome
cell. Where the summary table shows another arm with smaller bias, the *paired* difference
is indistinguishable (*p* = 0.11–0.34) and in two of three cases favours `bartMI`. Those
arms only look better because they start at **+1.6% to +1.8%** bias and a uniform −3 pp
shift carries them *through* zero — cancellation, not robustness, and it would reverse if
the curvature had the opposite sign.

**The V7 hypothesis is not supported.** All four arms degrade by 2.86–3.30 pp; the spread is
0.44 pp. The confound was real but inert, and the V6/V7 verdict stands on evidence that
could have overturned it.

**The finding that matters is pipeline-wide.** A non-linear outcome costs ≈3 pp of bias
*whichever imputer is used* — `bartMI` goes from 0.53% to 3.10% mean absolute bias, outside
the ≤1% / ≤2.7% the README quotes for the additive path. **Coverage is unchanged**
(0.953 → 0.953), because at `n` = 800 a 3% relative bias is only ≈0.28 SE. Third
consecutive track where calibration diagnostics missed a bias finding (V3, V8, V11).

Control was bit-identical to V10 on all three shared cells.

---

---

## Document map

| Document | What it is | Read it when |
|---|---|---|
| [`INTEGRATION_SUMMARY.md`](INTEGRATION_SUMMARY.md) | **Start here.** Requirement → design → implementation → evidence traceability; open blocks with theoretical resolutions; the claims ledger | Orienting, writing methods text, or checking what may be claimed |
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
