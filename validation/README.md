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

Last updated: **2026-08-29** · after V3 / R11 — the last substantive track. Mixture verdict final

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

### Validation plan — Tracks V0–V8 (V3 last to close)

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
| **07** | **Substantive-model-compatible imputation** — the only remaining fix for the mixture path, and now the sole substantive track | [`phase1/FINDINGS_v3.md`](phase1/FINDINGS_v3.md) |
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

### ✅ Last completed run — V3 / R11, 2026-08-28 14:16 → 2026-08-29 02:50 (12.6 h)

600 tasks, 200 reps in all 84 cells, zero errors. Written up in
[`phase1/FINDINGS_v3.md`](phase1/FINDINGS_v3.md).

**The mixture verdict is confirmed and is worse than the scaffold suggested.** On the
estimands a BKMR analysis actually reports, the shipped censored-exposure engine is the
**worst of the four arms** — mean 14.6% absolute paired excess over the oracle, against
6.4% for LOD/√2 substitution and 10.5% for complete case.

| | `curv_X1` (pure curvature) | `int_X1X2` (pure interaction) |
|---|---|---|
| `nd20` | **−11.7%** | +5.7% |
| `nd40` | **−56.7%** | +2.6% |
| `nd40_all` | **−58.4%** | **−12.5%** |

*(paired excess over the oracle, 200 reps, all p < 1e-4)*

Three things worth carrying forward:

1. **It is a curvature failure, not a blanket mixture failure.** Interaction estimands
   survive single-exposure censoring (+2.6% to +5.7%) and break only when the second
   exposure is censored too. Curvature is destroyed at 40% non-detects.
2. **Imputing with the wrong functional form is worse than not imputing.** On curvature,
   LOD/√2 substitution is statistically indistinguishable from the oracle (p = 0.42, 0.47)
   while the congenial engine loses 57%. The X block's conditional is *linear in the
   predictors*, so it imposes linearity on exactly the rows where curvature would show.
3. **Coverage was blind to it** — 0.945–0.960 alongside a −42% point bias, because the
   pooled intervals run 1.5–1.6× wider than the true sampling spread. Second track running
   to find this (see V8), and the lesson is now explicit: *coverage is not a substitute for
   a bias check against known truth.*

**The gate failed as written**, on `curv_X1` (oracle +15–17% against a 10% bar). The bar was
mis-set on a 44-rep calibration carrying ±5.3% MC error — reps 1–44 of this run reproduce
that calibration exactly, so nothing drifted; the calibration was simply underpowered.
The run remains interpretable because the headline is *paired excess over the oracle*, in
which the floor cancels — a design decision taken before the run. Full diagnosis in the
findings.

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
