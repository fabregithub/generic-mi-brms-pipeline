# Validation — status index

**Single source of truth for what has been validated, what is open, and what is next.**
Keep this file current: the validation study is the evidence base for the claims the
README makes about bias, coverage and scope, so a stale entry here silently becomes a
wrong claim in the documentation.

Last updated: **2026-08-26** · after the v1.4.0 imputation fix

---

## Current status

### Design plan — Phases 1–5

Tracks the *build* of the censored-exposure feature.
Document: [`PLAN_leftcensored_exposure_integration.md`](PLAN_leftcensored_exposure_integration.md)

| Phase | What | Status |
|---|---|---|
| 1 | Confirm the ERF bias; test the congenial fix | ✅ done 2026-08-13 |
| 2 | Decision gate — build or not | ✅ resolved: build linear for additive; mixtures need SMC |
| 2b | Estimand extension (BKMR estimands) | ⬜ **deferred** → roadmap 04 |
| 3 | Build the §4 component in `leftcens` | ✅ core shipped in `leftcens` 0.9.0 |
| 4 | Wire block-FCS into the pipeline | ✅ built + hardened |
| 5 | Pilot `m`, document FMI | ⬜ **deferred** → roadmap 03 |

### Validation plan — Tracks V0–V6

Tracks the *validation of the shipped code*.
Document: [`PLAN_pipeline_validation.md`](PLAN_pipeline_validation.md)

| Track | Question | Status | Headline |
|---|---|---|---|
| V0 | Do Steps 5–11 + Quarto work on the censored path? | ✅ passed | 78 s, 5/5 criteria |
| V1 | Is the *shipped* engine unbiased with nominal coverage? | ✅ passed | bias ≤0.8%, coverage 0.957 |
| V2 | Does that hold outside the tested corner? | ✅ 9/10 + 1 finding | max bias 2.71%; one under-coverage cell |
| V3 | Does the mixture verdict hold on real BKMR estimands? | ⬜ **not started** | needs `bkmr`; analytic estimand truth is the real work |
| V4 | Where does the V2 under-coverage come from? | ✅ resolved | improper MI in the Z block understates `B` by 8–59% |
| V5 | Does the candidate fix work? | ✅ partial | `proper_draw` default TRUE in v1.4.0; 66% of gap closed at 40% MCAR, 23% under heavy combined missingness |
| V6 | Forest+bootstrap or parametric+posterior? | ✅ resolved | `properBoot` stays: `mice pmm` bias −4.70% under non-linear covariates. Neither uniformly best |

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
| 06 | **Flexible + properly dispersed imputer** — carries the residual gap; BART hypothesis | [`INTEGRATION_SUMMARY.md`](INTEGRATION_SUMMARY.md) §2 |
| 03 | Pilot `m` / FMI *(= design-plan Phase 5)* | [`ROADMAP.md`](ROADMAP.md) |
| 04 | BKMR estimands *(= design-plan Phase 2b, = Track V3)* | [`ROADMAP.md`](ROADMAP.md) |

Plus loose ends — the unexplained `resume` failure, stale translated READMEs, unswept MAR
covariate missingness, external testers — all listed in [`ROADMAP.md`](ROADMAP.md).

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

**Runners** live in `phase1/`: `run_phase1.sh` (Phase 1), `run_v1_pipeline.sh` (V1),
`run_v2_robustness.sh` (V2), `run_v4_variance.sh` (V4 and V5). `results/` and `logs/` are
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
