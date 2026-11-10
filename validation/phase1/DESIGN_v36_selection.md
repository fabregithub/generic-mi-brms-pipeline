# V36 — outcome-dependent selection: does censored-exposure imputation amplify it?

**Status: RUN AND SCORED 2026-11-10 — see `FINDINGS_v36.md`.** P1 MET; **P2 NOT MET and
REVERSED** (imputation offsets selection bias, it does not amplify it); P4 confirmed in
magnitude and refuted in direction; P5 MET. Two design defects surfaced by running it — the
`sub_lod2` mechanism control is confounded by the realised non-detect fraction, and the
scorer's sign convention was wrong — both documented in `FINDINGS_v36.md` §7.

**Status at design time: BUILT; PILOT RUN 2026-11-09.** Registered predictions
below were fixed *before* any V36 code ran. Two things changed after the pilot and both are
dated and explained in place: the **strength lever** (§3) and the **measured cost** (§7).
Neither touches a prediction.

**MEASURED COST (full-wave pilot, 2026-11-09, tag `v36pilot`): 44 tasks on 22 cores in
0.6 min → 16.9 s/task → 1.36 min per 100 reps per cell.** The whole seven-cell design at
100 reps/cell is ~9.5 min; at 1000 reps/cell ~95 min. This is two orders of magnitude cheaper
than the BKMR tracks, which is what choosing an additive ERF bought.

**P1 IS ALREADY MET** (§5): oracle arm, n = 4000, 400 reps — the falsification control
`descnull` reads **−0.20% ± 0.30** (no bias) while `desc` reads −17.2%, `outcome` −13.4% and
`collider` −32.8%.

## 1. Why this track exists

It is the last gap blocking the manuscript's demonstration. `DRAFT.md` §6.3 prices three
deviations of Itoh et al. from Mezawa et al.'s DAG; two are derivable from the existing record
and the third is not:

| deviation | status |
|---|---|
| adjusted for a mediator (breastfeeding) | derivable |
| **excluded on a mediator** (low birth weight, 1,534 of 16,709 ≈ **9.2%**) | **no harness mechanism — this track** |
| MRL/2 substitution vs imputation | derived |

Nothing in the harness removes rows as a function of anything. Every track to date varies
*which covariates are conditioned on* (`z_role`: precision / fork / pipe / mixed / collider),
never *which rows survive*. Selection on a descendant of the exposure is a different defect
class and the project currently cannot say whether it is worse, better, or merely additive
when combined with censored-exposure imputation.

**The question is specifically the interaction, not selection bias on its own.** Collider
bias from selecting on a mediator is textbook and we would not publish a re-derivation. What
is unknown — and what this programme is uniquely placed to measure — is whether imputing a
censored exposure **amplifies** it. There is a mechanism to suspect it does: `THEORY.md` §3.1
splits the bias into attenuation (`T₁`) and **outcome borrowing** (`T₂`). The imputation draws
the censored exposure conditioning on `Y`. Under outcome-dependent selection the `Y`
distribution in the analysed sample is itself distorted, so `T₂` borrows from a corrupted
source. That predicts amplification, and predicts it is absent when selection is independent
of `Y`.

## 2. Estimand, and the target population — read this before arguing about bias

The estimand is `b1`, the focal censored exposure's coefficient on the log scale, **in the
unselected population**. Bias is measured against that target, not against the selected
sample's own coefficient.

This is deliberate and it is the opposite of `DRAFT.md` §6.2's C1/C2 rule, where changing the
adjustment set changes the *estimand* and must be reported separately. The difference: Itoh's
paper interprets its coefficient as the exposure–response in the cohort, not in the
LBW-excluded subpopulation. A reader is not told the estimand changed. So for this track
selection is a **bias against the advertised target**, and the design must say so out loud in
the findings, because a reviewer can reasonably take the other view.

## 3. Cells

Additive ERF (`erf_form = "additive"`, `y_form = "linear"`), so the analysis model is
correctly specified and every deviation is attributable to censoring or selection. **This also
keeps the track cheap** — no BKMR, unlike V3/V15/V33/V35.

| cell | mechanism | mirrors |
|---|---|---|
| **S0** | no selection — control | every track to date |
| **S1** | **selection on a mediator.** `M = α·logX1 + δ·U + ν`; drop rows with `M` below a threshold set to exclude ≈9.2%. `U` is an unmeasured common cause of `M` and `Y`, so selecting on `M` opens `X1 → M ← U → Y` | **Itoh's LBW exclusion** |
| **S2** | **selection on the outcome.** `P(S=1 \| Y)` logistic in `Y`, same ≈9.2% excluded | the literal outcome-dependent case; the sharpest test of the `T₂` mechanism |
| **S3** | **live-birth collider.** `S ← X1` and `S ← U_Y`, conditioned on by design (all arms run within `S = 1`) | **Mezawa's stated collider** |
| **S1-null** | S1's structure with `δ = 0`, so `M` is a pure descendant of `X1` with no link to `Y` | **falsification control** — selection here should NOT bias `b1`, and must not |

Each live cell runs at **two selection strengths**, because a single strength cannot
distinguish a real response from an offset.

> ⚠️ **Corrected 2026-11-09 after the pilot: the strength lever is the DROP FRACTION, not
> `δ`.** The design originally said "doubled `δ`". Measured on the oracle arm at n = 4000,
> 400 reps: `δ` = 0.80 → **−17.5%**, `δ` = 1.60 → **−12.8%**. Collider bias here is
> **non-monotone in `δ`** — at `δ` = 0 selection is on a pure descendant of the exposure and
> there is no bias (that is `descnull`), and as `δ` grows selection tends to selection on `U`
> alone, which is independent of the exposure, so the bias vanishes again. It peaks in
> between. Anyone "strengthening" a cell by raising `δ` would have measured a *smaller* bias
> and read it as a refutation. The drop fraction is monotone (0.092 → **−17.2%**, 0.20 →
> **−26.9%**) and is also the quantity an epidemiologist observes. Encoded in
> `R/selection.R`; the cells are `desc`/`desc2` and `outcome`/`outcome2`.

## 4. Arms

| arm | what it isolates |
|---|---|
| `oracle_full` | no censoring, **no selection** — defines the target `b1` |
| `oracle_sel` | no censoring, **with selection** — **pure selection bias** |
| `complete_case` | reference |
| `sub_lod2` | MRL/2 substitution — **Itoh's actual method**, and a built-in **mechanism control**: it imputes a constant and borrows nothing from `Y`, so if the amplification in P2 is the `T₂` borrowing channel, this arm must show **no** interaction. Pilot: **−0.56 ± 1.55** pp, consistent with zero |
| `pipeline` | the shipped block-FCS engine (`proc_pipeline_block_fcs`) |
| `smc_dec_true` | the V33 declared-basis fix |

The decomposition that matters, per cell:

```
selection bias          = oracle_sel      - target
imputation bias         = pipeline_S0     - oracle_full        (known from V17/V18)
imputation under sel.   = pipeline_sel    - oracle_sel
INTERACTION  (the point) = (pipeline_sel - oracle_sel) - (pipeline_S0 - oracle_full)
```

The interaction is the quantity this track exists to measure. Report it with MC error and with
direction (`PLAN` §11b), per cell and per arm.

## 5. Registered predictions

Scored MET / NOT MET / **UNDECIDED**, the third a real state. All comparisons carry MC error;
"significant" means the whole CI excludes the null value.

**P1 — instrument check (must pass or the track is void).**
Under S1, `oracle_sel` is biased against the target, and `S1-null` is not.
*Falsified if* `oracle_sel` under S1 is indistinguishable from the target (the selection
mechanism is not wired up), **or** if `S1-null` shows bias (selection is leaking through a
path it should not have). Either outcome invalidates every other reading and must be fixed
before the rest is interpreted.

**P2 — THE GATE. Censored-exposure imputation amplifies outcome-dependent selection bias.**
The interaction defined in §4 is **positive and significant** under S2 and S1, and
**absent under S1-null**.
*Mechanism:* the imputation borrows from `Y` (`T₂`); selection distorts `Y`; so the borrowed
information is systematically wrong. Under S1-null nothing distorts `Y`, so there is nothing
to amplify.
*Falsified if* the interaction is zero within MC error (the two defects are merely **additive**
— a perfectly publishable and simpler result), or negative (imputation **damps** selection
bias, which would be surprising and would need its own explanation before being believed).
*Direction registered:* **away from the null**, matching the pattern in 75 of 118 shipped-arm
cells (`PLAN` §11b). A toward-the-null interaction is a NOT MET even if significant.

**P3 — the declared-basis fix does not rescue selection.**
Under S1 and S2, `smc_dec_true` removes the *imputation* component but leaves the *selection*
component: it lands near `oracle_sel`, not near the target.
*Why registered:* it **bounds what the manuscript may claim for the fix.** V33 earns "declared
basis is the mixture fix"; nothing earns "it fixes a bad sampling frame", and the two will be
conflated by a reader unless we measure it.
*Falsified if* `smc_dec_true` is significantly closer to the target than `oracle_sel` is —
which would be a genuine and unexpected finding, not a failure.

**P4 — the amplification GROWS with the censored fraction.**
*Registered 2026-11-09, before any V36 run at a censored fraction other than 0.30.*
The interaction under `desc` and `outcome` is **larger at `nd` = 0.50 than at `nd` = 0.10**.
*Mechanism:* if the amplification is the `T₂` borrowing channel, then what matters is the
share of the censored exposure's information that comes from `Y` rather than from an observed
value — and that share grows with the censored fraction. A distorted `Y` then contaminates
more of the draw.
*Falsified if* the interaction is flat in `nd`, or decreases.
*Paired control:* `sub_lod2` borrows nothing from `Y` at any censored fraction, so its
interaction must stay at zero across the whole `nd` sweep. If it instead grows with `nd`, the
mechanism named here is wrong and the amplification is something else.

**P5 — the amplification does not vanish at cohort scale.**
*Registered 2026-11-09, before any V36 run at `n` ≠ 800.*
The interaction under `desc` and `outcome` is **still significantly positive at `n` = 12800**,
and its `bias/SE` is **larger** there than at `n` = 800.
*Why it matters:* `DRAFT.md` §1 claim 3 is that the cost does not shrink as the cohort grows.
JECS is `n` = 15,175. If the interaction is an `n` = 800 artefact the manuscript must not cite
it; if it persists while the SE falls, coverage collapses at exactly the scale the paper is
about — the V18 pattern, in a new mechanism.
*Falsified if* the interaction shrinks toward zero in `n` (a finite-sample effect of imputing
from few rows, not a structural one).

## 5b. Run plan and measured cost

Full-wave pilots, 2026-11-09, 22 cores, m = 30, all five arms, at the real configuration:

| `n` | s/task | min per wave of 22 |
|---|---|---|
| 800 | 16.5 | 0.30 |
| 3200 | 35.0 | 0.60 |
| 12800 | 112.8 | 1.90 |

| stage | `n` | `nd` | reps | cells | gate | measured cost |
|---|---|---|---|---|---|---|
| 1 | 800 | 0.30 | 3000 | all seven | P1, P2, P3 | 4.8 h |
| 2 | 3200 | 0.30 | 2000 | none, descnull, desc, outcome | P5 | 3.6 h |
| 3 | 12800 | 0.30 | 1500 | none, descnull, desc, outcome | P5 | 8.6 h |
| 4 | 800 | 0.10 | 2000 | none, desc, outcome | P4 | 0.9 h |
| 5 | 800 | 0.50 | 2000 | none, desc, outcome | P4 | 0.9 h |

**Total ≈ 18.8 h** against a 44 h window (2026-11-09 13:00 → 2026-11-11 09:00), i.e. ~2.3×
margin. Reps are sized from the pilot's paired SDs: at 3000 reps the interaction SE on
`pipeline_block_fcs` is ≈0.11 pp against a pilot point estimate of 1.11 pp. Going beyond that
buys resolution the question does not need, and the roadmap's rule — a cell must be justified
by a choice it would change — rules out filling the remaining window for its own sake.

## 6. Implementation

| file | what |
|---|---|
| `R/selection.R` | **new.** `apply_selection(bundle, spec)` → selected bundle + the retained target truth. Generates `M`/`U`/`S` from the existing bundle rather than editing `R/dgp.R`, so no existing track changes behaviour |
| `run_v36_selection.sh` | **new.** Must follow the current runner contract: `setsid(2)` detach (NOT `nohup &`), absolute `$SELF` after any `cd`, BLAS pinning, `logs/<tag>_<stamp>.pid` so `run_status.sh` can see it |
| `analyze_v36.R` | **new.** Three-state gates; reports the §4 decomposition with direction |

Reuse `proc_oracle`, `proc_complete_case`, `proc_pipeline_block_fcs` and the V33 SMC arm
unchanged — selection is applied to the bundle *before* any arm sees it, so no procedure code
needs to know selection exists.

## 7. Cost — measure it, do not estimate it

**Run a timed FULL-WAVE pilot first** (reps = worker count, 22) and record the measured
per-wave cost in this file's header before scheduling the run. This project's extrapolations
have been wrong by 3×, 8× and 13×, and a 3-rep pilot gave 1.8 h and 5.0 h for the same 3.1 h
run. Do not quote a duration from anything smaller than a full wave.

Expected to be **much cheaper than the BKMR tracks** — the analysis model is `lm`, so the cost
is the block-FCS imputation, as in V18 — but that expectation is exactly the kind of reasoning
the rule above exists to override.

## 8. What this track will NOT answer

- **The real JECS magnitude.** `α`, `δ` and the mediated share are simulated, not estimated
  from the cohort. The track prices the *mechanism* and its interaction, not Itoh's actual
  number. Any manuscript sentence attaching a percentage to Itoh's exclusion needs a separate
  argument.
- **MNAR**, which remains out of scope by the 2026-08-29 decision.
- **Whether LBW is in fact a mediator.** That is Mezawa's DAG's claim and we adopt it; this
  track asks what the exclusion costs *given* that structure.
- **Non-additive surfaces.** The additive ERF is chosen for cost and attributability; whether
  the interaction behaves the same on a mixture surface is a follow-on, not this track.
