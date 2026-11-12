# V36 — selection and censored-exposure imputation do NOT compound

**Run 2026-11-09 11:16 → 2026-11-10 05:19, five stages, `rc=0`. 18.05 h measured against
18.8 h projected from full-wave pilots (4% under).** Design and registered predictions:
`DESIGN_v36_selection.md`, fixed before any V36 code ran.

**The gate (P2) is REFUTED, and in the opposite direction.** Censored-exposure imputation does
not amplify outcome-dependent selection bias; where selection induces real collider bias it
slightly *offsets* it. The result that matters for the manuscript is not that, though — it is
that **selection bias alone destroys the analysis at cohort scale**, and the imputation
neither rescues nor materially worsens it.

## 1. P1 — instrument check: MET

Oracle arm, no censoring, vs the unselected-population target (3000 reps, `n` = 800):

| cell | selection bias | |
|---|---|---|
| `none` | +0.16% ± 0.18 | |
| **`descnull`** | **−0.08% ± 0.23** | falsification control — selection on a pure descendant of the exposure, no link to `Y`. **No bias, as required** |
| `desc` | −17.44% ± 0.23 | |
| `desc2` (20% dropped) | −27.40% ± 0.24 | monotone in the drop fraction |
| `outcome` | −13.10% ± 0.21 | |
| `outcome2` | −20.96% ± 0.22 | |
| `collider` | −32.32% ± 0.23 | |

The mechanism is wired to the intended path and to no other.

## 2. P2 — THE GATE: NOT MET, and reversed

Paired per replication against the **same cell's** oracle: `|arm − truth| − |oracle − truth|`,
pp of truth. Positive = the imputation makes it worse, which is what P2 predicted.

| arm | `descnull` | `desc` | `desc2` | `outcome` | `outcome2` | `collider` |
|---|---|---|---|---|---|---|
| `complete_case` | +4.82 **W** | −0.99 **B** | −4.95 **B** | +1.95 **W** | −0.60 **B** | −14.70 **B** |
| `sub_lod2` | +1.21 **W** | −4.29 **B** | −4.61 **B** | −3.23 **B** | −5.01 **B** | −8.91 **B** |
| **`pipeline_block_fcs`** | **+0.56 W** | **−0.68 B** | **−0.95 B** | **+0.29 W** | **+0.17 W** | **−3.34 B** |
| `smc_xgrid` | +0.49 **W** | −0.59 **B** | −0.81 **B** | +0.28 **W** | +0.21 **W** | −3.08 **B** |

SEs are 0.07–0.24 pp throughout, so every mark is ≥3σ. **In every cell carrying real collider
bias the shipped engine is significantly BETTER than its own oracle, not worse.**

**Where the offset appears is the informative part.** Signed shift, `pipeline − oracle`:

| cell | what selection depends on | shift |
|---|---|---|
| `none` | — | +0.01 ± 0.09 |
| `descnull` | the exposure only | −0.01 ± 0.07 |
| `outcome` | the outcome only | +0.13 ± 0.08 |
| `desc` | the exposure **and** a cause of the outcome | **+1.30 ± 0.08** |
| `collider` | the exposure **and** the outcome | **+3.46 ± 0.07** |

The offset appears **only when selection depends on both the exposure and the outcome (or a
cause of it)** — i.e. exactly when selection induces collider bias — and it is absent when
selection involves only one of them. P2's registered mechanism was that the imputation borrows
from a `Y` distorted by selection and so amplifies. The sign says otherwise; the mechanism as
registered is refuted. **A replacement mechanism is conjectured but NOT derived**: the censored
rows are the low-exposure rows, which are also the rows that selection on an increasing
function of the exposure preferentially removes, so imputing them restores part of the lost
mass. Treat that as a hypothesis, not a result.

## 3. P4 — the effect scales with the censored fraction (magnitude confirmed, direction refuted)

`pipeline − oracle`, `n` = 800:

| cell | `nd` = 0.10 | `nd` = 0.30 | `nd` = 0.50 |
|---|---|---|---|
| `desc` | +0.26 ± 0.04 | +1.30 ± 0.08 | **+3.04 ± 0.20** |
| `outcome` | −0.03 ± 0.04 | +0.13 ± 0.08 | −0.06 ± 0.21 |
| `none` | +0.06 ± 0.04 | +0.01 ± 0.09 | −0.46 ± 0.23 |

Strongly monotone in `nd` under `desc` and flat at zero elsewhere. P4's magnitude claim — the
effect is proportional to how much is imputed, so it really is the imputation doing it — is
**confirmed**; its direction claim shares P2's refutation.

## 4. P5 — the headline. Selection bias is flat in `n`, so coverage collapses

`desc` cell, 1500–3000 reps per level:

| `n` | oracle bias | oracle `bias/SE` | **oracle coverage** | pipeline coverage |
|---|---|---|---|---|
| 800 | −17.44% | −1.39 | 0.718 | 0.761 |
| 3200 | −17.06% | −2.70 | 0.221 | 0.332 |
| **12800** | **−17.12%** | **−5.45** | **0.001** | 0.005 |

The bias does not shrink; the SE does. **At `n` = 12800 a 9.2% exclusion on a descendant of the
exposure leaves coverage at 0.001** — one interval in a thousand contains the truth. The
`outcome` cell behaves the same way (0.800 → 0.370 → 0.003).

The offset itself is **flat in `n`** (+1.30 / +1.30 / +1.39), so it is structural rather than a
finite-sample artefact of imputing from few rows. P5's "does not vanish at cohort scale" is
**MET**; its framing as amplification is not.

This is V18's pattern in a new mechanism, and it is the fourth time this programme has found a
defect that hides at `n` = 800 and surfaces at cohort scale.

## 5. P3 — does the fix rescue selection bias?

`smc_xgrid` lands within 1.2 pp of `oracle_sel` under `desc`/`desc2` and within 0.1 pp under
`outcome`/`outcome2`, against selection biases of 13–27%. It does **not** rescue selection
bias, which is what P3 predicted; the small movement it does show is the same offset as §2, not
a correction. Scored cell by cell the gate reads NOT MET in 3 of 5 because the scorer tests
for *any* significant movement — at 3000 reps a 1.1 pp shift is significant. **Read §2, not the
gate's verdict: the fix addresses congeniality, not the sampling frame.**

## 5b. A by-product that matters more to the manuscript than the gate

The `none` cells carry an incidental measurement nobody had made: **what MRL/2 substitution
does on the additive path with a correctly specified analysis model.** Every prior substitution
result in this programme is from the mixture/BKMR path, where it attenuates curvature.

| non-detect fraction | `sub_lod2` bias | its coverage | shipped engine | its coverage |
|---|---|---|---|---|
| 0.10 | **+1.92%** ± 0.23 | 0.952 | +0.27% ± 0.23 | 0.955 |
| 0.30 | **+8.56%** ± 0.21 | 0.880 | +0.18% ± 0.21 | 0.954 |
| 0.50 | **+18.60%** ± 0.31 | **0.729** | −0.26% ± 0.32 | 0.955 |

**Substitution INFLATES the coefficient, monotonically**, log-log slope in the non-detect
fraction **1.48 [0.76, 2.20]**, while the shipped engine is unbiased across the whole range with
coverage at nominal. Unlike the point estimate, **substitution's coverage failure is
unambiguous**: 0.729 at a 50% non-detect rate.

⚠️ **Scope, and it is a tight one.** The LOD here is a **single** value set at a quantile of
the exposure, as in every cell this programme has run. Where reporting limits vary by
participant — as in the JECS assay — the substituted value carries variation unrelated to the
true exposure, which is classical measurement error and attenuates instead. **The direction
measured here must not be assumed to transfer to a subject-specific-limit study**; the
manuscript's §6.5.4 registers both regimes and names a discriminator computed from the limits
alone.

## 6. What may be asserted

- ✅ "Selection on a descendant of the exposure that shares a cause with the outcome attenuates
  the exposure–response by **~17%** at a 9.2% exclusion rate, **flat in sample size**, so
  coverage falls from 0.72 at `n` = 800 to **0.001** at `n` = 12800."
- ✅ "Censored-exposure imputation and outcome-dependent selection **do not compound**. Where
  selection induces collider bias the imputation partially offsets it (0.6–3.3 pp against
  13–32 pp); where it does not, the imputation adds only its own standalone cost
  (+0.2 to +0.6 pp)."
- ✅ "The offset scales with the censored fraction (+0.26 → +1.30 → +3.04 pp at `nd` = 0.10 /
  0.30 / 0.50) and is flat in `n`."
- ⚠️ **"Do not compound" is not "safe".** The offset is a partial cancellation of two errors,
  the same shape as V35's `f` = 0.2 cell. The residual at `n` = 12800 is still −15.7% with
  coverage 0.005.
- ❌ Not "imputation corrects selection bias" — it moves ≤3.5 pp of a 17–32 pp bias, and the
  mechanism is conjectured, not derived.
- ❌ Not "the `sub_lod2` arm is a clean mechanism control" — see §7.

## 7. Defects in this design, found by running it

**The mechanism control does not work.** `sub_lod2` was registered as a null control: it
imputes a constant, borrows nothing from `Y`, so under a `Y`-borrowing mechanism it should show
no interaction. It shows large significant ones (+1.21 under `descnull`, −4.29 under `desc`).
The reason is a confound the design missed: **selection on a descendant of the exposure changes
the realised non-detect fraction** (0.300 → 0.240 under `descnull`, 0.261 under `desc`),
because it removes low-exposure rows, and substitution's bias is a direct function of how many
rows are substituted. So the arm moves for a reason that has nothing to do with `Y`. Fixing it
means holding the *realised* non-detect fraction constant across cells by setting the LOD
within the selected sample — which is a different, and arguably less realistic, design. **Do
not cite `sub_lod2` here as evidence about the borrowing channel.**

**The gate's sign convention was wrong and was fixed after the fact.** `analyze_v36.R`
originally reported the signed interaction and labelled positive values "amplification".
Selection bias is negative in every live cell, so a positive interaction means *nearer* the
truth. The gate is now stated on `|bias|`, paired per replication, which is sign-safe. The
registered prediction is unchanged; only the scorer's arithmetic was wrong, and the corrected
reading is the one above.

**`descnull` is not a perfect null for the imputation arms** (+0.56 pp, ≥3σ) for the same
non-detect-fraction reason. It remains a valid null for the *oracle*, which is what P1 needs.

## 8. Files

`DESIGN_v36_selection.md` (registered predictions), `R/selection.R` (cells + the additive-path
LOD/√2 arm), `run_v36_selection.R` (driver), `run_v36_selection.sh` (five-stage runner),
`analyze_v36.R` (three-state gates), `results/v36_n{800,3200,12800}_nd{0.10,0.30,0.50}_*`.
