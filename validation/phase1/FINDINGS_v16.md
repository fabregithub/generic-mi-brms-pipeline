# V16 findings — the root claim, isolated at last: **−13.1 pp**, and the leak prediction fails

**Run:** 2026-09-03 10:42 → 12:28 JST (**106.0 min**, stage 1 of `run_v16_v17.sh`; 6 cells ×
4 arms + oracle × 500 reps = 3,000 tasks, `n_ok` = 500 in all 6 cells, zero task errors,
zero non-finite estimates). Results: `results/v16_latest.rds`. Verdict arithmetic:
[`analyze_v16_v17.R`](analyze_v16_v17.R). Criteria: `../PLAN_pipeline_validation.md` §8i.

**The question.** Every track from V1 to V15 tests a *refinement* of `THEORY.md`'s
congeniality condition — whether the imputation conditional has the right functional form,
shape or margin. **None varies the condition's own antecedent:** every arm in fifteen tracks
had `Y` in the imputation model. Phase 1's H1 is the cited evidence, but its no-`Y` arm also
drops the `Z` covariates *and* uses a different estimator, and two cells in that same
document have the no-`Y` arm **beating** the `Y`-aware one. So the root claim had never been
isolated inside a single estimator.

**Answer: it holds, and the size was predicted.** Removing `Y` from both blocks of the
shipped engine costs **−13.1 pp** against a derived **−14.6**; from the exposure block alone,
**−14.4** against **−15.8**. Both inside the registered ±3 pp. **The leak prediction is
refuted**, and one bar was registered against cells it was never derived for.

---

## Measured against predicted

Paired shift against `pipeline_bartMI` (`Y` in both blocks), within replication, in
percentage points of `b₁ = 0.40`. MC error ≈ 0.13–0.16 pp.

| cell | design | `noYx` | `noYz` | `noYboth` | leak |
|---|---|---|---|---|---|
| **`mcar_z40`** | X1 censored 40%, Z MCAR 40% | **−14.39** | +1.17 | **−13.14** | +0.09 |
| `nl_mcar_z40` | + non-linear Z1 | −14.68 | +1.08 | −13.52 | +0.08 |
| `mar_z40` | Z **MAR** on `Y` 40% | −14.09 | +0.93 | −13.05 | +0.11 |
| `nl_mar_z40` | both | −14.55 | +1.07 | −13.31 | +0.17 |
| `combined` | **all 3** censored, Z 20%, Y 20% | −10.53 | +0.52 | −10.15 | −0.14 |
| `nl_combined` | both | −10.45 | −0.46 | −10.47 | +0.45 |
| *predicted* | | *−15.8* | *≈ 0* | *−14.6* | *+1.6* |

**The root claim is confirmed in the configuration a real analysis produces.** `noYboth` is
the realistic arm — an analyst who omits `Y` omits it everywhere — and it lands at −13.1 pp
in the cell the prediction was derived for, inside the ±3 pp bar. `noYx` at −14.4 also
passes. The `Y`-omission penalty is real, is about 13% of the effect, and is now attributable
to `Y` alone rather than to `Y` plus a covariate set plus an estimator change.

**`noYz` ≈ 0 everywhere (+1.17 to −0.46).** Omitting `Y` from the *covariate* block costs
essentially nothing here. That is not a null result but the content of the finding, and V17
explains it: in all six of these cells `Z` is causally inert with respect to the exposure, so
there is no path along which a mis-imputed covariate can reach `β₁`. See
[FINDINGS_v17.md](FINDINGS_v17.md), where the same arm moves **+25 pp**.

**MAR makes no difference.** `mar_z40` (−14.09 / +0.93 / −13.05) is within MC error of
`mcar_z40` (−14.39 / +1.17 / −13.14), and the same holds for the non-linear pair. Missingness
in the MAR cells depends on `Y`, so omitting `Y` from the Z block genuinely invalidates that
imputation — and `β₁` still does not move. Consistent with V10, and now with the no-`Y` arms
V10 lacked.

---

## The leak prediction is refuted

Registered: `noYboth − (noYx + noYz)` should be **+1.4 to +1.6 pp and positive**, meaning the
*mixed* configuration — `Y` in the covariate block but not the exposure block — is **worse**
than omitting `Y` from both. The mechanism offered was that a `Y`-informed `Z̃` partly absorbs
`Y`, so adjusting for it alongside a `Y`-blind `X̃` over-adjusts.

**Measured: +0.09 ± 0.13 (`mcar_z40`), −0.14 ± 0.14 (`combined`).** Zero, at 500 reps. The
two blocks are **additive** in this structure, and *"partial compliance is worse than none"*
is **false** here.

**Why the derivation got it wrong.** `predict_roles.R` uses a correctly-specified *linear*
Gaussian draw for `Z`; the pipeline uses BART. The proposed over-adjustment depends on the
precise `Y`-loading of the covariate draw, and the two do not match on that. This is a
fidelity limit of the stand-in, not of the theory — and it is the first time the stand-in has
been caught being wrong about a *sign*, which is worth remembering the next time it is used
to register a prediction.

**The leak is real, just not here.** In V17's cells it is clearly non-zero and **negative**
(−2.9 pp in fork/pipe/mixed, −6.2 to −6.4 in the collider cells, all at ≈0.12–0.17 MC error).
So block non-additivity exists and matters; it simply vanishes when the covariate is off any
X–Y path. The general claim survives in a corrected form: **the blocks interact only when the
covariate block has a route to the estimand.**

---

## One bar was registered against the wrong cells — my error, not the theory's

`combined` and `nl_combined` missed the `noYx` and `noYboth` bars (−10.5 and −10.2 against
−15.8 and −14.6). They are not counter-evidence: the prediction was derived at **focal-only
censoring with 40% covariate missingness**, which is `mcar_z40`'s design. `combined` censors
**all three** exposures with 20% covariate and 20% outcome missingness — a cell shape the
derivation never modelled. I applied one number to six cells when it was computed for one.

The direction is also explicable: with 20% rather than 40% covariate missingness there is
less imputation to be uncongenial about, so a smaller penalty (−10.5) is what one should
expect. **Registering it as a failure of the prediction would be wrong; so would quietly
re-reading it against a softer bar.** The fix for next time is to derive per cell shape, or
to register the bar only for the cells it covers.

---

## What this revises in the theory

`THEORY.md` §3b recorded the root condition as *stated but never tested on its own*. It is
now tested: **the antecedent holds, at −13.1 pp, inside one estimator, with `Y`'s
contribution separated from the covariate set and the method.** Phase 1's H1 number
(−14.4% at 40% non-detects) turns out to have been dominated by the `Y` omission after all,
so its confounded design did not mislead — but that could only be known by isolating it.

Two things do **not** survive: the positive leak, withdrawn above; and the loose form of the
claim, *"include `Y` or you attenuate"*, which V16 leaves standing only for the **exposure**
block. For the covariate block, whether `Y` matters at all depends on the covariate's causal
role — which is V17.

---

## Caveats (scope of this run)

- **Six cells, one DGP, `n` = 800, `m` = 30, additive surface.** The −13.1 pp figure is for
  40% non-detects on the focal exposure with 40% covariate missingness.
- **`Y` is fully observed in four of six cells.** Where `Y` itself is missing (`combined`,
  `nl_combined`), MID deletes those rows before the fit, so the no-`Y` arms there are also
  a smaller-`m`-of-complete-cases comparison.
- **The prediction was derived at `m` = 10 and the run is `m` = 30.** Declared in §8i before
  the run. The passing cells passed comfortably, so this did not decide anything.
- **`noYx`/`noYz` are not clean marginal effects** — the blocks alternate, so `Y` reaches one
  block through the other. That was known and registered; the measured leak says the
  contamination is nil in this structure.
- **Nothing here licenses dropping the covariate block's `Y`.** `noYz` ≈ 0 in these cells and
  **+25 pp** in V17's. The safe reading of both is the one the pipeline already implements:
  `Y` in both blocks.
