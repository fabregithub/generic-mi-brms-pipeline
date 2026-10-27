# V35 — the `f²` exponent belongs to the ANALYSIS MODEL, not the imputation

**Run 2026-10-26/27, completed 2026-10-27 17:05.** Status: **all three registered predictions
MET, on the paired design as registered.** 4 arms × 6 censoring fractions × 100 reps = 2400
fits, no cell short.

**This was a RESUMED run** (`README.md` requires saying so). The first launch died at 88/600
tasks on 2026-10-26 (see §4); `RESUME=1` carried those 88 `nd20` rows forward and ran the
remaining 512 on 2026-10-27, 07:29–17:05. Same machine, same R and package versions — nothing
was installed between the two segments — and the per-task seeds are a function of
`(scenario, rep)`, so the resumed tasks drew what they would have drawn in one continuous run.

## 1. The question

Two independent routes in this project put the censored fraction into the bias as `f²`:
V15's `excess% = −323.4·f²` (2026-09-02) and V34's `f^1.85 ± 0.20` (2026-10-24). Four
mechanisms had been proposed and refuted — constant `Γ`, `Γ` growing in `f`, direct
propagation of the per-row imputation error, and V15's own evaluation-range conjecture
(retired without a run; see `../THEORY.md`).

P-B (`derive_declared_basis.R`, 2026-10-26) localised the residual: a replica matching V33's
imputation in every respect checked against the source, differing **only** in using OLS on the
true basis instead of BKMR, agreed with the asymptotic fixed point to ≤0.025 in all twelve
cells and still gave `f`-exponents of 0.83–1.33 against a measured 1.46–2.28. V35 tests that
on V15's own ground, with both analysis models run on **the same imputed datasets**.

**The OLS arm is an instrument, not a rival.** It is handed the generator's own basis, so it
is correctly specified and efficient by construction — which BKMR is never asked to be. Its
response to information loss is known exactly, so whatever BKMR adds shows up as the
difference. **Nothing here recommends OLS over BKMR for real data**, where the basis is not
known; that is the entire point of the declared-basis work.

## 2. The result

`curv_X1`, excess over each analysis model's **own** oracle, paired on identical datasets,
100 reps per cell:

| `f` | OLS | BKMR | difference | σ |
|---|---|---|---|---|
| 0.10 | −15.8 ± 2.7 | −1.2 ± 4.6 | +14.6 | 2.7 |
| 0.20 | −34.3 ± 2.9 | −11.9 ± 5.0 | +22.4 | 3.9 |
| 0.30 | −45.0 ± 2.6 | −25.5 ± 4.0 | +19.5 | 4.1 |
| 0.40 | −56.2 ± 3.1 | −49.3 ± 5.2 | +6.9 | 1.2 |
| 0.50 | −67.5 ± 2.7 | −78.5 ± 4.9 | −11.0 | 2.0 |
| 0.60 | −77.9 ± 2.7 | −101.7 ± 5.9 | −23.8 | 3.7 |

**Log-log slope in `f` (`f` ≥ 0.2, matching V15's fit):**

> **OLS 0.75 [0.69, 0.81]  ·  BKMR 2.01 [1.80, 2.21]  — disjoint, gap 1.26**

Over all six points: OLS **0.87 [0.75, 0.99]**, BKMR **2.47 [2.00, 2.95]** — still disjoint.
The OLS interval **contains `derive_declared_basis.R` STAGE 4's registered prediction of
0.85**, and the BKMR interval contains V15's 2.

**The curves cross at `f` = 0.44** and diverge in both directions, 2.7σ one way at `f` = 0.10
and 3.7σ the other at `f` = 0.60. Below the crossing BKMR reports far *less* damage than OLS;
above it, far *more* — at `f` = 0.60 BKMR reports −101.7%, i.e. the fitted curvature has
changed sign, which is V15's recorded inversion.

**So `f²` is not a property of censoring, of the imputation, or of the estimand.** With a
correctly specified parametric analysis, on the same imputed data, the curvature loss is close
to **linear** in `f`. The quadratic shape is manufactured by the analysis model: BKMR masks the
attenuation at low censoring and amplifies it at high censoring, and the two together bend a
near-linear loss into a quadratic.

## 3. Supporting checks

**P2 (the instrument) — MET.** `oracle_lm` is flat in `f` (trend p = 0.84 for `curv_X1`,
0.62 for `int_X1X2`) and consistent with zero (max |z| = 1.39). Coverage 0.90–0.97. The 22-rep
pilot's apparent drift (p = 0.017) was Monte-Carlo noise, as the gate's three-state scoring
expected. The instrument reads zero when it should.

**P3 (V15 reproduces on the current harness) — MET.** `pipeline_bkmr` now gives −11.9% ± 5.0
at `f` = 0.20 against V15's −13.9% ± 2.0 (0.4σ), and −49.3% ± 5.2 at `f` = 0.40 against
−55.6% ± 3.3 (1.2σ). V15 predates the `mice` default (v1.7.0) and other changes; its curve is
still the curve this harness produces.

**BKMR inflates curvature on complete data, and this is flat in `f`.** `oracle_bkmr` vs
`oracle_lm` on `curv_X1`: **+15.1% ± 1.5** pooled, per-`f` values +16.2 +15.1 +15.1 +11.9
+14.6 +17.6 — exactly reproducing V33's independently measured +15.1%, with no trend in `f`.
On `int_X1X2` the two agree (−0.4% ± 1.0). So the complete-data inflation is a constant
property of the GP prior; it cannot by itself produce an `f`-dependent exponent. What produces
the exponent is its **interaction** with a degraded tail.

**BKMR is not an offset or a rescaling of OLS.** The difference changes sign between `f` = 0.10
(+14.6, 2.7σ) and `f` = 0.60 (−23.8, 3.7σ), which rules out any additive offset outright; the
ratio `excess_bkmr / excess_lm` runs 0.07 → 0.35 → 0.57 → 0.88 → 1.16 → 1.31, an 18-fold
monotone climb, which rules out a multiplicative one. BKMR *recovers* curvature that OLS loses
while the imputed tail is still a minority, then loses more than OLS once it is not —
consistent with a GP whose smoothness prior is dominated by the uncensored majority, while OLS
uses the corrupted tail rows at full leverage.

**Nominal coverage hides this.** `pipeline_bkmr` keeps 0.92–0.99 coverage across the whole
sweep, including at `f` = 0.60 where its point estimate has the wrong sign — the intervals are
simply wide (width/se 1.2–1.3 even at the oracle). `pipeline_lm` coverage collapses
0.91 → 0.75 → 0.51 → 0.41 → 0.31 → 0.15 over the same sweep. **Coverage is not a usable
warning signal for this failure mode under BKMR**, which is worth stating in the manuscript
wherever coverage is reported as reassurance.

## 4. What is still open

**Still not derived: WHY BKMR does this.** The mechanism is now localised and its shape is
described quantitatively — masking below `f` ≈ 0.44, amplification above, sign inversion by
`f` = 0.60, with a constant +15% complete-data inflation underneath — but no closed form
predicts the exponent. That is a better position than before (the question has moved from "why
does censoring enter as `f²`" to "what does a GP do to a degraded curvature signal"), and it is
a different question from the one this project had been failing to answer.

**Methodology note, recorded because it nearly cost the result.** The first launch died at
88/600 tasks with no R error, no reboot and no sleep event: `nohup … &` only ignores SIGHUP and
leaves the job in the launching shell's process group, so a session teardown took the whole
tree. Separately, the per-chunk checkpoint was **write-only** — nothing read it back, so the
header's claim that "a kill costs one chunk, not the run" was false. Both are fixed: all 18
`run_v*.sh` wrappers now detach via `setsid(2)`, and `run_v3_bkmr.R` honours `RESUME=1`
(which saved 1.6 h here). See `README.md`.

## 5. What may be asserted

- ✅ "The quadratic scaling of curvature loss in the censored fraction is a property of the
  analysis model, not of the imputation: on the same imputed datasets, a correctly specified
  parametric analysis loses curvature close to linearly in `f` (slope 0.75 [0.69, 0.81])
  while BKMR loses it quadratically (2.01 [1.80, 2.21])."
- ✅ "BKMR inflates curvature on complete data by +15% (±1.5), independently of `f`, and this
  masks imputation-induced attenuation while the censored fraction is below ≈0.44."
- ✅ "Interval coverage under BKMR stays near nominal across the whole sweep, including where
  the point estimate inverts in sign; it is not a diagnostic for this failure."
- ❌ Not "OLS is the better analysis model" — it was given the true basis by construction.
- ❌ Not "`f²` is explained" — localised, not derived.

## 6. Files

`run_v35_analysis_model.sh` (the paired runner, registered predictions in its header),
`R/procedures_lm.R` (`oracle_lm` / `pipeline_lm`, mirroring `proc_bkmr_pipeline` step for
step), `analyze_v35.R` (three-state gates), `results/v35_summary.csv` + `v35_raw.csv` (the
paired 100-rep sweep), `results/v35lm_*` (the earlier OLS-only sweep, superseded by the paired
run and consistent with it).
