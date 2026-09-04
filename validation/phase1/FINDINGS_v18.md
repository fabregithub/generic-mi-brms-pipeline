# V18 findings — both hypotheses wrong: the bias decays as **n^−1/3**, and my V17 alarm was overstated

**Run:** 2026-09-03 14:26 → 17:11 JST (**164.6 min** against 2.6 h predicted; 5 stages ×
6 cells × 300 reps = 9,000 tasks, `n_ok` = 300 in all 30 cells, zero task errors, zero
non-finite estimates). Results: `results/v18_n{800,1600,3200,6400,12800}_latest.rds`.
Verdict arithmetic: [`analyze_v18.R`](analyze_v18.R). Criteria:
`../PLAN_pipeline_validation.md` §8k.

**The question.** V17 found the shipped default biased +5% to +10.4% under a confounder or a
mediator, with coverage still 0.93–0.97. I explained that as a sample-size accident — bias
constant in relative terms, SE falling as `1/√n`, so coverage must collapse — and projected
**0.32–0.73 at n = 12,800**. That was arithmetic on one cell, and I said so; V18 measures it.

**Answer: both registered hypotheses are rejected, and my warning was too loud.** The bias
decays as **n^−0.335** (95% CI [−0.378, −0.291]) — neither the 0 of constant bias nor the
−0.5 of a finite-sample artefact. Coverage at n = 12,800 is **0.867–0.923**, not 0.32–0.73.

---

## The primary outcome: relative bias of `pipeline_bartMI`

| cell | n=800 | 1,600 | 3,200 | 6,400 | 12,800 | slope | 95% CI |
|---|---|---|---|---|---|---|---|
| `zr_fork` | 4.62% | 4.28% | 3.50% | 2.43% | **2.07%** | −0.313 | [−0.434, −0.192] |
| `zr_pipe` | 4.78% | 3.01% | 2.70% | 2.03% | **2.00%** | −0.309 | [−0.499, −0.120] |
| `zrmar_fork` | 7.71% | 7.24% | 5.70% | 4.02% | **2.90%** | −0.367 | [−0.520, −0.214] |
| `zrmar_pipe` | 9.51% | 7.21% | 6.15% | 4.55% | **3.56%** | −0.350 | [−0.402, −0.298] |
| **pooled** | | | | | | **−0.335** | **[−0.378, −0.291]** |

**Constant bias is rejected** (CI excludes 0) **and finite-sample bias is rejected** (CI
excludes −0.5). Registering both as falsifiable was what made this readable: the run does not
merely fail to confirm one, it rules out both, and the four per-cell bars at n = 12,800 miss
in the same direction and by similar margins — 2.0–3.6% measured against 5.0–10.4% predicted
by constant bias and 1.2–2.6% by the finite-sample account. **The truth sits between them,
and neither account predicts anything there.**

---

## What this does to the V17 warning

I wrote that a defect hiding at n = 800 and surfacing at n = 8,000 is "the worst possible
shape for a pipeline built for large-N data". **On these numbers that overstates it.**

| cell | coverage @ 12,800 | I projected |
|---|---|---|
| `zr_fork` | **0.867** | 0.693 |
| `zr_pipe` | **0.923** | 0.731 |
| `zrmar_fork` | **0.890** | 0.358 |
| `zrmar_pipe` | **0.877** | 0.321 |

**The coverage model was not wrong — its input was.** Fed the *measured* bias and width/SE it
reproduces observed coverage to within **0.024** at every one of the 20 points. It was the
constant-relative-bias assumption fed into it that failed, and separating those two matters:
the machinery is reusable, the assumption was not.

**But the bias does not go away, and coverage does degrade.** Bias falls as n^−0.335 while
the SE falls as n^−0.50, so **bias/SE still grows — as n^0.17** (measured per-cell slopes
+0.10 to +0.17, consistent with the pooled fit). Coverage is unbounded below in the limit;
it just gets there **3× more slowly in the exponent** than I claimed.

**Extrapolating on the fitted exponent** (beyond the run, and labelled as such): coverage
≈ **0.81–0.89 at n = 100,000**, reaching 0.80 at roughly n = 130,000 (`zrmar_pipe`) to
1.6 million (`zr_pipe`). So the failure is real and permanent, but reaching *serious*
under-coverage takes a sample well beyond this pipeline's realistic range.

**The honest summary is narrower than V17's and still not comfortable:** at every sample size
tested, the shipped default carries 2–10% bias under a causally active covariate that
coverage does not reveal, and more data shrinks it only slowly.

---

## Controls

**The oracle is unbiased at every n** — worst |relative bias| **1.75%** (`zr_collider` at
n = 800, a single-cell fluctuation; every other value is ≤0.75%). The estimand and the DGP
are not n-dependent, so the rest of the run is readable.

**`zr_collider` stays flat** (+1.08% worst, no trend). **`mcar_z40`** — the precision-covariate
cell — runs −2.01% → −0.74% with slope **−0.418, CI [−0.685, −0.151]**. That CI contains both
−0.5 and the pooled −0.335, so **this run cannot say whether the inert cell's small bias has
the same origin as the fork/pipe cells' large one.** I flagged before the run that a
difference there would be a finding; the honest report is that it is underpowered to tell.

---

## A conjecture with a cheap discriminating test

**n^−1/3 is the classical nonparametric convergence rate** for estimating a Lipschitz
function in one dimension by local averaging — and the Z block here is BART, a nonparametric
regression. The conjecture is therefore: **the shipped default's bias under a causally active
covariate is the Z-block imputer's own convergence rate, showing up in the estimand.**

What makes it testable rather than decorative: **in the `fork` and `pipe` cells the true
`Z₁` conditional is linear-Gaussian**, so a *parametric* imputer is correctly specified there
and should converge at the parametric rate. The conjecture therefore predicts

| Z-block imputer | predicted slope |
|---|---|
| `bartMI` (nonparametric) | **≈ −1/3** ✓ measured |
| `micePmm` (parametric, correctly specified here) | **≈ −1/2** |
| a misspecified parametric imputer | **≈ 0** (asymptotic bias) |

Those are three visibly different curves, and both arms already exist
(`pipeline_micePmm`, `pipeline_properBoot`). The run is this same driver with a different
`ARMS` — no new code, and ~2.6 h. **If the exponent tracks the imputer, the theory gains its
first mechanism for a rate rather than a magnitude.**

**This is a conjecture, not a result.** The rate coincidence is suggestive and nothing more;
`−1/3` also sits close to several other plausible rates, and one DGP cannot distinguish them.

---

## What this revises

| claim | status |
|---|---|
| V17: "the concealment is a sample-size accident; coverage collapses" | **Overstated.** Coverage at n = 12,800 is 0.87–0.92, not 0.32–0.73. Direction right, magnitude wrong by ~3× in the exponent |
| "relative bias is ~constant in n" (V17 §8j, the projection's premise) | **Refuted** — slope −0.335, CI excludes 0 |
| "the bias is a finite-sample artefact that washes out" | **Refuted** — CI excludes −0.5; 2.0–3.6% remains at n = 12,800 |
| The coverage model itself | **Validated** — reproduces all 20 measured points to within 0.024 from measured inputs |
| V17: "the shipped default is biased under a confounder or mediator" | **Confirmed at every n tested**, and it does not wash out |
| `THEORY.md` | Gains a **third measured exponent with no derivation**, beside V15's `f²`. Two rates now sit in the record unexplained |

---

## Caveats (scope of this run)

- **One DGP, five n levels, 300 reps, `m` = 30 held fixed.** `m` is not scaled with `n`;
  whether the exponent survives `m` growing with `n` is untested and is a plausible
  confounder for the rate itself.
- **The extrapolation beyond n = 12,800 is arithmetic on the fitted exponent**, not a
  measurement, and it inherits the CI: at the lower end (−0.29) coverage falls faster than
  quoted, at the upper end (−0.38) slower.
- **`width/SE` drifts downward with n** (1.125 → 1.032 in `zrmar_pipe`; 1.001 → 0.900 in
  `zr_fork`) and `zr_fork` at n = 12,800 lands at **0.8998**, a hair under the 0.9 bound I
  registered. A technical miss on an arbitrary bound, reported rather than rounded away: it
  voids nothing, since the coverage model is validated against *measured* width, but the
  drift itself is unexplained and compounds the coverage decline.
- **The n^−1/3 conjecture is untested.** The discriminating run is specified above.
- **`m` = 30 with FMI 0.4–0.5** is adequate at every n here, but the pooling correction's
  behaviour at large `n` and fixed `m` was not separately checked.
