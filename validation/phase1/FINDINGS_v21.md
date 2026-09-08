# V21 findings — it is BART's **smoothing**, all of it; a correctly specified estimated draw is unbiased

**Run:** 2026-09-07 13:12 → 18:28 JST (**315.9 min** against 4.9 h predicted; 3 `n` levels ×
2 cells × 6 arms + oracle × 500 reps = 3,000 tasks, `n_ok` = 500 in all 6 cells, zero task
errors, zero non-finite estimates). Results: `results/v21_n*_latest.rds`. Verdict arithmetic:
[`analyze_v21.R`](analyze_v21.R). Registered in
[`run_v21_zdraw_ladder.sh`](run_v21_zdraw_ladder.sh)'s header, committed before the run;
criteria `../PLAN_pipeline_validation.md` §8n.

**The question.** V20 located the confounder/mediator bias in the covariate draw. V19 showed
"parametric and correctly specified in form" is not obviously the answer. **Which property of
a covariate draw has to be right?**

**Answer: only the flexibility. Nothing else on the ladder matters at all.** BART's
nonparametric fit of a linear conditional carries **100% of the bias**. A correctly specified
*estimated* linear draw carries **+0.26% / +0.58%** — essentially none. So does an
**improper** version of it, and so does a **donor-matched** one. **Every registered gate
passed.**

---

## The ladder

Relative bias of the exposure coefficient. X block held at the shipped `leftcens` draw in
every arm.

### `zr_fork` (confounder)

| arm | n = 800 | 3,200 | 12,800 | slope |
|---|---|---|---|---|
| `exact` *(anchor)* | −0.05% | +0.46% | −0.01% | — |
| `fit_proper` | **+0.26%** | +0.55% | −0.03% | — |
| `fit_improper` | +0.01% | +0.42% | −0.02% | — |
| `pmm` | +0.07% | +0.50% | −0.04% | — |
| **`bart`** | **+4.73%** | **+2.99%** | **+1.31%** | **−0.463** |
| `bart_inner3` | +4.68% | +2.87% | +1.25% | −0.476 |

### `zr_pipe` (mediator)

| arm | n = 800 | 3,200 | 12,800 | slope |
|---|---|---|---|---|
| `exact` *(anchor)* | +0.31% | −0.10% | +0.21% | — |
| `fit_proper` | **+0.58%** | −0.14% | +0.24% | — |
| `fit_improper` | +0.41% | −0.15% | +0.26% | — |
| `pmm` | +0.64% | −0.16% | +0.23% | — |
| **`bart`** | **+5.11%** | **+2.34%** | **+1.64%** | **−0.410** |
| `bart_inner3` | +5.13% | +2.30% | +1.67% | −0.404 |

**The anchor holds** (worst 0.46% from zero) and **the oracle is unbiased** (worst 0.37%), so
the ladder has a valid zero point. No slope is reported for the four arms whose bias is within
noise of zero — a slope fitted to ±0.5% is fitted to noise, and the naive fits confirm it
(pmm's CI is [−12.4, +12.0]).

**Share of the span the ladder attributes** (anchor → `bart`, ≈4.8 pp in both cells):

| arm | `zr_fork` | `zr_pipe` |
|---|---|---|
| `fit_proper` | 6.6% | 5.8% |
| `fit_improper` | 1.3% | 2.1% |
| `pmm` | 2.5% | 7.0% |
| **`bart`** | **100%** | **100%** |

---

## Three secondary results, each a confirmation

**Properness is an interval problem, not a bias one — exactly V4's signature.** Dropping the
parameter draw moves the bias by at most **0.25 pp**, while shrinking the between-imputation
variance `b` by **12–16%** and coverage from 0.946–0.956 to 0.938–0.954. That is what
improper MI looks like, measured here in a structure V4 could not see, and it means properness
is worth having for the interval and cannot explain the point-estimate defect.

**Donor matching costs nothing here.** `pmm` sits on the anchor. My registered prediction —
that pmm's bias would not decay with `n` — is **untestable rather than confirmed**: there is
no bias to decay. Worth saying plainly, because "the slope's CI includes 0" would be a
misleading way to report a null of that kind.

**The inner-FCS loop contributes nothing**: −0.12 to +0.03 pp. So V20's instrument omitting
it was harmless, and the 0.26–0.80 pp gap noted there has some other cause.

**The instruments agree.** V21's `z21_bart` and V20's `ef_bart_ship` match to **≤0.10 pp** at
all six cell × `n` combinations — independent code paths, same answer.

---

## The tension with V19, which was registered in advance

V19 measured `micePmm` at **+2.4 to +3.9%** and `properZ` at **+4.5 to +6.3%**, both
*asymptotic*. Both are parametric. Here every parametric variant sits at ±0.6%.

The driver's header registered this boundary before the run: *"if the ladder does not
reproduce micePmm's +4 to +6%, that is informative — it would mean the damage is in the
multi-target machinery rather than the draw."* **That is what happened**, so it is a located
discrepancy rather than a contradiction.

What differs between `z21_pmm` and V19's `micePmm`: this arm draws **`Z1` alone**, on the
**correct formula**, with `Z2` held at its exact conditional. `micePmm` runs the pipeline's
mice path — imputing `Z1`, `Z2` **and `Y`** together, with mice's own predictor matrix and its
own iteration. `properZ` is a harness instrument with its own construction.

**So "use a parametric covariate imputer" is still not the fix.** The *property* of being
parametric-and-correct is enough (V21); the *implementations* available in the pipeline are
biased for some other reason (V19). Identifying that reason is now the open question, and it
is a narrower one than before: multi-target joint imputation, the predictor matrix, or `Y`'s
presence as a target.

---

## What this means for a fix

**There is a shippable fix in principle, and V21 is what establishes that.** The most
consequential registered outcome was the opposite: *if `fit_proper` exceeded 2 pp, then
estimating a correctly specified conditional would itself be enough to cause the defect, and
there would be no fix — only the documented guidance.* It came in at **0.26% and 0.58%**.

**What the fix requires:** a covariate draw whose conditional is correctly *specified*, not
merely flexible. That is a real constraint — it means the analyst must get the covariate's
conditional form right, which is exactly the burden BART was adopted to remove. So the
honest framing is a **trade**, not an improvement:

| | BART (current default) | a correctly specified parametric draw |
|---|---|---|
| covariate conditional linear | +5% bias, decaying as n^−1/3 | **≈0%** |
| covariate conditional non-linear | flexible — handles it | **misspecified; untested here** |

**V21 does not test the second row**, and that is the decisive missing cell. The `fork` and
`pipe` cells have a linear-Gaussian `Z₁` conditional, so a parametric draw is correct *by
construction*. R8 adopted BART precisely because a parametric Z block is misspecified when the
covariate conditional is non-linear, and V6 measured `mice pmm` at −4.70% bias in that case.
**A fix that trades a 5% bias under linearity for a 5% bias under non-linearity is not a
fix**, and nothing here says which way that trade goes.

**The next step is therefore not to build the parametric draw.** It is to run this ladder in a
cell where the covariate conditional is genuinely non-linear *and* the covariate is on an X–Y
path. No scenario provides one — `z_form = "nonlinear"` makes `Z1` a descendant of the
exposures, which V17 showed is off-path. That cell has to be built.

---

## Caveats (scope of this run)

- **Two cells, both with a linear-Gaussian `Z₁` conditional.** This is the run's central
  limitation and it is the whole of the caveat above: every parametric arm is correct by
  construction here, so their success is not evidence they would succeed generally.
- **`Z2` is held at its exact conditional in every arm**, to keep the ladder about `Z1`. A
  real imputer has to draw both, and V19's multi-target discrepancy may live exactly there.
- **BART's slope here is −0.41 / −0.46**, steeper than the −0.31 to −0.36 measured for the
  whole pipeline in V18/V19/V20. The estimand is the same; the arm is not (X block shipped,
  `Z2` exact), so the exponents are not directly comparable.
- **`m` = 30 fixed at every `n`**, as in V18–V20.
- **Nothing here revisits `properBoot`**, which V19 measured at −23% asymptotic. It is a
  bootstrapped forest — nonparametric, so the smoothing account predicts a bias, but not one
  of that size or sign.
