# V19 findings — the conjecture is refuted, and **BART is the only Z-block imputer whose bias goes away**

**Run:** 2026-09-04 09:32 → 2026-09-07 01:47 JST (**64.2 h** against ~63 h projected; 5 `n`
levels × 3 cells × 4 arms + oracle × 300 reps = 4,500 tasks, `n_ok` = 300 in all 45 cells,
zero task errors, zero non-finite estimates). Results: `results/v19_n*_latest.rds`. Verdict
arithmetic: [`analyze_v19.R`](analyze_v19.R). Registered in the driver header of
[`run_v19_imputer_rate.sh`](run_v19_imputer_rate.sh), committed before the run; criteria
`../PLAN_pipeline_validation.md` §8l.

**The conjecture.** V18 measured the shipped default's bias under a causally active covariate
decaying as **n^−1/3** — the classical nonparametric convergence rate, and the Z block is
BART. So perhaps the bias simply *is* the Z-block imputer's own rate. In the `fork` and `pipe`
cells the true `Z₁` conditional is linear-Gaussian, so a **parametric** imputer is correctly
specified there and should reach −1/2. Four arms, two per family.

**Answer: refuted — family membership predicts nothing.** But the run answers a better
question than it asked. **Only one of four imputers has a bias that shrinks with sample size
at all**, and it is the shipped default.

---

## Pooled slopes: three of four are flat

log|relative bias| on log `n`, pooled over `zr_fork` and `zr_pipe`:

| arm | family | predicted | **measured slope** | verdict |
|---|---|---|---|---|
| `bartMI` | nonparametric | −1/3 | **−0.311 [−0.389, −0.234]** | **decays** |
| `properBoot` | nonparametric | −1/3 | **−0.029 [−0.087, +0.028]** | asymptotic |
| `micePmm` | parametric, correct | −1/2 | **+0.097 [−0.031, +0.225]** | asymptotic |
| `properZ` | parametric, correct | −1/2 | **+0.044 [−0.003, +0.091]** | asymptotic |

**Every registered gate on the conjecture failed.** The families are not separated — the
nonparametric pair spans 0.282 between its own members, the two parametric arms sit *above*
zero rather than at −1/2, and no arm but `bartMI` lands near the rate its family predicts.
The split is not between nonparametric and parametric. **It is between BART and everything
else.**

**`bartMI` is the only arm whose CI excludes zero.** The other three have bias that does not
vanish with data — the "constant bias" hypothesis V18 rejected for BART is exactly right for
them.

**`bartMI`'s −0.311 also replicates V18's −0.335** under a different arm set and therefore a
different RNG stream (PLAN §10b), which is the check that V18's exponent was not an artefact
of its arm configuration.

---

## The result that matters: bias and coverage at scale

| arm | bias @ n=800 | bias @ n=12,800 | coverage @ 800 | **coverage @ 12,800** |
|---|---|---|---|---|
| **`bartMI`** | +4.6 / +4.8% | **+2.1 / +2.0%** | 0.933 / 0.953 | **0.867 / 0.923** |
| `micePmm` | +2.4 / +2.9% | +3.5 / +3.9% | 0.953 / 0.960 | 0.780 / 0.837 |
| `properZ` | +4.5 / +5.4% | +5.0 / +6.3% | 0.943 / 0.950 | 0.623 / 0.607 |
| **`properBoot`** | **−22.8 / −26.4%** | **−22.0 / −23.3%** | 0.557 / 0.577 | **0.000 / 0.003** |

*(`zr_fork` / `zr_pipe`.)*

**The v1.5.0 decision to default to BART is vindicated on a question nobody had asked of it.**
R8 adopted BART for calibration and flexibility under a *causally inert* covariate. Under a
confounder or a mediator it is the only imputer here whose bias shrinks and whose coverage
survives to n = 12,800.

**`properBoot` fails catastrophically and does not hide it.** −23% bias that never decays, and
coverage reaching **zero** by n = 6,400. Unlike every other failure this project has recorded,
this one *is* visible in the diagnostics.

**`properZ` and `micePmm` are the quiet ones.** Both look fine at n = 800 (coverage 0.94–0.96)
and both degrade to 0.61–0.84 by n = 12,800 on a bias that never moves. They are the shape
V17 warned about and V18 found BART had escaped.

---

## What this does to `properBoot`, and the risk it carries today

`properBoot` is `proper_draw = TRUE` — the **v1.4.0 default**, still reachable two ways:

1. a config written against v1.4.0 that sets `proper_draw = TRUE`, honoured for
   back-compatibility (`resolve_z_imputer()`);
2. **`z_imputer = "bart"` when `dbarts` is not installed.** `run_row_level_imputation()` warns
   loudly and falls back to `forest_boot`, so the run completes.

Path 2 is the one to worry about. The warning is `immediate. = TRUE` and the imputer is
logged, so it is not silent — but a user who misses it, on data with a confounder, gets a
**−23% bias that more data will not fix**. The saving grace is that this is the one failure
mode coverage catches: it would be 0.0 at any reasonable `n`.

**This is not a call to change the fallback** — completing with a warning beats aborting, and
V19 measured one DGP. It is a call to say so in the documentation, which currently presents
`forest_boot` as merely the older default rather than as one that carries asymptotic bias
under a causally active covariate.

---

## Control cells

**The oracle is unbiased at every `n` and cell** (worst 0.59%), so the estimand and the DGP
are not `n`-dependent and the exponents are readable.

**`mcar_z40`, the causally inert cell**, was carried to see whether any of this appears
without a confounding path. It does not, in the way that matters: all four arms sit at
−0.7% to −4.2% with no arm exceeding it, and only `bartMI` shows a decaying trend
(−0.418 [−0.685, −0.151]) — the same ordering as the on-path cells at a twentieth of the
magnitude. **The imputer choice barely matters when the covariate is causally inert**, which
is precisely why fifteen tracks measured it and found little.

---

## What this revises

| claim | status |
|---|---|
| "the bias exponent tracks the Z-block imputer's family" (V18 conjecture) | **Refuted.** `properBoot` is nonparametric and flat; both parametric arms are flat and positive, not −1/2 |
| V18: `bartMI` decays as n^−1/3 | **Replicated** at −0.311 under a different arm set |
| "n^−1/3 is the nonparametric convergence rate showing up in the estimand" | **Unsupported.** A second nonparametric imputer shows no decay at all. `THEORY.md` keeps the exponent as measured-and-underived |
| R8 / v1.5.0: BART as the Z-block default | **Vindicated far beyond its original evidence** — the only imputer whose bias and coverage survive scale under a causally active covariate |
| `forest_boot` (v1.4.0 default, and the `dbarts`-missing fallback) | **New defect: −23% asymptotic bias under a confounder or mediator**, coverage → 0 |

---

## Caveats (scope of this run)

- **One DGP, two on-path cells, arrow strengths fixed at 0.60.** The ranking is measured;
  the magnitudes are specific to this surface.
- **`m` = 30 fixed at every `n`.** Untested whether the flat exponents move if `m` grows with
  `n` — a plausible confounder for all four arms, and the same caveat V18 carries.
- **The bias is still not attributed to a block.** All four arms share the `leftcens` X block,
  and three of four are flat while one decays, so a common X-block floor does not explain the
  pattern on its own. A V13-style decomposition under a fork remains the open step.
- **`properBoot`'s −23% is one configuration.** It was not investigated *why* — the arm is
  also 17× the cost of every other, and that too is unexplained (`ROADMAP.md` loose ends).
- **Nothing here tests a misspecified parametric imputer**, which the original design wanted
  as the ≈0 arm. Both parametric arms turned out to be ≈0 while *correctly specified*, so
  that contrast is now more interesting, not less, and needs a cell with a non-linear `Z₁`
  conditional on an X–Y path — which no scenario currently provides.
