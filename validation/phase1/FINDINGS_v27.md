# V27 — should the covariate imputer default to `mice` instead of `bart`?

**Run 2026-10-11/12.** 3 `n` levels × 11 cells × 500 reps = 16,500 tasks, 3 arms each,
**14 h 20 min against 15.9 h predicted**, `rc = 0`, zero task errors, all cells 500/500.
Driver [`run_v27_zimputer.sh`](run_v27_zimputer.sh), analysis [`analyze_v27.R`](analyze_v27.R).

**Outcome: the default switched to `mice` in v1.7.0.**

## The problem it addresses

The *recommended* configuration (`z_imputer = "bart"`, `Y` in both blocks) carried **+4.96%**
(confounder), **+5.39%** (mediator), **+8.52%** and **+10.40%** (both under MAR) in the exposure
coefficient. The oracle is unbiased in all four, so it is the imputation's bias and not the
design's — and **coverage does not detect it** (0.93–0.97). It decays only as `n^−1/3`, leaving
2.0–3.6% at n = 12,800. This was the last item where a user running the documented
configuration was exposed.

## Criteria

| # | criterion | 800 | 3200 | 12800 |
|---|---|---|---|---|
| C1 | mice within 0.3 `bias/SE` of bart on inert + non-linear cells | **MET** | **MET** | **MET** |
| C2 | mice improves `bias/SE` by ≥ 1.0 in causally-active cells | no | no | no |
| C3 | coverage not traded away | **MET** | **MET** | **MET** |
| C4 | the shipped option equals the harness one | **MET** | **MET** | **MET** |

## 1. C1 — the gate that decides, and it passed

A parametric covariate imputer is misspecified where the covariate's own conditional is
non-linear, and V6 measured `mice pmm` at **−4.70%** on exactly that DGP. That measurement is
why BART was adopted, and it is the only thing that could have stopped this change. The `nl_*`
cells are that ground. At n = 12,800:

| cell | bart | mice | mice's loss |
|---|---|---|---|
| `mcar_z40` | 0.07 | 0.01 | **−0.06** |
| `mar_z40` | 0.06 | 0.04 | **−0.02** |
| `nl_mcar_z40` | 0.13 | 0.12 | **−0.01** |
| `nl_mar_z40` | 0.16 | 0.08 | **−0.09** |

Every loss is **negative** — mice is slightly *better* even where it was supposed to lose, at
every `n`. V6's objection did not reproduce.

## 2. The gain, where the bias actually lives

Relative bias, bart → mice:

| cell | n = 800 | n = 3200 | n = 12,800 |
|---|---|---|---|
| `zrmar_pipe` | 10.21% → **2.95%** | 5.38% → 0.87% | 3.50% → **0.66%** |
| `zrmar_fork` | 8.02% → **0.84%** | 4.98% → 0.54% | 2.66% → **0.08%** |
| `zr_fork` | 5.23% → **0.79%** | 2.75% → 0.35% | 1.45% → **0.00%** |
| `zr_pipe` | 4.99% → **0.52%** | 2.42% → 0.09% | 1.87% → **0.36%** |
| `zr_pipenl` (non-linear arrow, on-path) | — | — | +0.46 → **−0.00** `bias/SE` |

**mice was better or equal in every cell at every `n`.** Coverage was never worse by more than
**0.008** (C3), and the shipped `z_imputer = "mice"` option reproduced the harness arm to
**0.00 Monte-Carlo errors** (C4), so the thing that shipped is the thing the evidence was
gathered on.

> **⚠️ C4 passed while the shipped function could not run at all** (found 2026-10-13, on the
> first end-to-end run after the default was switched). `run_row_level_imputation_mice()` was
> ported from the harness and kept a call to `.v4_is_binary()` — a harness helper that does not
> exist in shipped code. C4 executes the shipped option *inside the harness*, which defines that
> symbol, so the dependency was always satisfied and the agreement to 0.00 MC errors was real
> but said nothing about whether the code runs outside. `Rscript 03_impute.R` died immediately
> with `could not find function ".v4_is_binary"`.
>
> **The numbers stand; the test's scope did not cover what it was read as covering.** "The
> shipped option reproduces the harness arm" was true. "The shipped option works" was not
> tested, and is a different claim. A port can only be verified against the environment it will
> actually run in — the fix makes the function self-contained, and the `test/` end-to-end suite,
> which loads no harness code, is what caught it. Note also that the obvious substitution was
> wrong: shipped `is_binary_like()` accepts ANY two-level variable, while the harness predicate
> requires {0,1}, so swapping it in would have routed extra variables through `logreg` on
> evidence that never covered them.

## 3. ⚠️ C2 failed on a bar I set badly — the fifth such gate in this programme

C2 required a ≥ 1.0 improvement in `|bias/SE|`. But **bart's own `|bias/SE|` never exceeds 0.78
in any cell at any `n`**, so a 1.0 improvement was unreachable by construction. Four of six
cells came back VOID for precisely that reason (bart at 0.43, 0.48, 0.01, 0.01), and the two
informative ones gave +0.68 and +0.61.

This is the same failure as V22's, V23's and both of V26's: **a threshold registered against a
quantity whose achievable value is below it.** The guard I added after V26 — a floor on the
quantity being compared — worked exactly as intended here, correctly voiding the four
uninformative cells. What it cannot do is notice that the *bar itself* is unreachable. The
lesson to carry: when registering a threshold, check it against the largest value the quantity
has ever taken in a comparable cell.

**The decision therefore rests on C1, C3, C4 and the uniform improvement — not on C2.** C1 was
designated in the driver, before the run, as the gate that decides the default; that is why its
pass is load-bearing and C2's failure is not.

## 4. A registered prediction, partly held

`THEORY.md` §4a derived `n^{−1/3}` for a tree imputer's smoothing bias and registered **P1**
before this run: a correctly specified parametric imputer has no approximation bias, so `mice`
should sit at ≈ 0 at *every* `n` rather than decaying.

**Qualitatively held.** mice: 0.79% → 0.35% → 0.00% (`zr_fork`); bart decays from 5.23% to
1.45%.

**But bart's measured exponents are −0.35, −0.40, −0.39, −0.46** (mean −0.40) against the
derived −1/3 = −0.333. With three points per fit, −1/3 and −2/5 cannot be separated. −2/5 would
mean BART behaves like a *twice-differentiable* smoother rather than piecewise-constant on a
Lipschitz target — plausible for a shrunken ensemble. **The framework survives; the smoothness
exponent is unresolved**, and §4a's P2 (that adding continuous predictors flattens the exponent
to `1/(d+2)`) remains the sharper test.

## What this does not cover

- **Only the censored-exposure path's covariate block.** The row-level strategy uses the same
  `z_imputer` resolution, so it inherits the change, but no cell here exercises it directly.
- **`mice` is now a hard dependency of the default path.** If it is absent the run fails rather
  than falling back — deliberately, so a silent substitution cannot happen.
- **The inert cells show no change**, so an analysis whose covariates are all inert is unaffected
  and its previous results stand.
