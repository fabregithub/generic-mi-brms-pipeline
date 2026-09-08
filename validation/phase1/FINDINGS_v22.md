# V22 findings — **the trade does not exist**: a parametric covariate draw wins on both sides. And a new 20% defect in the exposure draw

**Run:** 2026-09-08 09:01 → 16:16 JST (**435.3 min** against 7.0 h predicted; 3 `n` levels ×
3 cells × 7 arms + oracle × 500 reps = 4,500 tasks, `n_ok` = 500 in all 9 cells, zero task
errors, zero non-finite estimates). Results: `results/v22_n*_latest.rds`. Verdict arithmetic:
[`analyze_v22.R`](analyze_v22.R). Registered in
[`run_v22_deciding_cell.sh`](run_v22_deciding_cell.sh)'s header, committed before the run.

**The question.** V21 traced the confounder/mediator bias entirely to BART's *flexibility* —
a correctly specified estimated draw was unbiased. But both V21 cells had a
**linear-Gaussian** covariate conditional, so every parametric arm was correct by
construction. R8 adopted BART precisely because a parametric Z block is misspecified when
that conditional is non-linear (V6: `mice pmm` at −4.70%). **Trading a 5% bias under
linearity for a 5% bias under non-linearity is not a fix**, and nothing could test both sides.

**Answer: there is no trade.** In the non-linear cell the parametric draw is **−0.14%** and
BART is **+2.84%**. The parametric draw wins in *both* cells, by 4.66 pp and 2.98 pp. **A fix
is real.**

And the run turned up something larger than the defect it was built to study: **with the
exposure censored under a non-linear covariate arrow, `leftcens` carries ~20% bias that does
not decay with sample size.**

---

## The primary result

Relative bias of the exposure coefficient. Exposure fully observed in both primary cells (see
*Why the exposure is observed*, below).

| arm | `zr_pipe_nc` (linear) 800 / 3,200 / 12,800 | `zr_pipenl` (**non-linear**) |
|---|---|---|
| `exact` *(anchor)* | −0.34 / +0.07 / +0.14% | −0.59 / +0.03 / +0.14% |
| **`fit_proper`** | **−0.22 / +0.27 / +0.16%** | **−0.14 / +0.29 / +0.17%** |
| `fit_improper` | −0.61 / +0.20 / +0.15% | −0.63 / +0.21 / +0.15% |
| `pmm` | −0.20 / +0.30 / +0.17% | −0.75 / +0.22 / +0.04% |
| **`bart`** | **+4.44 / +2.85 / +1.64%** | **+2.84 / +2.43 / +1.51%** |

**The anchor holds** — within 0.59% of zero in both primary cells — and the oracle is unbiased
(worst 0.78%), so the ladder has a valid zero point. That check is not a formality: it is what
caught the design error described below.

**`fit_proper` is unbiased even though it is now misspecified.** The mechanism was derived
before the run and it holds: the analysis conditions on `logX1`, so a linear draw's error in
approximating `g(logX1)` is absorbed by `logX1`'s own coefficient. What matters for `b₁` is
the part of `Z1` carrying information about `Y` beyond `X` — and `Y` enters the exact
conditional **linearly**, which a linear fit represents exactly.

**This is estimand-specific and that limit is the whole caveat.** Misspecifying the covariate
conditional is harmless *here* because the analysis adjusts for the variable that carries the
misspecification. For an estimand that depends on `Z1`'s own coefficient, or on a surface in
`Z1`, there is no such absorption and this result would not transfer.

---

## One gate missed, and its pre-registered reading was wrong

`bart` came in at **+2.84%** in `zr_pipenl`, under the registered ≥+3% bar. The registered
interpretation of that miss was *"flexibility stops costing anything once the truth is
non-linear, so there is nothing to fix."*

**The neighbouring numbers refute that.** BART's penalty **shrinks** under non-linearity
(+4.44% → +2.84%, exactly as one would expect when a flexible method finally has something to
learn) but it does **not** vanish — and `fit_proper` still beats it by **2.98 pp** in that
very cell.

A binary threshold on one arm was the wrong instrument; the registered **arm comparison** was
the right one, and it passed in both cells. The pre-written conclusion is withdrawn rather
than applied, and `analyze_v22.R` now prints the comparison instead of the sentence its own
numbers contradict.

---

## Why the exposure is observed in the primary cells — and the 20% defect that forced it

The first build of this cell censored `logX1` at 40%, like every other cell in the harness.
**Its exact-`Z` anchor came out at +17%.** Diagnosing rather than patching: with `logX1`
observed, the same anchor is +1.39% ± 0.89 — so the exact conditional is right and the X block
is the problem. `leftcens` draws `logX1` from a conditional **linear in `Z1`**, while the truth
has `Z1 = g(logX1)`.

Kept as a secondary cell, and it is worse than the +17% pilot suggested:

| arm | n = 800 | 3,200 | 12,800 |
|---|---|---|---|
| **`exact`** *(exact Z, shipped X)* | **+20.79%** | **+20.29%** | **+20.31%** |
| `fit_proper` | +22.12% | +22.03% | +22.40% |
| `bart` | +29.21% | +25.92% | +24.60% |

**The X block alone carries ~20%, and it is flat in `n`** — asymptotic, not a finite-sample
artefact. That is 2–4× the covariate-draw defect this whole line of work has been chasing, and
it lands on the censored-exposure draw, which is the strategy's entire purpose.

**Two things it is not.** It is not a decomposition — the cell has no valid anchor, so the
+20% is a *total* attributable to the X block only because the Z draw was made exact. And it
is not measured against a shipped baseline: `pipeline_bartMI` was not in this arm set.

**It needs its own track**, and it is now the largest open defect on the list. `leftcens`'s
conditional is linear in its predictors by construction (Phase 1 §7.5 chose the shash *margin*
for skew, not a non-linear mean), so a covariate with a non-linear relationship to the censored
exposure is outside what it can represent — and nothing before V22 had such a covariate.

---

## Two secondary results

**V19's discrepancy is not the `Z2` draw.** Drawing `Z2` by the same method as `Z1` moves the
estimate by at most **0.17 pp** in either primary cell. V19 measured `micePmm` at +2.4 to
+3.9% and `properZ` at +4.5 to +6.3% *asymptotic*, against ±0.6% on this ladder; the boundary
was registered in advance, and this eliminates one of the three candidates. What remains: `Y`
as an imputation target, or mice's own predictor matrix.

**Properness reproduces V4, weakly.** Dropping the parameter draw moves bias by ≤0.49 pp while
coverage falls 0.962 → 0.958 and 0.946 → 0.944. The direction is V4's signature; the magnitude
is smaller than V21's (0.95 → 0.94), and with `Z2`'s draw held exact in the base arms there is
less parameter uncertainty in play to lose. Worth noting that coverage in `zr_pipenl` runs
0.912–0.946 for *every* arm including the anchor — somewhat below nominal, and unexplained.

---

## What this changes

| claim | status |
|---|---|
| V21: "a correctly specified estimated covariate draw is unbiased" | **Extended** — it is unbiased even when *misspecified*, for this estimand, because the analysis absorbs the error |
| "the fix is a trade: flexibility for correct specification" | **Refuted.** The parametric draw wins on both sides, by 4.66 pp and 2.98 pp |
| "BART is doing its job once the truth is non-linear" | **Not supported** — its penalty shrinks 4.44 → 2.84% but does not vanish |
| V19's +2.4 to +6.3% is the `Z2` draw | **Ruled out** (≤0.17 pp) |
| The censored-exposure draw under a non-linear covariate arrow | **New defect: ~20%, asymptotic** — larger than anything previously measured here |

**For the fix:** a parametric covariate draw is now supported by evidence on both sides of the
one trade that stood against it. It is not yet a proposal — it needs the analyst to supply a
form, and V19's implementations remain unexplained — but the objection that killed it is gone.

**For the docs:** `docs/covariate-roles.md` says switching to a parametric covariate imputer
is "**not** a workaround", citing V19. That is still true of the *available* implementations
and false of the *approach*. It needs rewording, and the ~20% X-block finding needs a home.

---

## Caveats (scope of this run)

- **One non-linear form** (`1.2·tanh(1.8x) + 0.35·(x²−1)`), one DGP, one estimand. The
  absorption argument is general for a coefficient the analysis adjusts for; the magnitudes
  are not.
- **The primary cells have no censored exposure**, which is exactly why they are readable —
  but it means neither primary cell is the shipped configuration.
- **The ~20% X-block figure has no anchor and no shipped baseline.** It is a signal to
  investigate, not a decomposed result.
- **`Z2`'s draw is exact in the five base arms.** The two `z2same` arms probe that, and find
  it does not matter — but the multi-target question is only half-answered.
- **Coverage runs 0.912–0.946 in `zr_pipenl` for every arm**, anchor included. Below nominal,
  common to all arms, and not explained here.
- **`m` = 30 fixed at every `n`**, as in V18–V21.
