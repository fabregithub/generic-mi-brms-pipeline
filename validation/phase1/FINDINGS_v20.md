# V20 findings — it is the **covariate draw**, essentially all of it, and every registered gate passed

**Run:** 2026-09-07 09:27 → 12:31 JST (**183.9 min** against 2.9 h predicted; 3 `n` levels ×
2 cells × 5 arms + oracle × 500 reps = 3,000 tasks, `n_ok` = 500 in all 6 cells, zero task
errors, zero non-finite estimates). Results: `results/v20_n{800,3200,12800}_latest.rds`.
Verdict arithmetic: [`analyze_v20.R`](analyze_v20.R). Registered in
[`run_v20_attribute.sh`](run_v20_attribute.sh)'s header, committed before the run; criteria
`../PLAN_pipeline_validation.md` §8m.

**The question.** V17 found the shipped default biased +5% to +10.4% wherever the covariate
is a confounder or a mediator. V18 measured the scaling (n^−1/3, permanent). V19 ruled out
the imputer-family account. **None of them said which block produces it**, and no fix could
be written without that.

**Answer: the covariate draw, and almost nothing else.** Replacing the Z block's BART draw
with the true conditional — leaving the exposure draw exactly as shipped — removes **85% to
104%** of the bias. Replacing the exposure draw instead removes **0.1% to 11.7%**. The
interaction is **≤0.40 pp**. This is the first track since V15 where every registered gate
passed.

---

## The control, which gates everything else

`ef_bart_ship` is the instrument configured to do what the pipeline does. If it did not
reproduce `pipeline_bartMI`, nothing below would attribute anything.

| cell | `n` | `pipeline_bartMI` | instrument | paired difference |
|---|---|---|---|---|
| `zr_fork` | 800 | +5.01% | +4.63% | −0.38 ± 0.11 |
| `zr_fork` | 3,200 | +3.35% | +2.99% | −0.36 ± 0.06 |
| `zr_fork` | 12,800 | +2.13% | +1.33% | −0.80 ± 0.05 |
| `zr_pipe` | 800 | +5.44% | +5.14% | −0.31 ± 0.13 |
| `zr_pipe` | 3,200 | +2.64% | +2.38% | −0.26 ± 0.06 |
| `zr_pipe` | 12,800 | +2.00% | +1.67% | −0.33 ± 0.03 |

**Passes** — worst 0.80 pp against a registered 1.5 pp bar. The differences are small but
systematically negative, i.e. the instrument is very slightly *better* than the pipeline;
the likely reason is that it omits the pipeline's Z-block inner-FCS loop. Worth noting, not
worth chasing: it is a fifth of the effect being attributed.

The **oracle** is unbiased at every `n` and cell (worst 0.34%), so the estimand and the DGP
are not `n`-dependent.

---

## The decomposition

Paired within replication against `ef_bart_ship`, in percentage points of `b₁ = 0.40`.

### `zr_fork` (confounder)

| `n` | reference | fix-Z | fix-X | both | interaction |
|---|---|---|---|---|---|
| 800 | 4.63% | **−4.64** | −0.31 | −4.71 | +0.24 |
| 3,200 | 2.99% | **−2.55** | −0.35 | −2.67 | +0.22 |
| 12,800 | 1.33% | **−1.33** | −0.07 | −1.34 | +0.06 |

### `zr_pipe` (mediator)

| `n` | reference | fix-Z | fix-X | both | interaction |
|---|---|---|---|---|---|
| 800 | 5.14% | **−4.87** | −0.01 | −4.48 | +0.40 |
| 3,200 | 2.38% | **−2.47** | −0.06 | −2.32 | +0.22 |
| 12,800 | 1.67% | **−1.46** | −0.11 | −1.46 | +0.10 |

*(MC error ≈ 0.11–0.20 pp on each contrast.)*

**Share of the reference bias removed:**

| cell | `n` | fix-Z | fix-X |
|---|---|---|---|
| `zr_fork` | 800 / 3,200 / 12,800 | **100.3% / 85.1% / 100.0%** | 6.6% / 11.7% / 5.6% |
| `zr_pipe` | 800 / 3,200 / 12,800 | **94.8% / 103.7% / 87.5%** | 0.1% / 2.7% / 6.5% |

Registered bars: fix-Z ≥ 70%, fix-X ≤ 30%, |interaction| ≤ 2 pp, `ef_exact_exact` within
1.5 pp of zero. **All four pass**, the last at worst +0.66%.

**The blocks are additive here.** The interaction never exceeds 0.40 pp against a 2 pp bar,
so the alternation contributes essentially nothing — the opposite of the collider cells in
V17, where the leak was −6.4 pp. Whether the blocks interact is itself structure-dependent.

**`ef_exact_exact` landing on zero is what validates the instrument.** Drawing both blocks
from the true conditional gives −0.08% to +0.66% at every `n`. Had it not, the "exact" draws
would have been suspect and the attribution void.

---

## The sharper test dissolved rather than answered

Registered: if the bias is the Z draw's smoothing error, the n^−1/3 exponent should vanish
once Z is fixed and survive when only X is fixed.

| arm | `n` = 800 | 3,200 | 12,800 | slope |
|---|---|---|---|---|
| `pipeline_bartMI` (fork / pipe) | 5.01 / 5.44% | 3.35 / 2.64% | 2.13 / 2.00% | −0.309 / −0.361 |
| `ef_bart_ship` | 4.63 / 5.14% | 2.99 / 2.38% | 1.33 / 1.67% | −0.450 / −0.406 |
| `ef_bart_exact` (fix X) | 4.32 / 5.13% | 2.64 / 2.32% | 1.26 / 1.56% | −0.446 / −0.430 |
| **`ef_exact_ship`** (fix Z) | **−0.01 / +0.27%** | **+0.45 / −0.09%** | **−0.00 / +0.21%** | *no exponent* |
| **`ef_exact_exact`** | −0.08 / +0.66% | +0.32 / +0.06% | −0.01 / +0.21% | *no exponent* |

**Fixing the Z block drives the bias to ~0 at every `n`, so there is no exponent left to
locate.** The registered phrasing — "the exponent should flatten" — presupposed a residual
decay to flatten; there is none. A slope fitted to a bias of ±0.4% is fitted to noise, and
reporting one (the naive fit gives −1.97 in one cell) would be reading the instrument wrong.

`pipeline_bartMI`'s −0.309 / −0.361 **replicates V18 (−0.335) and V19 (−0.311)** a third
time, under a third arm set.

---

## What this means for the fix, and for item 07

**The fix belongs in the covariate block.** For weeks item 07 — a substantive-model-compatible
**exposure** draw — has been the only open build track and the presumed remedy for everything.
It is not the remedy for this: replacing the exposure draw removes **0.1% to 11.7%** of this
bias. Item 07 remains the right fix for the *mixture* failure (V3/V9) and the
non-linear-outcome penalty (V13/V14), which are X-block defects. This is a different defect
in a different block, and it had been quietly folded into the same queue.

**What the Z-block fix has to be is now sharply constrained, and it is not obvious.** The
exact conditional here is available only because the DGP is known. A shippable version needs
a covariate draw that is *correct* rather than merely flexible — and V19 already showed that
"parametric and correctly specified in form" is not enough: `micePmm` and `properZ` were
correctly specified for this `Z₁` and still carried +4 to +6% asymptotic bias. So the target
is not "use a parametric imputer"; something else about those draws is wrong, and V20 does
not say what.

**The practical guidance in `docs/covariate-roles.md` is unaffected** — it reports the
measured effect, which stands — except that "which part of the imputation produces it" can
now be answered.

---

## Caveats (scope of this run)

- **Two cells, one DGP, arrow strengths fixed at 0.60.** The attribution is measured under a
  fork and a pipe with a linear-Gaussian `Z₁` conditional. Where that conditional is
  genuinely non-linear, BART's flexibility is an asset rather than a liability and the split
  could differ — that cell does not exist yet (`ROADMAP.md` loose ends).
- **The instrument is 0.26–0.80 pp better than the pipeline**, probably the missing Z-block
  inner-FCS loop. Small against a ~5 pp effect, but it means the reference here is not
  exactly the shipped code.
- **"Fix the Z block" is not yet a shippable proposal.** The exact conditional uses knowledge
  of the DGP. V19 rules out the obvious substitutes.
- **`m` = 30 fixed at every `n`**, as in V18 and V19 — the same untested confounder for the
  exponent.
- **Nothing here revisits the collider cells**, where V17 measured a −6.4 pp leak. The
  additivity found here is a property of the fork and pipe structures, not a general one.
