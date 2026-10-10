# V26 — does the record transfer to a non-Gaussian estimand?

**Run 2026-09-11.** 3 `n` levels × 6 cells × 500 reps = 9,000 tasks, 8 arms each,
**14 h 3 min against 20.7 h predicted**, `rc = 0`, zero task errors, all cells at 500/500.
Driver [`run_v26_family.sh`](run_v26_family.sh), analysis [`analyze_v26.R`](analyze_v26.R).

## The question

Every track V0–V25 estimated a **Gaussian linear coefficient by OLS**. That is what made
500-replication runs affordable, and for a Gaussian coefficient it loses nothing. But the
pipeline's real targets are logistic, ordinal `mo()`, spline and mixed models, and the transfer
of the whole record to those classes had never been tested — it was assumed.

Three cells were re-run with `Y` drawn from a **logistic** model on the *same* linear
predictor, each paired with its Gaussian twin **in one invocation with an identical arm set**.
Same seed ⇒ byte-identical exposures; only `Y` differs. The estimand becomes the focal
exposure's log-odds-ratio, still 0.40.

> **Currency.** A log-odds-ratio and a mean difference are not on the same scale, so **relative
> bias is not comparable across families**. Everything below is `bias/SE` and coverage.

## Criteria

| # | criterion | 800 | 3200 | 12800 |
|---|---|---|---|---|
| C1 | arm ordering by `\|bias/SE\|` survives | void | void | **MET** |
| C2 | each arm's `\|bias/SE\|` within 2× of its twin | no | no | **NOT MET** |
| C3 | `pipeline_noYx` still biased toward the null | **MET** | **MET** | **MET** |
| C4 | the scale defect by family | descriptive | | |

## 1. Direction and ordering transfer; size does not

At `n` = 12800, where the defects are large enough to rank:

| arm | `bin_fork` | `zr_fork` | ratio |
|---|---|---|---|
| `leftcens_prestep` | −2.25 (cov 0.37) | −4.89 (cov 0.00) | 0.46 |
| `pipeline_noYx` | −2.26 (0.39) | −4.76 (0.00) | 0.48 |
| `pipeline_noYz` | +3.60 (0.03) | +7.52 (0.00) | 0.48 |
| `pipeline_bartMI` | +0.05 (0.97) | +0.43 (0.94) | — both null |
| `pipeline_properZ` | −0.04 (0.98) | +0.01 (0.95) | — both null |
| `pipeline_properZ_ds` | +1.75 (0.59) | +1.61 (0.65) | 1.09 |

Rank correlation of arms by `|bias/SE|`: **+0.89 / +0.96 / +0.96** across the three cell pairs.
**Ordering survives** — the family does not change which imputation strategy is better, which
was the outcome that would have invalidated the most prior guidance.

**Sizes do not transfer within 2×, and the miss is systematic rather than noisy.** Over the 30
cell×arm combinations where the Gaussian twin has a real defect (`|bias/SE|` > 0.5):

| arm | cells | ratio binomial / gaussian |
|---|---|---|
| `leftcens_prestep` | 9 | **0.444** [0.349, 0.474] |
| `pipeline_noYx` | 9 | **0.475** [0.440, 0.524] |
| `pipeline_noYz` | 6 | **0.482** [0.471, 0.513] |
| `pipeline_properZ_ds` | 6 | **0.969** [0.109, 1.088] |

Median over all 30: **0.474**, IQR 0.445–0.493.

## 2. The practical consequence: the Gaussian record is CONSERVATIVE

A defect that reads **1.0** `bias/SE` in a Gaussian validation reads about **0.47** in the
logistic twin. A binary outcome carries less information per observation, so the standard error
is larger and the same structural error consumes a smaller share of it.

**So existing guidance can be quoted for logistic analyses as an upper bound, not as an
estimate.** That is the useful form of the result: C2 failing is not "the record does not
transfer", it is "the record transfers and overstates the damage by about a factor of two".

Coverage behaves the same way — the defects still destroy it, less severely: `leftcens_prestep`
0.00 → 0.37, `noYx` 0.00 → 0.39, `noYz` 0.00 → 0.03 at `n` = 12800.

## 3. One defect does not attenuate, and the split is mechanistic

The three arms at ~0.47 are all defects **mediated by `Y`'s role in the imputation** — `Y`
dropped from a block, or the uncongenial pre-step. The pre-v1.6.0 **covariate-block scale
defect** sits at **0.969** in the causally-active cells (1.09, 0.89, 1.13, 0.94): *unchanged* by
the family.

That is what the mechanism predicts. A defect that acts through `Y`'s contribution to the
imputation is diluted when `Y` carries less information; a defect in the covariate–exposure
relationship does not involve `Y` at all and so is untouched.

**The exception to the exception:** in the inert-covariate cell (`bin_base`) the scale defect
*does* attenuate hard (−0.17 against −1.05, ratio 0.16). With no causal path, what little the
defect does is routed through `Y`, so it behaves like the others. This is consistent but it was
not predicted in advance and should be treated as an observation, not a law.

## 4. The root claim survives (C3)

`pipeline_noYx` is biased **toward the null** in all three binomial cells at all three sample
sizes (−0.41 to −2.35, against a Monte-Carlo error of 0.04). V16's root claim — that dropping
`Y` from the exposure draw attenuates the exposure–response coefficient — is a property of
congeniality, not of the Gaussian family.

## ⚠️ Two of these gates were badly constructed, and I am recording that

**C2's informativeness guard was far too weak.** It skipped a pair only when the Gaussian arm
read `|bias/SE|` < 0.1, so ratios between two *null* arms dominated the verdict: the headline
failure at `n` = 800 was **23.56×** between arms reading **0.00 and −0.11**. A ratio between two
near-zero numbers is undefined, not evidence. Restricting the comparison to arms whose Gaussian
twin exceeds 0.5 gives a worst case of 5.49×, which is the number the verdict now uses. **Both
are printed by `analyze_v26.R`** — the as-registered figure stays on the record, because a bar
must not move after the data are seen.

**C1's guard had the same flaw.** It asked only that the cell's Gaussian spread exceed 0.3,
which a cell can satisfy while every individual arm sits at ±0.1. It now requires at least three
arms to individually have a defect, which is why `n` = 800 and 3200 are marked **void** rather
than failed.

This is the third and fourth time in this programme that a gate has resolved for the wrong
reason (after V22's `bart ≥ +3%` and V23's uninformative CI). PLAN §11 exists because of the
first two, and it did not prevent these — the lesson is that an informativeness guard needs a
*floor on the quantity being compared*, not on a spread or a difference.

## What this does NOT close

- **Step 6's pooling machinery still never fires.** This track pools with Rubin's rules in the
  harness, not through Step 6, so the finite-`m` variance correction on a support-respecting
  transform and its bimodality gate remain unexercised. That needs real `brms` fits through
  Steps 3–6 and is separate, smaller work.
- **Only the binomial family.** Ordinal `mo()`, spline and mixed models are untested. The
  mechanism here (less information per observation ⇒ attenuated `bias/SE`) suggests the same
  direction, but that is a conjecture.
