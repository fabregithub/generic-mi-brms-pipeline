# V31 — the bias transfer localises to ONE channel, and it is not an expansion-order problem

> ## ⚠️ SUPERSEDED IN PART BY V32 (2026-10-19) — read that first
>
> **The localisation claimed below is withdrawn.** V31 fitted the channel scales on four arms
> whose surviving fractions topped out at 0.73, and across the 0.40 rung where its own §1
> measured the expansion broken (χ²/df 115.8). With a fifth arm at surviving fraction 1.000 and
> the dirty rung excluded, `a` moved 0.41 → 0.85 (binomial) and 0.63 → 0.94 (ordinal) while `b`
> moved *away* from 1 — the opposite assignment — and the non-Gaussian fits are **rejected**
> (χ²/df 17.8, 4.4), leaving neither channel distinguishable from 1.
>
> **`a` and `b` were not identified here.** With only ill-conditioned arms the two terms are
> nearly proportional and the scales trade off along a ridge. §3's table below is the location of
> that ridge, not a measurement.
>
> **What stands:** the Gaussian result (now χ²/df = 0.24 on five arms), the refutation of the
> expansion-order explanation (§1), and the raw bias-ratio range 0.26–2.2, which involves no
> model fit. `phase1/FINDINGS_v32.md`.


**Question.** V30 left the bias half of outcome-family transfer open: a defect's absolute bias
is not family-invariant (0.26×–2.2×) and no reduction to a function of the information ratio
survived. The harness was instrumented to record the imputation error `d` directly
(`R/instrument_d.R`), and V31 asks whether the bias it predicts is the bias that occurs.

**Outcome.**

1. **Under a Gaussian outcome the instrument is essentially exact** — allowing one scale per
   channel gives `a = 0.947 ± 0.016`, `b = 1.028 ± 0.020`, **χ²/df = 1.1** on 12 points.
2. ~~**The attenuation channel transfers across families untouched** (`b` = 0.94–1.03 everywhere).~~
   **WITHDRAWN — V32.** With a well-conditioned arm set `b` reads 1.05 / 1.44 / 1.25 and the
   non-Gaussian fits are rejected outright.
3. **The borrowing channel does not**: `a` falls **0.95 → 0.63 → 0.41** for gaussian → ordinal
   (K=3) → binomial, ordered by how much information the outcome carries.
4. **Registered P1 failed on its chosen statistic, but the mechanism it named is real.**
   `pred/meas` is flat in `|d|` — yet it is dominated by a fixed, arm-specific cancellation
   amplification (×1.4–3.8). Separating the channels shows the borrowing scale DOES drift with
   `|d|` (binomial 0.630 → 0.372), on top of a `|d|`-independent family effect. See §1.

## 1. P1 failed, but NOT because the mechanism was wrong — the STATISTIC was

Registered P1: the miss is an expansion-order effect, so `pred/meas` rises toward 1 as `|d|`
shrinks. Measured `pred/meas` by rung (0.10 / 0.20 / 0.40) is flat:

| arm | gaussian | binomial | ordinal3 |
|---|---|---|---|
| `properZ_ds` | 0.86 0.84 0.81 | 0.19 0.03 −0.04 | 0.29 0.35 0.22 |
| `bartMI` | 0.69 0.69 0.72 | · −3.11 · | · · −4.24 |
| `noYx` | 1.11 1.10 1.08 | 2.29 2.37 1.68 | 2.54 1.28 1.43 |

**On that statistic, P1 is refuted. But the statistic cannot see what P1 was about.** Because
`pred = borrow + atten` with the two of opposite sign, a per-channel error `ε` appears in
`pred/meas` multiplied by the cancellation factor `Σ|·| / |sum|`, which is an **arm property and
independent of `|d|`**. With the fitted Gaussian channel errors (5.8% and 5.5%):

| arm | `|sum|/Σ|·|` | implied `pred/meas` | observed | amplification |
|---|---|---|---|---|
| `noYx` | 0.73 | 1.086 | 1.083 | **×1.4** |
| `bartMI` | 0.24 | 0.812 | 0.724 | ×3.1 |
| `properZ_ds` | 0.19 | 0.770 | 0.811 | ×3.8 |

Six-percent channel errors produce 9–23% ratio errors, and the implied values reproduce the
observed ones (`noYx` to 0.003). **So `pred/meas` is dominated by a fixed amplification that the
registered hypothesis said nothing about, and a real order effect underneath it is invisible.**

**And the order effect IS there.** Fitting the channel scales PER RUNG, which separates them:

| family | d_rms 0.20 | 0.32 | 0.54 | χ²/df by rung |
|---|---|---|---|---|
| gaussian | a = 0.964 | 0.962 | 0.942 | 0.1 / 0.5 / **2.5** |
| binomial | a = **0.630** | 0.451 | **0.372** | 0.5 / 7.8 / **115.8** |
| ordinal3 | a = 0.559 | 0.817 | 0.562 | 1.2 / 1.3 / **79.4** |

The binomial borrowing scale drifts monotonically with `|d|`, and the two-scale model **fits at
the smallest rung in every family** (χ²/df 0.1 / 0.5 / 1.2) while failing badly at the largest.
That is an expansion-order signature. For a Gaussian outcome `g = y − η` is linear so
`g_ηη = 0` *identically* and no such term exists — which is exactly why the Gaussian barely
drifts (0.964 → 0.942) while the families drift hard.

**But order is not the whole story either.** Extrapolating the binomial `a` to `|d| → 0` lands
near 0.67–0.78, not 1, and at the smallest rung `a = 0.630 ± 0.103` is **3.6σ below 1**. So two
distinct effects were conflated in one number:

* a **`|d|`-dependent** part — expansion order, real, and what P1 was about;
* a **`|d|`-independent** part — a genuine family effect in the borrowing channel, which more
  terms will not remove.

**Where the wrong prediction came from, concretely.** The self-test that motivated P1 used
*single-channel* perturbations (`indep` ≈ pure attenuation, `withY` ≈ mostly borrowing). There
`pred ≈ one channel`, the amplification is ≈1, and `pred/actual` does track the channel error —
so the statistic was valid *in the setting where it was calibrated*. I carried it to the real
arms, where the channels cancel and the amplification is 1.4–3.8×, without re-checking that the
mapping still held. **The instrument was already recording `borrow` and `atten` separately, so
that check was free and I did not make it.** The lesson is not about expansions: *a statistic
validated on a special case must be re-validated when the case changes, especially when the
ingredients for the check are already in hand.*

## 2. Why the per-arm ratios differ: cancellation, not mechanism

`pred = borrow + atten`, and the two have opposite signs. Measuring how much survives
(`|sum| / Σ|·|`) on the Gaussian cells:

| arm | borrow | atten | pred | `|sum|/Σ|·|` | ratio |
|---|---|---|---|---|---|
| `noYx` | −0.0766 | +0.0121 | −0.0645 | **0.73** | **1.08** |
| `properZ_ds` | −0.0373 | +0.0546 | +0.0173 | 0.19 | 0.81 |
| `bartMI` | −0.0148 | +0.0243 | +0.0095 | 0.24 | 0.72 |

**The less the channels cancel, the better the prediction.** A ~5% error in either channel is
invisible where one dominates and large where they nearly cancel. That is conditioning, not a
wrong mechanism — and it is why the raw `pred/meas` ratio is a misleading summary of this
instrument.

## 3. The localisation: `measured = a·borrow + b·atten`

Fitting one scale per channel, inverse-variance weighted, 12 points per family (4 arms × 3
fractions), 2 parameters:

| family | `a` (borrow) | `b` (atten) | χ²/df free | χ²/df at `a=b=1` |
|---|---|---|---|---|
| gaussian | **0.947 ± 0.016** | **1.028 ± 0.020** | **1.1** | 5.8 |
| binomial | **0.409 ± 0.098** | 0.942 ± 0.152 | 27.9 | 136.8 |
| ordinal3 | **0.630 ± 0.091** | 0.980 ± 0.130 | 22.1 | 63.0 |

**Gaussian fits.** χ²/df = 1.1 means the instrument, with two ~5% corrections, accounts for the
Gaussian bias of four structurally different arms across three missingness fractions. The
decomposition is real.

**The attenuation channel is family-invariant.** `b` is within its error of 1 in all three. That
is the term `−E[w x d']θ₀ − E[w d d']θ₀`, and it transfers exactly, which is itself a result:
V26's "same absolute bias under both families" assumption is correct **for attenuation** and
wrong only for borrowing.

**The borrowing channel is over-predicted under a non-Gaussian outcome**, by 1.6× (ordinal) and
2.4× (binomial). A plausible mechanism, **not yet tested**: `borrow = E[d·g]` is the covariance
between the imputation error and the analysis score `g`. The imputation sees `Y` on its raw
numeric scale (0/1, or 1..K) while `g` is the GLM score — so the imputer borrows less of what
the likelihood actually weights, and the gap should widen as the two scales diverge. That is
**congeniality in scale, in `Y` rather than in the exposure** (`THEORY.md` §1a(ii)), and the
observed ordering (gaussian → ordinal → binomial) is the ordering of that divergence.

**But a single scale does not fully fix the families either** (χ²/df 27.9 and 22.1). So there is
residual structure in the borrowing channel beyond a constant factor, and the mechanism above
predicts a *specific* factor rather than a free one. Both are open.

## 4. What is settled, and what to do next

**Settled.** The bias is `borrow + atten`; `atten` transfers across outcome families; the
Gaussian instrument is accurate; the transfer problem is **entirely in the borrowing channel**.
That is a much smaller target than "the bias half has no theory".

**Not settled.** Why `a` takes the values it does. The congeniality-in-scale hypothesis above is
the obvious candidate and is directly testable: re-run the imputation feeding it a `Y` on the
likelihood's own scale (the working response, or the latent for an ordinal outcome) and `a`
should move toward 1. That changes the imputer, not the instrument, and it would also say
something shippable — that for a non-Gaussian outcome the imputation model should carry `Y` on
the scale the analysis model uses.

**On extending the expansion.** §1 shows order is *part* of the limitation — the binomial
borrowing scale drifts 0.630 → 0.372 over the `|d|` ladder, and for a Gaussian outcome the term
responsible (`g_ηη`) is identically zero, which is why only the families drift. A third-order
term would therefore buy something at large `|d|`. But it would **not** reach `a = 1`: at the
smallest rung the binomial scale is already 3.6σ below 1. Order and family are separate
deficits, and only the second one is interesting — so the next move is the
congeniality-in-scale test above, not more terms.

## 5. Reproduce

```bash
RECORD_D=1 N_REP=500 bash run_v31_instrument.sh   # 9 cells x 4 arms + oracle, n=3200; 2.9 h measured
Rscript analyze_v31.R
```

Instrument in `R/instrument_d.R` (self-test for it in this file's §1 commentary and in
`run_v31_instrument.sh`'s header); cells `fr{10,20}_{gau,bin,or3}` in `R/robustness.R` (axis
`"dfrac"`); results in `results/v31_n3200_*`.
