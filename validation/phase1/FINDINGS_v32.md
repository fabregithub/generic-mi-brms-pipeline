# V32 — there is no transfer law, and V31's localisation does not survive

**Question.** V31 fitted `measured = a·borrow + b·atten` and concluded the attenuation channel
transfers across outcome families untouched while the borrowing channel does not. Is `a` a
property of the **family** (so a transfer law exists) or of the **arm** (so there is nothing to
find)?

**Outcome: of the arm. Both registered gates fail, and that closes the line as registered.**
In the course of answering it, V31's localisation is **withdrawn**.

## 1. Result

Five instrumented arms × two missingness fractions (0.10, 0.20), `n` = 3200, 500 reps, 18,000
rows. Arms chosen to span the conditioning rather than repeat it — surviving fraction
`|T₁+T₂|/(|T₁|+|T₂|)` from **0.069** (`properZ_ds`) to **1.000** (`noYz`).

| family | `a` (borrow) | `b` (atten) | χ²/df | verdict |
|---|---|---|---|---|
| gaussian | **0.983 ± 0.008** | **1.046 ± 0.008** | **0.24** | fits |
| binomial | 0.854 ± 0.159 | 1.435 ± 0.212 | **17.80** | rejected |
| ordinal3 | 0.940 ± 0.066 | 1.247 ± 0.084 | **4.38** | rejected |

* **P1 NOT MET** — a single `(a, b)` per family does not fit for either non-Gaussian family.
* **P2 NOT MET** — an arm-specific borrowing scale significantly improves the fit for binomial
  (F = 25.95, p = 0.0040) and ordinal (F = 17.54, p = 0.0084), not for gaussian (p = 0.12).

**Under a Gaussian outcome the instrument is essentially exact**: χ²/df = **0.24** across five
structurally different arms and two fractions, both scales within 5% of 1. That part of V31
stands and is now better supported.

## 2. ⚠️ V31's localisation is withdrawn

V31 reported `b` = 0.94–1.03 across all three families and `a` = 0.41 / 0.63 for binomial and
ordinal, and concluded the transfer problem was **entirely** in the borrowing channel. With a
better-conditioned arm set:

| | `a`: V31 → V32 | `b`: V31 → V32 |
|---|---|---|
| binomial | 0.409 → **0.854** | 0.942 → **1.435** |
| ordinal3 | 0.630 → **0.940** | 0.980 → **1.247** |

`a` moved *toward* 1 and `b` moved *away* from it — the opposite assignment. And because the
non-Gaussian models are **rejected**, their nominal standard errors understate the true uncertainty;
inflating by `√(χ²/df)`, as a rejected fit requires, leaves **neither** channel distinguishable
from 1:

| family | `a` from 1 | `b` from 1 |
|---|---|---|
| binomial | 0.22σ | 0.49σ |
| ordinal3 | 0.43σ | 1.41σ |

**So there is no localisation to withdraw to.** The honest statement is that the two-scale model
fails for non-Gaussian outcomes, not that one of its two channels is at fault.

## 3. Why V31 got it wrong — two identifiable causes

**(i) It fitted across a regime where the expansion is broken.** V31 included the 0.40 rung, and
its own §1 measured χ²/df = **115.8** (binomial) and **79.4** (ordinal) there. Coefficients
estimated partly inside a regime the model does not describe are not interpretable. V32 uses the
clean rungs only.

**(ii) It had no well-conditioned arm, so `a` and `b` were collinear.** V31's best surviving
fraction was 0.73 (`noYx`). When every arm cancels heavily, `borrow` and `atten` are nearly
proportional across arms and the two scales trade off against each other — almost any `(a, b)`
on a ridge fits equally well. V32 added `noYz` at surviving fraction **1.000**, where
`pred = borrow + atten` with no cancellation at all, which pins `b` directly. The estimates then
moved by 0.3–0.5.

**The lesson is about identification, not about imputation.** A two-parameter fit needs arms that
differ in the *ratio* of the two terms, not merely in their size. V31 had five points on a ridge
and read the ridge's location as a finding. The diagnostic was available — the surviving fraction
was already being recorded — and was used to explain the `pred/meas` misses but not to ask
whether the fit was identified.

## 4. What this leaves

* **Gaussian:** the first-order-plus-second-order expansion with recorded `d` predicts the bias
  to ~5% per channel. Solid, replicated across V31 and V32 on different arm sets.
* **Non-Gaussian:** no transfer law, no channel localisation, and the raw bias ratio
  `|bias_fam|/|bias_gauss|` spans **0.26–2.2** (V31 §3, a direct measurement, unaffected by the
  above since it involves no model fit). Outcome-family bias transfer is **open**, with two
  candidate reductions and one localisation now refuted.
* **Still true and still the practical point:** only the standard-error half of the transfer is
  exactly computable (`SE_g/SE_fam = σ√E[w]`, 24/24 cells to ≤0.001, V30).

## 5. Reproduce

```bash
N_REP=500 bash run_v32_channels.sh    # 6 cells x 5 arms + oracle, n=3200; 2.9 h measured
Rscript analyze_v32.R
```
