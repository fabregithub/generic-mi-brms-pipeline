# V24 — the DAG-factorised exposure draw, all four cases at once

**Run 2026-09-09.** 3 `n` levels × 6 cells × 200 reps = 3,600 tasks, 7 arms each, **320 min against
8.2 h predicted**, `rc=0`, zero task errors. Driver
[`run_v24_dagdraw.sh`](run_v24_dagdraw.sh), analysis [`analyze_v24.R`](analyze_v24.R), algorithm
[`R/dag_draw.R`](R/dag_draw.R), 32 assertions in [`test_v24_dag.R`](test_v24_dag.R).

## Criteria

| # | criterion | n = 800 | n = 3200 | n = 12800 |
|---|---|---|---|---|
| C1 | `dag_true`: `\|bias/SE\|` < 0.3 in all six cells | **MET** | **MET** | **MET** |
| C2 | child factor gains ≥ 1.0 `\|bias/SE\|` in both `pipe_nl` cells; no-op under fork | *not met — my bar, see below* | **MET** | **MET** |
| C3 | `dag_lin` strictly worse than `dag_Yonly` in the `pipe_nl` cells | **HOLDS** | **HOLDS** | **HOLDS** |
| C4 | `*_tot` cells recover 0.700 / 0.830, not `b₁` = 0.400 | **MET** | **MET** | **MET** |
| C5 | direction of the incongenial draw's bias | descriptive — 6/6 toward null at every `n` |

**C2's n = 800 miss is a registration error, not a result.** The bar asks for a ≥ 1.0
improvement in `|bias/SE|`, but at `n` = 800 `dag_Yonly` only *has* 0.48 (direct) and 0.84
(total) of bias to remove — so the bar was arithmetically unreachable in the cell it was
registered against, whatever the algorithm did. `dag_true` took both cells to `|bias/SE|` ≤
0.08, which is everything available. This is the same class of mistake as V16's `combined` bar:
a threshold in one currency registered against a cell whose value in that currency is bounded
below the threshold. The gate is correct at `n` = 3200 and 12800, where the bias has grown past
1.0 SE.

## 1. The factorisation is right — in all four cases simultaneously, at every `n`

With the child model supplied, `|bias/SE|` never exceeds **0.13** anywhere:

| cell | `dag_true` `bias/SE` (800 / 3200 / 12800) | coverage |
|---|---|---|
| `dag_fork` | +0.01 / −0.07 / +0.05 | 0.980 / 0.925 / 0.935 |
| `dag_pipe_dir` | −0.13 / −0.06 / −0.12 | 0.970 / 0.950 / 0.945 |
| `dag_pipe_tot` | −0.11 / −0.09 / −0.05 | 0.985 / 0.940 / 0.925 |
| `dag_pnl_dir` | −0.08 / +0.05 / −0.06 | 0.960 / 0.940 / 0.945 |
| `dag_pnl_tot` | −0.08 / +0.01 / +0.02 | 0.945 / 0.940 / 0.930 |
| `dag_collider` | +0.06 / +0.07 / +0.03 | 0.905 / 0.950 / 0.975 |

Not a vacuous pass: the incongenial draw in the same cells runs **−1.15 to −8.24** at
`n` = 12800, with coverage 0.000 in five of six. `dag_true` is indistinguishable from `oracle`
(complete data) in every cell at every `n`.

## 2. Omitting a mediator from the exposure draw is **asymptotically biased** — and that is a reachable configuration, not the default

> **CORRECTION, added after the run (2026-09-09).** This section first claimed `dag_Yonly` was
> "the class the shipped `leftcens` path belongs to". **That was wrong**, and checking it is
> what produced §2b. `leftcens` regresses `x` on its whole predictor set — outcome, covariates
> and other exposures together — so where everything is jointly Gaussian that single regression
> **is** the exact conditional. `dag_Yonly` is a rung I built for the ladder: it conditions on
> parents only and multiplies by `p(Y | x, ·)`, which under a pipe drops `Z` *entirely*. The two
> are different draws and they measure different things. Verified directly at `n` = 3200, 100
> reps: under a linear pipe the shipped form is **−0.02%** (direct) and **+0.19%** (total),
> while `dag_Yonly` is −6.14% and −6.67%.
>
> So the numbers below are correct, but they answer **"what if a covariate that is informative
> about the exposure is absent from the exposure draw?"** — not "what does the pipeline do by
> default". That question still matters, because **the configuration reaches it**: see the end
> of this section.

`dag_Yonly` — parents plus `p(Y | x, Pa(Y))`, no child factor — has relative bias that is
**flat in `n`**, so its `bias/SE` **grows as √n**:

| cell | rel. bias (800 / 3200 / 12800) | `bias/SE` | coverage | exponent |
|---|---|---|---|---|
| `dag_pipe_dir` | −7.48 / −6.18 / −6.18% | 0.58 → 0.97 → **1.94** | 0.945 → 0.830 → **0.535** | 0.43 |
| `dag_pipe_tot` | −7.37 / −6.90 / −6.68% | 1.04 → 1.97 → **3.80** | 0.830 → 0.480 → **0.040** | 0.47 |
| `dag_pnl_dir` | +7.52 / +9.25 / +8.82% | 0.48 → 1.19 → **2.28** | 0.940 → 0.780 → **0.375** | 0.56 |
| `dag_pnl_tot` | +5.10 / +5.59 / +5.76% | 0.84 → 1.86 → **3.83** | 0.895 → 0.565 → **0.050** | 0.55 |
| `dag_fork` | +0.06 / −0.40 / +0.15% | ≤ 0.07 | 0.980 / 0.925 / 0.935 | — |
| `dag_collider` | +0.71 / +0.44 / +0.16% | ≤ 0.08 | 0.925 / 0.945 / 0.975 | −0.03 |

*(Exponent = slope of `log |bias/SE|` on `log n`. 0.5 is what a constant bias against an
`n^−1/2` SE produces. `dag_fork`'s 0.84 is fitted through `|bias/SE|` ≤ 0.07 — a slope through
noise, and reported only so it is not mistaken for a signal.)*

**Three things follow.**

**(a) The bias here comes from the covariate being *absent*, not from the arrow being curved.**
`dag_pipe_dir` and `dag_pipe_tot` have a **linear-Gaussian** arrow, so `u` = 0 by construction,
and they still carry **6–7%** flat across a 16× range in `n`. So "unrepresentable curvature"
cannot be what generates *this* bias — discarding an informative covariate is enough on its
own. (§2b shows the shipped draw, which does *not* discard it, is exact in these same cells.)

**(b) Coverage collapse is a large-`n` phenomenon and is invisible at typical simulation
sizes.** `dag_pipe_tot` reads 0.830 at `n` = 800 — poor but arguably survivable — and **0.040**
at `n` = 12800. A validation exercise run only at `n` ≈ 1000 would have called this cell
acceptable.

**(c) It does not touch the fork or the collider.** Both stay at `|bias/SE|` ≤ 0.08 with
coverage 0.93–0.98 across all three `n`. The defect is specific to the covariate being a
**child** of the exposure *and* used in the analysis or informative about it.

**(d) How a real configuration reaches this cell.** The pipeline's automatic X-block predictor
set was `outcome + use_in_model covariates + other exposures`. A mediator in a **total-effect**
analysis is marked `use_in_model = FALSE` — correctly, because adjusting for it would give the
direct effect — and it therefore **dropped out of the exposure draw too**. Same for anything
marked `use_as_auxiliary = TRUE`, which the vignette documented as a sub-1-percentage-point
curiosity. It is not: for a mediator it is this cell, at 6–7% flat in `n` with coverage
reaching **0.040**. Fixed in `00_censored_exposure.R` (auxiliaries and `role = "auxiliary"` now
enter `auto_preds`); the fix needs no assumption, because conditioning *linearly* on the
covariate is enough whenever the relationship is linear.

## 2b. What the shipped draw actually does — measured

`n` = 3200, 100 reps, 40% non-detects, arms paired on data and seed. "Shipped form" is an
interval-censored Gaussian regression of `x` on `(Y, other X, Z)` followed by a truncated draw
— what `leftcens` does, with its shash margin an identity on this Gaussian DGP.

| case | truth | **shipped form** | `dag_Yonly` | DAG (child known) |
|---|---|---|---|---|
| fork | 0.400 | **+0.20%** (0.04) | +0.21% (0.04) | +0.21% (0.04) |
| pipe, direct | 0.400 | **−0.02%** (0.00) | −6.14% (0.96) | −0.25% (0.04) |
| pipe, total | 0.700 | **+0.19%** (0.06) | −6.67% (1.89) | −0.05% (0.01) |
| `pipe_nl`, direct | 0.400 | **+22.03%** (2.44) | +9.07% (1.16) | +0.02% (0.00) |
| `pipe_nl`, total | 0.830 | **+17.92%** (6.13) | +5.71% (1.89) | −0.00% (0.00) |

`bias/SE` in parentheses. Three things worth stating plainly:

1. **The default draw is exact for a confounder and for a linear mediator**, both estimands. The
   single regression absorbs the covariate correctly when the joint is Gaussian.
2. **It is badly biased for a *non-linear* mediator: +22.0% and +17.9%.** This is V22's ~20%
   defect, reproduced with an instrument that isolates it, and the DAG child factor takes it to
   ~0.
3. **Under a non-linear mediator, conditioning on the covariate linearly is *worse* than
   dropping it** — +22.03% against +9.07%. That is counter-intuitive and it is the same
   mechanism as §3: a linear conditional forced onto a curved relationship does positive damage,
   where omission merely discards information. It also gives a **zero-assumption mitigation**:
   removing a non-linear mediator from `censored_exposure$predictors` roughly halves the bias.
   Not a fix — 9% is still 9% — but strictly better than the default and it costs nothing.

## 3. A misspecified child factor is worse than none — and *more flexible is worse*

C3 holds at all three `n`, and the gap widens with `n`:

| cell | `dag_lin` minus `dag_Yonly`, in `\|bias/SE\|` |
|---|---|
| `dag_pnl_dir` | +0.80 → +1.56 → **+3.13** |
| `dag_pnl_tot` | +2.72 → +5.43 → **+10.92** |

And it is **monotone in the child model's flexibility** — at `n` = 12800:

| cell | no child | linear | quadratic | cubic | true |
|---|---|---|---|---|---|
| `dag_pnl_dir` | +8.82% | +25.75% | +31.07% | +37.38% | **−0.21%** |
| `dag_pnl_tot` | +5.76% | +22.27% | +29.04% | +41.66% | **+0.04%** |

Coverage is 0.000 for every fitted form in both cells. The **linear `pipe` controls** confirm
this is misspecification and not the factor itself: there the same linear form is correct by
construction, and it *improves* `|bias/SE|` by 0.92–3.78 against `dag_Yonly`.

> **This is the binding constraint on the fix, and it is sharper than "the child model is an
> assumption".** The child factor must be evaluated *below* the LOD, where `THEORY.md` §0b
> measured that the arrow's curvature is not identified. V24 measures the consequence: a
> **more flexible basis extrapolates worse, monotonically**, so the usual defence — fit
> something flexible and let the data decide — is exactly backwards here. The algorithm cannot
> ship with a fitted-by-default child model. It ships as a **sensitivity procedure over a
> declared form**, and a run that declares the form badly is worse off than one that never
> used the algorithm at all.

## 3b. Post-run: the feature is salvageable, because the *basis* is what has to be right

**Added 2026-09-09, after the run.** §3 leaves a question the track did not ask: every child
model V24 tested was either a fitted polynomial (all worse than nothing) or the true arrow **with
its true coefficients** (exact). A shipped feature can never have the second, so §3 as it stood
implied the option could not help anyone. Tested directly — `n` = 3200, 200 reps, `pipe_nl`
direct effect, child model **fitted on the observed rows only**:

| declared basis | direct effect | `bias/SE` | total effect | `bias/SE` |
|---|---|---|---|---|
| `~ tanh(1.8·x) + I(x²)` — true basis, coefficients fitted | **−0.51%** | −0.07 | **−0.14%** | −0.05 |
| `~ tanh(x) + I(x²)` — right form, wrong internal constant | **+3.22%** | +0.44 | **+1.19%** | +0.39 |
| *none* (omit the child factor) | +8.47% | +1.09 | +5.59% | +1.86 |
| `~ tanh(1.8·x)` — right family, quadratic term dropped | +20.71% | +2.33 | +15.52% | +4.97 |
| `~ x` — a straight line (V24's `dag_lin`) | +25.43% | +2.67 | +22.16% | +7.34 |

The ordering is identical for both estimands, which is what makes it a property of the declared
basis rather than of one cell.

**The coefficients do not need to be known.** With the right basis and coefficients estimated
above the LOD, the arm is indistinguishable from the oracle. So §3's finding is not "the child
model must be known", it is **"the child model's *form* must be right"** — a materially weaker
requirement, and one an analyst with mechanistic knowledge can sometimes meet.

**And being approximately right is enough.** `tanh(x)` where the truth is `tanh(1.8x)` — a wrong
parameter *inside* a non-linear function, the sort of error a real declaration would make — costs
3.7 pp, not 25. The option degrades gracefully in the direction users will actually miss.

**But an incomplete basis is still worse than nothing** (+20.71% against +8.47%). So the failure
mode is *missing structure*, not *insufficient flexibility* — which is why every polynomial fails
and why "fit something flexible" is the wrong instinct. It also means `child_form` cannot be
chosen by in-range goodness of fit, since all five rows above fit the observed region acceptably
and differ by 26 percentage points below the LOD.

**Consequence for the pipeline.** The option takes a declared **formula**, not a flexibility
setting. `phase1/derive_dag_basis.R`.

## 4. The collider does not need the algorithm

`dag_Yonly` is already at +0.06 to +0.16% (`bias/SE` ≤ 0.08, coverage 0.925–0.975) at every
`n`. With `Z₁` off every `X`–`Y` path and absent from the correct analysis model, ignoring it
in the imputation is **congenial**. Every child form is also fine there (`dag_lin` +0.04,
`dag_quad` +0.06, `dag_cubic` +0.23 at `n` = 12800) precisely *because* the factor barely
matters. The design test measured the payoff at **1.4 ± 0.1%** of the SE (paired).

**This corrects `THEORY.md` §6b's four-cases table as first written**, which listed the
collider's child factor as needed. It is *legitimate* there, not *necessary* — and spending an
unidentifiable assumption on 1.4% of the SE is not a trade worth making.

## 5. The fork is a no-op, exactly

All four `dag_*` arms are bit-identical under the fork at all three `n` (`dag_spec()` gates the
child factor off because `Z₁` is a *parent*). Registered as C2's second clause and asserted in
`test_v24_dag.R`; it holds by construction, and the test exists so a future refactor cannot
quietly turn the factor on where it does not belong.

## 6. Direction: 6/6 toward null, at every `n`

C5 was descriptive and stayed descriptive. The incongenial draw is **conservative in every cell
at every `n`** — including the two total-effect cells, which no prior track measured. Extends
`PLAN` §11b's direction ledger by 18 cells. Note this is the *incongenial* draw; the
*misspecified-child* arms of §3 are strongly **anti-conservative** (+22% to +42%), so
"attenuation" is a property of omitting information, not of getting the imputation wrong in
general.

## 3c. Post-run: why the cheap implementation route does not work

**Added 2026-09-09.** §1 validates the *algorithm*; this records what happened when it was
implemented in the pipeline, because the failure is informative and the working version is not
obvious.

The grid sampler builds the product of the three factors explicitly. A cheaper route is to ask
`leftcens` for K candidate draws — which already cover the parent and outcome factors, with its
skew-aware margin intact — and reweight them by the child factor. That avoids reimplementing the
margin, which caused four separate defects in `FINDINGS_v13.md`. Measured against the grid
sampler on identical data, `n` = 3200, 80 reps:

| | linear pipe | `pipe_nl`, direct | `pipe_nl`, total |
|---|---|---|---|
| weight `p(C \| x)` | +7.61% | +27.37% | +13.54% |
| weight `p(C \| x, Y)` | **+0.46%** | +32.53% | +13.02% |
| grid sampler | −0.60% | −0.41% | −0.09% |

**Finding 1 — the weight is `p(C | x, Y)`.** The proposal already conditions on `Y`, and
`p(Y|x,C)·p(C|x) = p(Y|x)·p(C|Y,x)`, so the residual factor carries `Y`. This is a real
correction and it fixes the linear pipe outright (+7.61% → +0.46%, against the grid's −0.60%).
Note it applies to **every** child, not only to a collider: conditioning on `Y` couples `C` and
`Y` whenever `Y` is a descendant of `C`.

**Finding 2 — and it still fails where the feature is needed.** The child must be *excluded* from
the proposal's predictor set, or its information enters twice (measured at +24.48%). But with it
excluded, `leftcens` regresses `x` on `Y` **linearly**, while marginalising the child out makes
`Y = b₁x + γ₁·g(x) + …` — genuinely non-linear in `x`. So the **proposal itself** is
misspecified, and importance reweighting can only *add* a missing factor, never repair a wrong
proposal; that needs the full target/proposal ratio and `leftcens` does not expose its density.

> **The contrast that isolates it:** same code, same weight, linear proposal **works** (+0.46%)
> and non-linear proposal **does not** (+32.53%). An earlier version of this section claimed the
> limitation was structural to SIR — that was wrong, and the linear-pipe result is what refuted
> it. The limitation is specific to the proposal being misspecified, which under non-linearity it
> unavoidably is.

**What a correct implementation needs.** The grid sampler, whose outcome factor conditions on the
child — where `Y` *is* linear in `x`. Which means evaluating the analysis model's linear predictor
as a function of a candidate exposure value, per row. `phase1/derive_sir_weight.R`,
`phase1/derive_sir_proposal.R`. Tracked in `ROADMAP.md`.

**Process note worth keeping.** Three successive conclusions about this code were wrong — the
child left in the proposal, the weight without `Y`, and then the claim of structural
impossibility — and every one produced plausible numbers. What caught all three was running two
independent implementations of the same target against each other. Unit tests, smoke tests and
inspection caught none of them.

## What V24 does not settle

- **Both blocks missing.** Every cell here has fully observed covariates (`mcar_frac` = 0), to
  isolate the exposure draw and keep V21's covariate-draw ladder out as a confound.
  `dag_impute_datasets()` refuses missing covariates with a named error. Open in `ROADMAP.md`.
- **Where the declared form comes from.** V24 shows a wrong form is harmful and a right form is
  exact. It says nothing about how a user is supposed to choose, beyond §0b's finding that the
  data cannot tell them.
- **Still a Gaussian linear estimand.** Same gap as the other 23 tracks: `fit_lm_estimand()`
  throughout, and Step 6's pooling machinery has never fired.
