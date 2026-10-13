# V28 — the arrow sweep: is the "fork residual" real?

**Question.** `THEORY.md` records the scale-defect chain closing to **1.011** on `zr_pipe` and
missing by **1.147** on `zr_fork`, and calls the fork miss "a real residual rather than noise".
Is it?

**Outcome: no — not as recorded.** The fork/pipe contrast at the canonical arrow is an artefact
of measuring at `n = 800`. What survives is a different thing: a trend in the arrow that appears
in the **fork sweep only**, which acquits the exposure conditional's coefficient and implicates
`Var(logX1)`. V29 tests that directly.

## 1. Design

The two arms are `pipeline_properZ` and `pipeline_properZ_ds`, which differ **only** by the
pre-v1.6.0 scale defect and run on byte-identical data, so the wrapper effect is measured
**paired** and the anchor is hit by construction. Three earlier attempts to rebuild an
independent instrument all missed the pipeline's own level by 31%–4×; this design avoids the
problem instead of fighting it.

The sweep separates two things the fork confounds:

| | what the arrow moves | `Var(logX1)` |
|---|---|---|
| fork (`delta_zx`) | the exposure conditional's `Z1` coefficient **and** `Var(logX1)` | `sd_x² + delta_zx²` |
| pipe (`delta_xz`) | the coefficient only | `sd_x²`, fixed |

So a trend in the pipe sweep implicates the coefficient; a trend in the fork sweep only
implicates `Var(logX1)`.

Three sample sizes, 300/300/200 reps: `n` = 800, 3200, 12800.

## 2. The ratio is not stable in `n` — and that is the whole of the recorded contrast

`derived / measured`, canonical arrow `delta = 0.6`:

| cell | n=800 | n=3200 | n=12800 | derived |
|---|---|---|---|---|
| `fkd06` (fork) | **1.103** | 1.006 | **0.992** | +5.02 |
| `ppd06` (pipe) | **0.992** | 0.943 | **0.953** | +5.58 |

At `n = 800` this reproduces the recorded 1.147 / 1.011 (the small difference is seeding). At
`n = 12800` **the fork agrees better than the pipe.** The measured effect climbs with `n` toward
the asymptote the closed form describes — fork +4.55 → +4.99 → +5.06 against derived +5.02 —
because **the chain is an asymptotic statement and `n = 800` is not asymptotic.** Every prior
derived-vs-measured comparison in this programme was made at `n = 800`.

**The "fork residual" at the canonical cell does not exist.** It was a finite-sample gap read as
a structural one.

## 3. What does survive: a trend in the fork sweep, flat in the pipe

`derived / measured` at `n = 12800`, `delta >= 0.4`:

| `delta` | fork | pipe |
|---|---|---|
| 0.4 | 0.944 ± 0.018 | 0.942 ± 0.043 |
| 0.6 | 0.992 ± 0.017 | 0.953 ± 0.029 |
| 0.8 | **1.112 ± 0.023** | 0.952 ± 0.022 |

The fork rises monotonically and ends **~5σ from 1**; the pipe is flat at ~0.95 across the same
arrow range. By the design above, **the coefficient channel is acquitted and `Var(logX1)` is
implicated** — in the fork, `delta = 0.4/0.6/0.8` carries `s²` = 1.16/1.36/1.64.

A second, separate fact: the pipe sits at a **uniform ~0.95**. The ratio is `derived/measured`, so
that is the chain **under**-predicting the pipe's wrapper effect by ~5% at every arrow strength —
a level offset, not a trend, and invisible at `n = 800` (where the same cell reads 0.992).

Note the fork crosses unity rather than sitting off it: it **under**-predicts at `delta = 0.4`
(0.944), lands on the measurement at 0.6 (0.992), and **over**-predicts at 0.8 (1.112). A trend
through 1 is a different object from a constant offset, and the two roles show one each.

## 4. This re-opens a question THEORY records as refuted

`THEORY.md` reports the `exp`-linearisation hypothesis **refuted**: swept over `s²` = 0.52 → 2.05
the ratio moved *toward* 1 (0.841 → 0.908) rather than away. But that sweep used the two-block
simulation which **misses the pipeline's level by 31%** — a limitation `THEORY.md` states
explicitly in the caution immediately below it, noting it "can refute a hypothesis about the
closed form's *shape* but cannot adjudicate the pipeline's *level*".

V28 uses the real pipeline arms and finds the agreement **degrading** as `s²` grows within the
fork. The earlier refutation should be read as not established rather than as settled.

**It is not yet a confirmation either.** In the fork sweep `s²` and the coefficient move
together, so V28 cannot attribute to `s²` alone — it can only say the pipe's flatness rules the
coefficient out. **V29** (`run_v29_sdx.sh`) moves `s²` via `sd_x` with the arrow fixed at 0.60,
in both roles, which is the lever that decides it.

## 5. Two methodological notes, both of which cost real time

**(a) Evaluating an asymptotic closed form at too small an `N` manufactures a role contrast.**
The first pass computed the derived column at `N = 20000` and read the pipe 15% above the
chain's established +5.42 while the fork agreed — which looked exactly like a role-specific
defect. It was the evaluation: the closed form's own finite-`N` bias is **+0.87 pp for a pipe
against +0.27 pp for a fork** at `N = 20000`. At `N = 300000` it returns +5.44 / +5.04,
reproducing `derive_scale_propagation.R` exactly. *Two implementations of the same chain
disagreeing is a fact about the implementations until proven otherwise.*

**(b) A run with no `.pid` file is invisible to `run_status.sh`.** The `n = 800` and `n = 3200`
waves were launched as a bare background command. `run_status.sh` reports only runs that write
`logs/<tag>_<stamp>.pid`, so a live 24-core job showed as **"Nothing running"**, and a
part-written checkpoint read as a crash. It was ~3 minutes from being killed and restarted on
top of itself. Both later waves use runners (`run_v28_arrow.sh`, `run_v29_sdx.sh`).

**(c) The `delta = 0.2` cells carry no information** and are excluded from every reading above.
Their measured effect is 0.26–0.73 pp, so `derived/measured` is a ratio of near-zero arms and
reads 1.19–2.23 with enormous error. This is the **sixth** gate/statistic in this programme
spoiled by a denominator that can approach zero — the standing lesson is that an informativeness
guard needs **a floor on the quantity being compared**, not only on the comparison.

## 6. Reproduce

```bash
bash run_v28_arrow.sh      # n = 12800 wave (the n=800/3200 waves predate the runner)
Rscript analyze_v28.R      # all three waves + the convergence table
```

Cells `fkd02–08` / `ppd02–08` in `R/robustness.R` (axis `"arrow"`); raw results in
`results/v28_n{800,3200}_raw.csv` and `results/v28b_n12800_raw.csv`.
