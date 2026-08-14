# Demo DGP provenance — every parameter traced to a JECS publication

The demonstration cohort is **synthetic** but its structure is taken from three
JECS papers so the demo reproduces a realistic exposure–outcome problem. Nothing
here is real participant data; the "truth" (coefficients) is known by construction
so the pipeline's recovery can be checked.

**Sources**
- **Lai et al. 2025**, *Ecotox. Environ. Saf.* 294:118107 — *X baseline*: PFAS
  distributions, detection limits, correlations, and Z→X determinants.
- **Iwata et al. 2024**, *Environ. Int.* 183:108321 — *Y ~ X + Z*: Kawasaki disease
  outcome, the 7 eligible PFAS, single-logistic / WQS / BKMR effect sizes,
  covariate adjustment set.
- *(Atagi et al. 2024, Environ. Res. 240:117499 — alternative outcome, not used.)*

## Exposures X — 7 PFAS (log scale, ng/mL)

| PFAS | Group | Median (ng/mL) | Detection | Source |
|---|---|---|---|---|
| PFOA   | PFCA | 1.60 | 100%  | Iwata (median); Lai (100% detected) |
| PFOS   | PFSA | 2.90 | ~99%  | Iwata (median) |
| PFNA   | PFCA | 1.40 | >90%  | Iwata (median) |
| PFDA   | PFCA | 0.49 | >90%  | Iwata (median) |
| PFUnDA | PFCA | 1.10 | >90%  | Iwata (median) |
| PFHxS  | PFSA | 0.33 | >90%  | Iwata (median) |
| PFTrDA | PFCA | ~0.15 | 78.5% | Iwata (detection 78.5%) |

- **Left-censoring:** per-analyte MDL in **0.016–0.16 ng/mL** (Lai). MDLs are
  calibrated so realised detection matches the table (PFOA ~100%, six others
  ~90–99%, **PFTrDA ~78.5%** — the heavily-censored analyte that most stresses the
  imputation). Represented as `leftcens` interval columns `<pfas>_lo` / `<pfas>_hi`.
- **Marginal spread:** `sdlog ≈ 0.8` (geometric SD ≈ 2.2), typical for serum PFAS.
- **Correlation (log scale)** — Lai's Schober-classified R values:
  long-chain PFCAs PFNA–PFOA 0.74, PFNA–PFDA 0.72, PFUnDA–PFDA 0.81,
  PFUnDA–PFTrDA 0.83; PFOS–PFDA 0.79, PFOS–PFUnDA 0.70; **PFHxS weakly correlated**
  (~0.2–0.3). Encoded as a block correlation matrix.

## Covariates Z — distributions (Lai Table 1) and roles

| Z | Distribution | Role | Source |
|---|---|---|---|
| region (15 centres) | Hokkaido…Okinawa, probs per Lai | strongest X determinant; Y confounder | Lai (ALE: region top determinant) |
| maternal age (yr)   | median 31, SD 5 | X ↑ with age; Y confounder | Lai (median/SD); age↑PFAS |
| parity (0/1/2/3+)   | 0.39 / 0.39 / 0.16 / 0.06 | X ↓ with parity; age×parity interaction | Lai (parity↓PFAS, age×parity) |
| BMI (kg/m²)         | median 21, SD 3.3 | Y confounder (sensitivity set) | Lai; Iwata sensitivity |
| smoking (never/past/current) | 0.58 / 0.36 / 0.045 | Y confounder (↑ in KD group) | Lai; Iwata (KD↑smoking) |
| family income (7 cat) | per Lai | Y confounder (↑ in KD group) | Iwata (primary adjust; KD↑income) |
| folic acid (cont.)  | ~normal | Y confounder (primary adjust) | Iwata (DAG-supported) |
| eGFR (mL/min/1.73m²)| median 131, SD 24 | Y confounder (sensitivity) | Lai; Iwata sensitivity |
| child sex (M/F)     | ~0.51 M | Y confounder (M↑ in KD) | Iwata (KD↑male) |
| low birth weight (<2500g) | ~0.09 | Y confounder (↑ in KD) | Iwata (KD↑LBW) |

**Z→X determinant model (Lai):** `log(PFAS_j) = μ_j + region_effect_j − b_parity·parity
+ b_age·(age−31)/5 + b_ap·parity·(age−31)/5 + ε` (region strongest; parity negative;
age positive; age×parity interaction), then correlated residuals `ε`.

## Outcomes — two, by design

**Two outcomes are generated from the same PFAS mixture**, because the two demo
goals need different outcome properties:

| Outcome | Paper | Prevalence | Role |
|---|---|---|---|
| **`kd`** — Kawasaki disease | Iwata 2024 | ~1.1% (rare) | mixture/BKMR + imputation-vs-substitution; matches Iwata's 7 PFAS + BKMR |
| **`asthma`** — childhood wheeze/asthma | Atagi 2024 | ~22% (common) | the **Y-inclusion (congeniality)** demonstration; a rare outcome carries too little information about X |

Atagi (asthma) prevalences: "wheeze ever" 29%, "asthma ever" 11.9% (we target ~22%).
Atagi's real PFAS effects were near-null (aOR ~0.94/doubling); the demo uses
**moderate inverse effects** so the method's behaviour is demonstrable — a methods
illustration, not an asthma claim. Atagi is also a **direct methodological
comparator**: it already used `mice` (25 imputations) + a Gibbs sampler for the
below-LCMRL PFAS, but not an outcome-aware censored draw.

**Emerging finding (verified in code at small N):** imputation (either arm) clearly
beats half-MDL substitution in every scenario. The **Y-inclusion** advantage is
*conditional*: negligible when the PFAS are strongly correlated (real data — the
other analytes already pin down a censored one, as with Lai's copula), and it
appears when correlation is **low and/or censoring heavy**. This is exactly what the
`correlation` × `limit_scale` scenario grid is designed to expose, and it reconciles
with the §7.7 / copula rationale.

## Outcome Y — Kawasaki disease (binary, rare)

- **Prevalence ~1.1%** (Iwata: 271 / 25,256). Intercept set to this baseline.
- **Model:** `logit P(KD) = α + Σ_j β_j · log2(PFAS_j, standardised) + γ·Z_confounders (+ mixture nonlinearity)`.
- **Single-PFAS effect sizes (per log2 unit)** matched to Iwata Table 3
  (all inverse, OR < 1): PFOS OR 0.839 (β≈−0.176), PFOA OR 0.894 (β≈−0.112),
  others OR 0.83–0.91 (β≈−0.09…−0.19).
- **Mixture (WQS/BKMR):** the PFAS mixture is **inversely** associated with KD
  (Iwata: WQS negative-effect OR ≈ 0.86; BKMR inverse). Dominant weights on
  PFOA, PFTrDA, PFHxS, PFOS (Iwata WQS weights > 1/7). A mild **interaction /
  curvature** on the mixture surface is added so the BKMR-on-imputed analysis has a
  non-linear signal to recover (the Phase-1 §7.7 stress test).
- **Confounding directions (Iwata):** KD ↑ with maternal smoking, male child,
  low birth weight, higher family income.

## What the demo shows

1. The pipeline **imputes the censored PFAS congenially** (Y-aware, skew-robust)
   rather than the papers' **half-MRL substitution** (Iwata's stated method).
2. **Single-PFAS** brms logistic models recover the per-analyte ORs.
3. **BKMR-on-imputed** (bkmr::kmbayes per imputation, pooled) recovers the inverse
   mixture surface — and lets us compare congenial imputation vs half-MRL
   substitution feeding BKMR (the §7.7 question, on a realistic problem).

Truth (all β, γ, MDLs, correlations) is written to a companion `*_truth.rds` so
every recovered estimate can be scored against the value that generated it.
