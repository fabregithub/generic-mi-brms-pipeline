# Manuscript demonstration — congenial imputation of censored PFAS (JECS-grounded)

A synthetic-but-grounded demonstration that the MI + brms pipeline (with the
`leftcens` censored-exposure block-FCS) imputes left-/interval-censored chemical
exposures **congenially** and beats half-limit substitution, for both
single-pollutant and mixture (BKMR) health-effect models.

**All data are synthetic.** Their structure is taken from three JECS papers; every
DGP parameter is traced in [`PROVENANCE.md`](PROVENANCE.md). Not participant data.

## Layout

```
manuscript_demo/
├── PROVENANCE.md            paper -> DGP parameter mapping (Lai / Iwata / Atagi)
├── MANUSCRIPT_SKELETON.md   section-by-section manuscript skeleton
├── R/
│   ├── make_demo_data.R     the DGP: 7 interval-censored PFAS + Z + two outcomes
│   └── analysis_methods.R   the 4 methods + single-pollutant estimand + scoring
├── analyze_imputation.R     one scenario, single-pollutant (the Y-inclusion test)
├── run_scenarios.R          the censoring x correlation x outcome grid
├── analyze_bkmr.R           the mixture (BKMR) estimand — pluggable, bounded
└── results/                 generated CSVs (git-ignored)
```

## The design

- **Exposures X:** 7 prenatal PFAS, three-tier interval censoring (ND / DNQ /
  quantified), correlations and MDL/LCMRL from Lai 2025.
- **Outcomes Y:** Kawasaki disease (rare; Iwata 2024) for the mixture/BKMR +
  substitution story; childhood wheeze/asthma (common; Atagi 2024) for the
  outcome-inclusion (congeniality) story.
- **Methods vs oracle:** half-limit `substitution`; `impute_noY` (outcome-blind);
  `impute_congenial` (outcome-aware, the pipeline).
- **Robustness factors:** `limit_scale` (censoring rate) × `correlation` (low/real/
  high) × `outcome` — supply real MDL/LCMRL to `make_demo_data(mdl=, lcmrl=)`.

## Workflow (start small, then scale)

```bash
# 1. verify the code fast (small N, few reps)
N=2000 REPS=5 OUTCOME=asthma LIMIT_SCALE=1.5 Rscript analyze_imputation.R

# 2. the full single-pollutant grid (cheap estimand)
N=25000 REPS=200 Rscript run_scenarios.R          # scale N/REPS for the study

# 3. the mixture estimand (BKMR is slow — bounded/subsampled; needs the bkmr package)
OUTCOME=kd N=25000 N_SUB=2000 M=10 ITER=2000 REPS=1 Rscript analyze_bkmr.R
```

Environment variables select N, m, reps, outcome, and the scenario factors (see each
script's header). Results are written to `results/` (git-ignored).

## Notes

- **Start with small N to verify, then scale to ~25,000** (KD is rare — ~285 cases at
  n=25k). Verified at small N; the scientific signal (esp. the Y-inclusion contrast)
  needs the full n.
- **BKMR compute.** `bkmr::kmbayes` is MCMC and slow at large n; the mixture script
  subsamples (Iwata-style) and the engine is pluggable — drop in an accelerated
  variant (e.g. A-BKMR) via `fit_mixture_<engine>()` and `MIXTURE_ENGINE=`.
- **Key expected finding.** Imputation ≫ substitution everywhere; the outcome-aware
  (congeniality) advantage is *conditional* — largest under low analyte correlation
  and heavy censoring, negligible when analytes are strongly correlated (they rescue
  each other, as in the Gaussian-copula rationale).
