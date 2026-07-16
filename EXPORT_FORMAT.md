# Federated Meta-Analysis Export Format

This document describes the files produced by Step 12 (`12_export_draws.R`)
for use with the companion meta-analysis pipeline
[generic-mi-brms-meta](https://github.com/fabregithub/generic-mi-brms-meta).

---

## Variable name harmonisation

The most important prerequisite for federated meta-analysis is that all
cohorts use **identical variable names** for the same constructs. Because
parameter names in the export file are derived directly from the variable
names in `00_variable_dictionary.csv` (e.g. a variable named `smoke_z`
produces a parameter named `b_smoke_z`), harmonising the variable
dictionary is sufficient to harmonise the export — no post-hoc renaming
is needed.

**Recommended practice:** the coordinating site distributes a common
`00_variable_dictionary.csv` template with agreed `var` column values
before each cohort runs the pipeline. Cohorts may add cohort-specific
covariates freely, but the `var` names for shared exposures and
confounders should not be changed. When this is followed, the
`parameter_map` in `generic-mi-brms-meta/00_config.R` can be left empty.

---

## Files

Step 12 writes three files to `results/export/`:

### `cohort_draws.rds` (primary transfer file)

An R data frame in long format, one row per posterior draw:

| Column | Type | Description |
|---|---|---|
| `cohort_id` | character | Short unique identifier for this cohort, set in `analysis_spec$export$cohort_id` in `00_config.R` |
| `parameter` | character | brms parameter name (e.g. `b_exposure_z`, `b_exposure_z:time`) |
| `imputation` | integer | Imputation index (1 to m; m may differ across cohorts) |
| `draw_index` | integer | Within-imputation draw index from the MCMC posterior |
| `value` | double | Posterior draw value on the linear predictor scale (log-OR, log-HR, or unstandardised coefficient depending on model family) |

**Scale note:** values are always on the *link* scale (log for binomial/Poisson/Cox, identity for Gaussian). The meta-analysis model operates on this scale. Back-transformation (OR, HR, RR) is applied at the reporting stage.

**Scope:** by default only parameters for variables with `role == "exposure"` in `00_variable_dictionary.csv` are exported (`scope = "exposure_only"`). Interaction terms involving an exposure variable are also included. Set `scope = "all"` to export all parameters matched by `parameter_draw_regex`.

### `cohort_draws.csv`

Same content as `cohort_draws.rds` in CSV format for interoperability with non-R workflows. For large `m` or many draws, the `.rds` file is preferred for transfer (smaller, lossless).

### `cohort_metadata.json`

A lightweight sidecar file for validation at the coordinating site. Does not contain any draws.

```json
{
  "cohort_id":  "cohort_japan_2024",
  "parameters": ["b_exposure_z"],
  "m":          10,
  "n_draws":    40000,
  "family":     "bernoulli",
  "created_at": "2026-07-16T17:00:00"
}
```

| Field | Description |
|---|---|
| `cohort_id` | Must match the `cohort_id` column in the draws file |
| `parameters` | Sorted list of exported parameter names |
| `m` | Number of imputations used |
| `n_draws` | Total rows in the draws file |
| `family` | brms model family (for scale reference) |
| `created_at` | ISO 8601 timestamp of export |

---

## What to transfer

Each cohort sends **two files** to the coordinating site:

```
cohort_draws.rds
cohort_metadata.json
```

Individual-level data are never included in the export. Only posterior
draws for the pre-specified exposure parameter(s) are transferred,
which substantially reduces the data-sharing burden compared to
individual participant data (IPD) meta-analysis.

---

## Checklist for each cohort before transfer

- [ ] `cohort_id` is set in `00_config.R` and is unique across cohorts
- [ ] Variable names for shared exposures match the coordinating site's template `00_variable_dictionary.csv`
- [ ] `scope = "exposure_only"` (default) — confirm only exposure parameters are exported, not covariates
- [ ] `cohort_metadata.json` lists the expected parameters and family
- [ ] Model diagnostics passed (Step 5): no excessive divergences or treedepth hits
- [ ] Imputation stability confirmed (Step 11): pooled estimates stable at the chosen `m`
