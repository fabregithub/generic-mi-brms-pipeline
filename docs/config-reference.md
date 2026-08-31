# `00_config.R` reference

Option-by-option reference for `00_config.R` — the main file you edit for a new analysis.

> Reference material: look up the option you need. The commented `00_config.R` in the repository root is the authoritative copy.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

## Important `00_config.R` options

The example configs are the easiest reference, but the following settings are worth checking whenever you adapt the pipeline.

### Outcome family and link

Set the outcome family and link in `analysis_spec$outcome`:

```r
outcome = list(
  y_var = "low",
  family = "bernoulli",
  link = "logit",
  y_prefix = NULL,
  y_wide_regex = NULL,
  predict_missing_y = TRUE
)
```

Common link examples:

| Family | Common links |
|---|---|
| `gaussian` | `identity`, `log` |
| `bernoulli` | `logit`, `probit`, `cloglog` |
| `poisson` | `log`, `identity` |
| `negbinomial` | `log`, `identity` |
| `beta` | `logit`, `probit`, `cloglog`, `log` |
| `ordinal` | often `logit` or `probit`, depending on the ordinal family |
| `categorical` | often `logit` |
| `cox` | `log` (log hazard ratio scale) |

Choose a link that is valid for the selected family and appropriate for the scientific question.

### `y_prefix` and `y_wide_regex`

These are only needed for repeated-outcome wide imputation. For ordinary single-time analyses, keep both as `NULL`.

For example, if the long outcome is `ps` and `time` is `1, 2, ..., 6`, the subject-wide auxiliary outcome columns may be `ps_1`, `ps_2`, `ps_3`, `ps_4`, `ps_5`, `ps_6`.

Then use:

```r
y_prefix = "ps_"
y_wide_regex = "^ps_"
```

### Imputation targets and missingness assumptions

The imputation step is intended for ordinary missing data where a MICE-style imputation model is defensible. In the variable dictionary, set `impute_target = TRUE` only for variables to be imputed under that assumption.

If a variable is left-censored, below a detection limit, structurally missing, not applicable by design or likely not missing at random (MNAR), handle that issue before running the pipeline. Such values should not be treated as ordinary `NA` values unless that choice is explicitly justified.

### Skipping imputation when there is no missing data

If the analysis dataset has no missing covariate data, or if you deliberately want to run a complete-data analysis without multiple imputation, edit the `imputation = list(...)` block in `00_config.R` directly:

```r
imputation = list(
  enabled = FALSE,
  strategy = "none",

  # These are ignored when enabled = FALSE, but keeping m = 1
  # makes the intended one-dataset workflow explicit.
  m = 1,

  maxiter = 0,
  mean_match_k = 5,
  verbose = FALSE,
  impute_y = FALSE,
  extra_exclude_targets = character(0)
)
```


With this setting, the pipeline prepares one model dataset and fits one `brms` model, rather than creating multiple imputed datasets.

### Imputation seeding and extending an existing run

Two further `imputation = list(...)` settings control reproducibility and safe extension:

```r
imputation = list(
  ...
  # Base seed for miceRanger. Defaults to analysis_spec$model$seed if not set.
  # miceRanger has no native per-chain seed/resume support, so each batch of
  # new imputations (the initial m, or a later extension batch) is seeded
  # deterministically as seed + (number of imputations that already existed
  # before that batch). This reproduces the same imputed datasets across
  # repeated runs, though it is not bit-identical to a hypothetical
  # single-shot miceRanger(m = target_m) call.
  seed = 12345,

  # FALSE by default. Step 3 normally refuses to fit more imputations than
  # already exist on disk, to avoid silently treating a config edit as a
  # request to overwrite existing imputed data. Set this to TRUE to let
  # Step 3 generate only the additional imputations needed to reach the
  # current m, appending them to the existing manifest without touching
  # any existing imputed-data file. run_all.R's automatic m-increment loop
  # (analysis_spec$mi_stability$auto_increment) sets this internally.
  allow_extend = FALSE,
  ...
)
```

For a no-missing-data analysis, the cleanest dictionary setup is usually `impute_target = FALSE` for all variables.

Before disabling imputation, confirm that there are no missing predictor values that the model needs:

```r
source("00_config.R")
d <- readRDS(paths$raw_data)
colSums(is.na(d))
```

If missing predictor values remain and imputation is disabled, those rows may be dropped or model-data checks may fail, depending on the variables and model formula.

### Variable groups

The recommended default is:

```r
variables = NULL
```

This means variable roles, types, timing, scaling, reference categories, imputation targets, model inclusion and auxiliary-variable status are read from `00_variable_dictionary.csv`.

Advanced users can override selected derived groups in `00_config.R`:

```r
variables = list(
  exposure_vars = c("exposure"),
  covariate_vars = c("age", "sex", "income"),
  auxiliary_vars = c("baseline_score"),

  continuous_vars = c("age", "baseline_score"),
  categorical_vars = c("sex"),
  ordinal_vars = c("income"),

  subject_level_vars = c("age", "sex", "income", "baseline_score"),
  time_varying_vars = c("current_exposure"),

  scale_vars = c("age", "baseline_score")
)
```

For most analyses, it is simpler and safer to keep `variables = NULL` and control the analysis through the dictionary.

### Initial values

The default example setting is:

```r
init = "random"
```

This uses ordinary Stan-style random initialisation.

A useful alternative for some difficult logistic models is:

```r
init = 0
```

This starts parameters at zero and can sometimes avoid early warm-up overflow warnings. It is not required for ordinary analyses.

Advanced users may also supply a custom initialisation function:

```r
init = function() {
  list()
}
```

Use `init = "random"` unless there is a clear reason to use `init = 0`.

### Posterior draw extraction regex

Step 6 extracts posterior draws using `analysis_spec$model$parameter_draw_regex`.

For ordinary fixed, random-effect SD and residual parameters only:

```r
parameter_draw_regex = "^(b_|sd_|sigma)"
```

The recommended general default is:

```r
parameter_draw_regex = "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)"
```

This includes:

| Prefix | Meaning |
|---|---|
| `b_` | ordinary fixed effects |
| `bsp_` | monotonic-effect coefficients used by `brms::mo()` |
| `sd_` | group-level standard deviations |
| `sigma` | residual SD for Gaussian models |
| `sds_` | smooth-term standard deviations |
| `bs_` | smooth basis coefficients |
| `simo_` | simplex parameters used by `brms::mo()` |

Use the general default if the model contains `s()` or `mo()` terms.

### How Step 6 pools draws across imputations

Step 6 does not simply concatenate every imputation's posterior draws and summarise the pooled sample directly. Two corrections are applied per parameter:

1. Each draw is weighted by `1 / (m * K_i)`, where `K_i` is the number of finite draws from imputation `i`. This makes every imputation contribute exactly `1/m` to every summary statistic, regardless of how many draws a particular imputation happened to produce after filtering non-finite values (e.g. for monotonic simplex parameters).
2. A finite-`m` variance correction (the `B/m` term from Rubin's combining rule) is applied on top of the weighted pooled sample, but only when a per-parameter shape diagnostic (a bimodality coefficient) indicates the pooled posterior is unimodal enough for the correction's symmetry assumption to be safe, and only on a support-respecting transform (log for strictly positive parameters, logit for (0, 1)-bounded parameters). For visibly multimodal parameters, the correction is skipped and the uncorrected weighted pooled sample is reported instead, since a symmetric rescale would distort genuine between-imputation structure rather than represent it.

`results/parameter_summary.csv`/`.rds` record this per parameter via the `m_imputations`, `between_var`, `within_var`, `variance_corrected`, `transform_used` and `bimodality_coef` columns, so you can audit whether and how each parameter was corrected.

### Posterior summaries and ROPE

The `summary` block controls posterior summaries:

```r
summary = list(
  effects = "fixed",
  component = "conditional",
  centrality = "median",
  ci = 0.95,
  ci_method = "HDI",
  test = c("p_direction", "rope"),

  rope = list(
    method = "fixed",
    fixed_range = c(-1, 1),
    width_probability = 0.05
  )
)
```

`p_direction` is the posterior probability of the dominant sign. ROPE means region of practical equivalence and should ideally be defined using scientific or clinical relevance.

Common ROPE options:

```r
# No ROPE summary
rope = list(
  method = "none"
)
```

```r
# Fixed ROPE on the model coefficient scale
rope = list(
  method = "fixed",
  fixed_range = c(-1, 1)
)
```

```r
# Logistic/logit example: odds ratio between 0.95 and 1.05
rope = list(
  method = "fixed",
  fixed_range = log(c(0.95, 1.05))
)
```

```r
# Convenience option for workflows that implement family-specific defaults
rope = list(
  method = "auto",
  width_probability = 0.05
)
```

For Gaussian models, the ROPE is on the outcome/model scale. For Bernoulli/logit models, the ROPE is on the log-odds scale; using `log(c(0.95, 1.05))` corresponds to an odds-ratio range from 0.95 to 1.05.

If no meaningful ROPE has been defined, use:

```r
test = c("p_direction")
rope = list(method = "none")
```

---

---

*[← Back to the main README](../README.md)*
