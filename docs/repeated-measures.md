# Repeated-outcome data and subject-wide imputation

Handling data with several rows per subject: when to impute subject-wide rather than row-wise, and how `timing` in the variable dictionary controls it.

> Only relevant if your data has repeated measures. Row-wise imputation of a time-invariant covariate can give one subject different values at different visits.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

The generic pipeline can also handle repeated-outcome data where the outcome is measured multiple times per subject and most covariates are subject-level or baseline variables.

Use this imputation strategy in `00_config.R`:

```r
imputation = list(
  ...
  strategy = "subject_wide_with_repeated_y_auxiliary",
  ...
)
```

This strategy does the following:

1. Starts from long data, one row per subject-time observation.
2. Creates a subject-wide imputation dataset, one row per subject.
3. Converts the repeated outcome into wide auxiliary variables, for example `ps_1, ps_2, ps_3, ps_4, ps_5, ps_6`.
4. Imputes subject-level covariates once per subject.
5. Joins the imputed subject-level variables back to the original long outcome data.
6. Fits the `brms` model on the long data.

Required `00_config.R` fields:

- `analysis_spec$data$id_var`
- `analysis_spec$data$time_var`
- `analysis_spec$outcome$y_var`
- `analysis_spec$outcome$y_prefix`
- `analysis_spec$outcome$y_wide_regex`

Example:

```r
data = list(
  raw_data_file = "data/my_repeated_data.rds",
  id_var = "ID",
  row_id_var = "row_id",
  data_structure = "repeated_y_subject_covariates",
  time_var = "time"
)

outcome = list(
  y_var = "ps",
  family = "bernoulli",
  link = "logit",
  y_prefix = "ps_",
  y_wide_regex = "^ps_",
  predict_missing_y = TRUE
)

imputation = list(
  enabled = TRUE,
  strategy = "subject_wide_with_repeated_y_auxiliary",
  m = 24,
  maxiter = 5,
  mean_match_k = 5,
  verbose = FALSE,
  impute_y = FALSE,
  extra_exclude_targets = character(0)
)
```

In `00_variable_dictionary.csv`, subject-level covariates should have:

```r
timing = single
```

or:

```r
timing = baseline
```

Repeated or time-varying variables should have:

```r
timing = repeated
```

or:

```r
timing = time_varying
```

If `impute_y = FALSE`, repeated outcome values are not imputed as targets. They are used as wide auxiliary predictors when possible, and missing outcome rows can later be summarised using posterior prediction.

---

*[← Back to the main README](../README.md)*
