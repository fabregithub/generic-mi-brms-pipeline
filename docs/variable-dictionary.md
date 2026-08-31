# Variable dictionary

Field-by-field reference for `00_variable_dictionary.csv` — the file that tells the pipeline what each column of your data *is* and how it should be treated.

> This is reference material: look up the field you need rather than reading it through.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

The file `00_variable_dictionary.csv` is the main machine-readable description of the analysis variables.

Expected columns are:

- `var`
- `label`
- `role`
- `type`
- `timing`
- `scale`
- `reference`
- `impute_target`
- `use_in_model`
- `use_as_auxiliary`

Some fields are straightforward:

- `var`: exact variable name in the dataset
- `label`: human-readable label for tables and reports
- `impute_target`: whether this variable should be imputed when missing
- `use_in_model`: whether this variable should appear in the final brms model

Set `impute_target = TRUE` only for variables whose missing values should be treated as ordinary missing data in the MICE-style imputation model. Do not set `impute_target = TRUE` for left-censored values, below-detection-limit values, structural missingness, skip-pattern missingness or known MNAR variables unless those issues have already been handled appropriately before the pipeline.

The other fields need more explanation.

## `role`

`role` describes the analytical purpose of the variable.

Recommended values:

| Value | Meaning | Example |
|---|---|---|
| `outcome` | Continuous or general outcome | `Ozone` |
| `binary_outcome` | Binary outcome for Bernoulli/logistic models | `low` |
| `exposure` | Main exposure or predictor of scientific interest | treatment group, air pollution |
| `covariate` | Adjustment variable / confounder / predictor | age, sex, income |
| `auxiliary` | Used for imputation only, not included in final model | extra baseline score |
| `id` | Subject, cluster or row identifier | `ID`, `row_id` |
| `time` | Measurement occasion, wave, visit or follow-up time | `time`, `wave` |
| `cluster` | Grouping variable for random effects or clustering | school ID, hospital ID |
| `strata` | Stratification variable, if relevant | site, cohort |

For most standard regression analyses, the most common roles are:

- `outcome`
- `binary_outcome`
- `exposure`
- `covariate`
- `id`
- `time`
- `auxiliary`

Example:

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

## `type`

`type` describes the statistical data type.

Recommended values:

| Value | Meaning | Typical R class |
|---|---|---|
| `continuous` | Numeric variable on a continuous scale | numeric |
| `integer` | Count-like integer variable | integer/numeric |
| `binary` | Two-level variable, often 0/1 | factor or integer |
| `categorical` | Unordered categorical variable with 3+ levels | factor |
| `ordinal` | Ordered categorical variable | ordered factor |
| `date` | Calendar date | Date |
| `id` | Identifier, not modelled as numeric | character/factor/integer |

Examples:

```text
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
```

Important notes:

- Binary and categorical variables should usually be factors before modelling.
- Continuous variables with `scale = z` get a new `_z` version for modelling.
- ID variables should not be scaled or imputed unless there is a specific reason.

## `timing`

`timing` describes when the variable was measured.

Recommended values:

| Value | Meaning | Example |
|---|---|---|
| `single` | Measured once per analytic unit | sex, baseline income |
| `baseline` | Measured at baseline only | baseline age |
| `repeated` | Measured repeatedly across time/visits | repeated outcome |
| `time_varying` | Predictor changes over time | current exposure, current medication |
| `derived` | Created from other variables | standardised score |
| `id` | Identifier variable | subject ID |

Examples:

```text
ps,Psychological stress,binary_outcome,binary,repeated,no,0,FALSE,TRUE,TRUE
time,Measurement wave,time,integer,repeated,no,,FALSE,TRUE,FALSE
c_sex,Child sex,covariate,binary,baseline,no,0,TRUE,TRUE,FALSE
```

How this matters:

- For `single_time` datasets, most variables usually have `timing = single`.
- For repeated-outcome data with subject-level covariates, the outcome and time variable are usually `repeated`, while most covariates are `baseline` or `single`.
- For repeated data with time-varying predictors, mark those predictors as `time_varying`.

## `scale`

`scale` tells the pipeline whether and how to transform a variable before modelling.

Recommended values:

| Value | Meaning | Result |
|---|---|---|
| `no` | No scaling or transformation | original variable used |
| `z` | Standardise to mean 0 and SD 1 | creates `var_z` |
| `centre` | Mean-centre only | creates centred version if supported |
| `log` | Log-transform | creates log version if supported |
| `custom` | User-defined transform outside dictionary | handled in config/functions |

Currently, the most important supported value is `z`.

Example:

```text
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
```

This creates `age_z`, and the model should use `age_z`, not `age`.

Use `z` for continuous predictors when coefficients should represent the effect per 1 SD increase. This is usually helpful for Bayesian models because weakly informative priors like `normal(0, 1)` or `normal(0, 1.5)` are easier to interpret on standardised predictors.

Do not use `z` for:

- outcomes
- ID variables
- binary/categorical variables
- already standardised variables

## `reference`

`reference` defines the reference category for binary, categorical or ordinal variables.

Examples:

```text
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
```

For a binary variable coded `0/1`, setting:

```r
reference = 0
```

means the coefficient compares level `1` against level `0`.

For categorical variables, choose the scientifically meaningful or most common category as the reference.

Important notes:

- The value in `reference` must match an actual value/level in the data.
- Leave `reference` blank for continuous variables.
- Reference levels affect coefficient interpretation but not overall model fit.

## `use_as_auxiliary`

`use_as_auxiliary` controls whether a variable is used in imputation but excluded from the final analysis model.

This is useful for variables that help predict missingness or missing values, but are not part of the scientific model.

Examples:

```text
extra_score,Extra baseline score,auxiliary,continuous,baseline,z,,TRUE,FALSE,TRUE
hospital_id,Hospital ID,auxiliary,categorical,single,no,,FALSE,FALSE,TRUE
```

Interpretation:

| impute_target | use_in_model | use_as_auxiliary | Meaning |
|---|---|---|---|
| `TRUE` | `TRUE` | `FALSE` | Impute if missing and include in model |
| `TRUE` | `FALSE` | `TRUE` | Impute/use for imputation but exclude from model |
| `FALSE` | `FALSE` | `TRUE` | Use as imputation predictor only |
| `FALSE` | `TRUE` | `FALSE` | Include in model but do not impute |
| `FALSE` | `FALSE` | `FALSE` | Keep as metadata/ID or ignore analytically |

Common use cases:

1. Auxiliary predictor only:

```text
baseline_score,Baseline questionnaire score,auxiliary,continuous,baseline,z,,FALSE,FALSE,TRUE
```

2. Variable to impute and include:

```text
income,Household income,covariate,categorical,baseline,no,1,TRUE,TRUE,FALSE
```

3. Outcome not imputed:

```text
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
```

4. Repeated outcome as imputation auxiliary:

```text
ps,Psychological stress,binary_outcome,binary,repeated,no,0,FALSE,TRUE,TRUE
```

Whether repeated outcomes are used as auxiliary predictors depends on the imputation strategy in `00_config.R`.

## Example: Gaussian demo

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
Ozone,Ozone concentration,outcome,continuous,single,no,,FALSE,TRUE,FALSE
Solar.R,Solar radiation,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Wind,Wind speed,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Temp,Temperature,covariate,continuous,single,z,,TRUE,TRUE,FALSE
Month,Month,covariate,categorical,single,no,5,FALSE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

## Example: Logistic demo

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
lwt,Maternal weight at last menstrual period,covariate,continuous,single,z,,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ptl,Previous premature labours,covariate,continuous,single,z,,TRUE,TRUE,FALSE
ht,History of hypertension,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ui,Uterine irritability,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ftv,Physician visits during first trimester,covariate,continuous,single,z,,TRUE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

## Example: Spline + monotonic demo

The `birthwt_spline_monotonic` example is used to test custom `brms` formula terms: `s(age_z, k = 5)` and `mo(lwt_q)`.

The corresponding dictionary demonstrates two important ideas:

- `scale = z` creates `age_z`, `ptl_z` and `ftv_z`
- `type = ordinal` creates an ordered factor suitable for `mo()`

Example:

```text
var,label,role,type,timing,scale,reference,impute_target,use_in_model,use_as_auxiliary
low,Low birth weight,binary_outcome,binary,single,no,0,FALSE,TRUE,FALSE
age,Maternal age,covariate,continuous,single,z,,TRUE,TRUE,FALSE
lwt_q,Maternal weight quintile,covariate,ordinal,single,no,1,TRUE,TRUE,FALSE
race,Race,covariate,categorical,single,no,1,TRUE,TRUE,FALSE
smoke,Smoking during pregnancy,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ptl,Previous premature labours,covariate,continuous,single,z,,TRUE,TRUE,FALSE
ht,History of hypertension,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ui,Uterine irritability,covariate,binary,single,no,0,TRUE,TRUE,FALSE
ftv,Physician visits during first trimester,covariate,continuous,single,z,,TRUE,TRUE,FALSE
row_id,Row ID,id,integer,single,no,,FALSE,FALSE,FALSE
```

In `00_config.R`, this dictionary is paired with a custom formula:

```r
custom_formula = brms::bf(
  low ~ s(age_z, k = 5) + mo(lwt_q) + race + smoke + ptl_z + ht + ui + ftv_z
)
```

Notes:

- `lwt_q` is created by the example-data script before the pipeline runs.
- `lwt_q` is marked as `type = ordinal`, so the pipeline converts it to an ordered factor.
- `mo(lwt_q)` then models the ordinal effect as monotonic but not necessarily equally spaced.
- The posterior draw regex for `mo()` models should include both `bsp_` and `simo_` parameters:

```r
model = list(
  ...
  parameter_draw_regex = "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)",
  ...
)
```

---


---

## Variable roles: dictionary by default, config as optional override

The pipeline uses `00_variable_dictionary.csv` as the default source of truth for variable metadata.

This means the following information should normally be specified only in the dictionary:

- `role`
- `type`
- `timing`
- `scale`
- `reference`
- `impute_target`
- `use_in_model`
- `use_as_auxiliary`

Therefore, the standard setting in `00_config.R` is:

```r
variables = NULL
```

The pipeline automatically derives internal variable groups from the dictionary, including:

- `exposure_vars`
- `covariate_vars`
- `auxiliary_vars`
- `continuous_vars`
- `categorical_vars`
- `ordinal_vars`
- `subject_level_vars`
- `time_varying_vars`
- `scale_vars`

For example, if `00_variable_dictionary.csv` contains:

```text
var,role,type,timing,scale,use_in_model
Solar.R,covariate,continuous,single,z,TRUE
Wind,covariate,continuous,single,z,TRUE
Temp,covariate,continuous,single,z,TRUE
Month,covariate,categorical,single,no,TRUE
```

and `analysis_spec$model$fixed_effects` is:

```r
fixed_effects = "auto"
```

the model uses `Solar.R_z + Wind_z + Temp_z + Month`.

## Optional overrides

Advanced users can override selected derived groups in `00_config.R`.

Example:

```r
variables = list(
  scale_vars = c("age", "income"),
  auxiliary_vars = c("baseline_score")
)
```

If `variables = NULL`, no override is applied.

If `exposure_vars` or `covariate_vars` are explicitly supplied in `analysis_spec$variables`, they are used by `fixed_effects = "auto"` instead of the dictionary-derived `use_in_model` list. For most analyses, it is simpler and safer to leave `variables = NULL` and control the model through the dictionary.

---

---

*[← Back to the main README](../README.md)*
