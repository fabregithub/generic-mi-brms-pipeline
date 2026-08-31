# Reporting: publication outputs and manuscript writing

What the pipeline produces for a paper, and how to turn those files into Methods and Results text.

> The outputs are meant to be quoted directly rather than re-derived by hand. The templates are starting points, not finished prose.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

## Publication outputs and inference guidance

After successful completion, publication outputs are written to `results/publication/`.

Typical outputs include:

- `results/publication/tables/main_effect_table_display.csv`
- `results/publication/tables/main_effect_table_full.csv`
- `results/publication/tables/diagnostics_summary.csv`
- `results/publication/tables/analysis_metadata.csv`
- `results/publication/tables/analysis_metadata.rds`
- `results/publication/tables/parameter_template_sentences.csv`
- `results/publication/figures/forest_plot_odds_ratios.png`
- `results/publication/report/bayesian_mi_report_template.qmd`

The generated Quarto report includes posterior results, diagnostics, figures, an imputation-count stability chapter (see below), and a methods/settings table based on `analysis_metadata.csv`. This table records key analysis settings such as the imputation strategy, target number of imputations, fitted imputations used in posterior summaries, model formula, family/link, priors, MCMC settings, parallel settings, posterior-summary settings and predictive-draw settings.

`08_publication_results.R` renders this report itself, to both HTML and DOCX, as the last thing it does. Running `run_all.R` therefore always produces a finished, rendered report with no separate manual step. If you need to re-render it by hand (for example after manually editing the `.qmd`), run:

```bash
# Bash command block
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

Self-rendering can be turned off with `analysis_spec$publication$render_quarto <- FALSE` in `00_config.R`, if you only want the `.qmd` written and prefer to render it yourself.

---

### Publication-ready template sentences

`08_publication_results.R` auto-generates one result sentence per fixed-effect parameter (and mo() main coefficients), following the [bayestestR reporting guidelines](https://easystats.github.io/bayestestR/articles/guidelines.html). Sentences are written to `results/publication/tables/parameter_template_sentences.csv` and embedded as a chapter in the main report. They follow the neutral format:

> *The effect of X has a probability of pd% of being negative (Median = median, 95% CrI [low, high]; OR = or, 95% CrI [or_low, or_high]) and can be considered as significant (0% in ROPE).*

The OR/RR clause appears only for families with a log or logit link. The ROPE clause appears only when a ROPE was specified in `analysis_spec$summary`. Copy-paste the sentences as a starting point and adjust wording for your manuscript.

**Scope control** — by default all fixed effects receive a sentence (`"all"`). For causal inference, where only the exposure effect is of primary interest, set:

```r
publication = list(
  template_sentences_scope = "exposure_only"
)
```

in `00_config.R`. When set to `"exposure_only"`, the pipeline reads variables with `role == "exposure"` from `00_variable_dictionary.csv` and generates sentences only for those parameters and any interaction terms that involve an exposure variable (e.g. `time × exposure`). All other fixed effects are still estimated and reported in the posterior table and forest plot — only the template sentences are filtered.

---

### Federated meta-analysis draw export (Step 12)

When multiple cohorts run this pipeline with a common DAG, their per-imputation posterior draws can be combined in a hierarchical Bayesian meta-analysis. Step 12 packages the draws into a portable long-format file that the companion repo [`generic-mi-brms-meta`](https://github.com/fabregithub/generic-mi-brms-meta) consumes.

To enable, set in `00_config.R`:

```r
export = list(
  cohort_id = "cohort_japan_2024",   # short unique ID for this dataset
  scope     = "exposure_only"        # "exposure_only" (default) or "all"
)
```

Step 12 is skipped automatically when `cohort_id` is `NULL`.

**Output files** (written to `results/export/`):

| File | Description |
|---|---|
| `cohort_draws.rds` | Long-format draws: `cohort_id`, `parameter`, `imputation`, `draw_index`, `value` |
| `cohort_draws.csv` | Same content as CSV for interoperability |
| `cohort_metadata.json` | Sidecar: cohort ID, parameter list, m, family, timestamp |

The `scope` field controls which parameters are exported. `"exposure_only"` keeps only parameters whose name corresponds to a variable with `role == "exposure"` in `00_variable_dictionary.csv` (plus interaction terms involving an exposure), keeping the export file small for causal inference analyses. `"all"` exports every parameter matched by `parameter_draw_regex`.

Ship `cohort_draws.rds` and `cohort_metadata.json` to the coordinating site; each cohort may have a different `m`.

> **Variable name harmonisation:** parameter names in the export are derived directly from the `var` column in `00_variable_dictionary.csv`. If all cohorts use the same variable dictionary template (distributed by the coordinating site), exported parameter names will be identical across cohorts and no post-hoc renaming is needed. See [EXPORT_FORMAT.md](../EXPORT_FORMAT.md) for the full export specification and a pre-transfer checklist.

---

### How Step 8's report embeds Step 11's results

The main report's "Imputation-count stability" chapter is not a separate report -- it is built from the same tables and figures that `11_check_imputation_stability.R` writes to `results/publication/mi_stability/`, referenced from the main report's `.qmd` by relative path. There is one report, one `.qmd`, and one rendered HTML/DOCX pair; `results/publication/mi_stability/` still exists alongside it as the underlying data (and as a `.qmd` of its own, for anyone who wants the full, unabridged stepwise detail across every evaluated batch), but it is not rendered separately.

This is why `run_all.R` runs Step 11 before Step 8 (see [Pipeline scripts](../README.md#pipeline-scripts)): the embedded chunks are evaluated when the report is rendered, so Step 11's files need to already be on disk at that point. `11_check_imputation_stability.R`'s own `render_quarto` option now defaults to `FALSE` for this reason -- its standalone report would just duplicate what is already in the main report's chapter. Set `analysis_spec$mi_stability$render_quarto <- TRUE` if you specifically want that standalone, more detailed report rendered as well.

---

### Optional monotonic-effect post-processing

For models that use `brms::mo()`, the standard posterior parameter table is not always the most interpretable summary. Monotonic effects are parameterised using an overall monotonic coefficient and simplex parameters, so category-specific odds ratios should be derived from posterior draws.

Two scripts are provided for this purpose:

- `09_check_mo_parameter_columns.R`
- `10_publication_mo_results.R`

`run_all.R` runs both of these automatically, but only when the fitted model's formula actually contains `mo()` terms (detected directly from `model_spec$formula`, not from any hardcoded variable list). For ordinary Gaussian, logistic, spline-only or factor-coded models, they are skipped automatically. You can also run them manually after the main pipeline has completed and `results/parameter_draws.rds` has been created:

```bash
# Bash command block
Rscript 09_check_mo_parameter_columns.R
Rscript 10_publication_mo_results.R
```

`09_check_mo_parameter_columns.R` discovers every `mo()` variable from the model formula, reports the matching `bsp_`/`simo_` columns and the implied number of ordered categories, and prints a ready-to-paste `analysis_spec$mo_effects` config block.

By default, `10_publication_mo_results.R` runs with **zero configuration**: it discovers the same `mo()` variables from the formula and labels categories generically as `"Level 1"`, `"Level 2"`, and so on. For publication-ready labels and category names, paste the block 09 printed into `00_config.R` and edit it, for example:

```r
analysis_spec$mo_effects <- list(
  vars = list(
    lwt_q = list(
      label = "Maternal weight quintile",
      levels = c("Q1", "Q2", "Q3", "Q4", "Q5")
    )
  ),

  # Only needed if your formula contains a time * mo(variable) interaction.
  time_var = NULL,
  time_values = NULL
)
```

Any `mo()` variable not listed in `analysis_spec$mo_effects$vars` still gets generic "Level N" labels automatically, so this config is optional and additive -- you only need to add entries for the variables where you want nicer labels.

Optional Quarto rendering:

```bash
# Bash command block
quarto render results/publication/mo_effects/report/mo_effects_report.qmd
```

The main output table is `results/publication/mo_effects/tables/mo_cumulative_or_table.csv`.

Additional outputs include:

- `results/publication/mo_effects/tables/mo_adjacent_or_table.csv`
- `results/publication/mo_effects/tables/mo_average_or_table.csv`
- `results/publication/mo_effects/tables/mo_simplex_table.csv`
- `results/publication/mo_effects/figures/mo_cumulative_or_plot.png`
- `results/publication/mo_effects/figures/mo_adjacent_or_plot.png`
- `results/publication/mo_effects/report/mo_effects_report.qmd`

For models using `mo()`, make sure `00_config.R` includes `bsp_` and `simo_` in the posterior draw extraction regex:

```r
model = list(
  ...
  parameter_draw_regex = "^(b_|bsp_|sd_|sigma|sds_|bs_|simo_)",
  ...
)
```

The `10_publication_mo_results.R` script can handle monotonic effects with interactions such as `time * mo(ordinal_variable)`.

In that case, set `analysis_spec$mo_effects$time_var`/`time_values`, and it calculates category-specific monotonic-effect odds ratios at the configured time values.

These scripts are not required for ordinary Gaussian, logistic, spline-only or factor-coded models, and `run_all.R` skips them automatically in that case. They only run when derived monotonic-effect odds ratios are needed.


### Optional imputation-count stability outputs

If `11_check_imputation_stability.R` is run, publication-ready imputation-count stability outputs are written to `results/publication/mi_stability/`.

These outputs summarise how selected posterior summaries change as the number of completed imputations increases. They are intended to support transparent reporting when an adaptive number of imputations is used, especially for expensive Bayesian models where blindly running `m = 100` or more may be impractical.

Main outputs include:

- `results/publication/mi_stability/tables/imputation_stability_all_batches.csv`
- `results/publication/mi_stability/tables/imputation_stability_final_comparison_display.csv`
- `results/publication/mi_stability/tables/imputation_stability_settings.csv`
- `results/publication/mi_stability/tables/imputation_stability_stepwise_summary.csv`
- `results/publication/mi_stability/tables/imputation_stability_stepwise_comparison_full.csv`
- `results/publication/mi_stability/figures/imputation_stability_trajectories.png`
- `results/publication/mi_stability/figures/imputation_stability_stepwise_change.png`
- `results/publication/mi_stability/report/imputation_stability_report.qmd`

The stability check should focus on prespecified primary estimands rather than every nuisance parameter.

---

---

## Manuscript writing guide

The publication outputs described [above](#publication-outputs-and-inference-guidance) are designed to be quoted directly rather than re-derived by hand. This section maps each manuscript paragraph to the file that supplies its numbers, and gives fill-in-the-blank templates for Methods and Results text.

These templates are starting points, not finished prose. They still need scientific interpretation, comparison to prior literature, and journal-specific formatting, all of which require human judgement that this pipeline does not provide.

### Where each paragraph's numbers come from

| Manuscript paragraph | Source file | What it gives you |
|---|---|---|
| Missing data and imputation | `results/publication/tables/analysis_metadata.csv` | Imputation strategy, target/used `m`, imputation iterations, mean-matching candidates |
| Model specification | `results/publication/tables/analysis_metadata.csv` | Model formula, family/link, priors |
| MCMC settings | `results/publication/tables/analysis_metadata.csv` | Chains, iterations, warm-up, seed, `adapt_delta`, `max_treedepth` |
| Convergence/diagnostics | `results/publication/tables/diagnostics_summary.csv` | Divergences, treedepth hits, Rhat/ESS-type checks per fit |
| Adaptive imputation-count justification | `results/publication/mi_stability/tables/imputation_stability_stepwise_summary_display.csv` | Whether/when posterior summaries stopped changing as `m` increased |
| Main fixed-effect results | `results/publication/tables/main_effect_table_display.csv` | Estimate, 95% CrI, `pd`, ROPE %, and the transformed effect (e.g. odds ratio) per variable |
| Smooth/monotonic supplementary results | `results/publication/tables/special_parameter_table.csv` | Auxiliary `s()`/`mo()` parameters not shown in the main table |
| Monotonic-effect (`mo()`) results | `results/publication/mo_effects/tables/mo_cumulative_or_table.csv` | Category-specific odds ratios for monotonic ordinal predictors |
| Posterior prediction for missing outcomes | `results/missing_y_summary.csv` | Summaries of predicted values for rows with missing outcome data |

Open `analysis_metadata.csv` directly to read off exact values for the templates below, for example:

```r
metadata <- readr::read_csv("results/publication/tables/analysis_metadata.csv")
print(metadata, n = Inf)
```

### Methods text templates

**Missing data and imputation.** Adapt to the actual `Imputation strategy` and `m` values from `analysis_metadata.csv`:

> Missing covariate data were handled using multiple imputation by chained
equations implemented with random forests (miceRanger). Variables with
missing values were imputed under a missing-at-random assumption,
conditional on the variables included in the imputation model. We
generated m = XX imputed datasets, fitted a separate Bayesian model to
each, and combined posterior draws across imputations.

If you used the subject-wide repeated-outcome strategy, adapt instead:

> For the repeated outcome, subject-level covariates were imputed once per
subject using a subject-wide imputation dataset, then merged back onto
the long-format repeated-measures data before model fitting.

**Adaptive number of imputations.** If you used the automatic m-increment loop or the manual staged workflow described under [Adapting the pipeline to your data](../README.md#adapting-the-pipeline-to-private-study-data):

> We used an adaptive multiple-imputation strategy because each imputed-data
Bayesian model was computationally expensive. We first fitted an initial
set of imputed datasets and assessed the stability of prespecified
primary posterior summaries. We increased the number of imputations until
the pooled posterior medians, credible intervals, posterior direction
probabilities and substantive conclusions changed negligibly with
additional imputations. The final analysis used m = XX imputations.

If the stepwise stability table (`imputation_stability_stepwise_summary_display.csv`) shows early flattening, a more specific statement can be used instead:

> The stability assessment showed that posterior summaries changed
negligibly after m = XX. The largest subsequent changes in posterior
medians, credible-interval endpoints and odds-ratio summaries were below
prespecified practical thresholds. The final analysis therefore used
m = XX imputations.

**Model specification.** Read `Model formula`, `Model family/link` and `Priors` from `analysis_metadata.csv`:

> We fitted a Bayesian [family] regression model with a [link] link function
using brms (formula: [paste Model formula here]) and the cmdstanr backend
for Stan. Priors were [paste Priors here]. One model was fitted separately
to each imputed dataset rather than using brm_multiple(), so that fitted
models could be checkpointed and restarted independently.

**MCMC settings.** Read `MCMC chains`, `Total/Warm-up/Post-warm-up iterations`, `Seed`, `adapt_delta` and `max_treedepth` from `analysis_metadata.csv`:

>Each model was fitted using XX chains of XX total iterations, including
XX warm-up iterations, yielding XX post-warm-up draws per chain per
imputed dataset (XX total post-warm-up draws pooled across XX imputations).

**Posterior pooling across imputations.** Worth including given the corrected pooling method (see [Adapting the pipeline to your data](../README.md#adapting-the-pipeline-to-private-study-data), "How Step 6 pools draws across imputations"):

> Posterior draws were pooled across imputed datasets with equal weight per
imputation, and an additional finite-imputation variance correction
(following Rubin's combining rules) was applied where the pooled
posterior distribution was sufficiently unimodal for this correction to
be appropriate.

**Software and reproducibility.** Fill in package versions from your R session (`sessionInfo()`) and the `Seed` value from `analysis_metadata.csv`:

> Analyses were conducted in R version XX using brms version XX with the
cmdstanr backend (CmdStan version XX), and miceRanger version XX for
multiple imputation. A fixed random seed (XX) was used for model fitting;
imputation batches were seeded deterministically for reproducibility.

### Results text templates

**Reporting a main effect**, using a row from `main_effect_table_display.csv`:

> [Variable label] was associated with [outcome]: estimate = XX (95% credible interval [CrI]: XX to XX; posterior probability of direction = XX). [If a transformed effect column is present, e.g. odds ratio:] On the odds-ratio scale, this corresponds to an odds ratio of XX (95% CrI: XX to XX).


If a ROPE was defined and is reported in the `ROPE %` column:

> The posterior probability that this effect fell within the prespecified region of practical equivalence (ROPE: XX to XX) was XX%, [supporting / not supporting] practical equivalence to no effect.


**Reporting a monotonic (`mo()`) effect**, using a row from `mo_cumulative_or_table.csv` (see [Publication outputs](#publication-outputs-and-inference-guidance) for how these are derived):

> [Variable label] showed a monotonic ordinal association with [outcome]. Compared with the lowest category, the odds ratio for [highest category] was XX (95% CrI: XX to XX; Pr(OR > 1) = XX), with intermediate categories shown in [Table/Figure XX].


### What this guide does not write for you

- Scientific interpretation of effect sizes
- Clinical or substantive significance
- Comparison with prior literature
- Study limitations
- Sample size or power justification
- Ethics, consent and data-availability statements
- Citations for software and methods (cite R, brms, Stan, miceRanger and the relevant statistical methods papers directly)

---

---

*[← Back to the main README](../README.md)*
