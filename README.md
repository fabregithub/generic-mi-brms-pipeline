# Generic MICE + brms Pipeline Template

📖 **Documentation site: <https://fabregithub.github.io/generic-mi-brms-pipeline/>** —
the same guides as in this repository, with a sidebar and full-text search.

This is a reusable R pipeline template for Bayesian regression analyses with optional multiple imputation.

It supports:

- data validation
- optional multiple imputation with `miceRanger`
- Bayesian regression with `brms` + `cmdstanr`
- one-fit-per-imputation checkpointing
- parallel model fitting across imputations
- diagnostics
- posterior summaries
- posterior prediction for rows with missing outcomes
- publication-ready tables, figures, methods/settings metadata and report templates

The default example uses the built-in public dataset `datasets::airquality`, so the template can be tested and demonstrated without private data.

---

## Contents

- [Background and purpose](#background-and-purpose)
- [Structure of the pipeline](#structure-of-the-pipeline)
- [Quick start](#quick-start)
- [Adapting the pipeline to private study data](#adapting-the-pipeline-to-private-study-data)
- [Variable dictionary](#variable-dictionary) → [`docs/variable-dictionary.md`](docs/variable-dictionary.md)
- [Parallelisation and performance tuning](#parallelisation-and-performance-tuning) → [`docs/performance.md`](docs/performance.md)
- [Running, monitoring and recovering](#running-monitoring-and-recovering) → [`docs/operations.md`](docs/operations.md)
- [Publication outputs and manuscript writing](#publication-outputs-and-manuscript-writing) → [`docs/reporting.md`](docs/reporting.md)
- [Examples and tests](#examples-and-tests) → [`docs/examples.md`](docs/examples.md)
- [Computing environment setup](#computing-environment-setup) → [`docs/setup.md`](docs/setup.md)

---

## Documentation

The README is the entry point: what this is, whether it fits your analysis, and how to
start. Detailed material lives in focused guides under [`docs/`](docs/index.md).

Everything below is also published as a browsable site with search — useful for the longer
reference pages, where GitHub's in-repo search struggles:

**<https://fabregithub.github.io/generic-mi-brms-pipeline/>**

The site is built from these same Markdown files by
[`_quarto.yml`](_quarto.yml) on every push to `main`, so it cannot drift from the
repository.

| Guide | What it covers |
|---|---|
| [`docs/setup.md`](docs/setup.md) | Installing R, CmdStan, Quarto and packages on macOS, Windows and Linux |
| [`docs/variable-dictionary.md`](docs/variable-dictionary.md) | Field reference for `00_variable_dictionary.csv`: `role`, `type`, `timing`, `scale`, `reference`, `use_as_auxiliary` |
| [`docs/performance.md`](docs/performance.md) | Sizing `impute_workers`, `fit_workers` and `chains`; memory-vs-cores trade-offs |
| [`docs/examples.md`](docs/examples.md) | Four worked examples on public data, and the automated test scripts |
| [`docs/operations.md`](docs/operations.md) | Restarting an interrupted run, monitoring, CmdStan cache, debugging a fit |
| [`docs/config-reference.md`](docs/config-reference.md) | Option-by-option reference for `00_config.R` |
| [`docs/choosing-m.md`](docs/choosing-m.md) | How many imputations, the auto m-increment loop, and the three run modes |
| [`docs/covariate-roles.md`](docs/covariate-roles.md) | Which covariates to adjust for — confounder, mediator, collider or precision — and how much bias the imputation adds when one is causally active |
| [`docs/censored-exposures.md`](docs/censored-exposures.md) | Below-detection-limit exposures: the engine, its validation, and its scope limit |
| [`docs/repeated-measures.md`](docs/repeated-measures.md) | Several rows per subject: subject-wide vs row-wise imputation |
| [`docs/reporting.md`](docs/reporting.md) | Publication outputs, and Methods/Results text templates |

Validation evidence — what has been tested, what may be claimed, and the limits of each
claim — is indexed in [`validation/README.md`](validation/README.md).


---

## Background and purpose

This pipeline is intended for applied Bayesian regression analyses where the data may contain missing covariates, repeated outcomes, large models or models that need careful checkpointing. It combines multiple imputation, one-fit-per-imputation Bayesian modelling, diagnostics, posterior summaries, posterior prediction and publication-oriented outputs.

The design prioritises reproducibility and restartability over keeping all fitted models in memory. This is why the pipeline fits one `brms` model per imputed dataset, saves each fit immediately and reuses valid checkpoint files on rerun.

By default (`analysis_spec$mi_stability$auto_increment = TRUE` in `00_config.R`), the pipeline also avoids fitting more imputations than necessary: `run_all.R` fits imputations in batches, checks whether posterior summaries have stabilised, and stops increasing `m` automatically once they have, rather than always fitting `analysis_spec$imputation$m` up front. This is especially useful for large-N or otherwise expensive models, where blindly fitting `m = 100` imputations can cost far more compute than the analysis needs. See [Choosing the number of imputations adaptively](docs/choosing-m.md) for details. If you prefer a single fixed `m` with no automatic stopping, set `auto_increment <- FALSE`.

### Brief cautions and limitations

This repository is a workflow scaffold, not a substitute for statistical judgement. Before using it for a scientific analysis, check that the imputation strategy, model formula, priors, diagnostics and posterior summaries are appropriate for your study question.


The multiple-imputation steps are intended for variables where a standard MICE-style assumption is scientifically defensible, typically missing completely at random (MCAR) or missing at random (MAR) after conditioning on observed variables included in the imputation model. The pipeline does not automatically handle non-ignorable missingness, MNAR mechanisms, censoring, truncation, limit-of-detection problems or structural missingness.

If a variable has non-standard missingness, such as left-censored values below a detection limit, skip-pattern missingness or values missing for design reasons, process or model that missingness appropriately before using this pipeline. Do not simply code such values as ordinary `NA` and rely on the default MICE workflow unless that is justified for the study.

Some models can be computationally expensive. In particular, large mixed logistic models, spline terms, monotonic ordinal terms and many imputations can take substantial time. Always start with a small quick test before a production run.

**The covariate list is a modelling decision, and coverage will not check it for you.** The
variable dictionary has one `role` value, `covariate`, for confounders, mediators, colliders
and precision covariates alike — it controls handling, not identification. Validation measured
what that costs: adjusting for a **collider** moved an exposure coefficient from **+0.40 to
−0.10** on complete data, and missingness in a **confounder or mediator** adds **+5% to +10%**
bias that the credible interval's coverage does not reveal (0.93–0.95). Neither is repaired by
imputation. **→ [`docs/covariate-roles.md`](docs/covariate-roles.md)**

**Install `dbarts` before a real analysis.** With it absent, `z_imputer = "bart"` warns and
falls back to `forest_boot`, which under a confounder carries **−23% bias that does not shrink
with sample size**. Confirm which imputer ran:
`grep "Z-block imputer" run_all_stdout.log`.

---

### Important design notes

This template does **not** use `brm_multiple()` for model fitting.

Instead, it fits one `brms` model per imputed dataset and saves each fit immediately:

```text
fits/fit_imp_001.rds
fits/fit_imp_002.rds
...
fits/fit_imp_100.rds
```

This is intentionally safer for large datasets because:

- the main R session does not hold all fitted models in memory;
- completed fits are preserved if the run stops;
- failed or slow imputations can be re-run separately
- valid existing fits are skipped on re-run; and
- worker processes return only small status objects to the main session.

Parallelism happens across imputations using `future` / `furrr`, with dynamic scheduling. This is configured inside the R scripts:

```r
furrr::furrr_options(
  seed = TRUE,
  scheduling = Inf
)
```

This improves load balancing when some imputed datasets take longer than others.

---

### Supported analysis patterns

The template is designed to support:

- `single-time outcome, row-level covariates`
- `repeated outcome with subject-level covariates`
- `repeated outcome with time-varying covariates`
- `complete-case analysis without imputation`
- `row-level multiple imputation`
- `subject-level multiple imputation`
- `subject-wide imputation using repeated Y as auxiliary variables`
- `left-/interval-censored exposures below a detection limit (congenial, outcome-aware imputation)`

Supported model families include:

- `gaussian`
- `bernoulli`
- `poisson`
- `negbinomial`
- `beta`
- `ordinal`
- `categorical`
- `cox` (Cox proportional hazards)

Model families and links are set in `00_config.R`.

### Censored (below-detection-limit) exposures

Exposures reported only as "below the limit of detection" can be imputed jointly with the
outcome via the `censored_exposure_block_fcs` strategy, rather than substituted with LOD/2
or LOD/√2 before modelling. Substituting a constant attenuates the exposure–response
coefficient, and imputing *without* conditioning on the outcome attenuates it too — by up
to −14% at 40% non-detects.

**→ [Censored exposures: `docs/censored-exposures.md`](docs/censored-exposures.md)** — how
the two-block engine works, what has been validated and to what precision, the correction
history (v1.4.0, v1.5.0, v1.5.1), and the recommended `m`.

The scope limit is reproduced here in full rather than linked, because acting on it matters
more than finding it:

> ### ⚠️ Scope: additive exposure–response functions only. **Do not use this strategy for
> mixture / BKMR analyses.**
>
> The claim above holds for **additive / per-analyte** exposure–response functions.
>
> For **mixture / BKMR-style surfaces** with interactions and curvature, this strategy is
> **not merely inadequate — it is worse than doing nothing sophisticated.** Measured on the
> estimands a BKMR analysis actually reports, against analytic truth, at 200 replications
> per cell:
>
> | | pipeline | LOD/√2 substitution | complete case |
> |---|---|---|---|
> | mean abs. error vs the oracle, 21 cells | **14.6%** | **6.4%** | 10.5% |
> | pure curvature, 40% non-detects | **−57%** | −1.9% *(n.s.)* | −18% |
>
> It **destroys 57% of the curvature** at 40% non-detects, and on that estimand LOD/√2
> substitution is statistically indistinguishable from the oracle while this strategy is
> not. The reason is structural: the X block draws the censored exposure from a conditional
> that is **linear in the predictors**, so it imposes linearity on precisely the rows where
> curvature would appear. Imputing with the wrong functional form is worse than not
> imputing.
>
> Interaction estimands are more robust — they survive censoring of a single exposure
> (+2.6% to +5.7%) and break down (−12.5%) only when a second exposure is censored too.
>
> **Coverage will not warn you.** In these runs a −42% point bias sat behind 0.945–0.960
> coverage, because the pooled intervals ran 1.5–1.6× wider than the estimates' true
> sampling spread.
>
> Mixture analyses need **substantive-model-compatible imputation**, which this pipeline
> does not implement. Evidence:
> [`validation/phase1/FINDINGS_v3.md`](validation/phase1/FINDINGS_v3.md); background in
> [`validation/PLAN_leftcensored_exposure_integration.md`](validation/PLAN_leftcensored_exposure_integration.md) §7.7.

## Structure of the pipeline

The repository is organised around a small set of user-edited files and a sequence of numbered pipeline scripts. In most projects, users only need to edit:

- `00_config.R`
- `00_variable_dictionary.csv`

All other scripts should usually be treated as pipeline code.

---

### Pipeline scripts

The core pipeline is run by `run_all.R`, which calls these scripts in order:

1. `01_validate_config.R`
2. `02_prepare_data.R`
3. `03_impute.R`
4. `04_fit_models.R`
5. `05_diagnostics.R`
6. `06_posterior_summary.R`
7. `07_posterior_prediction.R`
8. `11_check_imputation_stability.R` -- always runs, producing the imputation-count stability tables/figures
9. `08_publication_results.R` -- writes and renders the combined main report

Step 11 runs *before* Step 8, not after, even though Step 11's checks are about whether `m` was large enough -- a question that only makes sense once fitting is done. The reason for the order is that Step 8's report template embeds Step 11's tables and figures directly, in an "Imputation-count stability" chapter (see [How Step 8's report embeds Step 11's results](docs/reporting.md#how-step-8s-report-embeds-step-11s-results) below), so Step 11's output files need to already exist on disk by the time Step 8 writes and renders the report.

`run_all.R` then automatically runs two more scripts:

- `09_check_mo_parameter_columns.R` -- runs automatically only if the fitted model's formula contains `mo()` terms
- `10_publication_mo_results.R` -- runs automatically only if the fitted model's formula contains `mo()` terms

`run_all.R` detects `mo()` terms directly from the fitted model's own formula (via `extract_special_term_vars()`), so steps 09-10 are skipped automatically for ordinary Gaussian, logistic, spline-only or factor-coded models, and run automatically whenever the formula contains `mo()`. You do not need to remember to run them manually, and you do not need to edit them per study; see [the reporting guide](docs/reporting.md) for how they discover monotonic-effect variables generically.

Step 11 runs on whichever final `m` the run settled on. It creates stepwise publication-ready stability summaries, including tables and figures that quantify how much estimates, credible intervals, posterior direction probabilities and transformed summaries change when `m` is increased.

If `analysis_spec$mi_stability$auto_increment = TRUE`, `run_all.R` also drives Steps 3/4/6 through an automatic imputation-count loop instead of fitting a single fixed `m` up front: it fits a small batch of imputations, checks stability against the previous batch, and stops increasing `m` as soon as the configured stability thresholds are met (or once the configured `m` is reached). See [Choosing the number of imputations adaptively](docs/choosing-m.md) for details.

All of `01`-`07`, `11`, `08`, and `09`-`10` (when applicable) can also be run manually and standalone with `Rscript <script>.R`, which is useful for debugging or for resuming a partially completed run. If you run Step 8 manually before Step 11 has ever been run, its report's "Imputation-count stability" chapter is simply omitted, rather than failing.

---

### Main files to edit for a new project

Usually edit only:

- `00_config.R`
- `00_variable_dictionary.csv`

`00_config.R` defines the analysis structure, including:

- outcome variable
- model family and link
- data structure
- imputation strategy
- brms formula settings
- priors
- MCMC settings
- parallel settings
- memory guard settings

The `examples/` folder contains ready-to-run public-data configurations. These are useful for testing the pipeline and for learning how to structure a new analysis.

`00_variable_dictionary.csv` defines:

- variable labels
- variable roles
- variable types
- reference categories
- scaling
- imputation targets
- model inclusion

---

## Quick start

This quick start uses the public Gaussian example based on `datasets::airquality`.

If R, CmdStan or Quarto are not yet installed, see the [setup guide](docs/setup.md) first.

First, obtain a local copy of the template repository.

Using Git:

```bash
# Bash command block
git clone https://github.com/fabregithub/generic-mi-brms-pipeline.git
cd generic-mi-brms-pipeline
```

Alternatively, download the repository as a ZIP file from GitHub, unzip it, and open Terminal in the unzipped project folder.

Then run:

```bash
# Bash command block
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv
Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

Outputs are written to:

- `objects/`
- `fits/`
- `results/`
- `results/publication/`

### Running the pipeline

The launchers (`launch.R` and `launch.sh`) are the intended entry points for most users.
`run_all.R` is the internal orchestrator called by the launchers — you do not need to
run it directly.

#### RStudio — Windows, Mac, Linux (no terminal needed)

Open the project folder in RStudio, open `launch.R`, and click **Source**
(`Ctrl+Shift+S` / `Cmd+Shift+S`). A numbered menu appears in the R console:

```
1.  Validate config only
2.  Run full pipeline
3.  Prepare data
...
11. Clean ALL outputs
12. Clean fits/posteriors only
 q. Quit
```

Type a number and press Enter. No terminal required.

#### Mac / Linux terminal

```bash
bash launch.sh
```

The same numbered menu appears in the terminal. Alternatively, make it executable
once and run it directly:

```bash
chmod +x launch.sh
./launch.sh
```

#### Advanced: run a single imputed dataset

Use this when debugging a specific imputation (e.g. imputation 51 failed):

```bash
Rscript fit_single_imputation.R 51
```

This is available as a menu option in both `launch.sh` and `launch.R`.

`run_all.R` (via `08_publication_results.R`) renders the main Quarto report automatically, so this step is not normally needed. If you want to re-render it manually, for example after editing the generated `.qmd` by hand, run in Terminal:

```bash
# Bash command block
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

For what to do if a run is interrupted, how to monitor a long run, or how to debug a problematic imputed dataset, see [the operations guide](docs/operations.md).

---

## Adapting the pipeline to private study data

For a new study, the recommended workflow is:

1. Prepare one clean analysis dataset and save it as an `.rds` file.
2. Edit `00_variable_dictionary.csv`.
3. Edit `00_config.R`.
4. Run validation.
5. Run a small quick test.
6. Run a modest parallel test.
7. Run the full production analysis.
8. Render and inspect publication outputs.

### Data preparation

Before running the pipeline, prepare one clean input dataset, for example `data/my_analysis_data.rds`.

The dataset should already contain consistent variable names, explicit ID variables if needed, explicit time or wave variables for repeated data, and any derived variables that are not created by the pipeline. Save the dataset in R with:

```r
saveRDS(my_data, "data/my_analysis_data.rds")
```

Then point `00_config.R` to it by editing the `data = list(...)` block:

```r
data = list(
  raw_data_file = "data/my_analysis_data.rds",
  ...
)
```

### Check the missing-data mechanism before imputation

Before using the imputation step, review why each variable is missing.

The pipeline uses MICE-style multiple imputation through `miceRanger`. This is suitable when treating the missing values as MCAR or MAR is reasonable after conditioning on observed variables in the imputation model.

It is not a general solution for non-ignorable missingness or censoring. In particular, do not pass the following directly to the pipeline as ordinary `NA` values without careful pre-processing:

- left-censored measurements
- values below a detection limit
- right-censored or interval-censored measurements
- structural missingness
- skip-pattern missingness
- not-applicable responses
- missingness caused by study design
- known MNAR variables

For these cases, first create an analysis-ready representation outside the pipeline. Depending on the scientific context, this might involve:

- using a censored-data model outside this pipeline
- creating a below-detection-limit indicator
- using an appropriate substitution or interval representation
- creating explicit "not applicable" categories
- separating structural missingness from true missingness
- conducting sensitivity analyses for MNAR assumptions

Only variables that can reasonably be treated as ordinary missing values under the chosen imputation model should be set as imputation targets in `00_variable_dictionary.csv`.

### Decision tree: where to start

```text
If your data are one row per person or analytic unit
  -> use data_structure = "single_time"

If your data are one row per subject-time observation
  -> use a repeated-data structure and set id_var and time_var

If your covariates are mostly measured once per subject
  -> consider subject_wide_with_repeated_y_auxiliary imputation

If your predictors change over time
  -> use a repeated/time-varying pattern and check that timing is set correctly

If there are no missing covariates to impute, or if you want a complete-data analysis
  -> edit the imputation = list(...) block so enabled = FALSE,
     strategy = "none", and m = 1

If you need ordinary linear/logistic effects
  -> use fixed_effects = "auto" and control variables through the dictionary

If you need nonlinear continuous effects
  -> use custom_formula with s()

If you need ordinal monotonic effects
  -> mark the variable as type = ordinal, use custom_formula with mo(),
     and run optional scripts 09 and 10 for derived odds-ratio summaries
```

### Important `00_config.R` options

`00_config.R` holds a single nested list, `analysis_spec`, plus a `paths` list. With
`00_variable_dictionary.csv` it is one of only two files a new analysis normally needs to
edit.

**→ [`00_config.R` reference: `docs/config-reference.md`](docs/config-reference.md)**

Data and outcome specification, the imputation block (strategy, `m`, `z_imputer`, censored
exposures), model formula and priors, sampler settings, pooling and ROPE, parallelism, and
publication options.

### Notes for adapting to private study data

For a new project:

1. Replace or edit `00_variable_dictionary.csv`.
2. Edit `00_config.R`.
3. Run in Terminal:

```bash
# Bash command block
Rscript 01_validate_config.R
```

4. If validation passes, run in Terminal:

```bash
# Bash command block
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```


### Choosing `m`, and running at scale

`m` — the number of imputed datasets — drives both accuracy and runtime, and it is the
setting people most often guess at. The pipeline can determine it adaptively: it fits
imputations in batches and stops once the pooled posterior summaries stabilise.

**→ [Choosing `m`, and running at scale: `docs/choosing-m.md`](docs/choosing-m.md)**

The automatic m-increment loop, fixed-`m` runs, the FMI-based rule of thumb, and the three
run modes — quick test, parallel test, and full production.

### Repeated-outcome data and subject-wide imputation

With several rows per subject, imputing row-wise can give the same subject different values
of a time-invariant covariate at different visits. The pipeline can impute subject-wide
instead, driven by the `timing` column of the variable dictionary.

**→ [Repeated measures: `docs/repeated-measures.md`](docs/repeated-measures.md)**

### Ordinal predictors and `mo()`

Ordinal variables should be marked in the dictionary:

```r
type = ordinal
```

The pipeline converts these variables to ordered factors before modelling. They can then be used in a custom `brms` formula:

```r
custom_formula = brms::bf(
  y ~ mo(education) + mo(income) + age_z + sex + (1 | id)
)
```

For large mixed logistic models, `mo()` can be much slower than ordinary factor coding. A practical approach is to use the factor-coded model as the main analysis and use `mo()` as a sensitivity analysis with fewer imputations.

You do not need to hardcode `mo()` variable names anywhere. `09_check_mo_parameter_columns.R` and `10_publication_mo_results.R` discover them directly from the fitted model's own formula, and `run_all.R` runs both scripts automatically only when `mo()` terms are present. See [the reporting guide](docs/reporting.md) for how to optionally supply nicer labels and category levels.

---

### Variable roles: dictionary or config?

By default the variable dictionary is authoritative: `role`, `type`, `timing`, `scale`,
`reference` and `use_as_auxiliary` are read from `00_variable_dictionary.csv`. A few
entries in `00_config.R` may override it, which is occasionally useful and easy to do by
accident.

**→ [Dictionary vs config, and the optional overrides](docs/variable-dictionary.md#variable-roles-dictionary-by-default-config-as-optional-override)**

## Variable dictionary

`00_variable_dictionary.csv` tells the pipeline what each column of your data *is* — its
role in the model, its type, whether it varies within a subject, how it should be scaled,
and whether it is an imputation auxiliary. With `00_config.R` it is one of only two files
a new analysis normally needs to edit.

**→ [Full field reference: `docs/variable-dictionary.md`](docs/variable-dictionary.md)**

Covers `role`, `type`, `timing`, `scale`, `reference`, `use_as_auxiliary`, when a config
entry may override the dictionary, and worked dictionaries for each bundled example.

## Parallelisation and performance tuning

The pipeline parallelises at two levels — across imputed datasets during imputation, and
across model fits afterwards — and the two compete for the same cores and memory.

**→ [Performance tuning: `docs/performance.md`](docs/performance.md)**

Suggested settings by computing environment (laptop, workstation, HPC), parallel
`miceRanger` imputation, memory-vs-cores trade-offs, and recommended high-performance
configurations.

## Running, monitoring and recovering

Long runs get interrupted. The pipeline checkpoints each fit to disk as it completes, so a
rerun skips work that already finished rather than starting over.

**→ [Running, monitoring and recovering: `docs/operations.md`](docs/operations.md)**

Restarting after an interruption, logging and progress monitoring, troubleshooting the
CmdStan cache, and debugging model fitting.

## Publication outputs and manuscript writing

Step 8 writes a self-rendering Quarto report, posterior tables, forest plots and
ready-to-quote template sentences. Those outputs are designed to be pasted into a
manuscript rather than re-derived by hand.

**→ [Reporting guide: `docs/reporting.md`](docs/reporting.md)**

Publication outputs and inference guidance (template sentences, federated draw export,
monotonic-effect post-processing, imputation-stability outputs) and the manuscript writing
guide (where each paragraph's numbers come from, Methods and Results templates, and what
the guide deliberately does not write for you).

## Examples and tests

Four end-to-end examples ship with the pipeline, each on a public dataset, plus automated
test scripts that run them in isolation.

**→ [Examples and tests: `docs/examples.md`](docs/examples.md)**

Gaussian (`airquality`), logistic (`birthwt`), spline + monotonic effects (`birthwt`) and
Cox proportional hazards (`lung`); how to clean outputs between examples; and the `test/`
scripts that run each example in a throwaway copy without touching your own `objects/`,
`fits/` or `results/`.

## Computing environment setup

The pipeline needs R, CmdStan (via `cmdstanr`), Quarto for report rendering, and a handful
of R packages. Setup is a once-per-machine task.

**→ [Full installation guide: `docs/setup.md`](docs/setup.md)**

Step-by-step instructions for **macOS**, **Windows** and **Linux**, plus verifying CmdStan
and `brms` actually work, clearing a stale CmdStanR cache, `renv` for reproducibility, and
the dependency list.
