# Birthwt custom-formula demo: spline + monotonic effect

This example tests the modified custom-formula reporting path using the public `MASS::birthwt` dataset.

It fits a Bernoulli/logit model with:

```
low ~ s(age_z, k = 5) + mo(lwt_q) + race + smoke + ptl_z + ht + ui + ftv_z
```

The example checks:

```
s() validation and reporting
mo() validation and reporting
special parameter extraction
conditional-effect plot creation
publication report updates
```

## Files

```
00_create_birthwt_spline_monotonic_example_data.R
00_config_birthwt_spline_monotonic.R
00_variable_dictionary_birthwt_spline_monotonic.csv
```

## Run from the project root folder

Run in Terminal:

```
cp examples/birthwt_spline_monotonic/00_config_birthwt_spline_monotonic.R 00_config.R
cp examples/birthwt_spline_monotonic/00_variable_dictionary_birthwt_spline_monotonic.csv 00_variable_dictionary.csv
Rscript examples/birthwt_spline_monotonic/00_create_birthwt_spline_monotonic_example_data.R
```

Clean previous outputs before switching examples:

```
rm -rf objects fits results
rm -f pipeline_error.flag pipeline_success.flag pipeline_progress.log pipeline_heartbeat.txt run_all_stdout.log
```

Then run:

```
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_birthwt_spline_monotonic_stdout.log
```

## Expected custom-formula outputs

If the modified reporting scripts are working, publication outputs should include files such as:

```
results/publication/tables/special_parameter_table.csv
results/publication/tables/conditional_effects_manifest.csv
results/publication/figures/conditional_effect_age_z.png
results/publication/figures/conditional_effect_lwt_q.png
```

The exact set of parameter names may vary by brms version.
