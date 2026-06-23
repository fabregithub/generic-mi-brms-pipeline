# Plantilla genérica de pipeline MICE + brms

Idioma: [English](../README.md) | [Français](README.fr.md) | Español | [Deutsch](README.de.md) | [日本語](README.ja.md)

> Esta es una traducción inicial del README en inglés. Si hay alguna diferencia o ambigüedad, la versión en inglés es la referencia oficial.

Este repositorio proporciona una plantilla reutilizable de pipeline en R para análisis de regresión bayesiana con imputación múltiple opcional.

Admite:

- validación de datos;
- imputación múltiple opcional con `miceRanger`;
- regresión bayesiana con `brms` + `cmdstanr`;
- guardado de un ajuste por cada conjunto de datos imputado;
- ajuste paralelo de modelos entre imputaciones;
- diagnósticos;
- resúmenes posteriores;
- predicción posterior para filas con resultados faltantes;
- tablas, figuras, metadatos de métodos/configuración y plantillas de informes listos para publicación.

El ejemplo predeterminado usa el conjunto de datos público integrado `datasets::airquality`, por lo que la plantilla puede probarse y demostrarse sin datos privados.

---

## Contenido

1. [Antecedentes y propósito](#1-antecedentes-y-propósito)
2. [Estructura del pipeline](#2-estructura-del-pipeline)
3. [Inicio rápido](#3-inicio-rápido)
4. [Adaptar el pipeline a datos privados de estudio](#4-adaptar-el-pipeline-a-datos-privados-de-estudio)
5. [Diccionario de variables](#5-diccionario-de-variables)
6. [Paralelización y ajuste de rendimiento](#6-paralelización-y-ajuste-de-rendimiento)
7. [Registros, monitoreo, reinicio y solución de problemas](#7-registros-monitoreo-reinicio-y-solución-de-problemas)
8. [Salidas de publicación y guía de inferencia](#8-salidas-de-publicación-y-guía-de-inferencia)
9. [Guía para redactar el manuscrito](#9-guía-para-redactar-el-manuscrito)
10. [Ejemplos y pruebas](#10-ejemplos-y-pruebas)
11. [Configuración del entorno de cómputo](#11-configuración-del-entorno-de-cómputo)

---

## 1. Antecedentes y propósito

Este pipeline está diseñado para análisis bayesianos aplicados cuando los datos pueden contener covariables faltantes, resultados repetidos, modelos grandes o modelos que requieren guardado cuidadoso. Combina imputación múltiple, modelado bayesiano con un ajuste por imputación, diagnósticos, resúmenes posteriores, predicción posterior y salidas orientadas a publicación.

El diseño prioriza la reproducibilidad y la capacidad de reinicio por encima de mantener todos los modelos ajustados en memoria. Por eso el pipeline ajusta un modelo `brms` por cada conjunto de datos imputado, guarda cada ajuste inmediatamente y reutiliza archivos de checkpoint válidos al volver a ejecutar.

### Advertencias y limitaciones

Este repositorio es una estructura de flujo de trabajo, no un sustituto del criterio estadístico. Antes de usarlo para un análisis científico, verifique que la estrategia de imputación, la fórmula del modelo, los priors, los diagnósticos y los resúmenes posteriores sean adecuados para su pregunta de investigación.

Los pasos de imputación múltiple están pensados para variables donde una suposición estilo MICE sea científicamente defendible, normalmente MCAR o MAR después de condicionar en variables observadas incluidas en el modelo de imputación. El pipeline no maneja automáticamente mecanismos MNAR, censura, truncamiento, límites de detección o faltantes estructurales.

Algunos modelos pueden ser computacionalmente costosos. Empiece siempre con una prueba rápida pequeña antes de una ejecución de producción.

### Notas importantes de diseño

Esta plantilla **no** usa `brm_multiple()` para ajustar modelos.

En su lugar, ajusta un modelo `brms` por conjunto imputado y guarda cada ajuste:

```text
fits/fit_imp_001.rds
fits/fit_imp_002.rds
...
fits/fit_imp_100.rds
```

Esto es más seguro para datos grandes porque:

- la sesión principal de R no mantiene todos los modelos en memoria;
- los ajustes completados se conservan si la ejecución se detiene;
- las imputaciones fallidas o lentas pueden repetirse por separado;
- los ajustes válidos existentes se omiten al volver a ejecutar;
- los procesos trabajadores devuelven solo objetos pequeños de estado.

La paralelización ocurre entre imputaciones usando `future` / `furrr`, con planificación dinámica:

```r
furrr::furrr_options(
  seed = TRUE,
  scheduling = Inf
)
```

### Patrones de análisis admitidos

```text
resultado de un solo tiempo con covariables por fila
resultado repetido con covariables a nivel de sujeto
resultado repetido con covariables variables en el tiempo
análisis de casos completos sin imputación
imputación múltiple por fila
imputación múltiple por sujeto
imputación en formato ancho por sujeto usando Y repetidos como auxiliares
```

Familias de modelos admitidas:

```text
gaussian
bernoulli
poisson
negbinomial
beta
ordinal
categorical
```

Las familias y enlaces se definen en `00_config.R`.

---

## 2. Estructura del pipeline

El repositorio se organiza alrededor de pocos archivos editados por el usuario y una secuencia de scripts numerados. En la mayoría de los proyectos, solo es necesario editar:

```text
00_config.R
00_variable_dictionary.csv
```

Los demás scripts normalmente deben tratarse como código del pipeline.

### Scripts del pipeline

`run_all.R` ejecuta estos scripts en orden:

```text
01_validate_config.R
02_prepare_data.R
03_impute.R
04_fit_models.R
05_diagnostics.R
06_posterior_summary.R
07_posterior_prediction.R
08_publication_results.R
```

`run_all.R` ejecuta además automáticamente dos scripts adicionales:

```text
09_check_mo_parameter_columns.R      se ejecuta automáticamente solo si la fórmula del modelo contiene términos mo()
10_publication_mo_results.R          se ejecuta automáticamente solo si la fórmula del modelo contiene términos mo()
11_check_imputation_stability.R      siempre se ejecuta, y produce el informe final de estabilidad del número de imputaciones
```

`run_all.R` detecta los términos `mo()` directamente en la fórmula del modelo ajustado; no es necesario mantener una lista de variables. Por eso `09`/`10` se omiten automáticamente para modelos gaussianos, logísticos, con splines o con factores ordinarios.

Si `analysis_spec$mi_stability$auto_increment <- TRUE`, `run_all.R` también aumenta `m` automáticamente por lotes y detiene el ajuste en cuanto los resúmenes posteriores son estables; vea la sección 4.

---

## 3. Inicio rápido

Este inicio rápido usa `datasets::airquality`.

```bash
git clone https://github.com/fabregithub/generic-mi-brms-pipeline.git
cd generic-mi-brms-pipeline
```

Luego ejecute:

```bash
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv
Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

Las salidas se escriben en:

```text
objects/
fits/
results/
results/publication/
```

---

## 4. Adaptar el pipeline a datos privados de estudio

Flujo recomendado:

1. Prepare un conjunto de datos limpio y guárdelo como `.rds`.
2. Edite `00_variable_dictionary.csv`.
3. Edite `00_config.R`.
4. Ejecute la validación.
5. Ejecute una prueba rápida pequeña.
6. Ejecute una prueba paralela moderada.
7. Ejecute el análisis completo de producción.
8. Renderice e inspeccione las salidas de publicación.

Prepare el archivo de datos, por ejemplo:

```text
data/my_analysis_data.rds
```

En R:

```r
saveRDS(my_data, "data/my_analysis_data.rds")
```

Y apúntelo en `00_config.R`:

```r
data = list(
  raw_data_file = "data/my_analysis_data.rds",
  ...
)
```

### Revisar el mecanismo de faltantes

Antes de imputar, revise por qué falta cada variable. El pipeline usa `miceRanger`, apropiado cuando MCAR o MAR es razonable después de condicionar en variables observadas.

No pase directamente como `NA` ordinarios:

```text
mediciones censuradas por la izquierda
valores bajo un límite de detección
mediciones censuradas por la derecha o por intervalo
faltantes estructurales
respuestas no aplicables
faltantes causados por el diseño del estudio
variables MNAR conocidas
```

Estos casos deben tratarse antes del pipeline o mediante análisis de sensibilidad.

### Desactivar imputación

```r
imputation = list(
  enabled = FALSE,
  strategy = "none",
  m = 1,
  maxiter = 0,
  mean_match_k = 5,
  verbose = FALSE,
  impute_y = FALSE,
  extra_exclude_targets = character(0)
)
```

Compruebe valores faltantes:

```r
source("00_config.R")
d <- readRDS(paths$raw_data)
colSums(is.na(d))
```

### Elegir el número de imputaciones

No existe un valor universal de `m`. Con conjuntos de datos grandes, los resultados suelen estabilizarse con un `m` menor que con conjuntos pequeños, mientras que cada ajuste tarda más; fijar un `m` grande desde el inicio puede desperdiciar mucho cómputo.

**Bucle automático (recomendado).** Defina `analysis_spec$imputation$m` como tope, y luego:

```r
analysis_spec$mi_stability <- list(
  auto_increment = TRUE,
  increment_size = NULL  # por defecto, igual a fit_workers
)
```

`run_all.R` ajusta entonces `m` por lotes (tamaño por defecto = `fit_workers`, para ocupar plenamente los workers paralelos), comprueba la estabilidad después de cada lote y deja de aumentar `m` en cuanto se cumplen los umbrales configurados, o al alcanzar el `m` configurado. Cada lote se siembra de forma determinista (`analysis_spec$imputation$seed`, por defecto basado en `analysis_spec$model$seed`).

**Estrategia manual por etapas.** Para controlar usted mismo los lotes, fije `analysis_spec$imputation$allow_extend <- TRUE` y aumente `m` progresivamente:

```text
20 -> 40 -> 60 -> 80 -> 100
```

Con `allow_extend = TRUE`, volver a ejecutar el pipeline con un `m` mayor solo agrega las imputaciones nuevas; las imputaciones existentes nunca se modifican.

Después del paso 6, ejecute:

```bash
Rscript 11_check_imputation_stability.R
```

---

## 5. Diccionario de variables

`00_variable_dictionary.csv` es la descripción legible por máquina de las variables.

Columnas esperadas:

```text
var
label
role
type
timing
scale
reference
impute_target
use_in_model
use_as_auxiliary
```

### Campos principales

- `var`: nombre exacto de la variable en los datos.
- `label`: etiqueta legible para tablas e informes.
- `impute_target`: si la variable debe imputarse cuando falta.
- `use_in_model`: si aparece en el modelo `brms` final.
- `use_as_auxiliary`: si se usa para imputación pero no para el modelo final.

`role` describe el propósito analítico: `outcome`, `binary_outcome`, `exposure`, `covariate`, `auxiliary`, `id`, `time`, `cluster`, `strata`.

`type` describe el tipo estadístico: `continuous`, `integer`, `binary`, `categorical`, `ordinal`, `date`, `id`.

`timing` describe cuándo se midió la variable: `single`, `baseline`, `repeated`, `time_varying`, `derived`, `id`.

`scale = z` crea una versión estandarizada `var_z`. No use `z` para resultados, IDs, variables binarias/categóricas o variables ya estandarizadas.

---

## 6. Paralelización y ajuste de rendimiento

Ajustes principales:

| Ajuste | Uso | Significado |
|---|---|---|
| `impute_workers` | Paso 3 | Número de workers de `miceRanger` |
| `num_impute_threads_per_worker` | Paso 3 | Threads por worker de imputación |
| `fit_workers` | Paso 4 | Imputaciones ajustadas en paralelo |
| `cores_per_fit` | Paso 4 | Cadenas/núcleos por ajuste `brms` |
| `summary_workers` | Paso 6 | Workers para extraer draws posteriores |
| `prediction_workers` | Paso 7 | Workers para predicción posterior |
| `future_globals_maxsize_gb` | Pasos 4, 6, 7 | Tamaño máximo de globals future |

Ejemplo:

```r
parallel = list(
  impute_workers = 2,
  num_impute_threads_per_worker = 2,
  num_impute_threads = 2,
  fit_workers = 4,
  cores_per_fit = 4,
  summary_workers = 2,
  prediction_workers = 2,
  future_globals_maxsize_gb = 80
)
```

Empiece con valores conservadores, especialmente en portátiles o con datos grandes.

---

## 7. Registros, monitoreo, reinicio y solución de problemas

Comandos comunes:

```bash
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_stdout.log
Rscript fit_single_imputation.R 1
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

Si se interrumpe una ejecución, vuelva a ejecutar:

```bash
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Los ajustes válidos existentes se omiten.

Archivos de seguimiento:

```text
pipeline_progress.log
pipeline_stdout.log
pipeline_heartbeat.txt
pipeline_success.flag
pipeline_error.flag
results/worker_logs/
```

Monitorear progreso:

```bash
tail -f pipeline_progress.log
```

Si aparece un problema de caché de CmdStanR:

```bash
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

---

## 8. Salidas de publicación y guía de inferencia

Después de una ejecución exitosa, las salidas están en:

```text
results/publication/
```

Salidas típicas:

```text
results/publication/tables/main_effect_table_display.csv
results/publication/tables/main_effect_table_full.csv
results/publication/tables/diagnostics_summary.csv
results/publication/tables/analysis_metadata.csv
results/publication/figures/forest_plot_odds_ratios.png
results/publication/report/bayesian_mi_report_template.qmd
```

Renderice el informe:

```bash
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

En el paso 6, los draws posteriores no se concatenan simplemente entre imputaciones: cada imputación contribuye con peso igual (1/m), y se aplica una corrección de varianza para `m` finito (regla de Rubin) cuando la forma de la distribución posterior combinada lo permite. Vea las columnas `m_imputations`, `variance_corrected`, etc. en `results/parameter_summary.csv`.

Para modelos con `brms::mo()`, `run_all.R` ejecuta automáticamente:

```bash
Rscript 09_check_mo_parameter_columns.R
Rscript 10_publication_mo_results.R
quarto render results/publication/mo_effects/report/mo_effects_report.qmd
```

Estos scripts detectan las variables `mo()` directamente en la fórmula del modelo; no se necesita ninguna lista de variables editada manualmente. Para etiquetas y nombres de categorías personalizados, agregue un bloque opcional `analysis_spec$mo_effects` en `00_config.R` (plantilla sugerida impresa por `09_check_mo_parameter_columns.R`).

---

## 9. Guía para redactar el manuscrito

Las salidas de publicación están diseñadas para citarse directamente en un manuscrito, en lugar de recalcularse a mano.

| Párrafo del manuscrito | Archivo fuente |
|---|---|
| Datos faltantes / imputación | `results/publication/tables/analysis_metadata.csv` |
| Especificación del modelo, ajustes MCMC | `results/publication/tables/analysis_metadata.csv` |
| Justificación del número adaptativo de imputaciones | `results/publication/mi_stability/tables/imputation_stability_stepwise_summary_display.csv` |
| Resultados de efectos principales | `results/publication/tables/main_effect_table_display.csv` |
| Resultados de efectos monotónicos (`mo()`) | `results/publication/mo_effects/tables/mo_cumulative_or_table.csv` |

Ejemplo de texto de "Métodos" (ajuste los valores desde `analysis_metadata.csv`):

```text
Los datos faltantes de covariables se trataron mediante imputación
múltiple por ecuaciones encadenadas con bosques aleatorios (miceRanger),
generando m = XX conjuntos de datos imputados. Se ajustó un modelo
bayesiano independiente a cada conjunto imputado usando brms con el
motor cmdstanr, y los draws posteriores se combinaron entre imputaciones
con peso igual por imputación.
```

Ejemplo de texto de "Resultados" para un efecto principal (desde `main_effect_table_display.csv`):

```text
[Variable] se asoció con [resultado]: estimación = XX (IC del 95 %: XX a
XX; probabilidad posterior de dirección = XX).
```

Esta guía no redacta por usted la interpretación científica, la relevancia clínica, las limitaciones del estudio ni las citas bibliográficas; eso sigue requiriendo criterio humano.

---

## 10. Ejemplos y pruebas

Ejemplos incluidos:

```text
examples/airquality_gaussian
examples/birthwt_logistic
examples/birthwt_spline_monotonic
```

Pruebas rápidas:

```bash
bash test/test_all_examples_quick.sh
```

Pruebas paralelas:

```bash
bash test/test_all_examples_parallel.sh
```

Limpiar antes de cambiar de ejemplo:

```bash
rm -rf objects fits results
rm -f pipeline_error.flag pipeline_success.flag
rm -f pipeline_progress.log pipeline_heartbeat.txt pipeline_stdout.log run_all_stdout.log
```

---

## 11. Configuración del entorno de cómputo

Instale R, las herramientas del sistema, los paquetes R, CmdStan/CmdStanR y Quarto.

Paquetes R principales:

```r
install.packages(c(
  "tidyverse", "miceRanger", "brms", "posterior", "bayestestR",
  "future", "furrr", "doParallel", "foreach", "gt", "flextable",
  "officer", "forcats", "glue", "readr", "tibble", "dplyr",
  "stringr", "purrr", "rlang"
))
```

`cmdstanr` y CmdStan:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

library(cmdstanr)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_path()
```

Verifique Quarto:

```bash
quarto --version
```

Para reproducibilidad, considere `renv`:

```r
install.packages("renv")
renv::init()
renv::snapshot()
renv::restore()
```
