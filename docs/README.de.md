# Generische MICE- + brms-Pipeline-Vorlage

Sprache: [English](../README.md) | [Français](README.fr.md) | [Español](README.es.md) | Deutsch | [日本語](README.ja.md)

> Dies ist eine erste deutsche Übersetzung des englischen README. Bei Unterschieden oder Unklarheiten ist die englische Version maßgeblich.

Dieses Repository stellt eine wiederverwendbare R-Pipeline-Vorlage für bayesianische Regressionsanalysen mit optionaler multipler Imputation bereit.

Sie unterstützt:

- Datenvalidierung;
- optionale multiple Imputation mit `miceRanger`;
- bayesianische Regression mit `brms` + `cmdstanr`;
- Checkpointing mit einem Modell-Fit pro Imputation;
- paralleles Modell-Fitting über Imputationen hinweg;
- Diagnostik;
- Posterior-Zusammenfassungen;
- Posterior-Prädiktion für Zeilen mit fehlendem Outcome;
- publikationsreife Tabellen, Abbildungen, Methoden-/Einstellungsmetadaten und Berichtsvorlagen.

Das Standardbeispiel verwendet den integrierten öffentlichen Datensatz `datasets::airquality`, sodass die Vorlage ohne private Daten getestet und demonstriert werden kann.

---

## Inhalt

1. [Hintergrund und Zweck](#1-hintergrund-und-zweck)
2. [Struktur der Pipeline](#2-struktur-der-pipeline)
3. [Schnellstart](#3-schnellstart)
4. [Anpassung an private Studiendaten](#4-anpassung-an-private-studiendaten)
5. [Variablenwörterbuch](#5-variablenwörterbuch)
6. [Parallelisierung und Performance-Tuning](#6-parallelisierung-und-performance-tuning)
7. [Logging, Monitoring, Neustart und Fehlerbehebung](#7-logging-monitoring-neustart-und-fehlerbehebung)
8. [Publikationsausgaben und Inferenzhinweise](#8-publikationsausgaben-und-inferenzhinweise)
9. [Beispiele und Tests](#9-beispiele-und-tests)
10. [Einrichtung der Rechenumgebung](#10-einrichtung-der-rechenumgebung)

---

## 1. Hintergrund und Zweck

Diese Pipeline ist für angewandte bayesianische Regressionsanalysen gedacht, bei denen Daten fehlende Kovariaten, wiederholte Outcomes, große Modelle oder Modelle mit sorgfältigem Checkpointing enthalten können. Sie kombiniert multiple Imputation, bayesianische Modellierung mit einem Fit pro imputiertem Datensatz, Diagnostik, Posterior-Zusammenfassungen, Posterior-Prädiktion und publikationsorientierte Ausgaben.

Das Design priorisiert Reproduzierbarkeit und Wiederanlaufbarkeit gegenüber dem Halten aller angepassten Modelle im Speicher. Deshalb passt die Pipeline ein `brms`-Modell pro imputiertem Datensatz an, speichert jeden Fit sofort und verwendet gültige Checkpoint-Dateien bei erneuter Ausführung wieder.

### Hinweise und Grenzen

Dieses Repository ist ein Workflow-Gerüst, kein Ersatz für statistisches Urteilsvermögen. Prüfen Sie vor einer wissenschaftlichen Analyse, ob Imputationsstrategie, Modellformel, Priors, Diagnostik und Posterior-Zusammenfassungen zu Ihrer Fragestellung passen.

Die Imputationsschritte sind für Variablen gedacht, bei denen eine MICE-artige Annahme wissenschaftlich vertretbar ist, typischerweise MCAR oder MAR nach Konditionierung auf beobachtete Variablen im Imputationsmodell. Die Pipeline behandelt nicht automatisch MNAR-Mechanismen, Zensierung, Trunkierung, Nachweisgrenzen oder strukturell fehlende Werte.

Einige Modelle können rechenintensiv sein. Starten Sie immer mit einem kleinen Schnelltest vor einem Produktionslauf.

### Wichtige Designhinweise

Diese Vorlage verwendet **nicht** `brm_multiple()` für das Modell-Fitting.

Stattdessen wird pro imputiertem Datensatz ein `brms`-Modell angepasst und gespeichert:

```text
fits/fit_imp_001.rds
fits/fit_imp_002.rds
...
fits/fit_imp_100.rds
```

Das ist für große Datensätze sicherer, weil:

- die Haupt-R-Sitzung nicht alle Fits im Speicher hält;
- abgeschlossene Fits erhalten bleiben, wenn der Lauf stoppt;
- fehlgeschlagene oder langsame Imputationen separat erneut ausgeführt werden können;
- gültige bestehende Fits bei erneuter Ausführung übersprungen werden;
- Worker-Prozesse nur kleine Statusobjekte zurückgeben.

Parallelisierung erfolgt über Imputationen hinweg mit `future` / `furrr` und dynamischer Planung:

```r
furrr::furrr_options(
  seed = TRUE,
  scheduling = Inf
)
```

### Unterstützte Analysemuster

```text
Outcome zu einem Zeitpunkt mit Kovariaten auf Zeilenebene
wiederholtes Outcome mit Kovariaten auf Subjektebene
wiederholtes Outcome mit zeitvariierenden Kovariaten
Complete-Case-Analyse ohne Imputation
multiple Imputation auf Zeilenebene
multiple Imputation auf Subjektebene
subjektweite Imputation mit wiederholten Y als Hilfsvariablen
```

Unterstützte Modellfamilien:

```text
gaussian
bernoulli
poisson
negbinomial
beta
ordinal
categorical
```

Modellfamilien und Links werden in `00_config.R` gesetzt.

---

## 2. Struktur der Pipeline

Das Repository besteht aus wenigen vom Benutzer bearbeiteten Dateien und einer Reihe nummerierter Pipeline-Skripte. In den meisten Projekten müssen nur diese Dateien bearbeitet werden:

```text
00_config.R
00_variable_dictionary.csv
```

Alle anderen Skripte sollten normalerweise als Pipeline-Code behandelt werden.

### Pipeline-Skripte

`run_all.R` ruft diese Skripte in Reihenfolge auf:

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

Optionale Skripte:

```text
09_check_mo_parameter_columns.R      optional; prüft extrahierte mo()-Parameter-Spalten
10_publication_mo_results.R         optional; erzeugt abgeleitete mo()-Odds-Ratio-Zusammenfassungen
11_check_imputation_stability.R      optional; prüft Stabilität posteriorer Zusammenfassungen bei steigendem m
```

`09` und `10` sind nur für monotone ordinale Effekte mit `brms::mo()` erforderlich. `11` ist hilfreich, wenn Modellfitting teuer ist und die gewählte Anzahl an Imputationen begründet werden soll.

---

## 3. Schnellstart

Dieses Beispiel nutzt `datasets::airquality`.

```bash
git clone https://github.com/fabregithub/generic-mi-brms-pipeline.git
cd generic-mi-brms-pipeline
```

Dann ausführen:

```bash
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv
Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

Ausgaben werden geschrieben nach:

```text
objects/
fits/
results/
results/publication/
```

---

## 4. Anpassung an private Studiendaten

Empfohlener Ablauf:

1. Einen sauberen Analysedatensatz vorbereiten und als `.rds` speichern.
2. `00_variable_dictionary.csv` bearbeiten.
3. `00_config.R` bearbeiten.
4. Validierung ausführen.
5. Kleinen Schnelltest ausführen.
6. Moderaten Paralleltest ausführen.
7. Vollständigen Produktionslauf ausführen.
8. Publikationsausgaben rendern und prüfen.

Beispiel für Eingabedaten:

```text
data/my_analysis_data.rds
```

In R speichern:

```r
saveRDS(my_data, "data/my_analysis_data.rds")
```

Dann in `00_config.R` referenzieren:

```r
data = list(
  raw_data_file = "data/my_analysis_data.rds",
  ...
)
```

### Mechanismus fehlender Daten prüfen

Vor der Imputation sollte geprüft werden, warum jede Variable fehlt. `miceRanger` ist geeignet, wenn MCAR oder MAR nach Konditionierung auf beobachtete Variablen plausibel ist.

Nicht einfach als gewöhnliche `NA` übergeben:

```text
links-zensierte Messungen
Werte unterhalb einer Nachweisgrenze
rechts- oder intervallzensierte Messungen
strukturell fehlende Werte
nicht zutreffende Antworten
fehlende Werte durch Studiendesign
bekannte MNAR-Variablen
```

Diese Fälle sollten vor der Pipeline geeignet verarbeitet oder in Sensitivitätsanalysen behandelt werden.

### Imputation deaktivieren

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

Fehlende Werte prüfen:

```r
source("00_config.R")
d <- readRDS(paths$raw_data)
colSums(is.na(d))
```

### Anzahl der Imputationen wählen

Es gibt keinen universellen Wert für `m`. Für teure Analysen empfiehlt sich ein stufenweiser Ablauf:

```text
20 -> 40 -> 60 -> 80 -> 100
```

Nach Schritt 6:

```bash
Rscript 11_check_imputation_stability.R
```

---

## 5. Variablenwörterbuch

`00_variable_dictionary.csv` ist die maschinenlesbare Beschreibung der Analysevariablen.

Erwartete Spalten:

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

Wichtige Felder:

- `var`: exakter Variablenname im Datensatz.
- `label`: lesbare Bezeichnung für Tabellen und Berichte.
- `impute_target`: ob die Variable bei fehlenden Werten imputiert werden soll.
- `use_in_model`: ob die Variable im finalen `brms`-Modell verwendet wird.
- `use_as_auxiliary`: ob die Variable nur für die Imputation genutzt wird.

`role` beschreibt den analytischen Zweck: `outcome`, `binary_outcome`, `exposure`, `covariate`, `auxiliary`, `id`, `time`, `cluster`, `strata`.

`type` beschreibt den Datentyp: `continuous`, `integer`, `binary`, `categorical`, `ordinal`, `date`, `id`.

`timing` beschreibt den Messzeitpunkt: `single`, `baseline`, `repeated`, `time_varying`, `derived`, `id`.

`scale = z` erzeugt eine standardisierte Version `var_z`. Verwenden Sie `z` nicht für Outcomes, IDs, binäre/kategoriale Variablen oder bereits standardisierte Variablen.

---

## 6. Parallelisierung und Performance-Tuning

Zentrale Einstellungen:

| Einstellung | Verwendung | Bedeutung |
|---|---|---|
| `impute_workers` | Schritt 3 | Anzahl paralleler `miceRanger`-Worker |
| `num_impute_threads_per_worker` | Schritt 3 | Threads pro Imputations-Worker |
| `fit_workers` | Schritt 4 | Imputierte Datensätze, die parallel gefittet werden |
| `cores_per_fit` | Schritt 4 | Ketten/Kerne pro `brms`-Fit |
| `summary_workers` | Schritt 6 | Worker für Posterior-Draw-Extraktion |
| `prediction_workers` | Schritt 7 | Worker für Posterior-Prädiktion |
| `future_globals_maxsize_gb` | Schritte 4, 6, 7 | Maximale Future-Globals-Größe |

Beispiel:

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

Starten Sie konservativ, besonders auf Laptops oder bei großen Daten.

---

## 7. Logging, Monitoring, Neustart und Fehlerbehebung

Häufige Befehle:

```bash
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_stdout.log
Rscript fit_single_imputation.R 1
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

Nach einer Unterbrechung:

```bash
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Gültige vorhandene Fits werden übersprungen.

Monitoring-Dateien:

```text
pipeline_progress.log
pipeline_stdout.log
pipeline_heartbeat.txt
pipeline_success.flag
pipeline_error.flag
results/worker_logs/
```

Fortschritt beobachten:

```bash
tail -f pipeline_progress.log
```

Bei CmdStanR-Cache-Problemen:

```bash
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

---

## 8. Publikationsausgaben und Inferenzhinweise

Nach erfolgreichem Abschluss liegen Ausgaben in:

```text
results/publication/
```

Typische Ausgaben:

```text
results/publication/tables/main_effect_table_display.csv
results/publication/tables/main_effect_table_full.csv
results/publication/tables/diagnostics_summary.csv
results/publication/tables/analysis_metadata.csv
results/publication/figures/forest_plot_odds_ratios.png
results/publication/report/bayesian_mi_report_template.qmd
```

Bericht rendern:

```bash
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

Für Modelle mit `brms::mo()`:

```bash
Rscript 09_check_mo_parameter_columns.R
Rscript 10_publication_mo_results.R
quarto render results/publication/mo_effects/report/mo_effects_report.qmd
```

---

## 9. Beispiele und Tests

Enthaltene öffentliche Beispiele:

```text
examples/airquality_gaussian
examples/birthwt_logistic
examples/birthwt_spline_monotonic
```

Schnelltests:

```bash
bash test/test_all_examples_quick.sh
```

Paralleltests:

```bash
bash test/test_all_examples_parallel.sh
```

Vor dem Wechsel des Beispiels bereinigen:

```bash
rm -rf objects fits results
rm -f pipeline_error.flag pipeline_success.flag
rm -f pipeline_progress.log pipeline_heartbeat.txt pipeline_stdout.log run_all_stdout.log
```

---

## 10. Einrichtung der Rechenumgebung

Benötigt werden R, Systemwerkzeuge, R-Pakete, CmdStan/CmdStanR und Quarto.

Wichtige R-Pakete:

```r
install.packages(c(
  "tidyverse", "miceRanger", "brms", "posterior", "bayestestR",
  "future", "furrr", "doParallel", "foreach", "gt", "flextable",
  "officer", "forcats", "glue", "readr", "tibble", "dplyr",
  "stringr", "purrr", "rlang"
))
```

`cmdstanr` und CmdStan:

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

Quarto prüfen:

```bash
quarto --version
```

Für Reproduzierbarkeit kann `renv` verwendet werden:

```r
install.packages("renv")
renv::init()
renv::snapshot()
renv::restore()
```
