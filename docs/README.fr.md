# Modèle générique de pipeline MICE + brms

Langue : [English](../README.md) | Français | [Español](README.es.md) | [Deutsch](README.de.md) | [日本語](README.ja.md)

> Ceci est une traduction de démarrage du README anglais. En cas de différence ou d’ambiguïté, la version anglaise fait autorité.

Ce dépôt fournit un modèle réutilisable de pipeline R pour des analyses de régression bayésienne avec imputation multiple optionnelle.

Il prend en charge :

- la validation des données ;
- l’imputation multiple optionnelle avec `miceRanger` ;
- la régression bayésienne avec `brms` + `cmdstanr` ;
- l’enregistrement d’un ajustement par jeu de données imputé ;
- l’ajustement parallèle des modèles entre imputations ;
- les diagnostics ;
- les résumés postérieurs ;
- la prédiction postérieure pour les lignes dont la variable de résultat est manquante ;
- des tableaux, figures, métadonnées de méthodes/réglages et modèles de rapports prêts pour publication.

L’exemple par défaut utilise le jeu de données public intégré `datasets::airquality`, ce qui permet de tester et démontrer le modèle sans données privées.

---

## Table des matières

1. [Contexte et objectif](#1-contexte-et-objectif)
2. [Structure du pipeline](#2-structure-du-pipeline)
3. [Démarrage rapide](#3-démarrage-rapide)
4. [Adapter le pipeline à des données d’étude privées](#4-adapter-le-pipeline-à-des-données-détude-privées)
5. [Dictionnaire des variables](#5-dictionnaire-des-variables)
6. [Parallélisation et réglage des performances](#6-parallélisation-et-réglage-des-performances)
7. [Journalisation, surveillance, reprise et dépannage](#7-journalisation-surveillance-reprise-et-dépannage)
8. [Sorties de publication et conseils d’inférence](#8-sorties-de-publication-et-conseils-dinférence)
9. [Guide de rédaction de manuscrit](#9-guide-de-rédaction-de-manuscrit)
10. [Exemples et tests](#10-exemples-et-tests)
11. [Configuration de l’environnement informatique](#11-configuration-de-lenvironnement-informatique)

---

## 1. Contexte et objectif

Ce pipeline est destiné aux analyses bayésiennes appliquées lorsque les données peuvent contenir des covariables manquantes, des résultats répétés, de grands modèles ou des modèles nécessitant une sauvegarde prudente. Il combine l’imputation multiple, la modélisation bayésienne avec un ajustement par jeu imputé, les diagnostics, les résumés postérieurs, la prédiction postérieure et des sorties orientées publication.

La conception privilégie la reproductibilité et la capacité de reprise plutôt que la conservation de tous les modèles ajustés en mémoire. Le pipeline ajuste donc un modèle `brms` par jeu de données imputé, sauvegarde immédiatement chaque ajustement et réutilise les fichiers de sauvegarde valides lors d’une nouvelle exécution.

### Précautions et limites

Ce dépôt est un canevas de flux de travail, pas un substitut au jugement statistique. Avant de l’utiliser pour une analyse scientifique, vérifiez que la stratégie d’imputation, la formule du modèle, les priors, les diagnostics et les résumés postérieurs sont adaptés à votre question d’étude.

Les étapes d’imputation multiple supposent qu’une hypothèse de type MICE est scientifiquement défendable, en général MCAR ou MAR après conditionnement sur les variables observées incluses dans le modèle d’imputation. Le pipeline ne traite pas automatiquement les mécanismes MNAR, la censure, la troncature, les limites de détection ou les valeurs structurellement manquantes.

Certains modèles peuvent être coûteux en calcul. Commencez toujours par un petit test rapide avant une exécution de production.

### Notes importantes de conception

Ce modèle n’utilise **pas** `brm_multiple()` pour l’ajustement des modèles.

Il ajuste plutôt un modèle `brms` par jeu imputé et sauvegarde chaque ajustement :

```text
fits/fit_imp_001.rds
fits/fit_imp_002.rds
...
fits/fit_imp_100.rds
```

Cette approche est plus sûre pour les grands jeux de données car :

- la session R principale ne garde pas tous les modèles en mémoire ;
- les ajustements terminés sont conservés si l’exécution s’arrête ;
- les imputations lentes ou en échec peuvent être relancées séparément ;
- les ajustements valides existants sont ignorés lors d’une relance ;
- les processus de travail ne renvoient que de petits objets de statut à la session principale.

La parallélisation se fait entre imputations avec `future` / `furrr`, avec planification dynamique :

```r
furrr::furrr_options(
  seed = TRUE,
  scheduling = Inf
)
```

### Schémas d’analyse pris en charge

Le modèle prend en charge :

```text
résultat à un seul temps, covariables au niveau ligne
résultat répété avec covariables au niveau sujet
résultat répété avec covariables variant dans le temps
analyse en cas complets sans imputation
imputation multiple au niveau ligne
imputation multiple au niveau sujet
imputation large par sujet utilisant les Y répétés comme variables auxiliaires
```

Familles de modèles prises en charge :

```text
gaussian
bernoulli
poisson
negbinomial
beta
ordinal
categorical
```

Les familles et liens sont définis dans `00_config.R`.

---

## 2. Structure du pipeline

Le dépôt est organisé autour de quelques fichiers modifiés par l’utilisateur et d’une séquence de scripts numérotés. Dans la plupart des projets, seuls ces fichiers doivent être modifiés :

```text
00_config.R
00_variable_dictionary.csv
```

Tous les autres scripts doivent généralement être considérés comme du code de pipeline.

### Scripts du pipeline

Le pipeline principal est lancé par `run_all.R`, qui appelle les scripts suivants dans l’ordre :

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

`run_all.R` exécute ensuite automatiquement deux scripts supplémentaires :

```text
09_check_mo_parameter_columns.R      s’exécute automatiquement seulement si la formule du modèle contient des termes mo()
10_publication_mo_results.R          s’exécute automatiquement seulement si la formule du modèle contient des termes mo()
11_check_imputation_stability.R      s’exécute toujours, et produit le rapport final de stabilité du nombre d’imputations
```

`run_all.R` détecte les termes `mo()` directement dans la formule du modèle ajusté ; vous n’avez pas besoin de tenir à jour une liste de variables. Les scripts `09`/`10` sont donc ignorés automatiquement pour les modèles gaussiens, logistiques, à spline ou à facteurs ordinaires.

Si `analysis_spec$mi_stability$auto_increment <- TRUE`, `run_all.R` augmente aussi automatiquement `m` par lots et arrête l’ajustement dès que les résumés postérieurs sont stables ; voir la section 4.

### Principaux fichiers à modifier

`00_config.R` définit la structure de l’analyse : variable de résultat, famille et lien, structure des données, stratégie d’imputation, formule `brms`, priors, réglages MCMC, parallélisation et garde-fous mémoire.

`00_variable_dictionary.csv` définit les libellés, rôles, types, catégories de référence, mise à l’échelle, cibles d’imputation et inclusion dans le modèle.

---

## 3. Démarrage rapide

Cet exemple utilise le jeu de données public `datasets::airquality`.

```bash
git clone https://github.com/fabregithub/generic-mi-brms-pipeline.git
cd generic-mi-brms-pipeline
```

Puis :

```bash
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv
Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

Les sorties sont écrites dans :

```text
objects/
fits/
results/
results/publication/
```

---

## 4. Adapter le pipeline à des données d’étude privées

Flux recommandé :

1. Préparer un jeu de données d’analyse propre et l’enregistrer en `.rds`.
2. Modifier `00_variable_dictionary.csv`.
3. Modifier `00_config.R`.
4. Lancer la validation.
5. Lancer un petit test rapide.
6. Lancer un test parallèle modeste.
7. Lancer l’analyse de production complète.
8. Générer et inspecter les sorties de publication.

### Préparation des données

Avant de lancer le pipeline, préparez un fichier comme :

```text
data/my_analysis_data.rds
```

En R :

```r
saveRDS(my_data, "data/my_analysis_data.rds")
```

Puis indiquez ce fichier dans `00_config.R` :

```r
data = list(
  raw_data_file = "data/my_analysis_data.rds",
  ...
)
```

### Vérifier le mécanisme de données manquantes

Avant l’imputation, examinez pourquoi chaque variable est manquante. Le pipeline utilise `miceRanger`, adapté lorsque MCAR ou MAR est raisonnable après conditionnement sur les variables observées.

Ne transmettez pas directement comme `NA` ordinaires :

```text
mesures censurées à gauche
valeurs sous une limite de détection
mesures censurées à droite ou par intervalle
valeurs structurellement manquantes
réponses non applicables
manquants dus au plan d’étude
variables MNAR connues
```

Ces cas doivent être traités avant le pipeline ou faire l’objet d’analyses de sensibilité.

### Désactiver l’imputation

Si aucune covariable nécessaire au modèle n’est manquante, ou si vous souhaitez une analyse sans imputation :

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

Vérifiez les valeurs manquantes :

```r
source("00_config.R")
d <- readRDS(paths$raw_data)
colSums(is.na(d))
```

### Choisir le nombre d’imputations

Il n’existe pas de valeur universelle pour `m`. Avec de grands jeux de données, les résultats se stabilisent souvent à un `m` plus petit qu’avec de petits jeux de données, alors que chaque ajustement prend plus de temps — ajuster un `m` fixe et élevé dès le départ peut donc gaspiller beaucoup de calcul.

**Boucle automatique (recommandée).** Réglez `analysis_spec$imputation$m` comme plafond, puis :

```r
analysis_spec$mi_stability <- list(
  auto_increment = TRUE,
  increment_size = NULL  # par défaut, égal à fit_workers
)
```

`run_all.R` ajuste alors `m` par lots (taille par défaut = `fit_workers`, pour occuper pleinement les workers parallèles), vérifie la stabilité après chaque lot, et arrête d’augmenter `m` dès que les seuils configurés sont atteints, ou lorsque le `m` configuré est atteint. Chaque lot est ensemencé de façon déterministe (`analysis_spec$imputation$seed`, par défaut basé sur `analysis_spec$model$seed`).

**Approche manuelle par étapes.** Pour piloter vous-même les lots, fixez `analysis_spec$imputation$allow_extend <- TRUE`, puis augmentez progressivement `m` :

```text
20 -> 40 -> 60 -> 80 -> 100
```

ou, si quatre modèles sont ajustés en parallèle :

```text
24 -> 40 -> 60 -> 80 -> 100
```

Avec `allow_extend = TRUE`, relancer le pipeline avec un `m` plus grand n’ajoute que les nouvelles imputations ; les imputations déjà existantes ne sont jamais modifiées.

Après l’étape 6, lancez :

```bash
Rscript 11_check_imputation_stability.R
```

Les sorties sont créées dans :

```text
results/publication/mi_stability/
```

---

## 5. Dictionnaire des variables

`00_variable_dictionary.csv` est la description lisible par machine des variables d’analyse.

Colonnes attendues :

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

### `role`

Valeurs recommandées :

| Valeur | Signification |
|---|---|
| `outcome` | Résultat continu ou général |
| `binary_outcome` | Résultat binaire pour modèle Bernoulli/logistique |
| `exposure` | Exposition principale ou prédicteur d’intérêt |
| `covariate` | Variable d’ajustement / confondeur / prédicteur |
| `auxiliary` | Utilisée uniquement pour l’imputation |
| `id` | Identifiant sujet, cluster ou ligne |
| `time` | Occasion de mesure, vague, visite ou suivi |
| `cluster` | Variable de regroupement |
| `strata` | Variable de stratification |

### `type`

Valeurs recommandées :

| Valeur | Signification |
|---|---|
| `continuous` | Variable numérique continue |
| `integer` | Variable de comptage |
| `binary` | Variable à deux niveaux |
| `categorical` | Catégorie non ordonnée |
| `ordinal` | Catégorie ordonnée |
| `date` | Date calendrier |
| `id` | Identifiant |

### `timing`

Valeurs recommandées : `single`, `baseline`, `repeated`, `time_varying`, `derived`, `id`.

### `scale`

Valeurs recommandées :

| Valeur | Résultat |
|---|---|
| `no` | Pas de transformation |
| `z` | Standardisation moyenne 0, écart-type 1 ; crée `var_z` |
| `centre` | Centrage |
| `log` | Transformation logarithmique |
| `custom` | Transformation définie par l’utilisateur |

N’utilisez pas `z` pour les résultats, identifiants, variables binaires/catégorielles ou variables déjà standardisées.

### `reference`

Définit la catégorie de référence pour les variables binaires, catégorielles ou ordinales. Elle doit correspondre à une valeur réellement présente dans les données.

### `use_as_auxiliary`

Indique qu’une variable est utilisée pour l’imputation mais exclue du modèle final. C’est utile pour des variables qui aident à prédire les valeurs manquantes sans faire partie du modèle scientifique.

---

## 6. Parallélisation et réglage des performances

Réglages principaux :

| Réglage | Utilisé dans | Signification |
|---|---|---|
| `impute_workers` | Étape 3 | Nombre de workers `miceRanger` |
| `num_impute_threads_per_worker` | Étape 3 | Threads par worker d’imputation |
| `fit_workers` | Étape 4 | Jeux imputés ajustés en parallèle |
| `cores_per_fit` | Étape 4 | Chaînes/cœurs par ajustement `brms` |
| `summary_workers` | Étape 6 | Workers pour l’extraction des tirages postérieurs |
| `prediction_workers` | Étape 7 | Workers pour la prédiction postérieure |
| `future_globals_maxsize_gb` | Étapes 4, 6, 7 | Taille maximale des globals future |

Exemple recommandé :

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

Commencez prudemment, surtout sur ordinateur portable ou avec de grandes données.

---

## 7. Journalisation, surveillance, reprise et dépannage

Commandes courantes :

```bash
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_stdout.log
Rscript fit_single_imputation.R 1
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

Le pipeline est checkpointé. Après interruption, relancez :

```bash
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

Les ajustements valides déjà terminés seront ignorés.

Fichiers de suivi :

```text
pipeline_progress.log
pipeline_stdout.log
pipeline_heartbeat.txt
pipeline_success.flag
pipeline_error.flag
results/worker_logs/
```

Surveiller la progression :

```bash
tail -f pipeline_progress.log
```

En cas de problème CmdStanR, essayez :

```bash
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

---

## 8. Sorties de publication et conseils d’inférence

Après une exécution réussie, les sorties sont dans :

```text
results/publication/
```

Sorties typiques :

```text
results/publication/tables/main_effect_table_display.csv
results/publication/tables/main_effect_table_full.csv
results/publication/tables/diagnostics_summary.csv
results/publication/tables/analysis_metadata.csv
results/publication/figures/forest_plot_odds_ratios.png
results/publication/report/bayesian_mi_report_template.qmd
```

Générer le rapport :

```bash
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

À l’étape 6, les tirages postérieurs ne sont pas simplement concaténés entre imputations : chaque imputation contribue avec un poids égal (1/m), et une correction de variance pour `m` fini (règle de Rubin) est appliquée lorsque la forme de la distribution postérieure poolée le permet. Voir les colonnes `m_imputations`, `variance_corrected`, etc. dans `results/parameter_summary.csv`.

Pour les modèles avec `brms::mo()`, `run_all.R` lance automatiquement :

```bash
Rscript 09_check_mo_parameter_columns.R
Rscript 10_publication_mo_results.R
quarto render results/publication/mo_effects/report/mo_effects_report.qmd
```

Ces scripts détectent les variables `mo()` directement dans la formule du modèle ; aucune liste de variables à modifier manuellement n’est nécessaire. Pour des libellés et noms de catégories personnalisés, ajoutez un bloc `analysis_spec$mo_effects` optionnel dans `00_config.R` (modèle suggéré affiché par `09_check_mo_parameter_columns.R`).

---

## 9. Guide de rédaction de manuscrit

Les sorties de publication sont conçues pour être citées directement dans un manuscrit plutôt que recalculées à la main.

| Paragraphe du manuscrit | Fichier source |
|---|---|
| Données manquantes / imputation | `results/publication/tables/analysis_metadata.csv` |
| Spécification du modèle, réglages MCMC | `results/publication/tables/analysis_metadata.csv` |
| Justification du nombre d’imputations adaptatif | `results/publication/mi_stability/tables/imputation_stability_stepwise_summary_display.csv` |
| Résultats des effets principaux | `results/publication/tables/main_effect_table_display.csv` |
| Résultats des effets monotones (`mo()`) | `results/publication/mo_effects/tables/mo_cumulative_or_table.csv` |

Exemple de texte « Méthodes » (à adapter avec les valeurs de `analysis_metadata.csv`) :

```text
Les données de covariables manquantes ont été traitées par imputation
multiple par équations chaînées avec forêts aléatoires (miceRanger),
générant m = XX jeux de données imputés. Un modèle bayésien distinct a
été ajusté à chaque jeu imputé avec brms et le moteur cmdstanr, et les
tirages postérieurs ont été poolés entre imputations avec un poids égal
par imputation.
```

Exemple de texte « Résultats » pour un effet principal (à partir de `main_effect_table_display.csv`) :

```text
[Variable] était associée à [résultat] : estimation = XX (IC à 95 % :
XX à XX ; probabilité postérieure de direction = XX).
```

Ce guide ne rédige pas pour vous l’interprétation scientifique, la pertinence clinique, les limites de l’étude ni les citations bibliographiques — cela reste du jugement humain.

---

## 10. Exemples et tests

Le dépôt contient trois exemples publics :

```text
examples/airquality_gaussian
examples/birthwt_logistic
examples/birthwt_spline_monotonic
```

Tests rapides :

```bash
bash test/test_all_examples_quick.sh
```

Tests parallèles :

```bash
bash test/test_all_examples_parallel.sh
```

Nettoyage avant de changer d’exemple :

```bash
rm -rf objects fits results
rm -f pipeline_error.flag pipeline_success.flag
rm -f pipeline_progress.log pipeline_heartbeat.txt pipeline_stdout.log run_all_stdout.log
```

ou :

```bash
bash 99_cleanall.sh
```

---

## 11. Configuration de l’environnement informatique

Le pipeline nécessite R, les packages R requis, CmdStan/CmdStanR et Quarto pour le rendu des rapports.

Installer les packages R principaux :

```r
install.packages(c(
  "tidyverse", "miceRanger", "brms", "posterior", "bayestestR",
  "future", "furrr", "doParallel", "foreach", "gt", "flextable",
  "officer", "forcats", "glue", "readr", "tibble", "dplyr",
  "stringr", "purrr", "rlang"
))
```

Installer `cmdstanr` :

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
```

Installer CmdStan :

```r
library(cmdstanr)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_path()
```

Vérifier Quarto :

```bash
quarto --version
```

Pour une reproductibilité à long terme, envisagez `renv` :

```r
install.packages("renv")
renv::init()
renv::snapshot()
renv::restore()
```
