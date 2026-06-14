# 汎用 MICE + brms パイプラインテンプレート

言語: [English](../README.md) | [Français](README.fr.md) | [Español](README.es.md) | [Deutsch](README.de.md) | 日本語

> これは英語版 README の日本語要約です。内容に差異や曖昧さがある場合は、英語版を正とします。

これは、多重代入法（multiple imputation by chained equation; MICE）の補完後に、ベイズ回帰分析を行うためのRパイプラインテンプレートです。

対応内容：

- データ検証
- `miceRanger` によるMICE
- `brms` + `cmdstanr` によるベイズ回帰
- 補完データセットごとに1つのモデルを保存するチェックポイント方式
- 補完データセット間での並列モデルフィッティング
- 診断
- 事後分布の要約
- アウトカムが欠測している行に対する事後予測
- 論文・報告に使いやすい表、図、方法・設定メタデータ、レポートテンプレート

実用例ではR組み込みの公開データセット `datasets::airquality` を使うため、私的データなしでテンプレートをテスト・試行できます。

---

## 目次

1. [背景と目的](#1-背景と目的)
2. [パイプラインの構成](#2-パイプラインの構成)
3. [クイックスタート](#3-クイックスタート)
4. [実際の研究データへの適用](#4-私的な研究データへの適用)
5. [変数辞書](#5-変数辞書)
6. [並列化と性能調整](#6-並列化と性能調整)
7. [ログ、監視、再開、トラブルシューティング](#7-ログ監視再開トラブルシューティング)
8. [論文用出力と推論上の注意](#8-出版用出力と推論上の注意)
9. [実例とテスト](#9-例とテスト)
10. [計算環境のセットアップ](#10-計算環境のセットアップ)

---

## 1. 背景と目的

このパイプラインは、欠測共変量、反復アウトカム、大規模モデル、または慎重なチェックポイント保存が必要なモデルを含む応用ベイズ回帰分析を想定しています。MICE補完、補完データセットごとのベイズモデリング、診断、事後要約、事後予測、出版向け出力を組み合わせて一つのパイプラインにしています。

特に大規模データの場合、すべてのフィット済みモデルをメモリに保持するとメモリオーバーにより計算機がストップこともあります。そこで、再現性と再開しやすさを優先し、補完データセットごとに 1 つの `brms` モデルをフィットし、各 fit をただちに保存します。モデル毎にチェックポイントを記録し、再実行時には有効なチェックポイントファイルを再利用します。

### 注意点と限界

このリポジトリはワークフローの雛形であり、統計的判断の代替ではありません。科学的分析に使う前に、代入戦略、モデル式、事前分布、診断、事後要約が研究課題に適しているか確認してください。

多重代入ステップは、標準的な MICE 型の仮定が科学的に妥当な変数を対象とします。通常は、代入モデルに含まれる観測変数で条件づけた後に MCAR（missing completely at random）または MAR（missing at random）とみなせる場合を対象としています。このパイプラインは、MNAR（missing not at random）メカニズム、打ち切り、切断、検出限界（左打切データなど）、構造的欠測を自動的には扱いません。

一部のモデルは計算コストが高くなります。必ず本番実行の前に小規模のテストから始めてください。

### 重要な設計メモ

このテンプレートはモデルフィッティングに **`brm_multiple()` を使いません**。

代わりに、補完データセットごとに 1 つの `brms` モデルをフィットし、すぐに保存します。

```text
fits/fit_imp_001.rds
fits/fit_imp_002.rds
...
fits/fit_imp_100.rds
```

この方式は大規模データに対して用いることを想定しており、メモリオーバーを防止できるため、より安全です。

- メインの R セッションが全 fit をメモリに保持しない。
- 実行が止まっても完了済み fit が残る。
- 失敗または遅いモデルフィッティングを個別に再実行できる。
- 再実行時に有効な既存 fit をスキップできる。
- worker プロセス（多重代入やモデルフィッティングなど）がメインセッションへ返すのは小さなステータスオブジェクトだけである。

並列化は `future` / `furrr` を使って行い、動的スケジューリングを使います。

```
furrr::furrr_options(
  seed = TRUE,
  scheduling = Inf
)
```

### 対応する分析パターン

```
単一時点アウトカム、行レベル共変量
反復アウトカム、被験者レベル共変量
反復アウトカム、時間変動共変量
代入なしの complete-case 解析
行レベルMICE
被験者レベルMICE
反復 Y を補助変数として用いる被験者ワイド代入
```

対応するモデルファミリー：

```
gaussian
bernoulli
poisson
negbinomial
beta
ordinal
categorical
```

ファミリーとリンクは `00_config.R` で設定します。

---

## 2. パイプラインの構成

リポジトリは、ユーザーが編集する少数のファイルと、番号付きのパイプラインスクリプトで構成されています。多くのプロジェクトでは、通常次の 2 つのファイルだけを編集します。

```
00_config.R
00_variable_dictionary.csv
```

その他のスクリプトは通常、パイプライン本体のコードとして使い、ユーザーが特に編集することはありません。

### パイプラインスクリプト

メインパイプラインは `run_all.R` で実行され、次のスクリプトを順に呼び出します。

```
01_validate_config.R
02_prepare_data.R
03_impute.R
04_fit_models.R
05_diagnostics.R
06_posterior_summary.R
07_posterior_prediction.R
08_publication_results.R
```

追加の任意スクリプト：

```
09_check_mo_parameter_columns.R      任意；抽出された mo() パラメータ列を確認
10_publication_mo_results.R         任意；mo() の派生 odds ratio 要約を作成
11_check_imputation_stability.R      任意；m を増やしたときの事後要約の安定性を確認
```

`09` と `10` は、`brms::mo()` による単調順序効果の論文向け要約が必要な場合だけ使います。`11` は、モデルフィッティングが高コストで、代入数`m`の妥当性を示したい場合に有用です。

---

## 3. クイックスタート

この例では `datasets::airquality` を使います。

まず、gitサイトから、最新のパイプラインを入手してください。任意のフォルダで以下を実行します。

```
git clone https://github.com/fabregithub/generic-mi-brms-pipeline.git
cd generic-mi-brms-pipeline
```

続いて以下を実行します。

```
cp examples/airquality_gaussian/00_config_airquality_gaussian.R 00_config.R
cp examples/airquality_gaussian/00_variable_dictionary_airquality_gaussian.csv 00_variable_dictionary.csv
Rscript examples/airquality_gaussian/00_create_airquality_example_data.R
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_airquality_stdout.log
```

出力先：

```
objects/
fits/
results/
results/publication/
```

---

## 4. 実際の研究データへの適用

推奨ワークフロー：

1. クリーンな解析用データセットを 1 つ準備し、`.rds` として保存する。
2. `00_variable_dictionary.csv` を編集する。
3. `00_config.R` を編集する。
4. 検証を実行する。
5. 小さなクイックテストを実行する。
6. 控えめな並列テストを実行する。
7. 本番解析を実行する。
8. 論文用出力をレンダリングして確認する。

入力データ例：

```
data/my_analysis_data.rds
```

R で保存：

```
saveRDS(my_data, "data/my_analysis_data.rds")
```

`00_config.R` で指定：

```
data = list(
  raw_data_file = "data/my_analysis_data.rds",
  ...
)
```

### 欠測メカニズムの確認

欠測補完前に、各変数がなぜ欠測しているのかを確認してください。このパイプラインは `miceRanger` を用いた MICE 型多重代入を使います。観測変数で条件づけた後に MCAR または MAR とみなすことが妥当な場合に適しています。

次のものを通常の `NA` としてそのまま渡さないでください。

```
左打ち切り測定値
検出限界未満の値
右打ち切りまたは区間打ち切り測定値
構造的欠測
該当なし回答
研究デザインに起因する欠測
既知の MNAR 変数
```

これらはパイプライン前に適切に処理するか、感度分析で扱う必要があります。

### 補完を使わない場合

```
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

欠測の確認：

```
source("00_config.R")
d <- readRDS(paths$raw_data)
colSums(is.na(d))
```

### 補完データセット数の選び方

補完データセット数`m` に万能の値はありません。計算コストが高い解析では、段階的に増やすのが安全です。

```
20 -> 40 -> 60 -> 80 -> 100
```

4 モデルを並列でフィッティングする場合は次も便利です。

```
24 -> 40 -> 60 -> 80 -> 100
```

ステップ 6 の後に実行：

```
Rscript 11_check_imputation_stability.R
```

出力先：

```
results/publication/mi_stability/
```

---

## 5. 変数辞書

`00_variable_dictionary.csv` は解析変数の機械可読な記述です。

期待される列：

```
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

主なフィールド：

- `var`: データセット内の正確な変数名。
- `label`: 表やレポート用の人間可読ラベル。
- `impute_target`: 欠測時に補完対象とするか。
- `use_in_model`: 最終 `brms` モデルに含めるか。
- `use_as_auxiliary`: 最終モデルには含めず、代入用の補助変数として使うか。

`role` は解析上の役割を表します。例：`outcome`, `binary_outcome`, `exposure`, `covariate`, `auxiliary`, `id`, `time`, `cluster`, `strata`。

`type` は統計的な型を表します。例：`continuous`, `integer`, `binary`, `categorical`, `ordinal`, `date`, `id`。

`timing` は測定タイミングを表します。例：`single`, `baseline`, `repeated`, `time_varying`, `derived`, `id`。

`scale = z` は標準化済み変数 `var_z` を作成します。アウトカム、ID、二値/カテゴリ変数、すでに標準化済みの変数には使わないでください。

`reference` は二値・カテゴリ・順序変数の基準カテゴリです。実データに存在する値と一致している必要があります。

---

## 6. 並列化と性能調整

主な設定：

| 設定 | 使用箇所 | 意味 |
|---|---|---|
| `impute_workers` | ステップ 3 | 並列 `miceRanger` worker 数 |
| `num_impute_threads_per_worker` | ステップ 3 | 補完 worker ごとのスレッド数 |
| `fit_workers` | ステップ 4 | 並列にフィッティングする補完データセット数 |
| `cores_per_fit` | ステップ 4 | `brms` fit ごとの chain/core 数 |
| `summary_workers` | ステップ 6 | 事後 draw 抽出用 worker 数 |
| `prediction_workers` | ステップ 7 | 事後予測用 worker 数 |
| `future_globals_maxsize_gb` | ステップ 4, 6, 7 | future globals の最大サイズ |

例：

```
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

ノート PC や大きなデータでは、保守的な設定から始めてください。

---

## 7. ログ、監視、再開、トラブルシューティング

通常使うコマンド：

```
Rscript 01_validate_config.R
Rscript run_all.R 2>&1 | tee run_all_stdout.log
Rscript fit_single_imputation.R 1
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

実行が中断された場合は、再度実行します。

```
Rscript run_all.R 2>&1 | tee run_all_stdout.log
```

既存の有効な fit ファイルはスキップされます。

ログとステータスファイル：

```
pipeline_progress.log
pipeline_stdout.log
pipeline_heartbeat.txt
pipeline_success.flag
pipeline_error.flag
results/worker_logs/
```

進捗監視：

```
tail -f pipeline_progress.log
```

CmdStanR キャッシュの問題が疑われる場合：

```
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

---

## 8. 出版用出力と推論上の注意

正常終了後、論文用出力は次に作成されます。

```
results/publication/
```

代表的な出力：

```
results/publication/tables/main_effect_table_display.csv
results/publication/tables/main_effect_table_full.csv
results/publication/tables/diagnostics_summary.csv
results/publication/tables/analysis_metadata.csv
results/publication/figures/forest_plot_odds_ratios.png
results/publication/report/bayesian_mi_report_template.qmd
```

レポートのレンダリング：

```
quarto render results/publication/report/bayesian_mi_report_template.qmd
```

`brms::mo()` を使うモデルでは：

```
Rscript 09_check_mo_parameter_columns.R
Rscript 10_publication_mo_results.R
quarto render results/publication/mo_effects/report/mo_effects_report.qmd
```

---

## 9. 実例とテスト

含まれる公開データ実例：

```
examples/airquality_gaussian
examples/birthwt_logistic
examples/birthwt_spline_monotonic
```

クイックテスト：

```
bash test/test_all_examples_quick.sh
```

並列テスト：

```
bash test/test_all_examples_parallel.sh
```

例を切り替える前のクリーンアップ：

```
rm -rf objects fits results
rm -f pipeline_error.flag pipeline_success.flag
rm -f pipeline_progress.log pipeline_heartbeat.txt pipeline_stdout.log run_all_stdout.log
```

または：

```
bash 99_cleanall.sh
```

---

## 10. 計算環境のセットアップ

R、システムツール、必要な R パッケージ、CmdStan/CmdStanR、レポートレンダリング用の Quarto が必要です。

主要 R パッケージ：

```
install.packages(c(
  "tidyverse", "miceRanger", "brms", "posterior", "bayestestR",
  "future", "furrr", "doParallel", "foreach", "gt", "flextable",
  "officer", "forcats", "glue", "readr", "tibble", "dplyr",
  "stringr", "purrr", "rlang"
))
```

`cmdstanr` と CmdStan：

```
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

library(cmdstanr)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_path()
```

Quarto の確認：

```
quarto --version
```

長期的な再現性のためには `renv` の利用も検討してください。

```
install.packages("renv")
renv::init()
renv::snapshot()
renv::restore()
```
