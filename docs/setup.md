# Computing environment setup

Installing R, CmdStan, Quarto and the required packages, on macOS, Windows and Linux.

> Do this once per machine. If the pipeline already runs, you do not need this page.

*Part of the [Generic MICE + brms Pipeline Template](../README.md) documentation.*

---

Before running the pipeline, prepare the R environment and CmdStan toolchain.

The exact setup depends on the operating system and computing environment. The sections below cover macOS, Windows and Linux. If you are using a managed workstation, shared server or high-performance computing cluster, some system tools may need to be installed by an administrator or loaded through environment modules.

---

## 11.1 macOS setup

Before running the pipeline, prepare the R environment and CmdStan toolchain. These instructions are for macOS. If running Windows, set up the Windows toolchain and R environment appropriately.

## 1. Install system tools on macOS

Run in Terminal to install the Apple command line tools:

```bash
# Bash command block
xcode-select --install
```

Run in Terminal to confirm that `make` and a C++ compiler are available:

```bash
# Bash command block
make --version
clang++ --version
```

If these commands fail, restart the Terminal and try again.

## 2. Install required R packages

Run in R or RStudio to install the packages required for this pipeline:

```r
install.packages(c(
  "tidyverse",
  "miceRanger",
  "brms",
  "posterior",
  "bayestestR",
  "future",
  "furrr",
  "doParallel",
  "foreach",
  "gt",
  "flextable",
  "officer",
  "forcats",
  "glue",
  "readr",
  "tibble",
  "dplyr",
  "stringr",
  "purrr",
  "rlang"
))
```

If you plan to use the Cox proportional hazards family (`family = "cox"`), the `survival` package is also required at model-fitting time. It ships with R as a recommended package and is normally already available, but you can install it explicitly if needed:

```r
install.packages("survival")
```

The default imputation path uses **BART** (`z_imputer = "bart"`, since v1.5.0), which needs the `dbarts` package:

```r
install.packages("dbarts")
```

Without it the pipeline still runs — it warns and falls back to the v1.4.0 imputer — but intervals will be slightly too narrow where much of the data is imputed. See [`validation/phase1/FINDINGS_v7.md`](../validation/phase1/FINDINGS_v7.md).

If you plan to use the censored-exposure block-FCS strategy (`strategy = "censored_exposure_block_fcs"`, for a left-/interval-censored focal exposure), the `leftcens` package is also required at imputation time. It is not on CRAN — install it from GitHub, pinned to a release tag:

```r
install.packages("remotes")

remotes::install_github("fabregithub/leftcens@v0.9.0")
```

Version **0.9.0 or later** is required: it exports `impute_censored_conditional()`, the outcome-aware, skew-aware conditional draw this strategy is built on. Earlier versions (≤ 0.8.0) do not, and `01_validate_config.R` will stop with an explanatory error.

Run in R or RStudio to install `cmdstanr` from the Stan R-universe repository:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
```

## 3. Install CmdStan

Run in R or RStudio:

```r
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

cmdstanr::install_cmdstan()
```

This may take several minutes because CmdStan is compiled locally.

Run in R or RStudio to check the CmdStan path:

```r
cmdstanr::cmdstan_path()
```

You should see something like:

```text
/Users/yourname/.cmdstan/cmdstan-2.xx.x
```

## 4. Verify CmdStan works

Run this small CmdStanR test in R or RStudio:

```r
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain(fix = TRUE)

stan_file <- file.path(
  cmdstanr::cmdstan_path(),
  "examples",
  "bernoulli",
  "bernoulli.stan"
)

mod <- cmdstanr::cmdstan_model(stan_file)

fit <- mod$sample(
  data = list(
    N = 10,
    y = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0)
  ),
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 10,
  iter_sampling = 10,
  refresh = 1
)

fit$summary()
```

If this runs successfully, CmdStan is ready.

## 5. Verify brms + cmdstanr works

Run this minimal `brms` test in R or RStudio:

```r
library(brms)
library(cmdstanr)

options(brms.backend = "cmdstanr")

fit_test <- brms::brm(
  mpg ~ wt,
  data = mtcars,
  family = gaussian(),
  chains = 1,
  iter = 500,
  warmup = 250,
  cores = 1,
  backend = "cmdstanr",
  refresh = 10,
  silent = 0
)

summary(fit_test)
```

If this succeeds, the local Bayesian modelling environment is ready for the pipeline.

## 6. Install Quarto for report rendering

The pipeline can generate a Quarto report template in `results/publication/report/bayesian_mi_report_template.qmd`.

To render this report to HTML or DOCX, install Quarto.

Run in Terminal on macOS with Homebrew:

```bash
# Bash command block
brew install --cask quarto
```

Alternatively, download and install Quarto from the Quarto website.

Run in Terminal to verify the installation:

```bash
# Bash command block
quarto --version
```

If this command prints a version number, Quarto is ready.

## 7. Optional: clear CmdStanR cache if fitting fails unexpectedly

If a simple model fails with an error such as:

```text
Fitting failed. Unable to retrieve the metadata.
No chains finished successfully. Unable to retrieve the fit.
```

and the data look valid, run the following in Terminal to clear the CmdStanR cache:

```bash
# Bash command block
rm -rf ~/.cmdstanr-cache
mkdir -p ~/.cmdstanr-cache
```

Then re-run the validation or pipeline.

## 8. Recommended reproducibility option: renv

For long-term reproducibility, consider using `renv`. Run in R or RStudio:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

This creates a project-specific package lockfile so the same package versions can be restored later. Run in R or RStudio:

```r
renv::restore()
```

---

## 11.2 Windows setup

Install:

- R
- RStudio
- Rtools
- Quarto

Use the version of Rtools that matches your R version. After installing Rtools, open R or RStudio and run:

```r
Sys.which("make")
Sys.which("g++")
```

Both should return valid paths.

Install the required R packages and `cmdstanr` using the same R commands shown in the macOS section. Then run:

```r
library(cmdstanr)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_path()
```

If toolchain problems persist, restart RStudio after installing Rtools. Then run the same CmdStan and `brms` verification examples shown in the macOS section.

To verify Quarto from Command Prompt or PowerShell:

```bash
# Command Prompt / PowerShell
quarto --version
```

---

## 11.3 Linux setup

Linux distributions differ, but the required system tools are generally:

- `make`
- `g++`
- `tar`
- `gzip`

For Ubuntu/Debian:

```bash
# Bash command block
sudo apt update
sudo apt install -y build-essential gfortran make
```

For Fedora:

```bash
# Bash command block
sudo dnf groupinstall -y "Development Tools"
sudo dnf install -y gcc-c++ gcc-gfortran make
```

For Red Hat Enterprise Linux or compatible systems, the exact commands may depend on the system configuration and repositories.

On a shared server or cluster, system tools may be provided through environment modules, for example:

```bash
# Bash command block
module avail gcc
module load gcc
```

Check with your system administrator or cluster documentation.

Install the required R packages and `cmdstanr` using the same R commands shown in the macOS section. Then run:

```r
library(cmdstanr)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_path()
```

Verify `brms` + `cmdstanr` using the same test shown in the macOS section.

Install Quarto using the installer appropriate for your distribution. On Ubuntu/Debian, for example:

```bash
# Bash command block
sudo dpkg -i quarto-*-linux-amd64.deb
```

Then check:

```bash
# Bash command block
quarto --version
```

On a shared server, Quarto may be available as a module or may need to be installed in a user directory.

---

## Dependencies

The pipeline requires both R packages and the Quarto command-line tool for report rendering.

Core R packages:

```r
install.packages(c(
  "tidyverse",
  "miceRanger",
  "brms",
  "posterior",
  "bayestestR",
  "future",
  "furrr",
  "doParallel",
  "foreach",
  "gt",
  "flextable",
  "officer",
  "forcats",
  "glue"
))
```

Quarto is also required if you want to render the generated `.qmd` report. Run in Terminal:

```bash
# Bash command block
# macOS with Homebrew
brew install --cask quarto

# Check installation
quarto --version
```

Install `cmdstanr` and CmdStan. Run in R or RStudio:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

cmdstanr::install_cmdstan()
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
```

If CmdStan is already installed, confirm the path in R or RStudio:

```r
cmdstanr::cmdstan_path()
```

---

---

*[← Back to the main README](../README.md)*
