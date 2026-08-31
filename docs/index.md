# Guides

Detailed documentation for the pipeline. The [main README](../README.md) is the entry
point — what this is, whether it fits your analysis, and how to start. These guides hold
the detail.

## Getting started

| Guide | What it covers |
|---|---|
| [Computing environment setup](setup.md) | Installing R, CmdStan, Quarto and the required packages on macOS, Windows and Linux |
| [Examples and tests](examples.md) | Four worked analyses on public data, and the automated test scripts |

## Configuring an analysis

Two files normally need editing for a new analysis — `00_config.R` and
`00_variable_dictionary.csv`.

| Guide | What it covers |
|---|---|
| [`00_config.R` reference](config-reference.md) | Option-by-option reference for the main configuration file |
| [Variable dictionary](variable-dictionary.md) | `role`, `type`, `timing`, `scale`, `reference`, `use_as_auxiliary`, and when config may override the dictionary |
| [Choosing `m`, and running at scale](choosing-m.md) | How many imputations, the automatic m-increment loop, and the three run modes |

## Data situations

| Guide | What it covers |
|---|---|
| [Censored (below-detection-limit) exposures](censored-exposures.md) | Imputing exposures reported only as "below the LOD", what has been validated, **and the restriction against mixture/BKMR analyses** |
| [Repeated measures](repeated-measures.md) | Several rows per subject: subject-wide versus row-wise imputation |

## Running and reporting

| Guide | What it covers |
|---|---|
| [Parallelisation and performance](performance.md) | Sizing `impute_workers`, `fit_workers` and `chains`; memory-versus-cores trade-offs |
| [Running, monitoring and recovering](operations.md) | Restarting an interrupted run, logging, the CmdStan cache, debugging a fit |
| [Publication outputs and manuscript writing](reporting.md) | What the pipeline produces for a paper, and Methods/Results text templates |

---

## Before you rely on a claim

Everything the documentation asserts about bias, coverage and scope is backed by a
Monte-Carlo validation study. What has been tested, what may legitimately be claimed, and
the limits of each claim are indexed in
[`validation/README.md`](../validation/README.md).

The one restriction worth knowing before you start: the censored-exposure strategy is
validated for **additive / per-analyte** exposure–response functions and
**must not be used for mixture / BKMR analyses** — see
[Censored exposures](censored-exposures.md).

---

*[← Back to the main README](../README.md)*
