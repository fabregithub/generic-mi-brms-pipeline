# Create the JECS PFAS -> Kawasaki disease demonstration input for the pipeline.
#
# SYNTHETIC data (grounded in Lai 2025 + Iwata 2024; see
# ../../manuscript_demo/PROVENANCE.md). Seven prenatal PFAS with three-tier
# interval censoring (leftcens columns <pfas>_lo / <pfas>_hi), covariates with
# some MAR, and a rare binary outcome (kd). The canonical generator lives in
# manuscript_demo/R/make_demo_data.R and is reused here.
#
# Env: N (default 3000 for a quick example; scale to ~25000 for a realistic run),
#      LIMIT_SCALE (1), CORRELATION (real).

dir.create("data", recursive = TRUE, showWarnings = FALSE)

# Locate the canonical DGP (relative to THIS script, with a repo-root fallback).
.args <- commandArgs(FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
.here <- if (length(.self)) dirname(normalizePath(.self)) else "."
gen <- file.path(.here, "..", "..", "manuscript_demo", "R", "make_demo_data.R")
if (!file.exists(gen)) gen <- "manuscript_demo/R/make_demo_data.R"
if (!file.exists(gen)) stop("Cannot find manuscript_demo/R/make_demo_data.R (needed to generate the example).")
source(gen)

N   <- if (nzchar(Sys.getenv("N"))) as.integer(Sys.getenv("N")) else 3000L
sc  <- if (nzchar(Sys.getenv("LIMIT_SCALE"))) as.numeric(Sys.getenv("LIMIT_SCALE")) else 1
cor <- if (nzchar(Sys.getenv("CORRELATION"))) Sys.getenv("CORRELATION") else "real"

out <- make_demo_data(n = N, seed = 2026L, limit_scale = sc, correlation = cor)
saveRDS(out$data, "data/jecs_pfas_kd_example.rds")

message("Saved: data/jecs_pfas_kd_example.rds")
message("Rows: ", nrow(out$data), ", KD cases: ", out$meta$n_kd,
        " (", sprintf("%.2f%%", 100 * out$truth$kd_prevalence), ")")
message("Censored (ND+DNQ) per PFAS: ",
        paste(sprintf("%s %.0f%%", out$truth$pfas, 100 * out$truth$censored_rate), collapse = ", "))
message("NOTE: rare outcome — use N ~ 25000 for a realistic number of cases.")
