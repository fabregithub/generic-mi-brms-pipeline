#!/usr/bin/env Rscript
# =============================================================================
# Scenario grid for the single-pollutant demonstration (cheap estimand): runs the
# four methods across the robustness factors and writes one master results table.
#
#   factors: outcome (kd, asthma) x limit_scale (censoring) x correlation
#   methods: oracle, substitution, impute_noY, impute_congenial   (scored vs oracle)
#
# Start SMALL (N ~ 2000, REPS ~ 5) to verify the grid; scale N/REPS for the study.
# The mixture (BKMR) grid is separate (analyze_bkmr.R) because it is far heavier.
#
# Env: N(2000) M(10) REPS(5) OUTCOMES(kd,asthma) SCALES(1,2) CORRS(low,real,high)
# Writes results/scenario_grid.csv (raw) and prints a headline summary.
# =============================================================================
suppressWarnings(suppressMessages(library(leftcens)))
.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "." })
source(file.path(.here, "R", "make_demo_data.R"))
source(file.path(.here, "R", "analysis_methods.R"))

geti <- function(k,d){v<-Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d}
gets <- function(k,d){v<-Sys.getenv(k); if (nzchar(v)) v else d}
N<-as.integer(geti("N",2000)); M<-as.integer(geti("M",10)); REPS<-as.integer(geti("REPS",5))
outcomes <- strsplit(gets("OUTCOMES","kd,asthma"),",")[[1]]
scales   <- as.numeric(strsplit(gets("SCALES","1,2"),",")[[1]])
corrs    <- strsplit(gets("CORRS","low,real,high"),",")[[1]]
methods  <- c("oracle","substitution","impute_noY","impute_congenial")

grid <- expand.grid(outcome=outcomes, scale=scales, corr=corrs, stringsAsFactors=FALSE)
cat(sprintf("Scenario grid: %d cells x %d reps x %d methods (N=%d, m=%d)\n",
            nrow(grid), REPS, length(methods), N, M))

acc <- list()
for (g in seq_len(nrow(grid))) {
  oc <- grid$outcome[g]; sc <- grid$scale[g]; cr <- grid$corr[g]
  cat(sprintf("[cell %d/%d] outcome=%s scale=%.2g cor=%s ", g, nrow(grid), oc, sc, cr))
  for (r in seq_len(REPS)) {
    dat <- make_demo_data(n=N, seed=3000+r, limit_scale=sc, correlation=cr)
    oref <- setNames(fit_single_pollutant(impute_pfas(dat$data, dat$truth, "oracle", outcome=oc), dat$truth, outcome=oc)$beta,
                     .PFAS)
    for (meth in methods) {
      fit <- fit_single_pollutant(impute_pfas(dat$data, dat$truth, meth, m=M, seed=13*r, outcome=oc), dat$truth, outcome=oc)
      fit$rep<-r; fit$method<-meth; fit$outcome<-oc; fit$scale<-sc; fit$corr<-cr
      fit$err_vs_oracle <- fit$beta - oref[fit$pfas]
      fit$covers_oracle <- fit$lo <= oref[fit$pfas] & oref[fit$pfas] <= fit$hi
      acc[[length(acc)+1]] <- fit
    }
    cat(".")
  }
  cat("\n")
}
raw <- do.call(rbind, acc)
dir.create(file.path(.here,"results"), showWarnings=FALSE)
utils::write.csv(raw, file.path(.here,"results","scenario_grid.csv"), row.names=FALSE)

# headline: mean |error vs oracle| by outcome x scale x corr x method
h <- aggregate(abs(err_vs_oracle) ~ outcome + scale + corr + method, raw, mean)
names(h)[5] <- "mean_abs_err"
h <- h[order(h$outcome, h$scale, h$corr, match(h$method, methods)), ]
cat("\n=== mean |error vs oracle| across PFAS (lower = better) ===\n")
key <- ""
for (i in seq_len(nrow(h))) {
  k <- sprintf("%s | scale %.2g | cor %s", h$outcome[i], h$scale[i], h$corr[i])
  if (k != key) { cat(sprintf("\n%s\n", k)); key <- k }
  cat(sprintf("  %-18s %.3f\n", h$method[i], h$mean_abs_err[i]))
}
cat(sprintf("\nWrote results/scenario_grid.csv (%d rows)\n", nrow(raw)))
