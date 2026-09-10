# Decompose properBoot's 24x. Three imputers, one Z-block problem, identical
# data, m = 20. `forest` is ordinary miceRanger; `forest_boot` is the same
# engine with a bootstrap resample per imputation (the shipped default before
# v1.5.0); `bart` is dbarts. So:
#   forest      vs bart        = the ENGINE cost
#   forest_boot vs forest      = the BOOTSTRAP wrapper's cost
# Measured 2026-09-10: bart 4.0 s, forest 29.3 s, forest_boot 34.0 s.
# In situ (zr_fork, n = 800, 22 reps, one invocation per arm): properZ 11 s,
# micePmm 16 s, bartMI 17 s, properBoot 263 s.
# NOTE `secs` in the runner's raw csv is the whole TASK's elapsed time copied
# onto every arm row, so a multi-arm run CANNOT attribute per-arm cost -- that
# is why this went unexplained. One invocation per arm is the only way.
paths <- list(cache = tempdir(), logs = tempdir(), objects = tempdir())
suppressMessages(suppressWarnings(source("00_common_functions.R")))
log_msg <- function(...) invisible(NULL)
options(mi.quiet_imputer_warning = TRUE)
set.seed(7)
n <- 800
d <- data.frame(Y = rnorm(n), logX1 = rnorm(n), logX2 = rnorm(n),
                logX3 = rnorm(n), Z1 = rnorm(n), Z2 = rbinom(n, 1, .5))
d$Z1[sample(n, .4*n)] <- NA
d$Z2[sample(n, .4*n)] <- NA
spec <- list(vars = list(Z1 = c("Y","logX1","logX2","logX3","Z2"),
                         Z2 = c("Y","logX1","logX2","logX3","Z1")),
             m = 20L, maxiter = 5L, mean_match_k = 5, seed = 1L)
for (imp in c("bart", "forest", "forest_boot")) {
  as_ <- list(imputation = list(z_imputer = imp), parallel = list(impute_workers = 1L))
  t <- tryCatch(system.time(invisible(run_row_level_imputation(d, spec, as_)))[["elapsed"]],
                error = function(e) NA_real_)
  cat(sprintf("  %-12s %7.1f s\n", imp, t))
}
