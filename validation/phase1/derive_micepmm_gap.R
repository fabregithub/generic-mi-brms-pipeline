# Item 3, last piece: why is micePmm biased (+2.4 to +3.9%, V19) when V21's
# z21_pmm -- the same PMM draw on the same correct formula -- is +/-0.6%?
#
# Both registered candidates are eliminated:
#   Y as an imputation target -- V19's cells have y_frac = 0, and
#     make_row_level_imputation_spec() filters targets to variables WITH
#     missingness, so Y is never a target there.
#   mice's predictor matrix -- measured identical to the correct set.
#
# Top surviving candidate: mice runs maxit = 5 FCS iterations over {Z1, Z2}
# JOINTLY per outer sweep; z21_pmm does inner_iter = 1. `inner_iter` is already
# a parameter, so this costs one run.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","exact_fork.R","robustness.R")) source(file.path("R",f))})
`%||%` <- function(a,b) if (is.null(a)) b else a
pool <- function(imps, fo) {
  e <- sapply(imps, function(w) coef(lm(fo, data = w))[["logX1"]])
  v <- sapply(imps, function(w) diag(vcov(lm(fo, data = w)))[["logX1"]])
  c(est = mean(e), se = sqrt(mean(v) + (1 + 1/length(e)) * var(e)))
}
sc <- Filter(function(s) s$name == "zr_fork", v2_scenarios())[[1]]
tr <- make_truth(p = 3L, erf_form = sc$erf_form, z_role = sc$z_role)
fo <- dgp_formula(tr); tv <- tr$estimand_true
R  <- as.integer(Sys.getenv("R", "200"))
res <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
  set.seed(31000 + r)
  b <- v2_make_bundle(sc, tr)
  out <- list()
  for (ii in c(1L, 3L, 5L)) {
    set.seed(32000 + r)
    p <- tryCatch(pool(ef21_impute_datasets(b$censored, tr, m = 20L, z_src = "pmm",
                                            inner_iter = ii), fo),
                  error = function(e) c(est = NA_real_, se = NA_real_))
    out[[paste0("e", ii)]] <- p[["est"]]; out[[paste0("s", ii)]] <- p[["se"]]
  }
  as.data.frame(out)
}, mc.cores = 12))
cat(sprintf("zr_fork, n = %d, %d reps, m = 20, z_src = pmm\n\n", sc$n, R))
cat(sprintf("%-14s %9s %9s %8s\n", "inner_iter", "bias%", "bias/SE", "mcse%"))
for (ii in c(1L, 3L, 5L)) {
  e <- res[[paste0("e", ii)]]; s <- res[[paste0("s", ii)]]
  ok <- is.finite(e)
  cat(sprintf("%-14d %+9.2f %+9.2f %8.2f\n", ii,
              100*(mean(e[ok])-tv)/tv, (mean(e[ok])-tv)/mean(s[ok]),
              100*(sd(e[ok])/sqrt(sum(ok)))/tv))
}
cat("\nV19 micePmm on this cell: +2.4 to +3.9%   |   V21 z21_pmm: +/-0.6%\n")
