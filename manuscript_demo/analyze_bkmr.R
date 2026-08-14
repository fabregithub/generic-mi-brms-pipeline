#!/usr/bin/env Rscript
# =============================================================================
# Mixture (BKMR) demonstration: the overall PFAS-mixture effect on Kawasaki
# disease, under the four censored-exposure methods, scored against the oracle.
#
# BKMR (bkmr::kmbayes) is MCMC and SLOW. Iwata ran it once on a ~10% subsample.
# This script is therefore BOUNDED and the mixture engine is PLUGGABLE:
#   * MIXTURE_ENGINE = "bkmr"  (default) — the published method.
#   * a faster variant (e.g. an accelerated / approximate BKMR such as "A-BKMR",
#     or a GP/gWQS surrogate) can be dropped in by defining fit_mixture_<engine>().
# Bound the cost with N_SUB (subsample), ITER (MCMC iters), M (imputations),
# REPS. Start tiny to verify, then scale.
#
# Env: OUTCOME(kd) N(25000) N_SUB(2000) M(10) ITER(2000) REPS(1)
#      LIMIT_SCALE(1) CORRELATION(real) METHODS(...) MIXTURE_ENGINE(bkmr)
# Writes results/mixture_<engine>_<tag>.csv
# =============================================================================
suppressWarnings(suppressMessages(library(leftcens)))
.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "." })
source(file.path(.here, "R", "make_demo_data.R"))
source(file.path(.here, "R", "analysis_methods.R"))

geti <- function(k,d) { v<-Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }
gets <- function(k,d) { v<-Sys.getenv(k); if (nzchar(v)) v else d }
OUTCOME<-gets("OUTCOME","kd"); N<-as.integer(geti("N",25000)); N_SUB<-as.integer(geti("N_SUB",2000))
M<-as.integer(geti("M",10)); ITER<-as.integer(geti("ITER",2000)); REPS<-as.integer(geti("REPS",1))
SCALE<-geti("LIMIT_SCALE",1); CORR<-gets("CORRELATION","real")
ENGINE<-gets("MIXTURE_ENGINE","bkmr")
METHODS<-strsplit(gets("METHODS","oracle,substitution,impute_noY,impute_congenial"),",")[[1]]
tag<-sprintf("%s_scale-%.2g_cor-%s_n%d", OUTCOME, SCALE, CORR, N)

# ---- pluggable mixture engine ----------------------------------------------
# Contract: fit_mixture(d, outcome) -> list(overall = estimate of the overall
# mixture effect (all exposures q25 -> q75, others at median), se = its SE).
# The exposures are the 7 log2 PFAS; adjust for the .Z_ADJUST confounders.
fit_mixture <- function(d, outcome, engine = ENGINE, n_sub = N_SUB, iter = ITER) {
  if (nrow(d) > n_sub) d <- d[sample(nrow(d), n_sub), , drop = FALSE]   # Iwata-style subsample
  Z <- as.matrix(sapply(.PFAS, function(j) log2(pmax(d[[j]], .Machine$double.eps))))
  colnames(Z) <- .PFAS
  zadj <- setdiff(intersect(.Z_ADJUST, names(d)), outcome)
  X <- stats::model.matrix(stats::as.formula(paste("~", paste(zadj, collapse = " + "))), d)[, -1, drop = FALSE]
  y <- d[[outcome]]
  fn <- get0(paste0("fit_mixture_", engine), ifnotfound = NULL)
  if (is.null(fn)) stop("No mixture engine '", engine, "'. Define fit_mixture_", engine, "().", call. = FALSE)
  fn(y, Z, X, iter)
}

# --- engine: bkmr (published) -----------------------------------------------
fit_mixture_bkmr <- function(y, Z, X, iter) {
  if (!requireNamespace("bkmr", quietly = TRUE))
    stop("Package 'bkmr' required for MIXTURE_ENGINE=bkmr. install.packages('bkmr').", call. = FALSE)
  fit <- bkmr::kmbayes(y = y, Z = Z, X = X, iter = iter, family = "binomial",
                       verbose = FALSE, varsel = TRUE)
  burn <- seq_len(floor(iter / 2))
  ov <- bkmr::OverallRiskSummaries(fit, qs = seq(0.25, 0.75, by = 0.5),
                                   q.fixed = 0.5, sel = setdiff(seq_len(iter), burn))
  r <- ov[ov$quantile == 0.75, ]
  list(overall = r$est, se = r$sd, pips = bkmr::ExtractPIPs(fit))
}

# --- engine: A-BKMR / accelerated (HOOK — plug the faster package here) ------
# fit_mixture_abkmr <- function(y, Z, X, iter) { ... same contract ... }

# ---- run --------------------------------------------------------------------
cat(sprintf("Mixture study | %s | engine=%s N=%d subsample=%d m=%d iter=%d reps=%d\n",
            tag, ENGINE, N, N_SUB, M, ITER, REPS))
.rubin_scalar <- function(est, se) { m<-length(est); q<-mean(est); u<-mean(se^2)
  b<-if(m>1) stats::var(est) else 0; se2<-sqrt(u+(1+1/m)*b); c(est=q, se=se2) }

acc <- list()
for (r in seq_len(REPS)) {
  dat <- make_demo_data(n = N, seed = 2000 + r, limit_scale = SCALE, correlation = CORR)
  oref <- fit_mixture(impute_pfas(dat$data, dat$truth, "oracle", outcome = OUTCOME)[[1]], OUTCOME)$overall
  for (meth in METHODS) {
    comp <- impute_pfas(dat$data, dat$truth, meth, m = M, seed = 11 * r, outcome = OUTCOME)
    fits <- lapply(comp, function(d) tryCatch(fit_mixture(d, OUTCOME), error = function(e) NULL))
    fits <- Filter(Negate(is.null), fits)
    if (!length(fits)) next
    p <- .rubin_scalar(vapply(fits, `[[`, 1, "overall"), vapply(fits, `[[`, 1, "se"))
    acc[[length(acc)+1]] <- data.frame(rep = r, method = meth, overall = p["est"],
      se = p["se"], overall_oracle = oref, err_vs_oracle = p["est"] - oref, row.names = NULL)
    cat(sprintf("  rep %d %-18s overall=%+.3f (oracle %+.3f, err %+.3f)\n",
                r, meth, p["est"], oref, p["est"] - oref))
  }
}
raw <- do.call(rbind, acc)
dir.create(file.path(.here, "results"), showWarnings = FALSE)
utils::write.csv(raw, file.path(.here, "results", sprintf("mixture_%s_%s.csv", ENGINE, tag)), row.names = FALSE)
cat(sprintf("\nWrote results/mixture_%s_%s.csv\n", ENGINE, tag))
