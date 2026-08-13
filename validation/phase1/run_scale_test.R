#!/usr/bin/env Rscript
# =============================================================================
# Scale test for the X-block draw at production size (PLAN §5 scale question)
# -----------------------------------------------------------------------------
# The block-FCS X block is a single call to leftcens::impute_censored_conditional()
# per (completed dataset x outer sweep): a shash-margin fit + one interval-censored
# survreg on n x k, per imputation. This times that call at n = 80k across a range
# of predictor widths k, so Phase 4 can pick a predictor-set size that keeps the
# per-sweep cost sane. Reports margin-fit time, per-draw time, and an extrapolation
# to a full run of  (m_datasets x T_outer)  draws.
#
# Env: N (default 80000), KS ("5,25,50,100,200,300"), M_DATASETS (20), T_OUTER (5).
# =============================================================================
suppressWarnings(suppressMessages(library(leftcens)))

geti <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.integer(v) else d }
n   <- geti("N", 80000L)
ks  <- { v <- Sys.getenv("KS"); if (nzchar(v)) as.integer(strsplit(v, ",")[[1]]) else c(5L,25L,50L,100L,200L,300L) }
m_datasets <- geti("M_DATASETS", 20L)
t_outer    <- geti("T_OUTER", 5L)
nd <- 0.30

cat(sprintf("Scale test: n=%d  widths k=%s  (extrapolate to m_datasets=%d x T_outer=%d draws)\n\n",
            n, paste(ks, collapse=","), m_datasets, t_outer))

simulate_wide <- function(n, k, skew = 0.75, seed = 1) {
  set.seed(seed)
  # k predictors: first 3 correlated with the focal exposure, rest noise.
  Z <- matrix(rnorm(n * k), n, k)
  g <- rowSums(Z[, seq_len(min(3L, k)), drop = FALSE]) / sqrt(min(3L, k))
  x1 <- sinh(asinh(0.5 * g + sqrt(1 - 0.25) * rnorm(n)) + skew)   # skewed, predictor-linked
  y  <- 0.4 * x1 + rowSums(Z) / sqrt(k) + rnorm(n)                # an outcome to condition on
  lod <- as.numeric(quantile(x1, nd)); cens <- x1 < lod
  yy <- x1; yy[cens] <- NA
  xpred <- as.data.frame(cbind(Y = y, Z))
  names(xpred) <- c("Y", paste0("v", seq_len(k)))
  list(y = yy, x = xpred, lower = ifelse(cens, -Inf, NA_real_),
       upper = ifelse(cens, lod, NA_real_), n_cens = sum(cens))
}

res <- data.frame()
# Unit cost = one complete impute_censored_conditional(m = 1): margin fit + one
# interval-censored survreg on n x k + truncated draw. Each block-FCS sweep does
# exactly this on the updated data, so a full run costs ~ (m_datasets x T_outer)
# of these. Warm up once (JIT/allocation), then take the median of 3 timings.
for (k in ks) {
  d <- simulate_wide(n, k)
  call1 <- function() leftcens::impute_censored_conditional(
    d$y, d$x, d$lower, d$upper, m = 1L, margin = "shash")
  imp1 <- call1()                                  # warmup (also a sanity draw)
  ts <- replicate(3, system.time(call1())["elapsed"])
  t_call <- stats::median(ts)
  full <- t_call * m_datasets * t_outer            # extrapolated X-block total
  ok <- all(is.finite(imp1[, 1]))
  res <- rbind(res, data.frame(k = k, n_cens = d$n_cens,
               sec_per_sweep = round(t_call, 2),
               full_xblock_min = round(full / 60, 1), finite = ok))
  cat(sprintf("k=%3d | %5.2fs per sweep-draw | extrapolated X-block %5.1f min\n",
              k, t_call, full / 60))
}

cat("\n=== summary ===\n"); print(res, row.names = FALSE)
results_dir <- file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, file.path(results_dir, "scale_test.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s/scale_test.csv\n", results_dir))
