#!/usr/bin/env Rscript
# =============================================================================
# Track V16 -- DERIVE the registered prediction (../PLAN_pipeline_validation.md §8i)
# -----------------------------------------------------------------------------
# Run this BEFORE the track and keep it: the prediction is part of the record,
# and a prediction nobody can reproduce is not one. V15's law was fitted to two
# measured points; this one is derived from the mechanism, which is the version
# of the cycle worth having.
#
# THE BLOCKS ARE NOT INDEPENDENT, AND THIS MUST MIRROR THAT.
# A first version of this script drew each block once, in isolation. That is not
# what the engine does and not what the arms measure. 00_censored_exposure.R
# ALTERNATES:
#
#   for (t in 1..sweeps) {
#     Z block:  draw Z (and Y) given the CURRENT X, Y
#     X block:  draw X        given the CURRENT Z, Y
#   }
#
# So in `noYx` the Z block still sees Y, and the Y-informed Z-hat it produces is
# then a PREDICTOR in the X draw: Y leaks into the exposure block through Z. The
# same leak runs the other way in `noYz`. The single-block arms therefore measure
# a DIRECT effect plus whatever the other block carries back, not a clean
# marginal effect -- and the size of that leak is itself a prediction, visible
# here as the gap between noYx + noYz and noYboth.
#
# Everything below is correctly specified for the linear cells (survreg for the
# censored exposure, lm for the covariate), so the only remaining difference from
# the shipped engine is BART-versus-linear in the Z block -- not a difference in
# kind where the true conditional is linear. `nl_*` cells are NOT derived here:
# there the true Z1 conditional is non-linear, so a linear stand-in would
# conflate the Y omission with V12's shape effect. Only a sign is registered for
# those.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
for (f in c("dgp.R", "censoring.R", "metrics.R", "procedures.R"))
  source(file.path(.here, "R", f))
suppressMessages({ library(MASS); library(survival) })

# --- one draw from each block, given the current state of the other -----------
draw_x <- function(d, is_c, lod, with_y) {
  pr <- c("logX2", "logX3", "Z1", "Z2", if (with_y) "Y")
  d$event <- as.integer(!is_c)
  d$resp  <- ifelse(is_c, lod, d$logX1)
  sr <- survival::survreg(
    as.formula(sprintf("Surv(resp, event, type='left') ~ %s", paste(pr, collapse = "+"))),
    data = d, dist = "gaussian")
  mm <- model.matrix(as.formula(paste("~", paste(pr, collapse = "+"))), data = d)
  k  <- length(coef(sr))
  dr <- MASS::mvrnorm(1, c(coef(sr), log(sr$scale)), vcov(sr))     # proper draw
  mu <- as.vector(mm %*% dr[seq_len(k)]); s <- exp(dr[k + 1L])
  x <- d$logX1
  x[is_c] <- leftcens::rnorm_trunc(sum(is_c), mu[is_c], s, -Inf, lod)
  x
}
draw_z <- function(d, zmiss, with_y) {
  pr <- c("logX1", "logX2", "logX3", "Z2", if (with_y) "Y")
  fit <- lm(as.formula(sprintf("Z1 ~ %s", paste(pr, collapse = "+"))), data = d[!zmiss, ])
  dr  <- MASS::mvrnorm(1, coef(fit), vcov(fit)); s <- summary(fit)$sigma
  mm  <- model.matrix(as.formula(paste("~", paste(pr, collapse = "+"))), data = d)
  z <- d$Z1
  z[zmiss] <- as.vector(mm[zmiss, , drop = FALSE] %*% dr) + rnorm(sum(zmiss), 0, s)
  z
}

# --- the engine's alternation, with Y switchable per block --------------------
block_fcs <- function(d0, is_c, lod, zmiss, y_in_x, y_in_z, sweeps = 3L) {
  d <- d0
  # Start the exposures at the LOD and the covariates at their observed mean, so
  # the first Z draw cannot see the true censored values.
  d$logX1[is_c] <- lod
  d$Z1[zmiss]   <- mean(d0$Z1[!zmiss])
  for (t in seq_len(sweeps)) {
    d$Z1    <- draw_z(d, zmiss, with_y = y_in_z)        # Z given current X, (Y)
    d$logX1 <- draw_x(d, is_c, lod, with_y = y_in_x)    # X given current Z, (Y)
  }
  d
}

ARMS <- list(bartMI  = c(x = TRUE,  z = TRUE),
             noYx    = c(x = FALSE, z = TRUE),
             noYz    = c(x = TRUE,  z = FALSE),
             noYboth = c(x = FALSE, z = FALSE))

set.seed(20260903)
tr <- make_truth(p = 3L, erf_form = "additive")
NR <- 300L; M <- 20L; N <- 800L; ND <- 0.40; ZM <- 0.40
acc <- setNames(lapply(names(ARMS), function(z) numeric(0)), names(ARMS))

for (r in seq_len(NR)) {
  d0 <- simulate_complete(n = N, truth = tr)$data
  cn <- inject_left_censoring(d0, nd_frac = ND, censor_which = 1L)
  is_c <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
  zmiss <- runif(N) < ZM
  for (nm in names(ARMS)) {
    a <- ARMS[[nm]]
    es <- vs <- numeric(M)
    for (i in seq_len(M)) {
      d <- block_fcs(d0, is_c, lod, zmiss, y_in_x = a[["x"]], y_in_z = a[["z"]])
      e <- fit_lm_estimand(d, tr); es[i] <- e[["est"]]; vs[i] <- e[["se"]]^2
    }
    p <- rubin_pool(es, vs)
    acc[[nm]] <- c(acc[[nm]], p[["est"]])
  }
}

tv <- tr$estimand_true
cat("\n  DERIVED PREDICTION -- with the engine's alternation\n")
cat(sprintf("  (n=%d, m=%d, %d reps, 3 sweeps, nd=%.2f, z_miss=%.2f, truth=%.2f)\n\n",
            N, M, NR, ND, ZM, tv))
cat(sprintf("  %-9s %9s %11s %8s\n", "arm", "mean est", "vs bartMI", "mcse"))
ref <- acc$bartMI; sh <- list()
for (nm in names(ARMS)) {
  d <- (acc[[nm]] - ref) / tv * 100                      # paired, per replication
  sh[[nm]] <- d
  cat(sprintf("  %-9s %9.4f %10.2f %8.2f\n", nm, mean(acc[[nm]]), mean(d),
              sd(d) / sqrt(length(d))))
}
gap <- mean(sh$noYboth) - (mean(sh$noYx) + mean(sh$noYz))
gse <- sd(sh$noYboth - sh$noYx - sh$noYz) / sqrt(NR)
cat(sprintf("\n  THE LEAK: noYboth - (noYx + noYz) = %+.2f pp  (se %.2f)\n", gap, gse))
cat("  Non-zero means the blocks are not separable: Y reaching one block through\n")
cat("  the other is doing measurable work, and each single-block arm understates\n")
cat("  its own block by that much.\n")
