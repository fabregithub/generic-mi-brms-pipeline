#!/usr/bin/env Rscript
# =============================================================================
# DERIVE the registered predictions for V16 and V17, from the mechanism
# -----------------------------------------------------------------------------
# Supersedes predict_v16.R, which this generalises in two ways:
#
#   1. ALL FIVE COVARIATE ROLES, not just the original "precision" design.
#   2. BOTH covariates are imputed. predict_v16.R drew only Z1 and left Z2
#      complete, but the harness makes Z1 AND Z2 MCAR together -- and Z2's draw
#      also loses Y in the noYz arms, so omitting it understated the Z block's
#      effect. The V16 numbers below therefore supersede the ones derived from
#      that script; both are on the record.
#
# The routine mirrors what 00_censored_exposure.R actually does -- alternate
# {Z block, X block} for three sweeps -- with correctly-specified conditionals
# (survreg for the left-censored exposure, lm for Z1, logistic for binary Z2).
# Where the true conditional is linear-Gaussian that makes this an exact stand-in
# for the engine up to BART-versus-linear; where it is not, no number is
# registered. The `descendant` role (z_form = "nonlinear") is therefore NOT
# derived: a linear stand-in there would conflate the Y omission with V12's
# shape effect.
#
# Output is the paired shift of each arm against the both-blocks reference, in
# percentage points of the true estimand -- the same quantity the runners report.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
for (f in c("dgp.R", "censoring.R", "metrics.R", "procedures.R", "robustness.R"))
  source(file.path(.here, "R", f))
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressMessages({ library(MASS); library(survival) })

XCOLS <- paste0("logX", 1:3)

draw_x1 <- function(d, is_c, lod, with_y) {
  pr <- c("logX2", "logX3", "Z1", "Z2", if (with_y) "Y")
  d$event <- as.integer(!is_c); d$resp <- ifelse(is_c, lod, d$logX1)
  sr <- survival::survreg(
    as.formula(sprintf("Surv(resp, event, type='left') ~ %s", paste(pr, collapse = "+"))),
    data = d, dist = "gaussian")
  mm <- model.matrix(as.formula(paste("~", paste(pr, collapse = "+"))), data = d)
  k  <- length(coef(sr))
  dr <- MASS::mvrnorm(1, c(coef(sr), log(sr$scale)), vcov(sr))
  mu <- as.vector(mm %*% dr[seq_len(k)]); s <- exp(dr[k + 1L])
  x <- d$logX1
  x[is_c] <- leftcens::rnorm_trunc(sum(is_c), mu[is_c], s, -Inf, lod)
  x
}
draw_z1 <- function(d, miss, with_y) {
  if (!any(miss)) return(d$Z1)
  pr  <- c(XCOLS, "Z2", if (with_y) "Y")
  fit <- lm(as.formula(sprintf("Z1 ~ %s", paste(pr, collapse = "+"))), data = d[!miss, ])
  dr  <- MASS::mvrnorm(1, coef(fit), vcov(fit)); s <- summary(fit)$sigma
  mm  <- model.matrix(as.formula(paste("~", paste(pr, collapse = "+"))), data = d)
  z <- d$Z1
  z[miss] <- as.vector(mm[miss, , drop = FALSE] %*% dr) + rnorm(sum(miss), 0, s)
  z
}
draw_z2 <- function(d, miss, with_y) {
  if (!any(miss)) return(d$Z2)
  pr  <- c(XCOLS, "Z1", if (with_y) "Y")
  fit <- suppressWarnings(glm(as.formula(sprintf("Z2 ~ %s", paste(pr, collapse = "+"))),
                              data = d[!miss, ], family = binomial()))
  dr  <- tryCatch(MASS::mvrnorm(1, coef(fit), vcov(fit)), error = function(e) coef(fit))
  mm  <- model.matrix(as.formula(paste("~", paste(pr, collapse = "+"))), data = d)
  p   <- stats::plogis(as.vector(mm[miss, , drop = FALSE] %*% dr))
  z <- d$Z2
  z[miss] <- rbinom(sum(miss), 1, p)
  z
}

block_fcs <- function(d0, is_c, lod, m1, m2, y_in_x, y_in_z, sweeps = 3L) {
  d <- d0
  d$logX1[is_c] <- lod                          # start at the bound, not the truth
  if (any(m1)) d$Z1[m1] <- mean(d0$Z1[!m1])
  if (any(m2)) d$Z2[m2] <- round(mean(d0$Z2[!m2]))
  for (t in seq_len(sweeps)) {
    d$Z1    <- draw_z1(d, m1, with_y = y_in_z)
    d$Z2    <- draw_z2(d, m2, with_y = y_in_z)
    d$logX1 <- draw_x1(d, is_c, lod, with_y = y_in_x)
  }
  d
}

ARMS <- list(bartMI  = c(x = TRUE,  z = TRUE),
             noYx    = c(x = FALSE, z = TRUE),
             noYz    = c(x = TRUE,  z = FALSE),
             noYboth = c(x = FALSE, z = FALSE))
ROLES <- c("precision", "fork", "pipe", "collider", "mixed")

NR <- as.integer(Sys.getenv("NR", "200"))
M  <- as.integer(Sys.getenv("M",  "10"))
# MECH: which covariate mechanism to derive under. MCAR costs information;
# MAR here is driven off Y (and logX2) by inject_mar_covariates(), so dropping Y
# from the Z block also drops the variable the mechanism depends on -- a
# validity failure on top of the congeniality one. The two are expected to differ
# and each is registered separately.
MECH <- Sys.getenv("MECH", "mcar")
N  <- 800L; ND <- 0.40; ZM <- 0.40

cat(sprintf("\n  DERIVED PREDICTIONS -- engine alternation, both covariates imputed\n"))
cat(sprintf("  (n=%d, m=%d, %d reps, 3 sweeps, nd=%.2f, Z missing %.2f %s)\n",
            N, M, NR, ND, ZM, toupper(MECH)))

out <- list()
for (role in ROLES) {
  tr <- make_truth(p = 3L, erf_form = "additive", z_role = role)
  set.seed(20260903)
  acc <- setNames(lapply(names(ARMS), function(z) numeric(0)), names(ARMS))
  for (r in seq_len(NR)) {
    d0 <- simulate_complete(n = N, truth = tr)$data
    cn <- inject_left_censoring(d0, nd_frac = ND, censor_which = 1L)
    is_c <- cn$logX1_cens == "left"; lod <- cn$logX1_lod[1]
    if (identical(MECH, "mar")) {
      # Reuse the harness's own mechanism rather than a copy of it, so the
      # derivation and the run share one definition of MAR.
      mm <- inject_mar_covariates(d0, mar_frac = ZM, strength = 1.0)
      m1 <- is.na(mm$Z1); m2 <- is.na(mm$Z2)
    } else {
      m1 <- runif(N) < ZM; m2 <- runif(N) < ZM
    }
    for (nm in names(ARMS)) {
      a <- ARMS[[nm]]; es <- vs <- numeric(M)
      for (i in seq_len(M)) {
        d <- block_fcs(d0, is_c, lod, m1, m2, y_in_x = a[["x"]], y_in_z = a[["z"]])
        e <- fit_lm_estimand(d, tr); es[i] <- e[["est"]]; vs[i] <- e[["se"]]^2
      }
      p <- rubin_pool(es, vs)
      acc[[nm]] <- c(acc[[nm]], p[["est"]])
    }
  }
  tv <- tr$estimand_true; ref <- acc$bartMI
  sh <- lapply(acc, function(v) (v - ref) / tv * 100)
  gap <- mean(sh$noYboth) - (mean(sh$noYx) + mean(sh$noYz))
  gse <- sd(sh$noYboth - sh$noYx - sh$noYz) / sqrt(NR)
  cat(sprintf("\n  --- %s ---   reference bias %+.2f%%\n", role, (mean(ref) - tv) / tv * 100))
  cat(sprintf("  %-9s %11s %8s\n", "arm", "vs bartMI", "mcse"))
  for (nm in names(ARMS))
    cat(sprintf("  %-9s %10.2f %8.2f\n", nm, mean(sh[[nm]]),
                sd(sh[[nm]]) / sqrt(NR)))
  cat(sprintf("  %-9s %10.2f %8.2f   (noYboth - noYx - noYz)\n", "leak", gap, gse))
  out[[role]] <- list(shift = vapply(sh, mean, 0),
                      mcse = vapply(sh, function(v) sd(v) / sqrt(NR), 0),
                      leak = c(est = gap, se = gse),
                      ref_bias = (mean(ref) - tv) / tv * 100)
}
f <- file.path(.here, "results", sprintf("predictions_v16_v17_%s.rds", MECH))
saveRDS(out, f)
cat(sprintf("\n  wrote %s\n", f))
