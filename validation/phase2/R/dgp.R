# =============================================================================
# Phase 2 harness -- data-generating process (birth-weight mixture)   DESIGN §1
# -----------------------------------------------------------------------------
# The known-truth cohort for the mixture-HBEB validation. 7 prenatal PFAS ->
# birth weight (-> LBW), grounded in `manuscript_demo/PROVENANCE.md`: correlated
# log-exposures from a Z + GFR determinant model, a KNOWN non-linear mixture
# surface (main effects + curvature + interactions), and a birth-weight marginal
# tail-matched to LBW 8.13%. Everything is fixed by construction, so the true
# mixture benchmark is computable (via benchmark.R) and every recovered estimate
# is scored against it.
#
# All parameters trace to PROVENANCE.md (tagged [grounded]/[by design] there).
# Base R + MASS (mvrnorm), matching the Phase-1 harness.
# =============================================================================

suppressWarnings(suppressMessages(library(MASS)))

# source benchmark.R for the true-benchmark computation (one inversion, applied
# to the true h -- DESIGN §0). Guarded so dgp.R is usable without it.
.dgp_here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.dgp_here) || !nzchar(.dgp_here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .dgp_here <- if (length(f)) dirname(normalizePath(f)) else "."
}
if (file.exists(file.path(.dgp_here, "benchmark.R")))
  source(file.path(.dgp_here, "benchmark.R"))

# ---- analyte specifications (PROVENANCE §1) ---------------------------------

#' Per-analyte medians (ng/mL) and log-scale spread. `sdlog` for PFOS/PFOA are
#' derived from Lai's IQRs; the others are plausible serum-PFAS values.
pfas_specs <- function() {
  data.frame(
    analyte = c("pfos","pfoa","pfna","pfda","pfunda","pfhxs","pftrda"),
    group   = c("PFSA","PFCA","PFCA","PFCA","PFCA","PFSA","PFCA"),
    median  = c(2.90, 1.60, 1.40, 0.49, 1.10, 0.33, 0.15),
    sdlog   = c(0.48, 0.61, 0.60, 0.65, 0.65, 0.55, 0.70),
    stringsAsFactors = FALSE)
}

#' Target log-scale correlation matrix (PROVENANCE §2). Grounded Schober pairs
#' from Lai; the remaining entries are plausible fills. `strength` scales the
#' off-diagonals for the low/real/high scenario axis (DESIGN §7).
pfas_correlation <- function(strength = c("real","low","high")) {
  strength <- match.arg(strength)
  an <- pfas_specs()$analyte; p <- length(an)
  R <- diag(p); dimnames(R) <- list(an, an)
  set <- function(a, b, r) { R[a,b] <<- R[b,a] <<- r }
  # grounded (Lai, Schober-classified)
  set("pfoa","pfna",0.74); set("pfna","pfda",0.72); set("pfda","pfunda",0.81)
  set("pfunda","pftrda",0.83); set("pfos","pfda",0.79); set("pfos","pfunda",0.70)
  # plausible fills [by design]
  set("pfos","pfoa",0.55); set("pfos","pfna",0.60); set("pfos","pftrda",0.55)
  set("pfoa","pfda",0.60); set("pfoa","pfunda",0.45); set("pfoa","pftrda",0.45)
  set("pfna","pfunda",0.65); set("pfna","pftrda",0.60); set("pfda","pftrda",0.70)
  # PFHxS weakly correlated with everything (~0.25)
  for (a in setdiff(an,"pfhxs")) set("pfhxs", a, 0.25)

  off <- R - diag(p)
  R <- switch(strength,
    real = R,
    low  = diag(p) + off * 0.40,                         # shrink toward independence
    high = diag(p) + pmin(off * 1.30 + 0.20*(off>0), 0.92) * (off > 0))
  make_pd(R)
}

#' Nearest positive-definite correlation matrix (eigenvalue clip + cov2cor).
make_pd <- function(R) {
  e <- eigen(R, symmetric = TRUE)
  R2 <- e$vectors %*% diag(pmax(e$values, 1e-6)) %*% t(e$vectors)
  R2 <- stats::cov2cor(R2); dimnames(R2) <- dimnames(R); R2
}

# ---- the truth (surface + confounding + benchmark spec) ---------------------

#' Assemble the ground truth. `erf_form="mixture"` includes curvature + two
#' interactions (the surface that stresses linear congenial imputation);
#' `"additive"` keeps main effects only.
make_truth <- function(erf_form = c("mixture","additive"),
                       corr = c("real","low","high")) {
  erf_form <- match.arg(erf_form); corr <- match.arg(corr)
  sp <- pfas_specs(); an <- sp$analyte
  med <- stats::setNames(sp$median, an); sdlog <- stats::setNames(sp$sdlog, an)

  # main effects, grams birth weight per DOUBLING (PROVENANCE §5)
  beta <- stats::setNames(c(pfos=-16.4, pfoa=-20.0, pfna=-10.0, pfda=-6.0,
                            pfunda=-6.0, pfhxs=-8.0, pftrda=-4.0)[an], an)
  mix <- erf_form == "mixture"
  list(
    analytes = an, median = med, sdlog = sdlog,
    R = pfas_correlation(corr), corr = corr, erf_form = erf_form,
    beta   = beta,
    b_quad = if (mix) c(pfos = 2.5) else numeric(0),               # u_PFOS^2 curvature
    inter  = if (mix) list(c("pfos","pfoa",-4.0), c("pfos","pfhxs",-3.0))
             else list(),                                          # synergies
    # covariate effects on birth weight (grams); renal factor R is the confounder
    gamma  = list(smoking_current = -180, smoking_past = -20,
                  ga_per_week = 180, ga_mean = 38.9, male = 100,
                  height_sd = 45, bmi_sd = 25, parity_per = 25),
    gfr    = list(c_R = 0.25, gamma_R = -95),   # R: +logPFAS, -BW (Verner); gamma_R tunes observed/causal ratio
    # Z -> X determinant model (log scale)
    det    = list(region_sd = 0.15, b_parity = -0.15, b_age = 0.08, b_ap = 0.03),
    # birth-weight marginal (effective SD tail-matched to LBW 8.13%)
    bw_mean = 3028.0, bw_sd = 378.2, lbw_cut = 2500, N = 720000,
    n_regions = 15L,
    estimand = "mixture_benchmark")
}

# ---- the surface h(X) (concentration matrix -> grams) -----------------------

#' Centred per-doubling exposure u_j = log2(X_j) - log2(median_j).
.u_of <- function(Xmat, med) sweep(log2(Xmat), 2, log2(med[colnames(Xmat)]), `-`)

#' The known surface predictor: matches PROVENANCE §5. Returns grams.
h_predict <- function(truth) {
  med <- truth$median; beta <- truth$beta; bq <- truth$b_quad; inter <- truth$inter
  function(Xmat) {
    Xmat <- Xmat[, names(beta), drop = FALSE]
    U <- .u_of(Xmat, med)
    h <- as.vector(U %*% beta)
    if (length(bq)) for (a in names(bq)) h <- h + bq[[a]] * U[, a]^2
    for (it in inter) h <- h + as.numeric(it[[3]]) * U[, it[[1]]] * U[, it[[2]]]
    h
  }
}

#' Wrap the true surface as a benchmark.R surface object (needs benchmark.R).
true_surface <- function(truth) make_surface(h_predict(truth), truth$analytes)

# ---- simulate one complete (uncensored, fully observed) cohort --------------

#' @return list(data, truth). `data` has concentration columns (ng/mL), the
#'   covariates, the latent renal factor R, birth weight `bw`, and `lbw`.
simulate_complete <- function(n, truth, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  an <- truth$analytes; p <- length(an); med <- truth$median; g <- truth$gamma

  # --- covariates Z + latent renal factor R ---------------------------------
  age    <- stats::rnorm(n, 31, 5)
  parity <- sample(0:3, n, replace = TRUE, prob = c(.39,.39,.16,.06))
  bmi_z  <- stats::rnorm(n)                                  # standardised
  height_z <- stats::rnorm(n)
  smoking <- sample(c("never","past","current"), n, TRUE, prob = c(.58,.36,.045))
  male    <- stats::rbinom(n, 1, 0.51)
  ga      <- stats::rnorm(n, g$ga_mean, 1.4)                 # gestational age (wk)
  region  <- sample.int(truth$n_regions, n, replace = TRUE)
  R       <- stats::rnorm(n)                                 # renal / plasma-volume

  # --- Z -> X determinant model on the log scale ----------------------------
  d <- truth$det
  region_eff <- matrix(stats::rnorm(truth$n_regions * p, 0, d$region_sd),
                       truth$n_regions, p)                   # fixed region loadings
  age_s <- (age - 31) / 5
  eff <- region_eff[region, ] + d$b_parity * parity + d$b_age * age_s +
    d$b_ap * parity * age_s + truth$gfr$c_R * R              # renal -> +log PFAS
  eff <- sweep(eff, 2, colMeans(eff), `-`)                   # centre: geomean = median
  det_lin <- sweep(eff, 2, log(med[an]), `+`)
  # correlated residual with per-analyte sd = sdlog
  Sigma <- diag(truth$sdlog[an]) %*% truth$R[an, an] %*% diag(truth$sdlog[an])
  eps   <- MASS::mvrnorm(n, rep(0, p), Sigma)
  logX  <- det_lin + eps
  X     <- exp(logX); colnames(X) <- an

  # --- birth weight ---------------------------------------------------------
  h  <- h_predict(truth)(X)
  zc <- ifelse(smoking == "current", g$smoking_current,
               ifelse(smoking == "past", g$smoking_past, 0)) +
        g$ga_per_week * (ga - g$ga_mean) + g$male * male +
        g$height_sd * height_z + g$bmi_sd * bmi_z + g$parity_per * parity
  rc <- truth$gfr$gamma_R * R
  lin <- (h - mean(h)) + (zc - mean(zc)) + (rc - mean(rc))
  # residual SD absorbs the remainder so the marginal SD = effective bw_sd
  var_expl <- stats::var(lin)
  sig_e <- sqrt(max(truth$bw_sd^2 - var_expl, (0.2 * truth$bw_sd)^2))
  if (truth$bw_sd^2 - var_expl <= 0)
    warning("explained BW variance exceeds target; SD floored -- reduce effects")
  bw <- truth$bw_mean + lin + stats::rnorm(n, 0, sig_e)
  lbw <- as.integer(bw < truth$lbw_cut)

  data <- data.frame(X, age = age, parity = parity, bmi_z = bmi_z,
                     height_z = height_z, smoking = smoking, male = male,
                     ga = ga, region = factor(region), R = R,
                     bw = bw, lbw = lbw, check.names = FALSE)
  attr(data, "analytes") <- an
  list(data = data, truth = truth)
}

# ---- the true mixture benchmark (needs benchmark.R) -------------------------

#' Compute the ground-truth mixture HBEB from the KNOWN surface, by simulating a
#' large complete population and inverting it (benchmark.R). One inversion,
#' applied to the true h (DESIGN §0). Directions and targets default to the
#' DESIGN §2 set.
#'
#' @return data.frame: direction, target_label, target_cases, theta, and the
#'   per-analyte benchmark concentrations.
true_benchmark <- function(truth, rel_targets = c(0.10, 0.25, 0.50),
                           directions = NULL, n_pop = 2e5, bg_p = 0.05,
                           seed = 20260814L) {
  stopifnot(exists("benchmark_mixture"))
  pop <- simulate_complete(n_pop, truth, seed = seed)$data
  Xmat <- as.matrix(pop[, truth$analytes]); surf <- true_surface(truth)
  tr <- default_translation(truth$bw_mean, truth$bw_sd, truth$lbw_cut, truth$N)
  if (is.null(directions))
    directions <- list(prop = dir_proportional(truth$analytes),
                       top  = dir_analyte(truth$analytes, "pfos"),
                       pc1  = dir_pc1(Xmat))
  curr <- lapply(directions, function(dv)
    residual_cases(Xmat, surf, 0, dvec = dv, bg_p = bg_p, tr = tr))
  rows <- list()
  for (dn in names(directions)) for (r in rel_targets) {
    b <- benchmark_mixture(Xmat, surf, target_cases = curr[[dn]] * (1 - r),
                           dvec = directions[[dn]], bg_p = bg_p, tr = tr)
    cc <- if (isTRUE(b$ok) && length(b$conc) == length(truth$analytes)) b$conc
          else stats::setNames(rep(NA_real_, length(truth$analytes)), truth$analytes)
    rows[[length(rows)+1]] <- data.frame(
      direction = dn, target_label = sprintf("%.0f%%", r*100),
      target_cases = curr[[dn]] * (1 - r), theta = b$theta,
      as.data.frame(as.list(cc)), check.names = FALSE)
  }
  do.call(rbind, rows)
}

# ---- smoke test -------------------------------------------------------------

#' Sanity-check the DGP against PROVENANCE: medians, sdlog, correlations,
#' realised birth-weight marginal / LBW, the GFR confounding ratio, and a true
#' benchmark. Not a pass/fail gate -- prints for eyeballing + tuning.
dgp_smoketest <- function(n = 40000, seed = 1) {
  truth <- make_truth("mixture", "real")
  sim <- simulate_complete(n, truth, seed = seed); d <- sim$data
  an <- truth$analytes; X <- as.matrix(d[, an])
  cat("=== marginals (target median | realised median | target sdlog | realised) ===\n")
  for (a in an) cat(sprintf("  %-7s %5.2f | %5.2f    %4.2f | %4.2f\n",
      a, truth$median[[a]], stats::median(X[,a]),
      truth$sdlog[[a]], stats::sd(log(X[,a]))))
  cat(sprintf("\nbirth weight: mean %.0f (target %.0f), sd %.0f (target %.0f), LBW %.2f%% (target 8.13%%)\n",
      mean(d$bw), truth$bw_mean, stats::sd(d$bw), truth$bw_sd, 100*mean(d$lbw)))
  cat("\nrealised inter-analyte correlation (log scale):\n")
  print(round(stats::cor(log(X)), 2))
  # GFR confounding: observed vs causal-adjusted PFOS slope (per doubling)
  u_pfos <- log2(X[,"pfos"]) - log2(truth$median[["pfos"]])
  obs <- stats::coef(stats::lm(d$bw ~ u_pfos))[2]
  adj <- stats::coef(stats::lm(d$bw ~ u_pfos + d$R))[2]
  cat(sprintf("\nPFOS slope (g/doubling): observed %.1f, R-adjusted %.1f, ratio %.2f (target ~2x)\n",
      obs, adj, obs/adj))
  if (exists("benchmark_mixture")) {
    cat("\ntrue mixture benchmark (proportional direction):\n")
    tb <- true_benchmark(truth, n_pop = 5e4)
    print(tb[tb$direction == "prop", c("target_label","theta","pfos","pfoa","pfhxs")],
          row.names = FALSE)
  }
  invisible(sim)
}

if (sys.nframe() == 0 && identical(Sys.getenv("DGP_SMOKETEST"), "1")) dgp_smoketest()
