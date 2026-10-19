# =============================================================================
# Phase 2 harness -- mixture HBEB benchmark: burden + inversion  (DESIGN.md §0,§2)
# -----------------------------------------------------------------------------
# THE shared inversion. One set of functions, applied to the TRUE surface h for
# the ground truth and to each procedure's estimated surface h-hat for the
# estimate (DESIGN §0). It is the mixture generalisation of the single-pollutant
# HBEB engine (`hbeb-framework-nakayamasf/R/burden_pipeline.R`):
#
#   single-pollutant :  Delta_mu = beta * ( E_obs[g(X)] - E_cf[g(X)] )
#   mixture (here)   :  Delta_mu = E[ h(X_obs) ] - E[ h(X_cf) ]
#
# Everything downstream -- the counterfactual mean birth weight `mu_cf = bw_mean
# - Delta_mu`, the normal lower-tail translation to attributable LBW cases, and
# the `uniroot` benchmark inversion -- is IDENTICAL. The single-pollutant engine
# is the special case h(X) = beta * g(X_1); `benchmark_selftest()` verifies that
# reduction against the shipped code, so this file is certified to invert the
# same way on real and simulated data.
#
# Base R only (matches burden_pipeline.R). Surface-agnostic: a "surface" is just
# a predictor contract (see make_surface), so A-BKMR / full bkmr / the known DGP
# h all plug in the same way (analysis.R supplies the wrappers later).
# =============================================================================

# ---- birth-weight -> LBW case translation (shared with HBEB) -----------------

#' Default continuous-to-LBW translation, tail-matched to JECS (HBEB CFG).
#'
#' `bw_sd` is the *effective* SD (378.2 g) that makes a normal reproduce the
#' observed 8.13% LBW at the 2500 g cut -- the singleton sample SD (416 g)
#' overpredicts LBW because the lower tail is thinner than normal. See
#' `manuscript_demo/PROVENANCE.md` §4 and burden_pipeline.R CFG.
default_translation <- function(bw_mean = 3028.0, bw_sd = 378.2,
                                lbw_cut = 2500, N = 720000) {
  list(bw_mean = bw_mean, bw_sd = bw_sd, lbw_cut = lbw_cut, N = N)
}

#' P(LBW) = P(birth weight < cut) under a normal with the given mean.
p_lbw <- function(mu, tr) stats::pnorm(tr$lbw_cut, mu, tr$bw_sd)

# ---- the surface contract ----------------------------------------------------

#' Wrap a mixture exposure-response surface as a predictor object.
#'
#' This is the single seam through which the truth and every procedure enter the
#' inversion (DESIGN §0). All exposures are CONCENTRATIONS (ng/mL); the surface
#' owns any internal rescaling (log2, centring at medians, etc.) and returns
#' birth-weight grams.
#'
#' Uncertainty is optional and, crucially, **predict-capable**: `predict_draw(Xmat,
#' s)` must evaluate the s-th posterior surface at *arbitrary* exposure points, so
#' the directional counterfactual can be re-predicted under each draw (BKMR/A-BKMR
#' expose exactly this). A fixed [n x S] matrix would NOT suffice, because the
#' benchmark shifts exposures to new points the draw must still score. The known
#' DGP surface is deterministic and leaves `predict_draw`/`n_draws` NULL/0.
#'
#' @param predict function(Xmat) -> numeric length n (grams of birth weight).
#' @param analytes character vector of the exposure column names, in order.
#' @param predict_draw optional function(Xmat, s) -> numeric length n: the s-th
#'   posterior draw's surface, evaluable at arbitrary Xmat.
#' @param n_draws number of posterior draws available (0 if deterministic).
make_surface <- function(predict, analytes, predict_draw = NULL, n_draws = 0L) {
  stopifnot(is.function(predict), length(analytes) >= 1L)
  if (!is.null(predict_draw)) stopifnot(is.function(predict_draw), n_draws >= 1L)
  structure(list(predict = predict, analytes = analytes,
                 predict_draw = predict_draw, n_draws = as.integer(n_draws)),
            class = "hbeb_surface")
}

#' Coerce an exposure input to a numeric concentration matrix in analyte order.
as_conc_matrix <- function(X, analytes) {
  M <- if (is.data.frame(X)) as.matrix(X[, analytes, drop = FALSE]) else X
  if (is.null(dim(M))) M <- matrix(M, ncol = length(analytes),
                                   dimnames = list(NULL, analytes))
  storage.mode(M) <- "double"
  M[, analytes, drop = FALSE]
}

# ---- counterfactuals: a directional log-shift of the whole distribution ------
# For a NON-LINEAR surface the counterfactual must be a distributional SHIFT
# (shape-preserving), not a collapse-to-a-point, because E[h(X)] depends on the
# whole distribution, not only its mean (unlike the linear single-pollutant case).
# All directions are one family: reduce log-concentration by `theta * dvec`.

#' Proportional reduction of every analyte (the overall-mixture direction, d_prop).
dir_proportional <- function(analytes) stats::setNames(rep(1, length(analytes)), analytes)

#' Reduce a single analyte only (d_top with j = the chosen analyte).
dir_analyte <- function(analytes, j) {
  d <- stats::setNames(rep(0, length(analytes)), analytes); d[[j]] <- 1; d
}

#' Reduce along the principal exposure axis (d_pc1): PC1 loadings of log-conc.
dir_pc1 <- function(Xmat) {
  L <- log(Xmat); pc <- stats::prcomp(L, center = TRUE, scale. = FALSE)$rotation[, 1]
  pc <- pc * sign(sum(pc))                       # orient toward "reduce exposure"
  stats::setNames(as.numeric(pc), colnames(Xmat))
}

#' Apply a counterfactual: shift log-concentrations by -theta*dvec (theta>=0).
apply_cf <- function(Xmat, theta, dvec) {
  dvec <- dvec[colnames(Xmat)]
  exp(sweep(log(Xmat), 2, theta * dvec, `-`))
}

#' theta that shifts the mean projected exposure index to its p-th percentile
#' (the near-background counterfactual, generalising HBEB's pct_floor). The index
#' is the projection of log-concentration onto dvec; shifting by theta lowers
#' every subject's index by theta*||dvec||^2, so this places the *mean* index at
#' the p-th percentile of the observed index distribution.
theta_background <- function(Xmat, dvec, p = 0.05) {
  dvec <- dvec[colnames(Xmat)]
  s <- as.numeric(log(Xmat) %*% dvec)            # per-subject index
  (mean(s) - as.numeric(stats::quantile(s, p))) / sum(dvec^2)
}

# ---- burden -----------------------------------------------------------------

#' Mean birth-weight shift Delta_mu between an observed and a counterfactual
#' exposure matrix under a surface (grams). Positive => reducing exposure raises
#' birth weight (for a harmful mixture), matching the single-pollutant sign.
mean_bw_shift <- function(Xobs, Xcf, surface) {
  mean(surface$predict(Xobs)) - mean(surface$predict(Xcf))
}

#' Attributable LBW cases between the OBSERVED population and a counterfactual.
#' `Xobs` is the observed exposure (mean birth weight = bw_mean by definition);
#' `Xcf` is a counterfactual reduction. Reduces EXACTLY to
#' burden_pipeline.R::attributable_cases_lbw when h(X)=beta*g(X_1) (verified in
#' benchmark_selftest A1). Do NOT pass a shifted scenario as `Xobs` -- that
#' re-anchors bw_mean; use `residual_cases()` for burden at a benchmark shift.
attributable_cases <- function(Xobs, Xcf, surface, tr = default_translation()) {
  d_mu  <- mean_bw_shift(Xobs, Xcf, surface)
  mu_cf <- tr$bw_mean - d_mu
  tr$N * (p_lbw(tr$bw_mean, tr) - p_lbw(mu_cf, tr))
}

#' Mean birth weight under a directional shift `theta`, anchored to the OBSERVED
#' distribution (not re-anchored): mu(theta) = bw_mean + (E[h(X_theta)] - E[h(X_obs)]).
mu_at_shift <- function(Xmat, surface, theta, dvec, tr) {
  tr$bw_mean + mean_bw_shift(apply_cf(Xmat, theta, dvec), Xmat, surface)
}

#' Attributable LBW REMAINING at a benchmark shift `theta`, relative to the
#' near-background scenario -- the quantity the benchmark inversion sets to a
#' target. Single-shot public version of the closure inside `benchmark_mixture`;
#' both share `mu_at_shift`, so they cannot drift. Matches
#' burden_pipeline.R::benchmark_central's residual (verified in selftest A2).
residual_cases <- function(X, surface, theta, dvec = NULL, bg_p = 0.05,
                           tr = default_translation()) {
  Xmat <- as_conc_matrix(X, surface$analytes)
  if (is.null(dvec)) dvec <- dir_proportional(surface$analytes)
  dvec <- dvec[surface$analytes]
  th_bg <- theta_background(Xmat, dvec, p = bg_p)
  mu_bg <- mu_at_shift(Xmat, surface, th_bg, dvec, tr)
  tr$N * (p_lbw(mu_at_shift(Xmat, surface, theta, dvec, tr), tr) - p_lbw(mu_bg, tr))
}

# ---- benchmark inversion ----------------------------------------------------

#' Invert the mixture burden to a benchmark along one reduction direction.
#'
#' Solves  AC_residual(theta*) = target_cases  for the directional shift theta*,
#' where AC_residual(theta) is the attributable LBW REMAINING at shift theta
#' relative to the near-background scenario -- the mixture analogue of
#' burden_pipeline.R::benchmark_central. Monotone decreasing in theta, so a
#' single uniroot on [0, theta_bg] suffices.
#'
#' @return list: theta (the shift), conc (named per-analyte benchmark
#'   concentrations = geometric mean * exp(-theta*dvec)), index (mean log-index
#'   reduction), current (current residual-vs-background cases), target, ok.
benchmark_mixture <- function(X, surface, target_cases,
                              dvec = NULL, bg_p = 0.05,
                              tr = default_translation()) {
  Xmat  <- as_conc_matrix(X, surface$analytes)
  if (is.null(dvec)) dvec <- dir_proportional(surface$analytes)
  dvec  <- dvec[surface$analytes]

  th_bg <- theta_background(Xmat, dvec, p = bg_p)
  if (!is.finite(th_bg) || th_bg <= 0)
    return(list(theta = NA_real_, conc = NA, index = NA, current = NA,
                target = target_cases, ok = FALSE, note = "no background range"))

  # attributable LBW remaining at shift theta, relative to near-background,
  # anchored to the observed distribution (see mu_at_shift / residual_cases).
  mu_bg    <- mu_at_shift(Xmat, surface, th_bg, dvec, tr)
  ac_resid <- function(theta)
    tr$N * (p_lbw(mu_at_shift(Xmat, surface, theta, dvec, tr), tr) - p_lbw(mu_bg, tr))
  current  <- ac_resid(0)
  if (target_cases >= current)             # target already met at observed exposure
    return(list(theta = 0, conc = gm(Xmat), index = 0, current = current,
                target = target_cases, ok = TRUE, note = "target >= current"))

  root <- tryCatch(
    stats::uniroot(function(th) ac_resid(th) - target_cases,
                   interval = c(0, th_bg))$root,
    error = function(e) NA_real_)
  if (!is.finite(root))
    return(list(theta = NA_real_, conc = NA, index = NA, current = current,
                target = target_cases, ok = FALSE, note = "uniroot failed"))

  conc <- gm(Xmat) * exp(-root * dvec)     # benchmark concentrations (per analyte)
  list(theta = root, conc = conc, index = root * sum(dvec^2),
       current = current, target = target_cases, ok = TRUE, note = "")
}

#' Geometric mean per column (the reported central concentration for a log-scale
#' mixture; matches burden_pipeline.R `central(scale="ln")`).
gm <- function(Xmat) exp(colMeans(log(Xmat)))

# ---- Monte Carlo uncertainty ------------------------------------------------

#' Monte-Carlo interval for a mixture benchmark.
#'
#' Two independent sources of benchmark uncertainty are propagated:
#'  * exposure sampling  -- bootstrap the exposure rows;
#'  * surface posterior  -- if the surface is predict-capable (`predict_draw`),
#'    fix one draw `s` per MC iteration and invert under it. Because `predict_draw`
#'    re-scores the *shifted* exposures, the directional counterfactual is honest
#'    under each draw (no fixed-row hack).
#' For the KNOWN DGP surface (no draws) only the exposure bootstrap contributes --
#' the truth benchmark carries no surface uncertainty, which is exactly the point.
#'
#' @return list(theta = c(median,lo,hi), conc = per-analyte quantile matrix).
mc_benchmark_mixture <- function(X, surface, target_cases, dvec = NULL,
                                 bg_p = 0.05, tr = default_translation(),
                                 n_mc = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Xmat <- as_conc_matrix(X, surface$analytes)
  n <- nrow(Xmat)
  has_draws <- !is.null(surface$predict_draw) && surface$n_draws >= 1L

  th <- numeric(n_mc)
  cc <- matrix(NA_real_, n_mc, length(surface$analytes),
               dimnames = list(NULL, surface$analytes))
  for (i in seq_len(n_mc)) {
    Xb <- Xmat[sample.int(n, n, replace = TRUE), , drop = FALSE]
    surf_i <- surface
    if (has_draws) {                                   # fix one predict-capable draw
      s <- sample.int(surface$n_draws, 1L)
      surf_i <- make_surface(function(Z) surface$predict_draw(Z, s),
                             surface$analytes)
    }
    b <- benchmark_mixture(Xb, surf_i, target_cases, dvec, bg_p, tr)
    th[i] <- b$theta
    if (isTRUE(b$ok) && length(b$conc) == ncol(cc)) cc[i, ] <- b$conc
  }
  list(theta = stats::quantile(th, c(.5, .025, .975), na.rm = TRUE),
       conc  = apply(cc, 2, stats::quantile, probs = c(.5, .025, .975), na.rm = TRUE))
}

# =============================================================================
# SELF-TEST -- the build-order step-1 gate (DESIGN §8): truth-in = truth-out
# =============================================================================

#' Run the correctness gate. Two checks:
#'  (A) single-pollutant REDUCTION parity: with h(X)=beta*log(X_1), the mixture
#'      burden + benchmark must match the shipped burden_pipeline.R to tolerance.
#'  (B) truth-in = truth-out round trip: for a known mixture h, the benchmark
#'      theta* must reproduce the target burden, AC_residual(theta*) == target.
#'
#' @param hbeb_burden_path path to the shipped burden_pipeline.R for the (A)
#'   parity check. If missing/not found, (A) uses an inline reference of the same
#'   formulas and says so.
#' @return invisibly TRUE if all checks pass; stops on failure.
benchmark_selftest <- function(
    hbeb_burden_path = Sys.getenv(
      "HBEB_BURDEN_R",
      "/Users/fabre/Library/CloudStorage/Dropbox-NIES-JECS/00-NIES/00-Publication/FY2026/04-Working/Nakayama-PFAS-EBD/hbeb-framework-nakayamasf/R/burden_pipeline.R"),
    tol = 1e-4, verbose = TRUE) {

  say <- function(...) if (verbose) cat(sprintf(...))
  tr  <- default_translation()
  N   <- tr$N

  # --- shared toy exposure sample (ln-scale analyte, like PFOS) --------------
  set.seed(20260814)
  x <- rlnorm(24000, meanlog = log(2.90), sdlog = 0.478)   # PFOS-like (Lai)
  beta <- -23.6                                            # Wright timing (g/ln)
  Xmat <- matrix(x, ncol = 1, dimnames = list(NULL, "pfos"))
  surf <- make_surface(function(M) beta * log(M[, "pfos"]), "pfos")

  # ===== (A) single-pollutant parity vs the shipped engine ===================
  have_ship <- nzchar(hbeb_burden_path) && file.exists(hbeb_burden_path)
  if (have_ship) {
    src_env <- new.env(); sys.source(hbeb_burden_path, envir = src_env)
    ref_ac  <- src_env$attributable_cases_lbw
    ref_bm  <- src_env$benchmark_central
    say("(A) parity vs shipped burden_pipeline.R\n")
  } else {
    # inline reference (identical formulas) if the HBEB repo is not on this machine
    ref_ac <- function(x, beta, scale, N, cf, bw_mean = tr$bw_mean,
                       bw_sd = tr$bw_sd, cut = tr$lbw_cut) {
      Eo <- mean(log(x))
      Ec <- if (cf$type == "pct_floor") log(as.numeric(quantile(x, cf$p)))
            else stop("test cf")
      d_mu <- beta * (Eo - Ec); mu_cf <- bw_mean - d_mu
      N * (pnorm(cut, bw_mean, bw_sd) - pnorm(cut, mu_cf, bw_sd))
    }
    ref_bm <- NULL
    say("(A) parity vs INLINE reference (shipped file not found at %s)\n",
        hbeb_burden_path)
  }

  # background counterfactual matched to HBEB's default pct_floor p=0.05
  bg_p <- 0.05
  dvec <- dir_proportional("pfos")

  # A1: attributable burden (observed vs 5th-percentile background)
  th_bg <- theta_background(Xmat, dvec, p = bg_p)
  Xbg   <- apply_cf(Xmat, th_bg, dvec)
  mine_ac <- attributable_cases(Xmat, Xbg, surf, tr)
  ref_ac_v <- ref_ac(x, beta, "ln", N, list(type = "pct_floor", p = bg_p))
  rel_ac <- abs(mine_ac - ref_ac_v) / abs(ref_ac_v)
  say("    A1 attributable cases: mine=%.2f ref=%.2f  rel.diff=%.2e\n",
      mine_ac, ref_ac_v, rel_ac)
  stopifnot(rel_ac < tol)

  # A2: 25%-reduction benchmark concentration
  curr <- mine_ac
  b <- benchmark_mixture(Xmat, surf, target_cases = curr * 0.75,
                         dvec = dvec, bg_p = bg_p, tr = tr)
  mine_bench <- as.numeric(b$conc[["pfos"]])
  if (!is.null(ref_bm)) {
    dat <- data.frame(pfos = x, pfoa = x)   # ref_bm indexes by tolower(chem)
    erf <- list(beta = beta, se = 7.68, tau = 3.0, scale = "ln")
    ref_bench <- ref_bm(dat, "PFOS", erf, N, target_cases = curr * 0.75,
                        cf_bg = list(type = "pct_floor", p = bg_p))
    rel_b <- abs(mine_bench - ref_bench) / abs(ref_bench)
    say("    A2 25%%-benchmark: mine=%.4f ref=%.4f  rel.diff=%.2e\n",
        mine_bench, ref_bench, rel_b)
    stopifnot(rel_b < tol)
  } else {
    say("    A2 25%%-benchmark: mine=%.4f (no shipped ref to compare)\n", mine_bench)
  }

  # ===== (B) truth-in = truth-out on a known MIXTURE surface ==================
  say("(B) round trip on a known non-linear mixture surface\n")
  set.seed(4)
  p <- 3L; an <- c("pfos", "pfoa", "pfhxs")
  med <- c(pfos = 2.90, pfoa = 1.60, pfhxs = 0.33)
  sdl <- c(pfos = 0.478, pfoa = 0.609, pfhxs = 0.55)
  n <- 8000
  L <- sapply(an, function(a) rnorm(n, log(med[[a]]), sdl[[a]]))
  Xmix <- exp(L); colnames(Xmix) <- an
  # known surface: main effects + curvature + interaction (grams), centred at medians
  u <- function(M) sweep(log2(M), 2, log2(med[an]), `-`)
  h_true <- function(M) {
    U <- u(M)
    -16.4*U[,"pfos"] -20.0*U[,"pfoa"] -8.0*U[,"pfhxs"] +
      2.5*U[,"pfos"]^2 - 4.0*U[,"pfos"]*U[,"pfoa"] - 3.0*U[,"pfos"]*U[,"pfhxs"]
  }
  surf_mix <- make_surface(h_true, an)

  dvec <- dir_proportional(an)
  th_bg <- theta_background(Xmix, dvec, p = bg_p)
  curr  <- residual_cases(Xmix, surf_mix, theta = 0, dvec = dvec, bg_p = bg_p, tr = tr)
  say("    current residual-vs-background: %.1f cases\n", curr)

  ok_all <- TRUE
  for (r in c(0.10, 0.25, 0.50)) {
    target <- curr * (1 - r)
    b <- benchmark_mixture(Xmix, surf_mix, target_cases = target,
                           dvec = dvec, bg_p = bg_p, tr = tr)
    # verify AC_residual(theta*) == target (the round-trip invariant), via the
    # public residual (which benchmark_mixture shares through mu_at_shift).
    got <- residual_cases(Xmix, surf_mix, theta = b$theta, dvec = dvec, bg_p = bg_p, tr = tr)
    rel <- abs(got - target) / abs(target)
    say("    %2.0f%% reduction: theta*=%.4f  GM(pfos)=%.3f  AC(theta*)=%.1f target=%.1f  rel=%.2e\n",
        r*100, b$theta, b$conc[["pfos"]], got, target, rel)
    ok_all <- ok_all && is.finite(rel) && rel < tol
  }
  stopifnot(ok_all)

  # B2: reducing to background must zero out the residual burden by construction
  b0 <- residual_cases(Xmix, surf_mix, theta = th_bg, dvec = dvec, bg_p = bg_p, tr = tr)
  say("    residual at background: %.3e (should be ~0)\n", b0)
  stopifnot(abs(b0) < 1e-6)

  say("ALL CHECKS PASSED.\n")
  invisible(TRUE)
}

if (sys.nframe() == 0 && identical(Sys.getenv("BENCHMARK_SELFTEST"), "1")) {
  benchmark_selftest()
}
