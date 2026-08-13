# =============================================================================
# Phase 1 harness -- the four "no new code" procedures (PLAN §7.3, §10 Phase 1)
# -----------------------------------------------------------------------------
# Each procedure fits the SAME matched analysis model (dgp_formula) and returns a
# one-row data.frame with the focal-exposure estimand and a 95% interval:
#   procedure, estimate, se, ci_lo, ci_hi, note
#
#   1. oracle            -- complete data, no missingness (best-case reference)
#   2. complete_case     -- listwise-drop censored/missing rows
#   3. leftcens_prestep  -- leftcens::gsimp_mi imputes X WITHOUT Y, then pool
#                           (the current architecture; HYPOTHESISED BIASED)
#   4. brms_joint        -- congenial joint model bf(Y ~ mi(X)) + bf(X|mi()+cens())
#                           (gold standard; needs no new package code)
#
# Nothing here calls a component that does not already exist -- that is the whole
# point of Phase 1 (confirm the problem and the joint-model fix before building
# the bespoke §4 engine).
# =============================================================================

# ---- shared helpers ---------------------------------------------------------

#' Fit the matched model and pull out the focal estimand coefficient.
#' @return named numeric c(est, se) or c(NA, NA) on failure.
fit_lm_estimand <- function(data, truth) {
  fm <- dgp_formula(truth)
  fit <- tryCatch(stats::lm(fm, data = data), error = function(e) NULL)
  if (is.null(fit)) return(c(est = NA_real_, se = NA_real_))
  cf <- stats::coef(fit); vc <- stats::vcov(fit)
  nm <- "logX1"
  if (!nm %in% names(cf)) return(c(est = NA_real_, se = NA_real_))
  c(est = unname(cf[nm]), se = unname(sqrt(vc[nm, nm])))
}

#' Rubin's rules for one scalar estimand across m completed-data fits.
rubin_pool <- function(est, var, conf = 0.95) {
  ok <- is.finite(est) & is.finite(var)
  est <- est[ok]; var <- var[ok]
  m <- length(est)
  if (m == 0L) return(c(est = NA, se = NA, ci_lo = NA, ci_hi = NA, df = NA))
  qbar <- mean(est); ubar <- mean(var)
  b <- if (m > 1L) stats::var(est) else 0
  Tt <- ubar + (1 + 1 / m) * b
  df <- if (b > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2 else Inf
  se <- sqrt(Tt); tcrit <- stats::qt(1 - (1 - conf) / 2, df)
  c(est = qbar, se = se, ci_lo = qbar - tcrit * se, ci_hi = qbar + tcrit * se, df = df)
}

one_row <- function(procedure, est, se, ci_lo, ci_hi, note = "") {
  data.frame(procedure = procedure, estimate = est, se = se,
             ci_lo = ci_lo, ci_hi = ci_hi, note = note,
             stringsAsFactors = FALSE)
}

# ---- procedure 1: oracle ----------------------------------------------------

proc_oracle <- function(bundle) {
  e <- fit_lm_estimand(bundle$complete, bundle$truth)
  ci <- e["est"] + c(-1, 1) * stats::qnorm(0.975) * e["se"]
  one_row("oracle", e["est"], e["se"], ci[1], ci[2])
}

# ---- procedure 2: complete-case --------------------------------------------

proc_complete_case <- function(bundle) {
  # `logX1` is NA on censored rows, so lm() listwise-drops them = complete case.
  e <- fit_lm_estimand(bundle$censored, bundle$truth)
  ci <- e["est"] + c(-1, 1) * stats::qnorm(0.975) * e["se"]
  one_row("complete_case", e["est"], e["se"], ci[1], ci[2])
}

# ---- procedure 3: independent leftcens pre-step (NO Y) ----------------------

proc_leftcens_prestep <- function(bundle, m = 10L, seed = NULL, n_cores = 1L) {
  d <- bundle$censored; truth <- bundle$truth
  xcols <- grep("^logX[0-9]+$", names(d), value = TRUE)
  censored_cols <- attr(d, "censored_cols")

  # Build interval-censored bounds on the CONCENTRATION scale (log_transform
  # handled internally by leftcens); imputations come back on the LOG scale,
  # matching mc_validation.R's convention. Y is deliberately NOT a column here.
  cens_list <- stats::setNames(lapply(xcols, function(xc) {
    val <- d[[xc]]                                  # NA where censored
    if (xc %in% censored_cols) {
      is_cens <- d[[paste0(xc, "_cens")]] == "left"
      lod_conc <- exp(d[[paste0(xc, "_lod")]])
      left  <- ifelse(is_cens, 0, exp(val))
      right <- ifelse(is_cens, lod_conc, exp(val))
    } else {
      left <- right <- exp(val)                     # fully observed analyte
    }
    leftcens::as_interval_data(left, right, log_transform = TRUE)
  }), xcols)

  bnds <- leftcens::build_bounds(cens_list)
  mi <- tryCatch(
    leftcens::gsimp_mi(bnds, m = as.integer(m), return_imputations = TRUE,
                       seed = seed, n_cores = as.integer(n_cores)),
    error = function(e) NULL
  )
  if (is.null(mi) || is.null(mi$imputations)) {
    return(one_row("leftcens_prestep", NA, NA, NA, NA, "gsimp_mi failed"))
  }

  # Fit the matched model on each completed dataset; pool with Rubin's rules.
  ests <- vars <- numeric(length(mi$imputations))
  for (i in seq_along(mi$imputations)) {
    filled <- mi$imputations[[i]]                   # n x p, LOG scale
    dfit <- d
    for (xc in xcols) dfit[[xc]] <- filled[, xc]
    e <- fit_lm_estimand(dfit, truth)
    ests[i] <- e["est"]; vars[i] <- e["se"]^2
  }
  p <- rubin_pool(ests, vars)
  one_row("leftcens_prestep", p["est"], p["se"], p["ci_lo"], p["ci_hi"],
          sprintf("m=%d", length(mi$imputations)))
}

# ---- procedure 4: congenial censored MI including Y (the gold standard) ------

#' Congenial multiple imputation of the left-censored focal exposure, INCLUDING
#' the outcome Y in the imputation model.
#'
#' This is the honest Phase-1 reference and a preview of the §4 component. It is
#' the congenial construction the pre-step (proc_leftcens_prestep) deliberately
#' omits: it draws the missing X1 from p(X1 | Y, Z, other X) respecting the LOD.
#'
#' Draw (proper MI):
#'  1. Fit an interval-censored (left-censored) Gaussian regression of X1 on
#'     (Y, other X, Z) by MLE  ->  survival::survreg(type = "left").
#'  2. Per imputation, draw the coefficients AND log-scale from their asymptotic
#'     posterior (MVN at the MLE), so between-imputation variance is honest.
#'  3. Impute each censored X1 from N(mu_i*, sigma*) TRUNCATED below its LOD
#'     (leftcens::rnorm_trunc) -- respects the censoring bound.
#'  4. Fit the matched analysis model on each completed dataset; pool by Rubin.
#'
#' Correctly specified for the additive (log-linear-Gaussian) DGP; for a
#' non-linear mixture surface the linear imputation model is only approximate
#' (the §7.7 congeniality question), so it is reported with a note there.
proc_censored_mi_y <- function(bundle, m = 20L, seed = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    return(one_row("cens_mi_y", NA, NA, NA, NA, "survival not installed"))
  }
  if (!is.null(seed)) set.seed(seed)
  d <- bundle$censored; truth <- bundle$truth
  p <- length(truth$b); q <- length(truth$gamma)
  other_x <- if (p >= 2L) paste0("logX", 2:p) else character(0)
  preds <- c("Y", other_x, c("Z1", "Z2")[seq_len(q)])
  rhs <- paste(preds, collapse = " + ")

  # Response for the censored analyte: the value where observed, the LOD where
  # left-censored; `event` flags exact (1) vs left-censored (0) observations.
  d$event <- as.integer(d$logX1_cens == "none")
  d$resp  <- ifelse(d$logX1_cens == "none", d$logX1, d$logX1_lod)

  sr <- tryCatch(
    survival::survreg(
      stats::as.formula(sprintf('survival::Surv(resp, event, type = "left") ~ %s', rhs)),
      data = d, dist = "gaussian"),
    error = function(e) NULL
  )
  if (is.null(sr)) return(one_row("cens_mi_y", NA, NA, NA, NA, "survreg failed"))

  mm <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = d)
  k <- length(stats::coef(sr))
  mu_hat <- c(stats::coef(sr), log(sr$scale))     # coefficients + log(scale)
  V <- stats::vcov(sr)                            # (k+1) x (k+1)
  cidx <- which(d$logX1_cens == "left")

  ests <- vars <- numeric(m)
  for (i in seq_len(m)) {
    dr <- MASS::mvrnorm(1, mu_hat, V)             # proper-MI parameter draw
    mu_i <- as.vector(mm %*% dr[seq_len(k)]); sigma_star <- exp(dr[k + 1L])
    x1 <- d$logX1
    x1[cidx] <- leftcens::rnorm_trunc(length(cidx), mu_i[cidx], sigma_star,
                                      lower = -Inf, upper = d$logX1_lod[cidx])
    dfit <- d; dfit$logX1 <- x1
    e <- fit_lm_estimand(dfit, truth); ests[i] <- e["est"]; vars[i] <- e["se"]^2
  }
  pooled <- rubin_pool(ests, vars)
  note <- if (truth$erf_form == "mixture") sprintf("m=%d; linear-approx (§7.7)", m) else sprintf("m=%d", m)
  one_row("cens_mi_y", pooled["est"], pooled["se"], pooled["ci_lo"], pooled["ci_hi"], note)
}

# ---- procedure 4-shash: skew-aware congenial censored MI including Y ---------

#' Skew-aware variant of proc_censored_mi_y: the §4 recipe in full.
#'
#' `cens_mi_y` uses a Gaussian (tobit) conditional, which is misspecified when the
#' log-exposure margin is skewed and then breaks (bias + coverage collapse under
#' right-skew; see validation/phase1/FINDINGS.md). This version fits X1's margin
#' with a **sinh-arcsinh** interval-censored MLE (`leftcens`'s `fit_shash_margin`),
#' maps X1 and the LOD to a **latent normal scale**, runs the `Y`-aware censored
#' regression *there* (correctly specified), draws truncated below the latent LOD,
#' and back-transforms. Proper-MI draws of BOTH the margin params and the
#' regression params. Degrades gracefully to Gaussian when eps≈0, so it is a valid
#' gold standard at any skew. This is the §4 component prototyped.
#'
#' Uses three `leftcens` internals (`:::`) — acceptable in this validation harness;
#' they become exported API when §4 lands.
proc_censored_mi_y_shash <- function(bundle, m = 20L, seed = NULL) {
  need <- all(vapply(c("fit_shash_margin", "draw_margin", "x_to_z", "z_to_x"),
                     function(fn) exists(fn, where = asNamespace("leftcens"), inherits = FALSE),
                     logical(1)))
  if (!requireNamespace("survival", quietly = TRUE) || !need) {
    return(one_row("cens_mi_y_shash", NA, NA, NA, NA, "leftcens shash internals unavailable"))
  }
  fsm <- get("fit_shash_margin", asNamespace("leftcens"))
  dm  <- get("draw_margin",      asNamespace("leftcens"))
  x2z <- get("x_to_z",           asNamespace("leftcens"))
  z2x <- get("z_to_x",           asNamespace("leftcens"))

  if (!is.null(seed)) set.seed(seed)
  d <- bundle$censored; truth <- bundle$truth
  p <- length(truth$b); q <- length(truth$gamma)
  preds <- c("Y", if (p >= 2L) paste0("logX", 2:p) else character(0), c("Z1", "Z2")[seq_len(q)])
  rhs <- paste(preds, collapse = " + ")
  cidx <- which(d$logX1_cens == "left"); oidx <- which(d$logX1_cens == "none")
  lod <- d$logX1_lod

  # 1. marginal sinh-arcsinh fit of X1 (observed exact + censored intervals).
  mfit <- fsm(d$logX1[oidx], lo_c = rep(-Inf, length(cidx)), hi_c = lod[cidx])
  mm <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = d)

  ests <- vars <- numeric(m)
  for (i in seq_len(m)) {
    mp <- dm(mfit)                                       # proper-MI margin draw
    z_obs <- x2z(d$logX1, mp$mu, mp$sigma, mp$eps)       # latent normal scores
    z_lod <- x2z(lod,     mp$mu, mp$sigma, mp$eps)
    dd <- d; dd$zresp <- ifelse(d$logX1_cens == "none", z_obs, z_lod)
    dd$event <- as.integer(d$logX1_cens == "none")
    sr <- tryCatch(survival::survreg(
      stats::as.formula(sprintf('survival::Surv(zresp, event, type = "left") ~ %s', rhs)),
      data = dd, dist = "gaussian"), error = function(e) NULL)
    if (is.null(sr)) { ests[i] <- NA; vars[i] <- NA; next }
    k <- length(stats::coef(sr))
    dr <- MASS::mvrnorm(1, c(stats::coef(sr), log(sr$scale)), stats::vcov(sr))
    mu_i <- as.vector(mm %*% dr[seq_len(k)]); s_reg <- exp(dr[k + 1L])
    z1 <- z_obs
    z1[cidx] <- leftcens::rnorm_trunc(length(cidx), mu_i[cidx], s_reg,
                                      lower = -Inf, upper = z_lod[cidx])
    x1 <- d$logX1; x1[cidx] <- z2x(z1[cidx], mp$mu, mp$sigma, mp$eps)
    dfit <- d; dfit$logX1 <- x1
    e <- fit_lm_estimand(dfit, truth); ests[i] <- e["est"]; vars[i] <- e["se"]^2
  }
  pooled <- rubin_pool(ests, vars)
  one_row("cens_mi_y_shash", pooled["est"], pooled["se"], pooled["ci_lo"], pooled["ci_hi"],
          sprintf("m=%d; shash eps=%.2f", m, mfit$eps))
}

# ---- procedure 4b: brms joint model (DEPRECATED in scaffold) -----------------
# NOTE: in the installed brms, `x1 | mi() + cens()` with the censored value set
# to the LOD creates NO latent per-observation parameters, so `mi(x1)` in the Y
# model silently uses the LOD-substituted values -- i.e. this collapses to
# LOD-substitution and is biased (verified: it equals `lm` on LOD-filled X).
# Kept for reference only; NOT in the default procedure set. Use proc_censored_mi_y.

#' @param cache An environment used to reuse the compiled Stan model across reps
#'   (compile once, then `update(recompile = FALSE)`); pass a fresh env per grid
#'   cell. If NULL, compiles every call.
proc_brms_joint <- function(bundle, control = list(), cache = NULL) {
  truth <- bundle$truth
  if (truth$erf_form != "additive") {
    # Non-linear mi() transforms (interaction, I(mi(x)^2)) are exactly the §7.7
    # congeniality question -- deliberately out of scope for the Phase-1 scaffold.
    return(one_row("brms_joint", NA, NA, NA, NA, "mixture: Phase-1 extension (§7.7)"))
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(one_row("brms_joint", NA, NA, NA, NA, "brms not installed"))
  }

  d <- bundle$censored
  p <- length(truth$b); q <- length(truth$gamma)
  # Response for the censored analyte: observed value, or the LOD on censored rows.
  d$x1_resp <- ifelse(d$logX1_cens == "left", d$logX1_lod, d$logX1)
  d$x1_cens <- d$logX1_cens

  other_x <- if (p >= 2L) paste0("logX", 2:p) else character(0)
  zt <- c("Z1", "Z2")[seq_len(q)]
  rhs_pred <- paste(c(other_x, zt), collapse = " + ")

  bf_y <- brms::bf(stats::as.formula(paste("Y ~ mi(x1_resp) +", rhs_pred)))
  bf_x <- brms::bf(stats::as.formula(paste("x1_resp | mi() + cens(x1_cens) ~", rhs_pred)))

  ctl <- utils::modifyList(
    list(chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0,
         seed = 1, backend = getOption("phase1.brms_backend", "rstan")),
    control
  )

  fit <- NULL
  if (!is.null(cache) && !is.null(cache$fit)) {
    fit <- tryCatch(
      stats::update(cache$fit, newdata = d, recompile = FALSE,
                    chains = ctl$chains, iter = ctl$iter, warmup = ctl$warmup,
                    cores = ctl$cores, refresh = ctl$refresh, seed = ctl$seed),
      error = function(e) NULL
    )
  }
  fit_err <- NULL
  if (is.null(fit)) {
    fit <- tryCatch(
      brms::brm(bf_y + bf_x + brms::set_rescor(FALSE), data = d,
                chains = ctl$chains, iter = ctl$iter, warmup = ctl$warmup,
                cores = ctl$cores, refresh = ctl$refresh, seed = ctl$seed,
                backend = ctl$backend),
      error = function(e) { fit_err <<- conditionMessage(e); NULL }
    )
    if (!is.null(cache) && !is.null(fit)) cache$fit <- fit
  }
  if (is.null(fit)) {
    msg <- if (!is.null(fit_err)) paste("brms fit failed:", fit_err) else "brms fit failed"
    return(one_row("brms_joint", NA, NA, NA, NA, substr(msg, 1, 120)))
  }

  # Estimand = coefficient of mi(x1_resp) in the Y model. Grep the draws
  # robustly (brms names it like `bsp_Y_mix1_resp`).
  draws <- as.data.frame(brms::as_draws_df(fit))
  vn <- grep("^bsp_Y_.*x1", names(draws), value = TRUE, ignore.case = TRUE)
  if (length(vn) == 0L) vn <- grep("x1_resp", names(draws), value = TRUE)[1]
  if (length(vn) == 0L || is.na(vn[1])) {
    return(one_row("brms_joint", NA, NA, NA, NA, "estimand param not found"))
  }
  post <- draws[[vn[1]]]
  ci <- stats::quantile(post, c(0.025, 0.975), names = FALSE)
  one_row("brms_joint", mean(post), stats::sd(post), ci[1], ci[2], vn[1])
}

# ---- registry ---------------------------------------------------------------

#' Run a set of procedures on one data bundle; return a stacked data.frame.
#'
#' Default set uses `cens_mi_y_shash` (the skew-aware congenial gold standard).
#' `cens_mi_y` (Gaussian/tobit) is available but breaks under skew (see FINDINGS);
#' `brms_joint` is deprecated in the scaffold (see proc_brms_joint).
run_procedures <- function(bundle, which = c("oracle", "complete_case",
                                             "leftcens_prestep", "cens_mi_y_shash"),
                           m = 10L, seed = NULL, n_cores = 1L,
                           brms_control = list(), brms_cache = NULL) {
  out <- list()
  if ("oracle" %in% which)           out$oracle <- proc_oracle(bundle)
  if ("complete_case" %in% which)    out$cc     <- proc_complete_case(bundle)
  if ("leftcens_prestep" %in% which) out$lc     <- proc_leftcens_prestep(bundle, m = m, seed = seed, n_cores = n_cores)
  if ("cens_mi_y" %in% which)        out$cmy    <- proc_censored_mi_y(bundle, m = m, seed = seed)
  if ("cens_mi_y_shash" %in% which)  out$cmys   <- proc_censored_mi_y_shash(bundle, m = m, seed = seed)
  if ("brms_joint" %in% which)       out$brms   <- proc_brms_joint(bundle, control = brms_control, cache = brms_cache)
  do.call(rbind, out)
}
