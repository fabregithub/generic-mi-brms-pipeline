# =============================================================================
# Analysis methods for the JECS PFAS -> Kawasaki disease demonstration.
# -----------------------------------------------------------------------------
# Four ways to handle the censored PFAS exposures, all scored against the
# complete-data ORACLE:
#   * oracle            — uncensored true concentrations (truth$conc_true)
#   * substitution      — half-limit fill (Iwata's stated method): ND/DNQ -> hi/2
#   * impute_noY        — congenial-style draw WITHOUT the outcome  (the wrong arm)
#   * impute_congenial  — congenial draw WITH the outcome           (the pipeline)
#
# Imputation uses leftcens::impute_censored_conditional() (the shipped Phase-4
# core) on the LOG scale, per PFAS, in a short block-FCS loop so each analyte
# conditions on the current values of the others (+ Z, + Y for the congenial arm).
#
# Estimand implemented here: SINGLE-POLLUTANT log-odds (per log2 ng/mL), one PFAS
# at a time adjusted for confounders, pooled across imputations by Rubin's rules.
# The BKMR mixture estimand lives in analyze_bkmr.R (heavier; pluggable engine).
# =============================================================================

# Covariate roles (columns of the demo data).
.PFAS      <- c("PFOA","PFOS","PFNA","PFDA","PFUnDA","PFHxS","PFTrDA")
.Z_IMPUTE  <- c("region","mat_age","parity","bmi","egfr","folate")      # X-block predictors (determinants)
.Z_ADJUST  <- c("mat_age","parity","bmi","smoking","income","child_male","low_bw","mat_asthma","region")  # Y-model confounders

# ---- shared covariate completion (MAR Z) -----------------------------------
# Simple, method-agnostic single fill so all four methods share identical Z; the
# production pipeline instead uses the miceRanger Z-block (00_censored_exposure.R).
complete_covariates <- function(data, cols = unique(c(.Z_IMPUTE, .Z_ADJUST))) {
  for (c in intersect(cols, names(data))) {
    x <- data[[c]]
    if (anyNA(x)) {
      if (is.numeric(x)) x[is.na(x)] <- stats::median(x, na.rm = TRUE)
      else { lv <- names(which.max(table(x))); x[is.na(x)] <- lv }
      data[[c]] <- x
    }
  }
  data
}

# ---- log-scale bounds for one PFAS from its _lo/_hi columns -----------------
.pfas_log_bounds <- function(data, j) {
  lo <- data[[paste0(j, "_lo")]]; hi <- data[[paste0(j, "_hi")]]
  obs <- data[[j]]                                   # value where quantified, else NA
  list(y     = ifelse(is.na(obs), NA_real_, log(obs)),
       lower = ifelse(is.na(obs), ifelse(lo <= 0, -Inf, log(lo)), NA_real_),
       upper = ifelse(is.na(obs), log(hi), NA_real_))
}

#' Produce `m` completed datasets (PFAS columns filled) under one method.
#' @return list of `m` data.frames (1 for oracle/substitution).
impute_pfas <- function(data, truth, method = c("oracle","substitution","impute_noY","impute_congenial"),
                        m = 10L, seed = 1L, sweeps = 2L, outcome = "kd") {
  method <- match.arg(method)
  d0 <- complete_covariates(data)

  if (method == "oracle") {
    d <- d0; for (j in .PFAS) d[[j]] <- truth$conc_true[, j]
    return(list(d))
  }
  if (method == "substitution") {
    d <- d0
    for (j in .PFAS) {
      obs <- data[[j]]; hi <- data[[paste0(j, "_hi")]]
      d[[j]] <- ifelse(is.na(obs), hi / 2, obs)      # half the upper limit (ND: MDL/2; DNQ: LCMRL/2)
    }
    return(list(d))
  }

  # --- congenial / no-Y multiple imputation ---------------------------------
  if (!requireNamespace("leftcens", quietly = TRUE) ||
      !"impute_censored_conditional" %in% getNamespaceExports("leftcens")) {
    stop("needs leftcens >= 0.9.0", call. = FALSE)
  }
  include_y <- method == "impute_congenial"
  bounds <- lapply(.PFAS, function(j) .pfas_log_bounds(data, j)); names(bounds) <- .PFAS
  z_pred <- d0[, intersect(.Z_IMPUTE, names(d0)), drop = FALSE]

  out <- vector("list", m)
  for (i in seq_len(m)) {
    set.seed(seed + i)
    d <- d0
    # initialise each PFAS on the log scale (observed, else interval midpoint).
    logx <- lapply(.PFAS, function(j) {
      b <- bounds[[j]]; v <- b$y
      cw <- is.na(v)
      mid <- ifelse(is.finite(b$lower[cw]) & is.finite(b$upper[cw]),
                    (b$lower[cw] + b$upper[cw]) / 2, b$upper[cw] - 0.5)
      v[cw] <- mid; v
    }); names(logx) <- .PFAS

    for (t in seq_len(sweeps)) for (j in .PFAS) {
      others <- setdiff(.PFAS, j)
      xpred <- z_pred
      for (o in others) xpred[[paste0("log", o)]] <- logx[[o]]     # current other-PFAS
      if (include_y) xpred[[outcome]] <- d[[outcome]]              # Y-aware (congeniality)
      b <- bounds[[j]]
      imp <- leftcens::impute_censored_conditional(
        y = b$y, x = xpred, lower = b$lower, upper = b$upper, m = 1L, margin = "shash")[, 1]
      logx[[j]] <- imp
    }
    for (j in .PFAS) d[[j]] <- exp(logx[[j]])                      # back to ng/mL
    out[[i]] <- d
  }
  out
}

# ---- single-pollutant estimand ---------------------------------------------
.rubin <- function(est, var, conf = 0.95) {
  ok <- is.finite(est) & is.finite(var); est <- est[ok]; var <- var[ok]; m <- length(est)
  if (m == 0) return(c(est = NA, se = NA, lo = NA, hi = NA))
  qbar <- mean(est); ubar <- mean(var); b <- if (m > 1) stats::var(est) else 0
  Tt <- ubar + (1 + 1/m) * b; df <- if (b > 0) (m-1)*(1 + ubar/((1+1/m)*b))^2 else Inf
  se <- sqrt(Tt); tc <- stats::qt(1 - (1-conf)/2, df)
  c(est = qbar, se = se, lo = qbar - tc*se, hi = qbar + tc*se)
}

#' Fit the single-pollutant logistic estimand (beta per log2 ng/mL) for every
#' PFAS, pooled across the completed datasets. Returns a data.frame keyed by pfas.
fit_single_pollutant <- function(completed, truth, outcome = "kd") {
  zadj <- setdiff(intersect(.Z_ADJUST, names(completed[[1]])), outcome)
  rhs <- paste(zadj, collapse = " + ")
  rows <- lapply(.PFAS, function(j) {
    est <- vars <- numeric(length(completed))
    for (i in seq_along(completed)) {
      d <- completed[[i]]; d$.l2 <- log2(pmax(d[[j]], .Machine$double.eps))
      fit <- tryCatch(stats::glm(stats::as.formula(paste(outcome, "~ .l2 +", rhs)),
                                 data = d, family = stats::binomial), error = function(e) NULL)
      if (is.null(fit) || !".l2" %in% rownames(summary(fit)$coefficients)) { est[i] <- NA; vars[i] <- NA; next }
      co <- summary(fit)$coefficients[".l2", ]; est[i] <- co[1]; vars[i] <- co[2]^2
    }
    p <- .rubin(est, vars)
    data.frame(pfas = j, beta = p["est"], se = p["se"], lo = p["lo"], hi = p["hi"],
               or = exp(p["est"]), row.names = NULL)
  })
  do.call(rbind, rows)
}
