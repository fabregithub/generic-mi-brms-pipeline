# STEP 7: does Rubin's finite-m correction apply to a POSTERIOR PREDICTIVE?
#
# DERIVED ABOVE: yes, and with the same formula. The m-mixture estimates the
# predictive variance V_inf = Ubar + B, but the interval is centred on Qbar,
# which estimates Q_inf with variance B/m. Covering a future y therefore needs
#     Var(y - Qbar) = V_inf + B/m = Ubar + (1 + 1/m) B.
# The magnitude differs from the parameter case because U_pred contains the
# outcome's RESIDUAL variance, so B/(Ubar+B) is much smaller.
#
# REGISTERED PREDICTIONS, before looking:
#   P1 the UNCORRECTED predictive interval under-covers, and more at small m
#   P2 the CORRECTED one is at nominal
#   P3 the effect is far smaller than for a parameter, because the residual
#      variance dominates U_pred
# P1 failing (no under-coverage at any m) would mean the derivation is wrong
# and the mixture needs no correction.
# MEASURED 2026-10-11, 300 reps (coverage of a 95% predictive interval against
# the KNOWN true value):
#     m      uncorrected        corrected
#     3    0.9446 +/- 0.0014   0.9469 +/- 0.0014
#     5    0.9467             0.9484
#    10    0.9497             0.9505
#    20    0.9508             0.9513
# All three predictions held: under-coverage without the correction, shrinking
# with m, improved by it, and small throughout. At m = 3 the corrected arm is
# still ~2 SE below nominal -- B is estimated on 2 df there and is itself very
# noisy, so the correction cannot fully close it.
paths <- list(cache = tempdir(), logs = tempdir(), objects = tempdir())
suppressMessages(suppressWarnings(source("../../00_common_functions.R")))

one_rep <- function(n = 400L, m = 5L, K = 400L, p_zmis = 0.5, p_ymis = 0.25) {
  Z <- rnorm(n); X <- rnorm(n)
  Y <- 1 + 0.8 * Z + 0.5 * X + rnorm(n)            # residual SD 1
  ymis <- runif(n) < p_ymis                         # outcome to be predicted
  # COVARIATE MISSINGNESS ONLY ON ROWS WHOSE OUTCOME IS OBSERVED.
  # The first version let Z be missing on prediction rows too and imputed it
  # from Y -- using the TRUE Y of a row whose Y is supposedly unknown. That
  # leaks the outcome into its own predictor, puts the predictive mean too
  # close to the truth while the interval keeps full residual width, and
  # OVER-covers in both arms (0.964-0.973 against nominal 0.95), which is what
  # exposed it. Restricting Z's missingness to observed-Y rows removes the leak
  # and leaves exactly the quantity under test: between-imputation variation in
  # the FITTED PARAMETERS, which is what B measures in the real Step 7.
  zmis <- (runif(n) < p_zmis) & !ymis
  Ztrue <- Z; Z[zmis] <- NA
  obs <- !ymis
  out_lo <- out_hi <- out_lo_c <- out_hi_c <- rep(NA_real_, sum(ymis))
  draws <- vector("list", m)
  for (i in seq_len(m)) {
    # impute Z from (Y, X) on rows where Y is observed; proper-ish draw
    zi <- Z
    fz <- lm(Z ~ Y + X, data = data.frame(Z = Z, Y = Y, X = X)[obs, ])
    nd <- data.frame(Y = Y, X = X)
    zi[zmis] <- predict(fz, newdata = nd[zmis, ]) + rnorm(sum(zmis), 0, sigma(fz))
    # fit on observed-Y rows, predict the missing-Y rows
    fit <- lm(Y ~ zi + X, data = data.frame(Y = Y, zi = zi, X = X)[obs, ])
    s2 <- sum(residuals(fit)^2) / rchisq(1, fit$df.residual)
    bh <- as.vector(coef(fit) + t(chol(summary(fit)$cov.unscaled)) %*%
                      rnorm(length(coef(fit))) * sqrt(s2))
    mm <- model.matrix(~ zi + X, data = data.frame(zi = zi, X = X))[ymis, , drop = FALSE]
    mu <- as.vector(mm %*% bh)
    draws[[i]] <- vapply(seq_len(K), function(k) mu + rnorm(length(mu), 0, sqrt(s2)),
                         numeric(length(mu)))
  }
  # per predicted row: pool across imputations, uncorrected vs corrected
  cov_u <- cov_c <- logical(sum(ymis)); wid_u <- wid_c <- numeric(sum(ymis))
  for (r in seq_len(sum(ymis))) {
    v <- unlist(lapply(draws, function(D) D[r, ]))
    imp <- rep(seq_len(m), each = K)
    q_u <- quantile(v, c(0.025, 0.975))                       # Step 7 today
    p <- mi_pool_draws(v, imp)
    q_c <- mi_weighted_quantile(p$value, p$weight, c(0.025, 0.975))
    yt <- Y[ymis][r]
    cov_u[r] <- yt >= q_u[1] && yt <= q_u[2]
    cov_c[r] <- yt >= q_c[1] && yt <= q_c[2]
    wid_u[r] <- diff(q_u); wid_c[r] <- diff(q_c)
  }
  c(cov_u = mean(cov_u), cov_c = mean(cov_c),
    wid_u = mean(wid_u), wid_c = mean(wid_c))
}

R <- as.integer(Sys.getenv("R", "200"))
cat(sprintf("%4s | %13s %13s | %9s %9s\n", "m",
            "coverage unc.", "coverage corr.", "width unc", "width corr"))
for (m in c(3L, 5L, 10L, 20L)) {
  set.seed(100 + m)
  res <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) one_rep(m = m),
                                           mc.cores = 10))
  se <- function(x) sd(x)/sqrt(length(x))
  cat(sprintf("%4d | %7.4f+/-%.4f %7.4f+/-%.4f | %9.4f %9.4f\n", m,
              mean(res[,"cov_u"]), se(res[,"cov_u"]),
              mean(res[,"cov_c"]), se(res[,"cov_c"]),
              mean(res[,"wid_u"]), mean(res[,"wid_c"])))
}
cat("\nnominal coverage is 0.95\n")
