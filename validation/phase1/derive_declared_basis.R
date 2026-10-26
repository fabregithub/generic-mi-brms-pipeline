#!/usr/bin/env Rscript
# =============================================================================
# P-B -- the QUANTITATIVE test of the projection account of V33/V34.
# -----------------------------------------------------------------------------
# `DERIVATION_declared_basis.md` says that declaring an INCOMPLETE outcome basis
# makes the plugin fit an omitted-variable regression,
#       theta_tilde_ret = theta_ret + Gamma theta_omit,
# so the omitted term is attenuated and the retained ones are INFLATED by
# Gamma*theta_omit. P-A (signs and targets) has survived every test. P-B is the
# magnitude: Gamma is computable from the design with no free parameters, so the
# account must predict V33's +29.6% curvature inflation (omit x1x2) and +10.1%
# interaction inflation (omit x1^2) -- and V34's f-scaling -- or be wrong.
#
# WHY THIS IS AN INDEPENDENT IMPLEMENTATION. The DGP, the censored draw, the
# plugin fit and the analysis model are all re-written here from the generator's
# definition rather than sourced from `R/`. Nothing but `MASS` is shared with the
# harness. Two implementations agreeing is then evidence; one implementation
# agreeing with itself is not. (Same discipline as `derive_scale_retention.R`.)
#
# WHAT READING THE CODE CHANGED. `R/smc_impute.R:238` fits the declared surface
# with `lm(fm, data = w)` -- on the CURRENT COMPLETED data, imputed rows
# included, not on the observed rows alone. So `h_dec` is not a one-shot
# projection over the observed region as section 3 of the derivation assumed: it
# is a FIXED POINT of
#       h_dec  ->  impute censored x1 given h_dec  ->  refit h_dec on all rows.
# That feedback is absent at f = 0 and strengthens as f grows, because the
# imputed rows increasingly determine the very surface used to impute them. It
# is therefore a candidate source of super-linear f-scaling that was never in
# the account. STAGE 3 computes the fixed point; STAGE 2 computes the one-shot
# projection the written account actually predicts. Comparing them separates
# "the projection argument" from "the feedback the code adds to it".
#
# STAGES
#   1  Gamma over the observed region, and theta_tilde_ret. Checks section 3.
#   2  One-shot projection prediction (x1 at its TRUE values on censored rows).
#      This is the literal content of section 4. No free parameters.
#   3  Fixed-point prediction WITH the x1 redraw, iterated to convergence.
#      This is what the imputer actually does, asymptotically. No free
#      parameters.
#   4  Harness-faithful finite-sample check: n = 800, m = 10, 3 sweeps, proper
#      beta draw, cold start -- but an OLS analysis model on the TRUE basis
#      instead of BKMR. Isolates whether the residual gap is BKMR's or the
#      imputation's.
#
# The estimands are pure: `curv_X1 = b_quad*((c-b)^2+(a-b)^2)` and
# `int_X1X2 = b_int*(c-a)^2` (R/estimands_bkmr.R), so relative bias in each
# estimand IS relative bias in the corresponding coefficient, and a multiplicative
# shrinkage common to the arm and the oracle cancels in the vs-oracle ratio.
# That is the one assumption carried into the comparison with BKMR numbers.
#
# Usage:  Rscript derive_declared_basis.R [N] [REPS]
# Cost:   ~3 min at the defaults on one core.
# =============================================================================

suppressPackageStartupMessages(library(MASS))

args  <- commandArgs(trailingOnly = TRUE)
NBIG  <- if (length(args) >= 1) as.integer(args[1]) else 300000L
REPS  <- if (length(args) >= 2) as.integer(args[2]) else 3L

# ---- the generator, re-written from its definition --------------------------
# make_truth(p=3, q=2, erf_form="mixture"); simulate_complete(rho=0.4, sd_x=1,
# mu_x=0, sigma_y=1); z_form="linear" so Z is independent of X.
B1 <- 0.40; B2 <- 0.20; B3 <- 0.10
G1 <- 0.50; G2 <- -0.30
BINT <- 0.25; BQUAD <- 0.15
RHO <- 0.4; SIG <- 1.0

h_true <- function(d) {
  B1*d$x1 + B2*d$x2 + B3*d$x3 + G1*d$z1 + G2*d$z2 +
    BINT*d$x1*d$x2 + BQUAD*d$x1^2
}

dgp <- function(n, seed) {
  set.seed(seed)
  S <- matrix(RHO, 3, 3); diag(S) <- 1
  X <- MASS::mvrnorm(n, rep(0, 3), S)
  d <- data.frame(x1 = X[,1], x2 = X[,2], x3 = X[,3],
                  z1 = rnorm(n), z2 = rbinom(n, 1, 0.5))
  d$Y <- h_true(d) + rnorm(n, 0, SIG)
  d
}

# conditional of x1 given (x2,x3) under the exchangeable Sigma -- closed form
S  <- matrix(RHO, 3, 3); diag(S) <- 1
W  <- solve(S[2:3, 2:3], S[2:3, 1])          # (0.285714, 0.285714)
VC <- as.numeric(1 - S[1, 2:3] %*% W)        # 0.771429
SC <- sqrt(VC)

# ---- bases ------------------------------------------------------------------
TERMS <- list(
  true   = c("x1","x2","x3","z1","z2","x1x2","x1sq"),
  noint  = c("x1","x2","x3","z1","z2","x1sq"),
  noquad = c("x1","x2","x3","z1","z2","x1x2"),
  lin    = c("x1","x2","x3","z1","z2")
)
OMIT <- list(noint = "x1x2", noquad = "x1sq", lin = c("x1x2","x1sq"))

mm <- function(d, terms) {
  cols <- lapply(terms, function(t) switch(t,
    x1 = d$x1, x2 = d$x2, x3 = d$x3, z1 = d$z1, z2 = d$z2,
    x1x2 = d$x1*d$x2, x1sq = d$x1^2))
  M <- cbind(1, do.call(cbind, cols))
  colnames(M) <- c("(Int)", terms)
  M
}

# OLS via normal equations, with optional row weights.
ols <- function(M, y, w = NULL) {
  if (is.null(w)) { A <- crossprod(M); b <- crossprod(M, y) }
  else            { Mw <- M * w; A <- crossprod(Mw, M); b <- crossprod(Mw, y) }
  drop(solve(A, b))
}

FS <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

# =============================================================================
# STAGE 1 -- Gamma over the observed region, and theta_tilde_ret
# =============================================================================
cat("\n================ STAGE 1: Gamma over the OBSERVED region ================\n")
cat("Gamma = E_obs[phi_ret phi_ret']^-1 E_obs[phi_ret phi_omit']; theta_omit =",
    sprintf("b_int %.2f, b_quad %.2f\n", BINT, BQUAD))
d0 <- dgp(NBIG, 1001L)
Mfull0 <- mm(d0, TERMS$true)

st1 <- NULL
for (f in FS) {
  L   <- qnorm(f)
  obs <- d0$x1 > L
  for (dn in c("noint","noquad")) {
    ret <- TERMS[[dn]]; om <- OMIT[[dn]]
    Mr <- mm(d0, ret)[obs, , drop = FALSE]
    Mo <- Mfull0[obs, om, drop = FALSE]
    Gam <- solve(crossprod(Mr), crossprod(Mr, Mo))        # (ret+1) x (omit)
    th_om <- if (om == "x1x2") BINT else BQUAD
    # the rerouting onto the OTHER non-linear term retained in this basis
    tgt <- if (dn == "noint") "x1sq" else "x1x2"
    st1 <- rbind(st1, data.frame(
      f = f, declared = dn, omitted = om, rerouted_onto = tgt,
      Gamma = Gam[tgt, 1], shift = Gam[tgt, 1] * th_om,
      rel_shift = Gam[tgt, 1] * th_om / (if (tgt == "x1sq") BQUAD else BINT)))
  }
}
print(st1, row.names = FALSE, digits = 4)

# =============================================================================
# STAGE 2 -- the one-shot projection: section 4 taken literally
# =============================================================================
# On observed rows the analysis model sees the true surface; on censored rows it
# sees h_dec, the declared-basis fit. x1 is NOT redrawn here -- this isolates the
# projection argument from the imputation's effect on the x1 distribution.
cat("\n================ STAGE 2: one-shot projection (no x1 redraw) ================\n")
st2 <- NULL
for (f in FS) {
  L   <- qnorm(f); obs <- d0$x1 > L
  for (dn in c("noint","noquad","lin")) {
    Mr <- mm(d0, TERMS[[dn]])
    bdec <- ols(Mr[obs, , drop = FALSE], d0$Y[obs])        # declared fit, obs rows
    heff <- ifelse(obs, h_true(d0), drop(Mr %*% bdec))
    th <- ols(Mfull0, heff)
    st2 <- rbind(st2, data.frame(
      f = f, declared = dn,
      curv = unname(th["x1sq"]/BQUAD - 1),
      int  = unname(th["x1x2"]/BINT  - 1)))
  }
}
print(st2, row.names = FALSE, digits = 3)

# =============================================================================
# STAGE 3 -- the fixed point WITH the x1 redraw (what the code actually does)
# =============================================================================
# Censored draw:  p(x1 | rest, Y) propto N(x1; mu_i, s^2) * N(Y; h_dec(x1), sig^2)
#                                   * 1{x1 <= L}
# h_dec is linear in its columns, each of which is const / linear / quadratic in
# x1, so h_dec(x1) = A_i + B_i x1 + C_i x1^2 and the log posterior is a quartic.
# Gridded, log-sum-exp, inverse-CDF. h_dec is refitted on ALL rows each sweep.
NGRID <- 600L

redraw_x1 <- function(d, idx, L, bdec, terms, sdraw, chunk = 20000L) {
  if (length(idx) > chunk) {
    sp <- split(idx, ceiling(seq_along(idx)/chunk))
    rs <- lapply(sp, function(ii) redraw_x1(d, ii, L, bdec, terms, sdraw, chunk))
    return(list(x = unlist(lapply(rs, `[[`, "x"), use.names = FALSE),
                edge = max(vapply(rs, `[[`, 0, "edge"))))
  }
  nb <- length(idx)
  g <- function(t) switch(t, x1 = NULL, x1x2 = NULL, x1sq = NULL,
                          x2 = d$x2[idx], x3 = d$x3[idx],
                          z1 = d$z1[idx], z2 = d$z2[idx])
  A <- rep(bdec[1], nb); B <- numeric(nb); C <- numeric(nb)
  for (t in terms) {
    cf <- bdec[t]
    if (t == "x1")        B <- B + cf
    else if (t == "x1x2") B <- B + cf * d$x2[idx]
    else if (t == "x1sq") C <- C + cf
    else                  A <- A + cf * g(t)
  }
  mu <- drop(cbind(d$x2[idx], d$x3[idx]) %*% W)
  lo <- pmin(mu, L) - 8 * SC
  gr <- lo + outer(L - lo, seq(0, 1, length.out = NGRID))  # nb x NGRID
  resid <- d$Y[idx] - (A + B*gr + C*gr^2)
  lp <- -0.5*((gr - mu)/SC)^2 - 0.5*(resid/sdraw)^2
  lp <- lp - apply(lp, 1, max)
  p  <- exp(lp); p <- p / rowSums(p)
  edge <- mean(p[, 1])
  cp <- t(apply(p, 1, cumsum))
  u  <- runif(nb)
  j  <- max.col(cp >= u, ties.method = "first")
  list(x = gr[cbind(seq_len(nb), j)], edge = edge)
}

fixed_point <- function(d, f, dn, iters = 25L, seed = 7L, proper = FALSE) {
  set.seed(seed)
  terms <- TERMS[[dn]]; L <- qnorm(f)
  cens <- which(d$x1 <= L); w <- d
  w$x1[cens] <- L - 0.5                       # the harness's cold start
  Mfull <- NULL; edge <- 0; path <- numeric(0)
  for (t in seq_len(iters)) {
    Md <- mm(w, terms)
    fit <- ols(Md, w$Y)
    r   <- w$Y - drop(Md %*% fit)
    sd_ <- sqrt(sum(r^2)/(nrow(Md) - ncol(Md)))
    if (proper) {
      sd_ <- sqrt(sum(r^2)/rchisq(1, nrow(Md) - ncol(Md)))
      V   <- chol2inv(chol(crossprod(Md)))
      fit <- fit + drop(t(chol(V)) %*% rnorm(length(fit))) * sd_
    }
    names(fit) <- colnames(Md)
    dr <- redraw_x1(w, cens, L, fit, terms, sd_)
    w$x1[cens] <- dr$x; edge <- max(edge, dr$edge)
    th <- ols(mm(w, TERMS$true), w$Y)
    path <- c(path, unname(th["x1sq"]))
  }
  th <- ols(mm(w, TERMS$true), w$Y)
  list(curv = unname(th["x1sq"]/BQUAD - 1), int = unname(th["x1x2"]/BINT - 1),
       edge = edge, path = path)
}

cat("\n================ STAGE 3: fixed point WITH redraw ================\n")
cat(sprintf("n = %d, reps = %d, grid = %d, 8 sweeps (converged by sweep 2)\n",
            NBIG, REPS, NGRID))
acc3 <- list()
for (r in seq_len(REPS)) {
  dr <- dgp(NBIG, 2000L + r)                 # one dataset per rep, reused
  for (f in FS) for (dn in c("noint","noquad","lin")) {
    fp <- fixed_point(dr, f, dn, iters = 8L, seed = 70L + r)
    k <- paste(f, dn)
    acc3[[k]] <- rbind(acc3[[k]],
                       c(f = f, curv = fp$curv, int = fp$int, edge = fp$edge))
  }
  cat(sprintf("  rep %d/%d done\n", r, REPS))
}
st3 <- do.call(rbind, lapply(names(acc3), function(k) {
  A <- acc3[[k]]; pr <- strsplit(k, " ")[[1]]
  data.frame(f = as.numeric(pr[1]), declared = pr[2],
             curv = mean(A[,"curv"]), curv_se = sd(A[,"curv"])/sqrt(nrow(A)),
             int  = mean(A[,"int"]),  int_se  = sd(A[,"int"])/sqrt(nrow(A)),
             grid_edge = max(A[,"edge"]))
}))
st3 <- st3[order(st3$f, match(st3$declared, c("noint","noquad","lin"))), ]
print(st3, row.names = FALSE, digits = 3)
if (max(st3$grid_edge) > 1e-6)
  cat("WARNING: grid lower edge carries mass", max(st3$grid_edge), "\n")

# ---- how many sweeps does the fixed point need? ----------------------------
# The harness runs THREE sweeps from a cold start at L-0.5, not to convergence.
# If three sweeps is far from the fixed point, stage 3 is answering a different
# question from the one V33 measured.
cat("\n---- sweeps to convergence (n =", NBIG, ", declared = noint) ----\n")
for (f in c(0.20, 0.40)) {
  fp <- fixed_point(dgp(NBIG, 3001L), f, "noint", iters = 12L, seed = 11L)
  cat(sprintf("f = %.2f  theta_quad by sweep: %s\n", f,
              paste(sprintf("%.4f", fp$path[c(1:6, 12)]), collapse = " ")))
}

# =============================================================================
# STAGE 4 -- harness-faithful: n = 800, m = 10, 3 sweeps, proper draw
# =============================================================================
# Everything V33 did EXCEPT the analysis model, which is OLS on the true basis
# rather than BKMR. If this reproduces V33's numbers, the whole phenomenon is
# linear algebra and BKMR is irrelevant to it; if it does not, the gap is BKMR's.
# Normalised the way V33's table is: (mean_arm / mean_oracle) - 1, so a
# multiplicative shrinkage common to arm and oracle cancels.
cat("\n================ STAGE 4: harness-faithful (n = 800, m = 10, 3 sweeps) ================\n")
N4 <- 800L; M4 <- 10L; SW4 <- 3L
R4 <- if (length(args) >= 3) as.integer(args[3]) else 400L

one_rep4 <- function(f, dn, seed) {
  d <- dgp(N4, seed)
  L <- qnorm(f); cens <- which(d$x1 <= L)
  th_or <- ols(mm(d, TERMS$true), d$Y)
  terms <- TERMS[[dn]]
  acc <- matrix(0, M4, 2)
  for (i in seq_len(M4)) {
    w <- d; w$x1[cens] <- L - 0.5
    for (t in seq_len(SW4)) {
      Md <- mm(w, terms)
      fit <- ols(Md, w$Y)
      r <- w$Y - drop(Md %*% fit)
      sd_ <- sqrt(sum(r^2)/rchisq(1, nrow(Md) - ncol(Md)))
      V <- chol2inv(chol(crossprod(Md)))
      fit <- fit + drop(t(chol(V)) %*% rnorm(length(fit))) * sd_
      names(fit) <- colnames(Md)
      w$x1[cens] <- redraw_x1(w, cens, L, fit, terms, sd_)$x
    }
    th <- ols(mm(w, TERMS$true), w$Y)
    acc[i, ] <- c(th["x1sq"], th["x1x2"])
  }
  c(colMeans(acc), th_or["x1sq"], th_or["x1x2"])
}

set.seed(4242)
st4 <- NULL
for (f in c(0.20, 0.40)) {
  for (dn in c("noint","noquad","lin")) {
    A <- matrix(NA_real_, R4, 4)
    for (r in seq_len(R4)) A[r, ] <- one_rep4(f, dn, 50000L + r)
    cm <- colMeans(A)
    # MC error on the ratio, by the delta method on independent-rep means
    se <- apply(A, 2, sd)/sqrt(R4)
    st4 <- rbind(st4, data.frame(
      f = f, declared = dn,
      curv = cm[1]/cm[3] - 1, curv_mcse = abs(cm[1]/cm[3]) *
        sqrt((se[1]/cm[1])^2 + (se[3]/cm[3])^2),
      int  = cm[2]/cm[4] - 1, int_mcse  = abs(cm[2]/cm[4]) *
        sqrt((se[2]/cm[2])^2 + (se[4]/cm[4])^2),
      oracle_curv_rel = cm[3]/BQUAD - 1))
  }
}
print(st4, row.names = FALSE, digits = 3)

# ---- the comparison table ---------------------------------------------------
meas <- data.frame(
  f = rep(c(0.20, 0.40), each = 3),
  declared = rep(c("noint","noquad","lin"), 2),
  m_curv = c( 0.0610, -0.1443, -0.1023,  0.2960, -0.3972, -0.4601),
  m_int  = c(-0.0487,  0.0553,  0.0556, -0.1317,  0.1009,  0.0271))
cmp <- merge(merge(meas, st3[st3$f %in% c(0.20,0.40), c("f","declared","curv","int")],
                   by = c("f","declared")),
             st4[, c("f","declared","curv","int")],
             by = c("f","declared"), suffixes = c("_fp","_n800"))
names(cmp)[names(cmp)=="curv"] <- "curv_n800"
names(cmp)[names(cmp)=="int"]  <- "int_n800"
cat("\n================ PREDICTED vs MEASURED (vs-oracle relative bias) ================\n")
cat("m_* = V33/V34 measured (BKMR analysis). *_fp = stage 3 asymptotic fixed point.\n")
cat("*_n800 = stage 4, harness-faithful except OLS analysis model.\n")
print(cmp[order(cmp$declared, cmp$f),
          c("declared","f","m_curv","curv_fp","curv_n800","m_int","int_fp","int_n800")],
      row.names = FALSE, digits = 3)

ex <- function(a, b) log(abs(b/a))/log(2)
cat("\n---- f-scaling exponent, f = 0.20 -> 0.40 (measured 1.85 +/- 0.20 pooled) ----\n")
for (dn in c("noint","noquad","lin")) {
  r2 <- cmp[cmp$declared==dn & cmp$f==0.20, ]; r4 <- cmp[cmp$declared==dn & cmp$f==0.40, ]
  cat(sprintf("%-7s curv: measured %.2f | fixed point %.2f | n800 %.2f\n", dn,
              ex(r2$m_curv, r4$m_curv), ex(r2$curv_fp, r4$curv_fp),
              ex(r2$curv_n800, r4$curv_n800)))
}

saveRDS(list(stage1 = st1, stage2 = st2, stage3 = st3, stage4 = st4, cmp = cmp),
        "results/derive_declared_basis.rds")
cat("\nwrote results/derive_declared_basis.rds\n")
