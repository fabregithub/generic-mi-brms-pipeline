# =============================================================================
# Demonstration cohort generator — JECS prenatal PFAS -> Kawasaki disease.
# -----------------------------------------------------------------------------
# SYNTHETIC data whose structure is grounded in three JECS papers (see
# ../PROVENANCE.md for the parameter-by-parameter mapping):
#   * Lai et al. 2025  — 7 PFAS distributions, MDLs, correlations, Z->X determinants
#   * Iwata et al. 2024 — Kawasaki disease outcome, single/WQS/BKMR effect sizes
#
# Produces INTERVAL-censored PFAS in the three-tier analytical structure:
#   * non-detect (ND):        x < MDL             -> interval (0, MDL]
#   * detected-not-quantified: MDL <= x < LCMRL   -> interval [MDL, LCMRL)
#   * quantified:             x >= LCMRL          -> exact value
# encoded as leftcens columns <pfas>_lo / <pfas>_hi.
#
# Scenario factors for the robustness study:
#   * censoring   = "real" | "increased"   (raise MDL/LCMRL to censor more)
#   * correlation = "real" | "low" | "high" (scale the off-diagonal PFAS correlations)
#
# Generating coefficients + uncensored concentrations are returned as `truth`, so
# every censored-data estimate can be scored against the complete-data oracle.
# Not real participant data. `make_demo_data()` returns list(data, truth, meta).
# =============================================================================

.make_pd_cor <- function(R, eps = 1e-3) {
  e <- eigen((R + t(R)) / 2, symmetric = TRUE)
  A <- e$vectors %*% diag(pmax(e$values, eps)) %*% t(e$vectors)
  d <- 1 / sqrt(diag(A)); A <- diag(d) %*% A %*% diag(d)
  (A + t(A)) / 2
}
.rmvn <- function(n, Sigma) matrix(stats::rnorm(n * ncol(Sigma)), n) %*% chol(Sigma)

#' @param n Sample size. Default 2000 for fast code verification; scale to ~25000
#'   (~Iwata, ~285 KD cases at ~1.1% prevalence) for the production study.
#' @param seed RNG seed.
#' @param mdl,lcmrl Optional named vectors (over the 7 PFAS) of real method
#'   detection / lowest-concentration-quantification limits (ng/mL). When NULL,
#'   they are derived from target detection/quantification rates. Supplying real
#'   values (and scaling them via `limit_scale`) is the intended way to build
#'   censoring-rate scenarios.
#' @param limit_scale Multiplies MDL and LCMRL. 1 = base limits; >1 censors more
#'   (the "increased-censoring" scenario), <1 censors less.
#' @param correlation "real" (Lai's R), "low" (halved), or "high" (inflated).
#' @return list(data, truth, meta).
make_demo_data <- function(n = 2000L, seed = 2026L,
                           mdl = NULL, lcmrl = NULL, limit_scale = 1,
                           correlation = c("real", "low", "high")) {
  correlation <- match.arg(correlation)
  set.seed(seed)

  # --- 7 PFAS: medians (ng/mL, Iwata), marginal log-SD, detection/quantification --
  pfas  <- c("PFOA", "PFOS", "PFNA", "PFDA", "PFUnDA", "PFHxS", "PFTrDA")
  group <- c(PFOA="PFCA", PFOS="PFSA", PFNA="PFCA", PFDA="PFCA", PFUnDA="PFCA", PFHxS="PFSA", PFTrDA="PFCA")
  med   <- c(PFOA=1.60, PFOS=2.90, PFNA=1.40, PFDA=0.49, PFUnDA=1.10, PFHxS=0.33, PFTrDA=0.15)
  sdlog <- c(PFOA=0.80, PFOS=0.85, PFNA=0.80, PFDA=0.80, PFUnDA=0.85, PFHxS=0.85, PFTrDA=0.85)
  # Limits: use supplied real MDL/LCMRL if given, else derive from target
  # detection/quantification rates. P(x > MDL) = detection; P(x > LCMRL) =
  # quantified (< detection => LCMRL > MDL, a detected-not-quantified band).
  det_real  <- c(PFOA=0.999, PFOS=0.99, PFNA=0.95, PFDA=0.93, PFUnDA=0.94, PFHxS=0.92, PFTrDA=0.785)
  quant_real<- c(PFOA=0.99,  PFOS=0.95, PFNA=0.88, PFDA=0.85, PFUnDA=0.86, PFHxS=0.85, PFTrDA=0.65)
  if (is.null(mdl))   mdl   <- med * exp(sdlog * stats::qnorm(1 - det_real))   # Lai range ~0.016-0.16
  if (is.null(lcmrl)) lcmrl <- med * exp(sdlog * stats::qnorm(1 - quant_real))
  mdl   <- mdl[pfas]   * limit_scale                        # scale for censoring scenarios
  lcmrl <- lcmrl[pfas] * limit_scale
  p <- length(pfas)

  # --- inter-PFAS correlation (log scale), Lai's Schober R values ---------------
  R <- diag(p); dimnames(R) <- list(pfas, pfas)
  sr <- function(a, b, r) { R[a, b] <<- r; R[b, a] <<- r }
  sr("PFNA","PFOA",0.74); sr("PFNA","PFDA",0.72); sr("PFUnDA","PFDA",0.81)
  sr("PFUnDA","PFTrDA",0.83); sr("PFOS","PFDA",0.79); sr("PFOS","PFUnDA",0.70)
  sr("PFOA","PFOS",0.55); sr("PFOA","PFDA",0.70); sr("PFOA","PFUnDA",0.65); sr("PFOA","PFTrDA",0.60)
  sr("PFOS","PFNA",0.60); sr("PFOS","PFTrDA",0.60); sr("PFNA","PFUnDA",0.75)
  sr("PFNA","PFTrDA",0.70); sr("PFDA","PFTrDA",0.78)
  for (a in pfas) sr("PFHxS", a, if (group[[a]] == "PFSA") 0.45 else 0.28)
  R["PFHxS","PFHxS"] <- 1
  cor_scale <- switch(correlation, real = 1, low = 0.5, high = 1.25)
  Roff <- R; diag(Roff) <- 0
  R <- diag(p) + pmin(Roff * cor_scale, 0.95); dimnames(R) <- list(pfas, pfas)
  R <- .make_pd_cor(R)

  # --- covariates Z (Lai Table 1) ---------------------------------------------
  regions <- c("Hokkaido","Miyagi","Fukushima","Chiba","Kanagawa","Koshin","Toyama",
               "Aichi","Kyoto","Osaka","Hyogo","Tottori","Kochi","Fukuoka","Okinawa")
  region_p <- c(8.4,8.7,12,5.7,5.8,7.1,5.8,5.9,3.8,7.9,4.9,3.0,6.9,7.6,6.3)
  region <- factor(sample(regions, n, TRUE, region_p / sum(region_p)), levels = regions)
  mat_age <- round(pmax(16, pmin(48, stats::rnorm(n, 31, 5)))); age_z <- (mat_age - 31) / 5
  parity  <- sample(0:3, n, TRUE, c(0.39, 0.39, 0.16, 0.06))
  bmi     <- pmax(14, stats::rnorm(n, 21, 3.3))
  smoking <- factor(sample(c("never","past","current"), n, TRUE, c(0.58, 0.375, 0.045)),
                    levels = c("never","past","current"))
  income  <- factor(sample(1:7, n, TRUE, c(0.06,0.34,0.33,0.16,0.07,0.03,0.01)), levels = 1:7)
  folate  <- stats::rnorm(n, 0, 1)
  egfr    <- stats::rnorm(n, 131, 24)
  child_male <- stats::rbinom(n, 1, 0.512); low_bw <- stats::rbinom(n, 1, 0.09)

  # --- Z -> X determinants (Lai): region strongest; parity down; age up; age x parity --
  region_eff <- matrix(stats::rnorm(length(regions) * p, 0, 0.25),
                       length(regions), p, dimnames = list(regions, pfas))
  b_parity <- 0.10; b_age <- 0.10; b_ap <- 0.04
  eps <- .rmvn(n, R)
  logPFAS <- matrix(0, n, p, dimnames = list(NULL, pfas))
  for (j in seq_len(p)) {
    mu <- log(med[j]) + region_eff[as.integer(region), j] -
      b_parity * parity + b_age * age_z + b_ap * parity * age_z
    logPFAS[, j] <- mu + sdlog[j] * eps[, j]
  }
  conc <- exp(logPFAS); colnames(conc) <- pfas

  # --- outcome Y: Kawasaki disease (rare), inverse PFAS mixture -----------------
  beta <- c(PFOA=log(0.894), PFOS=log(0.839), PFNA=-0.10, PFDA=-0.13,
            PFUnDA=-0.11, PFHxS=-0.15, PFTrDA=-0.14)
  log2c <- sweep(log2(conc), 2, log2(med), "-")
  lp_x  <- as.vector(log2c %*% beta)
  zbar  <- rowMeans(scale(log2c))
  lp_mix <- -0.06 * zbar^2 - 0.05 * as.vector(scale(log2c[,"PFOA"]) * scale(log2c[,"PFOS"]))
  lp_z  <- 0.35 * (smoking == "current") + 0.20 * child_male + 0.30 * low_bw + 0.08 * as.integer(income)
  lp    <- lp_x + lp_mix + lp_z
  alpha <- stats::uniroot(function(a) mean(stats::plogis(a + lp)) - 0.011, c(-12, 2))$root
  kd    <- stats::rbinom(n, 1, stats::plogis(alpha + lp))

  # --- second outcome Y2: childhood wheeze/asthma (Atagi) — COMMON & informative,
  #     so the Y-inclusion (congeniality) effect is demonstrable (KD is too rare).
  #     Atagi's real effects were near-null (aOR ~0.94/doubling); we use moderate
  #     inverse effects for a clean methods demonstration, not an asthma claim.
  mat_asthma <- stats::rbinom(n, 1, 0.10)                 # maternal asthma history (Atagi Model 3)
  beta_asthma <- c(PFOA=log(0.90), PFOS=log(0.93), PFNA=-0.06, PFDA=-0.08,
                   PFUnDA=-0.07, PFHxS=log(0.90), PFTrDA=-0.03)
  lp_a <- as.vector(log2c %*% beta_asthma) - 0.03 * zbar^2 +
    0.90 * mat_asthma + 0.30 * (smoking == "current") + 0.15 * child_male +
    0.20 * low_bw + 0.10 * (mat_age > 35)
  alpha_a <- stats::uniroot(function(a) mean(stats::plogis(a + lp_a)) - 0.22, c(-6, 4))$root
  asthma  <- stats::rbinom(n, 1, stats::plogis(alpha_a + lp_a))

  # --- three-tier interval censoring -> leftcens columns -----------------------
  df <- data.frame(row_id = seq_len(n), kd = kd, asthma = asthma)
  tier <- matrix("Q", n, p, dimnames = list(NULL, pfas))
  for (j in seq_len(p)) {
    x <- conc[, j]
    is_nd  <- x < mdl[j]
    is_dnq <- !is_nd & x < lcmrl[j]
    tier[is_nd, j]  <- "ND"; tier[is_dnq, j] <- "DNQ"
    df[[pfas[j]]]                <- ifelse(is_nd | is_dnq, NA, x)       # observed if quantified
    df[[paste0(pfas[j], "_lo")]] <- ifelse(is_nd, 0, ifelse(is_dnq, mdl[j],   x))
    df[[paste0(pfas[j], "_hi")]] <- ifelse(is_nd, mdl[j], ifelse(is_dnq, lcmrl[j], x))
  }
  df$region <- region; df$mat_age <- mat_age; df$parity <- parity; df$bmi <- bmi
  df$smoking <- smoking; df$income <- income; df$folate <- folate; df$egfr <- egfr
  df$child_male <- child_male; df$low_bw <- low_bw; df$mat_asthma <- mat_asthma
  df$income[sample(n, round(0.08 * n))] <- NA      # MAR (Lai) -> exercises the Z-block
  df$folate[sample(n, round(0.05 * n))] <- NA
  df$bmi[sample(n,    round(0.02 * n))] <- NA

  nd_rate  <- colMeans(tier == "ND")
  dnq_rate <- colMeans(tier == "DNQ")
  truth <- list(
    pfas = pfas, group = group, median = med, sdlog = sdlog,
    mdl = mdl, lcmrl = lcmrl, conc_true = conc, tier = tier,
    nd_rate = nd_rate, dnq_rate = dnq_rate, censored_rate = nd_rate + dnq_rate,
    cor_logPFAS = R,
    beta_per_log2 = list(kd = beta, asthma = beta_asthma),
    mixture = c(curvature = -0.06, PFOA_PFOS = -0.05),
    confounders = c(smoke_current = 0.35, child_male = 0.20, low_bw = 0.30, income = 0.08),
    intercept = c(kd = alpha, asthma = alpha_a),
    kd_prevalence = mean(kd), asthma_prevalence = mean(asthma),
    determinant = c(parity = -b_parity, age = b_age, age_parity = b_ap),
    scenario = c(limit_scale = limit_scale, correlation = correlation))
  meta <- list(n = n, seed = seed, n_kd = sum(kd), limit_scale = limit_scale,
               correlation = correlation, sources = "Lai 2025 (X); Iwata 2024 (Y~X+Z)")
  list(data = df, truth = truth, meta = meta)
}

# ---- run as a script: write one scenario's data + truth --------------------
if (sys.nframe() == 0L) {
  n   <- if (nzchar(Sys.getenv("N"))) as.integer(Sys.getenv("N")) else 2000L
  sc  <- if (nzchar(Sys.getenv("LIMIT_SCALE"))) as.numeric(Sys.getenv("LIMIT_SCALE")) else 1
  cor <- if (nzchar(Sys.getenv("CORRELATION"))) Sys.getenv("CORRELATION") else "real"
  out <- make_demo_data(n = n, limit_scale = sc, correlation = cor)
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  tag <- sprintf("scale-%.2g_cor-%s", sc, cor)
  saveRDS(out$data,  sprintf("data/jecs_pfas_kd_%s.rds", tag), compress = FALSE)
  saveRDS(out$truth, sprintf("data/jecs_pfas_kd_%s_truth.rds", tag), compress = FALSE)
  message(sprintf("Saved scenario [%s]: n=%d, KD=%d (%.2f%%)", tag, nrow(out$data),
                  out$meta$n_kd, 100 * out$truth$kd_prevalence))
  message("Censored (ND+DNQ) rate: ",
          paste(sprintf("%s %.0f%%(ND %.0f/DNQ %.0f)", out$truth$pfas,
                        100*out$truth$censored_rate, 100*out$truth$nd_rate, 100*out$truth$dnq_rate),
                collapse=", "))
}
