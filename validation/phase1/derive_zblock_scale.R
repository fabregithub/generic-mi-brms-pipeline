# THE HYPOTHESIS. Both pipeline arms sit 3-5 pp above their harness twins, with
# two different Z engines -- so the cause is the wrapper, not the draw. Located
# candidate: the pipeline's working column holds the exposure on the DATA scale
# between sweeps (00_censored_exposure.R:442 exponentiates the draw back), so
# its Z block regresses Z1 on exp(logX1). The harness regresses on logX1. The
# truth is linear in logX1, so the pipeline's Z model is MISSPECIFIED in scale.
#
# TEST: give the harness's Z1 draw exp(logX1) instead of logX1, change nothing
# else, and see whether the 3-5 pp appears.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","exact_fork.R","robustness.R")) source(file.path("R",f))})
`%||%` <- function(a,b) if (is.null(a)) b else a

draw_z1 <- function(w, idx, obs, raw_scale) {
  xc <- grep("^logX[0-9]+$", names(w), value = TRUE)
  ww <- w
  if (raw_scale) ww$logX1 <- exp(w$logX1)          # the pipeline's scale
  fo <- stats::as.formula(paste("Z1 ~", paste(c(xc, "Z2", "Y"), collapse = "+")))
  f <- lm(fo, data = ww[obs, , drop = FALSE])
  mm <- model.matrix(fo, data = ww)
  bh <- coef(f); V <- vcov(f)
  bh <- as.vector(bh + t(chol(V)) %*% rnorm(length(bh)))
  s2 <- sum(residuals(f)^2)/rchisq(1, f$df.residual)
  out <- w$Z1
  out[idx] <- as.vector(mm[idx, , drop=FALSE] %*% bh) + rnorm(length(idx), 0, sqrt(s2))
  out
}
imp <- function(d, truth, raw_scale, m = 20L, sweeps = 3L) {
  xc <- paste0("logX", seq_len(length(truth$b)))
  is_c <- is.na(d$logX1); lod <- d$logX1_lod
  z1na <- is.na(d$Z1); z2na <- is.na(d$Z2)
  idx1 <- which(z1na); obs1 <- which(!z1na)
  J <- .ef_joint(truth)
  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d; w$logX1[is_c] <- lod[is_c] - 0.5
    if (any(z1na)) w$Z1[z1na] <- mean(d$Z1, na.rm = TRUE)
    if (any(z2na)) w$Z2[z2na] <- round(mean(d$Z2, na.rm = TRUE))
    for (t in seq_len(sweeps)) {
      if (length(idx1)) w$Z1 <- draw_z1(w, idx1, obs1, raw_scale)
      if (any(z2na)) w$Z2 <- ef_draw_z2(w, truth, which(z2na), J)
      pr <- c("Y", setdiff(xc, "logX1"), "Z1", "Z2")
      yv <- w$logX1; yv[is_c] <- NA_real_
      lo <- rep(NA_real_, nrow(w)); hi <- rep(NA_real_, nrow(w))
      lo[is_c] <- -Inf; hi[is_c] <- lod[is_c]
      w$logX1 <- leftcens::impute_censored_conditional(
        y = yv, x = w[, pr, drop = FALSE], lower = lo, upper = hi,
        m = 1L, margin = "shash")[, 1]
    }
    out[[i]] <- w
  }
  out
}
pool <- function(imps, fo) {
  e <- sapply(imps, function(w) coef(lm(fo, data=w))[["logX1"]])
  v <- sapply(imps, function(w) diag(vcov(lm(fo, data=w)))[["logX1"]])
  c(est = mean(e), se = sqrt(mean(v) + (1+1/length(e))*var(e)))
}
R <- as.integer(Sys.getenv("R","150"))
for (cell in c("zr_fork","zr_pipe")) {
  sc <- Filter(function(s) s$name == cell, v2_scenarios())[[1]]
  tr <- make_truth(p=3L, erf_form=sc$erf_form, z_role=sc$z_role)
  fo <- dgp_formula(tr); tv <- tr$estimand_true
  res <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    set.seed(41000+r); b <- v2_make_bundle(sc, tr)
    set.seed(42000+r); a <- pool(imp(b$censored, tr, raw_scale = FALSE), fo)
    set.seed(42000+r); c2 <- pool(imp(b$censored, tr, raw_scale = TRUE), fo)
    data.frame(lg=a[["est"]], lgs=a[["se"]], rw=c2[["est"]], rws=c2[["se"]])
  }, mc.cores = 12))
  f <- function(e,s) sprintf("%+6.2f%% (bias/SE %+5.2f)", 100*(mean(e)-tv)/tv, (mean(e)-tv)/mean(s))
  cat(sprintf("%-8s Z1 on logX1 %s | Z1 on exp(logX1) %s | gap %+.2f pp\n",
              cell, f(res$lg,res$lgs), f(res$rw,res$rws),
              100*(mean(res$rw)-mean(res$lg))/tv))
}
cat("\nmeasured pipeline-vs-harness gap to explain: +2.96 / +3.39 pp (micePmm),\n")
cat("+4.25 / +5.26 pp (properZ)\n")
