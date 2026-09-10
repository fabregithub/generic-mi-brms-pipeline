# THE DECIDING QUESTION FOR THE PIPELINE FEATURE.
# V24 tested child models that were either a fitted POLYNOMIAL (all worse than
# omitting the factor) or the TRUE arrow with its true coefficients (exact). A
# shipped feature can never have the second. So: is declaring the correct BASIS
# and ESTIMATING its coefficients from the observed rows enough?
#
# If yes, the interface is "declare the shape, we fit it" and the feature can
# help. If no, the feature cannot help anyone and must not ship.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","dag_draw.R")) source(file.path("R",f))})

# g(x) = nl_a*tanh(nl_b*x) + nl_c*(x^2-1), with nl_b = 1.8 INSIDE the tanh.
# Four declared bases, in decreasing order of how much the analyst knows:
BASES <- list(
  exact_basis = ~ tanh(1.8 * .xm) + I(.xm^2),   # right basis, coefficients fitted
  tanh_scale1 = ~ tanh(.xm) + I(.xm^2),         # right FORM, wrong inner scale
  tanh_only   = ~ tanh(1.8 * .xm),              # right basis, one term missing
  linear      = ~ .xm                            # V24's dag_lin, for reference
)

fit_declared <- function(child, xm, form) {
  d <- data.frame(.child = child, .xm = xm)
  d <- d[!is.na(d$.xm) & !is.na(d$.child), ]
  f <- stats::lm(stats::update(form, .child ~ .), data = d)
  list(mu = function(x, yv = 0) {          # dag_draw_x1 calls mu(G, sub$Y)
         nd <- data.frame(.xm = as.vector(x))
         matrix(as.vector(stats::predict(f, newdata = nd)), nrow = nrow(as.matrix(x)))
       },
       sd = stats::summary.lm(f)$sigma)
}

# Same grid draw as dag_draw.R, but with a DECLARED-basis child model.
draw_declared <- function(d, truth, form, m = 20L, sweeps = 3L, ngrid = 400L) {
  spec <- dag_spec(truth)
  xc <- paste0("logX", seq_len(length(truth$b))); zc <- intersect(c("Z1","Z2"), names(d))
  isc <- is.na(d$logX1); lod <- d$logX1_lod
  zy <- if (spec$z1_in_y) zc else setdiff(zc, "Z1")
  fo_y <- as.formula(paste("Y ~", paste(c(xc, zy), collapse="+")))
  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d; w$logX1[isc] <- lod[isc] - 0.5
    for (t in seq_len(sweeps)) {
      fy <- lm(fo_y, data = w); dfr <- fy$df.residual
      s2 <- sum(residuals(fy)^2)/rchisq(1, dfr)
      cf <- as.vector(coef(fy) + t(chol(summary(fy)$cov.unscaled)) %*%
                        rnorm(length(coef(fy))) * sqrt(s2))
      names(cf) <- names(coef(fy))
      xm <- w$logX1; xm[isc] <- NA                    # observed-only for the child fit
      ch <- fit_declared(w$Z1, xm, form)
      w$logX1 <- as.vector(dag_draw_x1(w, truth, which(isc), lod, spec, ch, cf,
                                       sqrt(s2), obs_na = isc, ngrid = ngrid))
    }
    out[[i]] <- w
  }
  out
}
pool <- function(imps, fo) {
  e <- sapply(imps, function(w) coef(lm(fo, data=w))[["logX1"]])
  v <- sapply(imps, function(w) diag(vcov(lm(fo, data=w)))[["logX1"]])
  c(est=mean(e), se=sqrt(mean(v) + (1+1/length(e))*var(e)))
}
N <- as.integer(Sys.getenv("N","3200")); R <- as.integer(Sys.getenv("R","100"))
for (tg in c("direct","total")) {
  tr <- make_truth(z_role="pipe_nl", target=tg); fo <- dgp_formula(tr); tv <- tr$estimand_true
  res <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    set.seed(4000+r)
    cm <- simulate_complete(N, tr)$data
    cn <- inject_left_censoring(cm, nd_frac=0.4, censor_which=1L)
    row <- list()
    for (nm in names(BASES)) {
      set.seed(6000+r); p <- pool(draw_declared(cn, tr, BASES[[nm]]), fo)
      row[[paste0(nm,"_e")]] <- p[["est"]]; row[[paste0(nm,"_s")]] <- p[["se"]]
    }
    set.seed(6000+r); p <- pool(dag_impute_datasets(cn, tr, m=20L, use_child=FALSE), fo)
    row$none_e <- p[["est"]]; row$none_s <- p[["se"]]
    as.data.frame(row)
  }, mc.cores=20))
  cat(sprintf("\n=== pipe_nl, %s effect (truth %.4f, n=%d, %d reps) ===\n", tg, tv, N, R))
  for (nm in c("none", names(BASES))) {
    e <- res[[paste0(nm,"_e")]]; s <- res[[paste0(nm,"_s")]]
    cat(sprintf("  %-12s %+7.2f%%  bias/SE %+6.2f\n", nm,
                100*(mean(e)-tv)/tv, (mean(e)-tv)/mean(s)))
  }
}
