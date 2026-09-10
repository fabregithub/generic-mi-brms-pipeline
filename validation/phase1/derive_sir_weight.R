# THE WEIGHT MUST CONDITION ON Y.
# leftcens's proposal already conditions on Y, so
#   target   ~ p(x|Pa) . p(Y|x,C) . p(C|x) . 1{x<=L}
#   proposal ~ p(x|Pa) . p(Y|x)   . 1{x<=L}
# and since p(Y|x,C).p(C|x) = p(Y|x).p(C|Y,x), the ratio is p(C | x, Y) -- NOT
# p(C | x). Using p(C|x) double-counts nothing but MISSES Y, which is why the
# linear pipe (where the proposal is exactly right) still came out at +7.63%.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","dag_draw.R")) source(file.path("R",f))})

sir <- function(d, truth, m=20L, sweeps=3L, K=100L, with_y=TRUE) {
  xc <- paste0("logX", seq_len(length(truth$b)))
  isc <- is.na(d$logX1); lod <- d$logX1_lod
  nl <- identical(truth$z_role, "pipe_nl")
  # The declared basis, as the pipeline would take it (coefficients fitted).
  rhs <- if (nl) "tanh(1.8 * .x) + I(.x^2)" else ".x"
  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d; w$logX1[isc] <- lod[isc] - 0.5
    for (t in seq_len(sweeps)) {
      b <- list(y = ifelse(isc, NA_real_, d$logX1),
                lower = ifelse(isc, -Inf, NA_real_), upper = ifelse(isc, lod, NA_real_))
      xmat <- w[, c("Y", xc[-1], "Z2"), drop=FALSE]           # Z1 (the child) out
      cand <- as.matrix(leftcens::impute_censored_conditional(
        y=b$y, x=xmat, lower=b$lower, upper=b$upper, m=K,
        margin="gaussian"))[isc, , drop=FALSE]
      # Child model fitted on OBSERVED-exposure rows, with or without Y.
      df <- data.frame(.child=w$Z1, .x=b$y, .y=w$Y)
      df <- df[!is.na(df$.x), ]
      fo <- as.formula(paste(".child ~", rhs, if (with_y) "+ .y" else ""))
      f <- lm(fo, data=df)
      nd <- data.frame(.x=as.vector(cand))
      if (with_y) nd$.y <- rep_len(w$Y[isc], nrow(nd))   # column-major recycling
      mu <- matrix(predict(f, newdata=nd), nrow=nrow(cand))
      lw <- dnorm(w$Z1[isc], mu, sigma(f), log=TRUE)
      lw <- lw - apply(lw,1,max); W <- exp(lw)
      pick <- vapply(seq_len(nrow(W)), function(r) sample.int(ncol(W),1L,prob=W[r,]), 1L)
      w$logX1[isc] <- cand[cbind(seq_len(nrow(cand)), pick)]
    }
    out[[i]] <- w
  }
  out
}
pool <- function(imps, fo) {
  e <- sapply(imps, function(w) coef(lm(fo,data=w))[["logX1"]])
  v <- sapply(imps, function(w) diag(vcov(lm(fo,data=w)))[["logX1"]])
  c(est=mean(e), se=sqrt(mean(v)+(1+1/length(e))*var(e)))
}
N <- as.integer(Sys.getenv("N","3200")); R <- as.integer(Sys.getenv("R","80"))
for (cs in list(c("pipe","direct"), c("pipe_nl","direct"), c("pipe_nl","total"))) {
  tr <- make_truth(z_role=cs[1], target=cs[2]); fo <- dgp_formula(tr); tv <- tr$estimand_true
  res <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    set.seed(4000+r)
    cn <- inject_left_censoring(simulate_complete(N,tr)$data, nd_frac=0.4, censor_which=1L)
    set.seed(8500+r); a <- pool(sir(cn,tr,with_y=FALSE), fo)
    set.seed(8500+r); b <- pool(sir(cn,tr,with_y=TRUE),  fo)
    set.seed(8500+r); g <- pool(dag_impute_datasets(cn,tr,m=20L,child_form="true"), fo)
    data.frame(n=a[["est"]],ns=a[["se"]], y=b[["est"]],ys=b[["se"]], g=g[["est"]],gs=g[["se"]])
  }, mc.cores=20))
  f <- function(e,s) sprintf("%+7.2f%% (%+5.2f)", 100*(mean(e)-tv)/tv, (mean(e)-tv)/mean(s))
  cat(sprintf("%-8s %-6s | weight p(C|x) %s | weight p(C|x,Y) %s | grid %s\n",
              cs[1], cs[2], f(res$n,res$ns), f(res$y,res$ys), f(res$g,res$gs)))
}
