# Is SIR-over-leftcens salvageable at all?
#
# HYPOTHESIS. SIR can only ADD a factor to a proposal; it cannot repair a
# proposal that is itself misspecified, because that would need the full
# target/proposal ratio and leftcens does not expose its density.
#
# Under pipe_nl the proposal is wrong EITHER WAY:
#   child IN  -- (x, Y, Z1) is not jointly Gaussian, since Z1 = g(x) + noise
#   child OUT -- marginalising Z1 makes Y NON-LINEAR in x (Y = b1x + g1*g(x) +...),
#                so leftcens's linear regression of x on Y is misspecified
# The harness grid draw avoids both because its outcome factor conditions on Z1,
# where Y IS linear in x.
#
# THE TEST. Run the same SIR on the LINEAR pipe, where the proposal (child out)
# IS correctly specified. If SIR lands at ~0 there and only fails under pipe_nl,
# the hypothesis holds and the SIR composition cannot implement this feature.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","dag_draw.R")) source(file.path("R",f))})
sir <- function(d, truth, m=20L, sweeps=3L, K=100L, DROP=character(0)) {
  xc <- paste0("logX", seq_len(length(truth$b))); zc <- c("Z1","Z2")
  isc <- is.na(d$logX1); lod <- d$logX1_lod
  nl <- identical(truth$z_role, "pipe_nl")
  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d; w$logX1[isc] <- lod[isc] - 0.5
    for (t in seq_len(sweeps)) {
      b <- list(y = ifelse(isc, NA_real_, d$logX1),
                lower = ifelse(isc, -Inf, NA_real_), upper = ifelse(isc, lod, NA_real_))
      xmat <- w[, c("Y", xc[-1], setdiff(zc, DROP)), drop=FALSE]
      cand <- as.matrix(leftcens::impute_censored_conditional(
        y=b$y, x=xmat, lower=b$lower, upper=b$upper, m=K,
        margin="gaussian"))[isc, , drop=FALSE]
      mu <- if (nl) ef_nl_g(cand, truth) else truth$delta_xz * cand
      lw <- dnorm(w$Z1[isc], mu, if (nl) truth$nl_sd else 0.8, log=TRUE)
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
for (role in c("pipe","pipe_nl")) {
  tr <- make_truth(z_role=role); fo <- dgp_formula(tr); tv <- tr$estimand_true
  res <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    set.seed(4000+r)
    cn <- inject_left_censoring(simulate_complete(N,tr)$data, nd_frac=0.4, censor_which=1L)
    set.seed(8000+r); a <- pool(sir(cn,tr,DROP="Z1"), fo)
    set.seed(8000+r); g <- pool(dag_impute_datasets(cn,tr,m=20L,child_form="true"), fo)
    data.frame(s=a[["est"]], ss=a[["se"]], g=g[["est"]], gs=g[["se"]])
  }, mc.cores=20))
  f <- function(e,s) sprintf("%+7.2f%% (bias/SE %+5.2f)", 100*(mean(e)-tv)/tv, (mean(e)-tv)/mean(s))
  cat(sprintf("%-8s | SIR (child out of proposal) %s | grid draw %s\n",
              role, f(res$s,res$ss), f(res$g,res$gs)))
}
