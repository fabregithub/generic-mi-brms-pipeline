# Is the SHIPPED draw class the same as V24's dag_Yonly? It should not be:
# leftcens regresses x on (Y, Z, other X) -- a single conditional that USES Z --
# whereas dag_Yonly conditions on parents only and multiplies by p(Y|x), which
# under a pipe drops Z entirely. Under a LINEAR pipe (x, Y, Z) are jointly
# Gaussian, so the shipped single regression may be the EXACT conditional.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","dag_draw.R")) source(file.path("R",f))})

# The shipped form: interval-censored Gaussian regression of x on EVERYTHING,
# then a truncated draw. This is what leftcens does (its extra machinery is the
# shash margin, which is the identity here since the DGP is Gaussian).
draw_shipped <- function(d, truth, m = 20L, sweeps = 3L) {
  xc <- paste0("logX", seq_len(length(truth$b)))
  zc <- intersect(c("Z1","Z2"), names(d))
  isc <- is.na(d$logX1); lod <- d$logX1_lod
  pr <- c("Y", xc[-1], zc)
  out <- vector("list", m)
  for (i in seq_len(m)) {
    w <- d; w$logX1[isc] <- lod[isc] - 0.5
    for (t in seq_len(sweeps)) {
      lo <- w$logX1; hi <- w$logX1; lo[isc] <- -Inf; hi[isc] <- lod[isc]
      f <- survreg(as.formula(paste("Surv(lo,hi,type='interval2') ~",
                                    paste(pr, collapse="+"))),
                   data = cbind(w, lo=lo, hi=hi), dist = "gaussian")
      mu <- predict(f, newdata = w[isc, , drop=FALSE], type="response"); s <- f$scale
      u  <- runif(sum(isc)) * pnorm(lod[isc], mu, s)
      w$logX1[isc] <- qnorm(pmax(u, 1e-300), mu, s)
    }
    out[[i]] <- w
  }
  out
}
pool <- function(imps, fo) {
  e <- sapply(imps, function(w) coef(lm(fo, data=w))[["logX1"]])
  v <- sapply(imps, function(w) diag(vcov(lm(fo, data=w)))[["logX1"]])
  mm <- length(e); c(est=mean(e), se=sqrt(mean(v) + (1+1/mm)*var(e)))
}
N <- as.integer(Sys.getenv("N","3200")); R <- as.integer(Sys.getenv("R","100"))
cases <- list(c("pipe","direct"), c("pipe","total"),
              c("pipe_nl","direct"), c("pipe_nl","total"), c("fork","direct"))
for (cs in cases) {
  tr <- make_truth(z_role=cs[1], target=cs[2]); fo <- dgp_formula(tr)
  res <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    set.seed(4000+r)
    cm <- simulate_complete(N, tr)$data
    cn <- inject_left_censoring(cm, nd_frac=0.4, censor_which=1L)
    set.seed(5000+r); a <- pool(draw_shipped(cn, tr), fo)
    set.seed(5000+r); b <- pool(dag_impute_datasets(cn, tr, m=20L, use_child=FALSE), fo)
    set.seed(5000+r); c3 <- pool(dag_impute_datasets(cn, tr, m=20L, child_form="true"), fo)
    data.frame(sh=a[["est"]], sh_se=a[["se"]], yo=b[["est"]], yo_se=b[["se"]],
               dg=c3[["est"]], dg_se=c3[["se"]])
  }, mc.cores=20))
  tv <- tr$estimand_true
  f <- function(e, se) sprintf("%+7.2f%% (bias/SE %+5.2f)", 100*(mean(e)-tv)/tv,
                               (mean(e)-tv)/mean(se))
  cat(sprintf("%-8s %-6s truth %.4f | SHIPPED-form %s | dag_Yonly %s | DAG %s\n",
              cs[1], cs[2], tv, f(res$sh,res$sh_se), f(res$yo,res$yo_se), f(res$dg,res$dg_se)))
}
