# V29 -- does s^2 drive the chain's residual, with the arrow held fixed?
#
# See run_v29_sdx.sh and R/robustness.R for the registered prediction (P1/P2/P3).
# Measured effect is PAIRED: pipeline_properZ vs pipeline_properZ_ds differ only
# by the pre-v1.6.0 scale defect and run on byte-identical data.
#
# DERIVED IS EVALUATED AT N = 300000, NOT SMALLER. The closed form is an
# asymptotic statement evaluated on a finite sample, and its finite-N bias is
# role-dependent: at N = 20000 it reads +0.87 pp high for a pipe against
# +0.27 pp for a fork, which manufactures a fork/pipe contrast out of the
# evaluation rather than the mechanism. At N = 300000 it reproduces
# derive_scale_propagation.R's +5.44 / +5.04. (Cost me a wrong reading of V28.)
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f))})
local({ ex <- parse("derive_fork_residual.R")
  for (e in ex) if (is.call(e) && identical(as.character(e[[1]]),"<-") && is.name(e[[2]]) &&
                    as.character(e[[2]]) %in% c("closed_form","oth","%||%")) eval(e, envir=globalenv()) })
stopifnot(is.function(closed_form))

cells <- list(
  sxf07=list(role="fork", sdx=0.7), sxf10=list(role="fork", sdx=1.0),
  sxf13=list(role="fork", sdx=1.3), sxp07=list(role="pipe", sdx=0.7),
  sxp10=list(role="pipe", sdx=1.0), sxp13=list(role="pipe", sdx=1.3))

raw <- read.csv("results/v29_n12800_raw.csv")
meas <- do.call(rbind, lapply(names(cells), function(cn) {
  a <- raw[raw$scenario==cn & raw$procedure=="pipeline_properZ", ]
  b <- raw[raw$scenario==cn & raw$procedure=="pipeline_properZ_ds", ]
  k <- intersect(a$rep, b$rep); a <- a[match(k,a$rep),]; b <- b[match(k,b$rep),]
  d <- 100*(b$estimate - a$estimate)/a$estimand_true[1]
  data.frame(cell=cn, role=cells[[cn]]$role, sdx=cells[[cn]]$sdx, n_rep=length(k),
             meas=mean(d), meas_se=sd(d)/sqrt(length(k)))
}))

der <- do.call(rbind, lapply(names(cells), function(cn) {
  cf <- cells[[cn]]
  v <- vapply(1:4, function(i) { set.seed(4000L+97L*i)
    tr <- make_truth(p=3L, erf_form="additive", y_form="linear", z_role=cf$role, target="direct",
                     delta_xz=if (cf$role=="pipe") 0.6 else NULL,
                     delta_zx=if (cf$role=="fork") 0.6 else NULL)
    s <- simulate_complete(300000L, tr, rho=0.4, skew=0.0, sd_x=cf$sdx, z_form="linear")
    closed_form(s$data, s$truth, mcar=0.4, ndf=0.4)$pred }, numeric(1))
  s2 <- { set.seed(11)
    tr <- make_truth(p=3L, erf_form="additive", y_form="linear", z_role=cf$role, target="direct",
                     delta_xz=if (cf$role=="pipe") 0.6 else NULL,
                     delta_zx=if (cf$role=="fork") 0.6 else NULL)
    var(simulate_complete(60000L, tr, rho=0.4, skew=0.0, sd_x=cf$sdx)$data$logX1) }
  data.frame(cell=cn, s2=s2, der=mean(v), der_se=sd(v)/sqrt(4))
}))

x <- merge(meas, der, by="cell")
x$ratio <- x$der/x$meas
x$ratio_se <- with(x, abs(ratio)*sqrt((der_se/der)^2 + (meas_se/meas)^2))
x <- x[order(x$role, x$sdx), ]
cat("\n================ V29  n=12800  (arrow fixed at 0.60) ================\n")
print(format(x[,c("cell","role","sdx","s2","n_rep","meas","meas_se","der","der_se","ratio","ratio_se")],
             digits=3), row.names=FALSE)
cat("\n")
for (rr in c("fork","pipe")) {
  s <- x[x$role==rr,]; s <- s[order(s$sdx),]
  fit <- lm(ratio ~ s2, data=s, weights=1/s$ratio_se^2)
  sl <- summary(fit)$coefficients["s2", ]
  cat(sprintf("%-5s ratio %.3f -> %.3f over s^2 %.2f -> %.2f | slope %+.4f per unit s^2 (SE %.4f, t=%.2f)\n",
              rr, s$ratio[1], s$ratio[nrow(s)], s$s2[1], s$s2[nrow(s)], sl[1], sl[2], sl[1]/sl[2]))
}
cat("\nVERDICT: P1 = both roles trend up; P2 = fork only; P3 = neither.\n")
