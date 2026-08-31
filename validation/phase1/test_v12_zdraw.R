#!/usr/bin/env Rscript
# =============================================================================
# V12 -- self-test for the Z-block samplers (no model fitting, ~2 s)
# -----------------------------------------------------------------------------
# The V12 contrast only means anything if:
#   * "exact" really samples p(Z | Y, rest), shape included, and
#   * "gaussian" really has the SAME MEAN but homoscedastic Gaussian noise.
# If "gaussian" also got the mean wrong, a difference between the two would
# confound shape with mean and the experiment would answer nothing.
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error=function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly=FALSE); f <- sub("^--file=","",a[grep("^--file=",a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.here,"R","dgp.R")); source(file.path(.here,"R","smc_impute.R"))
fails <- 0L
ok <- function(l,p,d="") { cat(sprintf("  [%s] %s%s\n", if(p)"PASS" else "FAIL", l,
                            if(nzchar(d)) paste0("  --  ",d) else "")); if(!p) fails <<- fails+1L }

ref <- function(Y,A,B,C,sig,m,s,N=400001L){
  g <- seq(m-10*s, m+10*s, length.out=N)
  lp <- -(Y-(A+B*g+C*g^2))^2/(2*sig^2) - (g-m)^2/(2*s^2)
  w <- exp(lp-max(lp)); w <- w/sum(w)
  mu <- sum(w*g); v <- sum(w*(g-mu)^2)
  list(mean=mu, sd=sqrt(v), skew=sum(w*(g-mu)^3)/v^1.5)
}

cat("\n=== 1. \"exact\" reproduces the true conditional (mean, SD and SKEW) ===\n\n")
set.seed(4); n <- 200000L
cases <- list(c(Y=-1,C=0.4), c(Y=0,C=0.4), c(Y=1,C=0.4), c(Y=2,C=0.4), c(Y=0.5,C=0.0))
for (cs in cases) {
  Y <- unname(cs["Y"]); C <- unname(cs["C"]); A <- -C; B <- 0.5
  d <- .smc_draw_z(rep(Y,n), rep(A,n), rep(B,n), rep(C,n), 1, rep(0,n), 1, mode="exact")
  r <- ref(Y,A,B,C,1,0,1)
  ok(sprintf("Y=%+.1f C=%.1f", Y, C),
     abs(mean(d)-r$mean)<0.01 && abs(sd(d)-r$sd)<0.01 && abs(e<-(mean((d-mean(d))^3)/sd(d)^3)-r$skew)<0.05,
     sprintf("mean %+.3f/%+.3f  sd %.3f/%.3f  skew %+.3f/%+.3f",
             mean(d),r$mean, sd(d),r$sd, mean((d-mean(d))^3)/sd(d)^3, r$skew))
}

cat("\n=== 2. \"gaussian\" has the SAME MEAN but loses the shape ===\n")
cat("    (this is what makes the contrast attributable to shape alone)\n\n")
for (Y in c(-1,0,1,2)) {
  C <- 0.4; A <- -C; B <- 0.5
  de <- .smc_draw_z(rep(Y,n), rep(A,n), rep(B,n), rep(C,n), 1, rep(0,n), 1, mode="exact")
  dg <- .smc_draw_z(rep(Y,n), rep(A,n), rep(B,n), rep(C,n), 1, rep(0,n), 1, mode="gaussian")
  ok(sprintf("Y=%+.1f  means agree", Y), abs(mean(de)-mean(dg))<0.01,
     sprintf("exact %+.3f vs gaussian %+.3f", mean(de), mean(dg)))
  ok(sprintf("Y=%+.1f  gaussian has ~zero skew", Y),
     abs(mean((dg-mean(dg))^3)/sd(dg)^3) < 0.05,
     sprintf("exact skew %+.3f vs gaussian %+.3f",
             mean((de-mean(de))^3)/sd(de)^3, mean((dg-mean(dg))^3)/sd(dg)^3))
}

cat("\n=== 3. under a LINEAR outcome the two must coincide ===\n")
cat("    (C = 0 makes the true conditional Gaussian, so there is nothing to lose)\n\n")
Y <- 0.7; A <- 0; B <- 0.5; C <- 0
de <- .smc_draw_z(rep(Y,n), rep(A,n), rep(B,n), rep(C,n), 1, rep(0,n), 1, mode="exact")
dg <- .smc_draw_z(rep(Y,n), rep(A,n), rep(B,n), rep(C,n), 1, rep(0,n), 1, mode="gaussian")
ok("means agree",  abs(mean(de)-mean(dg))<0.01, sprintf("%+.4f vs %+.4f", mean(de), mean(dg)))
ok("SDs agree",    abs(sd(de)-sd(dg))<0.01,     sprintf("%.4f vs %.4f", sd(de), sd(dg)))
ok("both ~symmetric", abs(mean((de-mean(de))^3)/sd(de)^3)<0.05)

cat("\n=== 4. binary conditional matches a direct Bayes calculation ===\n\n")
set.seed(9)
Y <- rnorm(n, 0, 1.5); e0 <- rep(0.2,n); e1 <- rep(0.9,n)
z <- .smc_draw_z_binary(Y, e0, e1, 1)
p <- 1/(1+exp((-(Y-e0)^2/2 + log(0.5)) - (-(Y-e1)^2/2 + log(0.5))))
ok("P(Z2=1) matches", abs(mean(z)-mean(p))<0.005, sprintf("drawn %.4f vs analytic %.4f", mean(z), mean(p)))

cat(sprintf("\n%s  (%d failed)\n\n", if(fails==0L)"ALL CHECKS PASSED" else "FAILURES", fails))
quit(status=if(fails==0L) 0L else 1L)
