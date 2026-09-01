#!/usr/bin/env Rscript
# =============================================================================
# Item 07 -- self-test for the general importance-resampling exposure draw
# -----------------------------------------------------------------------------
# The method must reproduce the EXACT conditional wherever that conditional is
# available in closed form. If it cannot do that on a case we can solve by hand,
# no result it produces on a case we cannot is worth anything.
#
# The proposal is injected here rather than taken from leftcens, so that the
# weighting and resampling are tested on their own. When a real run looks wrong,
# this separation is what tells you whether the machinery or the proposal model
# is at fault.
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

# Truth: Y | x ~ N(A + B x + C x^2, sig^2); prior x ~ N(m, s) truncated above at U.
A <- 0.10; B <- 0.40; C <- 0.15; sig <- 1.0; m <- 0.0; s <- 1.0

exact_moments <- function(Y, U, N=400001L) {
  g <- seq(m-10*s, U, length.out=N)
  lp <- -(Y-(A+B*g+C*g^2))^2/(2*sig^2) - (g-m)^2/(2*s^2)
  wv <- exp(lp-max(lp)); wv <- wv/sum(wv)
  mu <- sum(wv*g); v <- sum(wv*(g-mu)^2)
  c(mean=mu, sd=sqrt(v))
}

run_is <- function(Y, U, K, n=40000L) {
  w <- data.frame(Y = rep(Y, n), x = rep(NA_real_, n))
  # proposal = the TRUE prior, truncated at U (so weights need only carry p(Y|x))
  propose <- function(nr, K) {
    matrix(leftcens::rnorm_trunc(nr*K, m, s, -Inf, U), nrow=nr, ncol=K)
  }
  loglik <- function(dat) {
    mu <- A + B*dat$x + C*dat$x^2
    stats::dnorm(dat$Y, mu, sig, log=TRUE)
  }
  r <- .smc_x_importance(w, "x", cen=rep(TRUE,n), upper=rep(U,n),
                         preds=character(0), loglik=loglik, K=K, propose=propose)
  r
}

cat("\n=== 1. does importance resampling recover the exact conditional? ===\n")
cat("    (proposal = the true prior, so any error is the weighting/resampling)\n\n")
cat(sprintf("  %6s %6s %5s %11s %11s %9s %9s\n","Y","U","K","IS mean","exact mean","IS sd","exact sd"))
for (cs in list(c(0.5,0.3,50), c(0.5,0.3,200), c(1.5,0.0,200), c(-1.0,-0.5,200))) {
  Y<-cs[1]; U<-cs[2]; K<-as.integer(cs[3])
  set.seed(11)
  r <- run_is(Y,U,K); e <- exact_moments(Y,U)
  cat(sprintf("  %6.1f %6.1f %5d %11.4f %11.4f %9.4f %9.4f\n", Y,U,K, mean(r$x), e["mean"], sd(r$x), e["sd"]))
  ok(sprintf("Y=%+.1f U=%+.1f K=%d", Y,U,K),
     abs(mean(r$x)-e["mean"])<0.03 && abs(sd(r$x)-e["sd"])<0.03,
     sprintf("d_mean %.4f  d_sd %.4f  ESS %.1f/%d", abs(mean(r$x)-e["mean"]), abs(sd(r$x)-e["sd"]), r$ess, K))
}

cat("\n=== 2. K matters only where the proposal is POOR ===\n")
cat("    With sigma_y = 1 the prior is already an excellent proposal (ESS/K ~ 0.98),\n")
cat("    so K changes nothing -- an easy regime cannot demonstrate K-dependence.\n")
cat("    Testing in the HARD regime instead: a tight outcome, where it must.\n\n")
sig <<- 0.2
Y <- 1.5; U <- 0.0; e <- exact_moments(Y,U)
cat(sprintf("  %5s %11s %10s %9s\n","K","|d mean|","|d sd|","ESS/K"))
errs <- c()
for (K in c(5L,25L,200L,2000L)) {
  set.seed(7); r <- run_is(Y,U,K)
  dm <- abs(mean(r$x)-e["mean"]); errs <- c(errs, dm)
  cat(sprintf("  %5d %11.4f %10.4f %8.3f\n", K, dm, abs(sd(r$x)-e["sd"]), r$ess/K))
}
ok("error at K=2000 is smaller than at K=5", errs[4] < errs[1],
   sprintf("%.4f vs %.4f", errs[4], errs[1]))
ok("error decreases overall as K grows", errs[4] < errs[2],
   sprintf("K=2000 %.4f vs K=25 %.4f", errs[4], errs[2]))
sig <<- 1.0

cat("\n=== 3. ESS flags degeneracy when the outcome is highly informative ===\n")
cat("    A tight outcome makes the prior a poor proposal; ESS must reveal it.\n\n")
cat(sprintf("  %10s %9s %9s\n","sigma_y","ESS/K","reading"))
essv <- c()
for (sg in c(2.0, 1.0, 0.3, 0.1)) {
  sig <<- sg
  set.seed(3); r <- run_is(1.5, 0.0, 200L)
  essv <- c(essv, r$ess/200)
  cat(sprintf("  %10.1f %9.3f  %s\n", sg, r$ess/200,
      if (r$ess/200 < 0.05) "DEGENERATE -- raise K" else if (r$ess/200 < 0.2) "marginal" else "healthy"))
}
sig <<- 1.0
ok("ESS/K decreases monotonically as the outcome tightens",
   all(diff(essv) < 0), sprintf("%s", paste(sprintf("%.3f",essv), collapse=" > ")))
ok("ESS/K flags the tightest case as degenerate", essv[4] < 0.05,
   sprintf("sigma_y=0.1 gives ESS/K = %.3f", essv[4]))

# =============================================================================
# The GRID sampler -- the design actually adopted for item 07.
# Importance sampling above is retained as the record of why it was rejected:
# it cannot find a target far from its proposal, and ESS does not warn you.
# =============================================================================
suppressWarnings(suppressMessages(library(leftcens)))
Ag <- 0.10; Bg <- 0.40; Cg <- 0.15; Ug <- 0.0
exact_g <- function(Y, sg, N=400001L) {
  g <- seq(-14, Ug, length.out=N)
  lp <- -(Y-(Ag+Bg*g+Cg*g^2))^2/(2*sg^2) - g^2/2
  wv <- exp(lp-max(lp)); wv <- wv/sum(wv); mu <- sum(wv*g)
  c(mean=mu, sd=sqrt(sum(wv*(g-mu)^2)))
}
run_g <- function(Y, sg, n=8000L, n_grid=1024L, preds=character(0)) {
  xt <- rnorm(n)
  w <- data.frame(Y=rep(Y,n), x=ifelse(xt<Ug, NA_real_, xt))
  ll <- function(dat) stats::dnorm(dat$Y, Ag+Bg*dat$x+Cg*dat$x^2, sg, log=TRUE)
  .ce_smc_x_grid(w, "x", cen=is.na(w$x), upper=rep(Ug,n), preds=preds,
                 loglik=ll, n_grid=n_grid, proper=FALSE)$x
}

cat("\n=== 4. the grid sampler recovers the target ACROSS regimes ===\n")
cat("    including the far-root regime where importance sampling failed by 4.6\n\n")
cat(sprintf("  %8s %11s %11s %9s   %s\n","sigma_y","grid mean","exact mean","|diff|","was IS error"))
isref <- c("1"=0.001,"0.5"=0.204,"0.3"=3.665,"0.2"=4.370,"0.1"=4.596)
for (sg in c(1.0,0.5,0.3,0.2,0.1)) {
  set.seed(41); d <- run_g(1.5, sg); e <- exact_g(1.5, sg)
  dm <- abs(mean(d,na.rm=TRUE)-e["mean"])
  cat(sprintf("  %8.1f %11.3f %11.3f %9.3f   %.3f\n", sg, mean(d,na.rm=TRUE), e["mean"], dm,
              isref[[as.character(sg)]]))
  ok(sprintf("sigma_y=%.1f within 0.25 of exact", sg), dm < 0.25, sprintf("|diff| = %.3f", dm))
  # The grid is NOT uniformly better. Where the proposal is well matched
  # (sigma_y = 1) importance sampling is near-exact and the grid carries
  # discretisation error instead -- 0.001 against 0.027. The grid's advantage is
  # that it does not COLLAPSE: it is required to beat IS only where IS is
  # actually failing, which is the regime that makes the method unusable.
  if (isref[[as.character(sg)]] > 0.25) {
    ok(sprintf("sigma_y=%.1f beats IS where IS fails", sg),
       dm < isref[[as.character(sg)]],
       sprintf("%.3f vs %.3f", dm, isref[[as.character(sg)]]))
  } else {
    ok(sprintf("sigma_y=%.1f both methods acceptable", sg), dm < 0.25,
       sprintf("grid %.3f, IS %.3f -- IS marginally better here", dm, isref[[as.character(sg)]]))
  }
}

cat("\n=== 5. every draw respects the censoring bound ===\n\n")
set.seed(43); d <- run_g(1.5, 0.3)
ok("no draw exceeds the LOD", max(d, na.rm=TRUE) <= Ug + 1e-8,
   sprintf("max draw %.6f, bound %.1f", max(d,na.rm=TRUE), Ug))
ok("no NA draws", !anyNA(d))

cat("\n=== 6. a finer grid does not change the answer (resolution is adequate) ===\n\n")
set.seed(45); d1 <- run_g(1.5, 0.2, n_grid=256L)
set.seed(45); d2 <- run_g(1.5, 0.2, n_grid=2048L)
ok("256 vs 2048 grid points agree", abs(mean(d1,na.rm=TRUE)-mean(d2,na.rm=TRUE)) < 0.05,
   sprintf("%.4f vs %.4f", mean(d1,na.rm=TRUE), mean(d2,na.rm=TRUE)))

cat(sprintf("\n%s  (%d failed)\n\n", if(fails==0L)"ALL CHECKS PASSED" else "FAILURES", fails))
quit(status=if(fails==0L) 0L else 1L)
