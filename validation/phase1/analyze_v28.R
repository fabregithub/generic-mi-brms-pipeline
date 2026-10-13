# V28 -- THE ARROW SWEEP. Does the fork residual (derived/measured = 1.147,
# against the pipe's 1.011) track the arrow, and in WHICH sweep?
#
# The design separates two things the fork confounds. Under a fork, delta_zx
# moves BOTH Var(logX1) = 1 + delta^2 AND the exposure conditional's Z1
# coefficient. Under a pipe, delta_xz moves the Z1-X1 relationship while
# Var(logX1) stays 1. So:
#   ratio tracks the arrow in the PIPE sweep too -> implicates the COEFFICIENT
#   ratio tracks it only in the FORK sweep        -> implicates Var(logX1)
#
# Measured effect is PAIRED: the two arms differ only by the pre-v1.6.0 scale
# defect and run on byte-identical data, so the paired SE is the right one.
suppressMessages({library(survival)
  for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f))})
# Pull closed_form() out of the derivation WITHOUT running its own sweep:
# sourcing the file executes the s^2 study at the bottom (minutes, and noise
# in this log). Evaluate only the definitions this analysis needs.
local({
  ex <- parse("derive_fork_residual.R")
  keep <- c("closed_form", "oth", "%||%")
  for (e in ex) if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
                    is.name(e[[2]]) && as.character(e[[2]]) %in% keep)
    eval(e, envir = globalenv())
})
stopifnot(is.function(closed_form), exists("oth"))
`%||%` <- function(a,b) if (is.null(a)) b else a

cells <- list(
  fkd02=list(role="fork", d=0.2), fkd04=list(role="fork", d=0.4),
  fkd06=list(role="fork", d=0.6), fkd08=list(role="fork", d=0.8),
  ppd02=list(role="pipe", d=0.2), ppd04=list(role="pipe", d=0.4),
  ppd06=list(role="pipe", d=0.6), ppd08=list(role="pipe", d=0.8))

## ---- measured, paired ------------------------------------------------------
measured <- function(tag) {
  f <- sprintf("results/v28_%s_raw.csv", tag)
  if (!file.exists(f)) f <- sprintf("results/v28b_%s_raw.csv", tag)
  raw <- read.csv(f)
  do.call(rbind, lapply(names(cells), function(cn) {
    a <- raw[raw$scenario==cn & raw$procedure=="pipeline_properZ", ]
    b <- raw[raw$scenario==cn & raw$procedure=="pipeline_properZ_ds", ]
    k <- intersect(a$rep, b$rep)
    a <- a[match(k,a$rep),]; b <- b[match(k,b$rep),]
    tr <- a$estimand_true[1]
    dif <- 100*(b$estimate - a$estimate)/tr          # per-rep wrapper effect, pp
    data.frame(cell=cn, role=cells[[cn]]$role, delta=cells[[cn]]$d, n_rep=length(k),
               meas=mean(dif), meas_se=sd(dif)/sqrt(length(k)))
  }))
}

## ---- derived ---------------------------------------------------------------
# N MATTERS, AND NOT EQUALLY ACROSS ROLES. A first pass ran this at N=20000 and
# read the pipe 15% above the chain's established +5.42 while the fork agreed --
# which looked like a role-specific defect and was nothing of the kind. The
# closed form is an ASYMPTOTIC statement evaluated on a finite sample, and its
# finite-N bias is +0.87 pp for the pipe against +0.27 pp for the fork at
# N=20000. At N=300000 it returns +5.44 / +5.04, reproducing
# derive_scale_propagation.R's +5.42 / +5.05. Anything smaller manufactures a
# role contrast out of the evaluation, not the mechanism.
derived <- function(nrep = 4L, n = 300000L) {
  do.call(rbind, lapply(names(cells), function(cn) {
    cf <- cells[[cn]]
    v <- vapply(seq_len(nrep), function(i) {
      set.seed(4000L + 97L*i)
      tr <- make_truth(p=3L, erf_form="additive", y_form="linear",
                       z_role=cf$role, target="direct",
                       delta_xz = if (cf$role=="pipe") cf$d else NULL,
                       delta_zx = if (cf$role=="fork") cf$d else NULL,
                       y_family="gaussian")
      sim <- simulate_complete(n, tr, rho=0.4, skew=0.0, z_form="linear")
      closed_form(sim$data, sim$truth, mcar=0.4, ndf=0.4)$pred
    }, numeric(1))
    data.frame(cell=cn, der=mean(v), der_se=sd(v)/sqrt(nrep))
  }))
}

der <- derived()
waves <- c("n800","n3200","n12800")
waves <- waves[vapply(waves, function(w)
  file.exists(sprintf("results/v28_%s_raw.csv", w)) ||
  file.exists(sprintf("results/v28b_%s_raw.csv", w)), logical(1))]

for (tag in waves) {
  m <- measured(tag); x <- merge(m, der, by="cell")
  x$ratio <- x$der/x$meas
  x$ratio_se <- with(x, abs(ratio)*sqrt((der_se/der)^2 + (meas_se/meas)^2))
  cat("\n================ V28  ", tag, " ================\n")
  print(format(x[order(x$role, x$delta),
       c("cell","role","delta","n_rep","meas","meas_se","der","der_se","ratio","ratio_se")],
       digits=3), row.names=FALSE)
  for (rr in c("fork","pipe")) {
    s <- x[x$role==rr,]; s <- s[order(s$delta),]
    ct <- cor.test(s$delta, s$ratio)
    cat(sprintf("  %s: ratio %.3f -> %.3f over delta 0.2->0.8 | slope cor=%.3f p=%.3f\n",
                rr, s$ratio[1], s$ratio[nrow(s)], unname(ct$estimate), ct$p.value))
  }
}

## ---- convergence in n -------------------------------------------------------
# THE POINT OF THE THIRD WAVE. The closed form is asymptotic; the measurement is
# not. If the chain is right, measured -> derived as n grows, for BOTH roles.
cat("\n\n================ convergence of measured -> derived ================\n")
mm <- do.call(rbind, lapply(waves, function(w) cbind(measured(w), n = w)))
tab <- merge(mm, der, by = "cell")
tab <- tab[order(tab$role, tab$delta), ]
for (rr in c("fork", "pipe")) {
  cat("\n--", rr, "--\n")
  s2 <- tab[tab$role == rr, ]
  for (dd in sort(unique(s2$delta))) {
    r <- s2[s2$delta == dd, ]
    r <- r[match(waves, r$n), ]
    cat(sprintf("  delta=%.1f  derived %+5.2f |", dd, r$der[1]))
    for (i in seq_len(nrow(r)))
      if (!is.na(r$meas[i]))
        cat(sprintf("  %s: %+5.2f+/-%.2f (ratio %.3f)", r$n[i], r$meas[i], r$meas_se[i],
                    r$der[i]/r$meas[i]))
    cat("\n")
  }
}
cat("\nFLOOR: delta=0.2 cells are reported but carry no information --\n")
cat("their measured effect is ~0.4-0.7 pp, so derived/measured is a ratio of\n")
cat("near-zero arms. Read delta >= 0.4 only.\n")
