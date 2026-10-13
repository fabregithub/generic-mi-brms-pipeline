# Stronger test: the derivation predicts a DIFFERENT ratio per cell, because
# E[p(1-p)] depends on Var(eta). Matching a pooled median could be luck;
# tracking the cell-to-cell ORDERING cannot be.
suppressMessages(for (f in c("dgp.R","censoring.R","robustness.R")) source(file.path("R",f)))
`%||%` <- function(a,b) if (is.null(a)) b else a
PAIR <- list(c("bin_base","mcar_z40"), c("bin_fork","zr_fork"), c("bin_pipe","zr_pipe"))
ARMS <- c("leftcens_prestep","pipeline_noYx","pipeline_noYz")   # the real defects
pred <- sapply(PAIR, function(p) {
  tr <- make_truth(p = 3L, z_role = Filter(function(s) s$name == p[1],
                     v2_scenarios())[[1]]$z_role %||% "precision", y_family = "binomial")
  set.seed(3); d <- simulate_complete(300000L, tr)$data
  eta <- qlogis(pmin(pmax(fitted(glm(dgp_formula(tr), data=d, family=binomial())),1e-9),1-1e-9))
  sqrt(mean(plogis(eta) * (1 - plogis(eta))))
})
names(pred) <- sapply(PAIR, `[`, 1)
rd <- function(n) { d <- as.data.frame(readRDS(sprintf("results/v26_n%d_latest.rds", n))$summary)
  scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
  d$truth <- vapply(as.character(d$scenario), function(nm) { s <- scs[[nm]]
    make_truth(p=3L, erf_form=s$erf_form, z_role=s$z_role %||% "precision",
               y_family=s$y_family %||% "gaussian")$estimand_true }, 0)
  d$bias_se <- (d$mean_est - d$truth)/d$claimed_se; d }
cat(sprintf("%-10s %10s | %s\n", "cell", "derived", "measured ratio by arm (n=3200, 12800)"))
for (k in seq_along(PAIR)) {
  p <- PAIR[[k]]; rr <- c()
  for (n in c(3200, 12800)) { d <- rd(n)
    for (a in ARMS) {
      bb <- abs(d$bias_se[d$scenario==p[1] & d$procedure==a])
      bg <- abs(d$bias_se[d$scenario==p[2] & d$procedure==a])
      if (length(bb)==1 && length(bg)==1 && bg > 0.5) rr <- c(rr, bb/bg) } }
  cat(sprintf("%-10s %10.4f | median %.4f  [%.3f, %.3f]  n=%d\n",
              p[1], pred[[p[1]]], median(rr), min(rr), max(rr), length(rr)))
}
cat("\nDoes the derivation track the cell ORDERING?\n")
m <- sapply(PAIR, function(p) { rr <- c()
  for (n in c(3200,12800)) { d <- rd(n)
    for (a in ARMS) { bb <- abs(d$bias_se[d$scenario==p[1] & d$procedure==a])
      bg <- abs(d$bias_se[d$scenario==p[2] & d$procedure==a])
      if (length(bb)==1 && length(bg)==1 && bg > 0.5) rr <- c(rr, bb/bg) } }
  median(rr) })
cat(sprintf("  derived  : %s\n", paste(sprintf("%.4f", pred), collapse="  ")))
cat(sprintf("  measured : %s\n", paste(sprintf("%.4f", m), collapse="  ")))
cat(sprintf("  same ordering: %s | max |miss|: %.4f\n",
            identical(order(pred), order(m)), max(abs(pred - m))))
