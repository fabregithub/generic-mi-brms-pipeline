# How much of the validation record did the Z-block scale defect contaminate?
# Identical arm set, same cells, 200 reps, pre- vs post-fix.
# ANCHOR: micePmm and properZ were independently measured at -3.16 and -4.40 pp.
# If they do not reproduce, distrust the rest of this table.
`%||%` <- function(a,b) if (is.null(a)) b else a
source("R/dgp.R"); source("R/robustness.R")
scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
rd <- function(tag) {
  d <- as.data.frame(readRDS(sprintf("results/%s_latest.rds", tag))$summary)
  d$truth <- vapply(as.character(d$scenario), function(nm)
    make_truth(p=3L, erf_form=scs[[nm]]$erf_form,
               z_role=scs[[nm]]$z_role %||% "precision")$estimand_true, 0)
  d$bias_pct <- 100*(d$mean_est - d$truth)/d$truth
  d$mcse <- 100*(d$emp_se/sqrt(d$n_ok))/d$truth
  d
}
a <- rd("v25_scope_pre"); b <- rd("v25_scope_post")
g <- function(d,cl,ar,f) { v <- d[[f]][d$scenario==cl & d$procedure==ar]
                           if (length(v)!=1) NA_real_ else v }
ARMS <- c("oracle","pipeline_block_fcs","pipeline_bartMI",
          "pipeline_micePmm","pipeline_properZ")
cat("Z-block scale fix: effect on the PIPELINE arms (200 reps, identical arm set)\n\n")
cat(sprintf("%-9s %-22s %9s %9s %9s %7s\n","cell","arm","pre","post","change","mcse"))
for (cl in c("zr_fork","zr_pipe")) {
  for (ar in ARMS) {
    x <- g(a,cl,ar,"bias_pct"); y <- g(b,cl,ar,"bias_pct")
    if (is.na(x)||is.na(y)) next
    cat(sprintf("%-9s %-22s %+9.2f%% %+9.2f%% %+9.2f pp %7.2f\n",
                cl, ar, x, y, y-x, g(b,cl,ar,"mcse")))
  }
  cat("\n")
}
cat("--- ANCHOR CHECK (independently measured: micePmm -3.16, properZ -4.40 on zr_fork) ---\n")
cat(sprintf("  micePmm %+.2f pp | properZ %+.2f pp  -> %s\n",
  g(b,"zr_fork","pipeline_micePmm","bias_pct")-g(a,"zr_fork","pipeline_micePmm","bias_pct"),
  g(b,"zr_fork","pipeline_properZ","bias_pct")-g(a,"zr_fork","pipeline_properZ","bias_pct"),
  "compare to -3.16 / -4.40"))
cat("\n--- V17's arm: how much of its headline was the scale defect? ---\n")
for (cl in c("zr_fork","zr_pipe")) {
  x <- g(a,cl,"pipeline_bartMI","bias_pct"); y <- g(b,cl,"pipeline_bartMI","bias_pct")
  cat(sprintf("  %-8s bartMI %+.2f%% -> %+.2f%%  (scale defect accounted for %.0f%% of it)\n",
              cl, x, y, 100*(x-y)/x))
}
cat("\n--- the SHIPPED DEFAULT, i.e. what a user gets ---\n")
for (cl in c("zr_fork","zr_pipe")) {
  x <- g(a,cl,"pipeline_block_fcs","bias_pct"); y <- g(b,cl,"pipeline_block_fcs","bias_pct")
  cat(sprintf("  %-8s block_fcs %+.2f%% -> %+.2f%%  (%+.2f pp)\n", cl, x, y, y-x))
}
