# Item 1: does discarding an INFORMATIVE covariate generate the bias, or does
# unrepresentable curvature? Three cells, byte-identical exposures, one arm set.
#   dag_pipe_null  arrow ABSENT   (u = 0, Z1 carries NO information about x)
#   dag_pipe_dir   arrow LINEAR   (u = 0, Z1 IS informative)
#   dag_pnl_dir    arrow CURVED   (u > 0, Z1 informative and unrepresentable)
`%||%` <- function(a,b) if (is.null(a)) b else a
source("R/dgp.R"); source("R/robustness.R")
scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
tv <- function(nm) { s <- scs[[nm]]
  make_truth(p=3L, erf_form=s$erf_form, z_role=s$z_role,
             target=s$target %||% "direct", nl_a=s$nl_a, nl_c=s$nl_c,
             delta_xz=s$delta_xz)$estimand_true }
CELLS <- c("dag_pipe_null","dag_pipe_dir","dag_pnl_dir")
LAB <- c(dag_pipe_null="arrow ABSENT", dag_pipe_dir="arrow LINEAR",
         dag_pnl_dir="arrow CURVED")
for (n in as.integer(strsplit(Sys.getenv("NS","800"), ",")[[1]])) {
  f <- sprintf("results/v25_arrow_n%d_latest.rds", n)
  if (!file.exists(f)) { cat("(n =", n, "not yet available)\n"); next }
  d <- as.data.frame(readRDS(f)$summary)
  d$truth <- vapply(as.character(d$scenario), tv, 0)
  d$bias_pct <- 100*(d$mean_est - d$truth)/d$truth
  d$bias_se  <- (d$mean_est - d$truth)/d$claimed_se
  d$mcse     <- (d$emp_se/sqrt(d$n_ok))/d$claimed_se
  cat(sprintf("\n=== n = %d, %d reps ===\n", n, max(d$n_ok, na.rm=TRUE)))
  cat(sprintf("%-14s %-13s %-10s %8s %9s %7s %7s\n","cell","arrow","arm",
              "bias%","bias/SE","mcse","cover"))
  for (cl in CELLS) for (ar in c("oracle","dag_noY","dag_Yonly","dag_lin","dag_true")) {
    r <- d[d$scenario==cl & d$procedure==ar, ]
    if (!nrow(r)) next
    cat(sprintf("%-14s %-13s %-10s %+8.2f %+9.2f %7.2f %7.3f\n",
                cl, LAB[[cl]], ar, r$bias_pct, r$bias_se, r$mcse, r$coverage))
  }
  g <- function(cl, ar) d$bias_se[d$scenario==cl & d$procedure==ar]
  cat("\n--- the contrast (dag_Yonly: the draw that discards Z1) ---\n")
  for (cl in CELLS)
    cat(sprintf("  %-13s bias/SE %+6.2f +/- %.2f\n", LAB[[cl]], g(cl,"dag_Yonly"),
                d$mcse[d$scenario==cl & d$procedure=="dag_Yonly"]))
  a <- abs(g("dag_pipe_null","dag_Yonly")); b <- abs(g("dag_pipe_dir","dag_Yonly"))
  cat(sprintf("\n  -> %s\n", if (a < 0.3 && b > 0.3)
    "CONFIRMED within one instrument: u = 0 in BOTH the absent and linear cells, yet only the LINEAR one is biased. So DISCARDING INFORMATION generates the bias; u only modifies it." else
    sprintf("NOT the expected pattern (absent %.2f, linear %.2f) -- read the table.", a, b)))
}
