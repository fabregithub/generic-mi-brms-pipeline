#!/usr/bin/env Rscript
# =============================================================================
# CAN THE FOUR CASES BE TESTED AT ONCE?  A design test, before any real run.
# -----------------------------------------------------------------------------
# THEORY.md §6b concluded the DAG-factorised draw is needed in FOUR cases, and
# they differ along TWO axes that the existing runner does not currently cross:
#
#     case            z_role      target    Z1 in analysis?   child factor?
#     fork            fork        direct    yes (backdoor)    no  (parent)
#     pipe, direct    pipe_nl     direct    yes               yes  p(Z1|x)
#     pipe, total     pipe_nl     total     NO                yes  p(Z1|x)
#     collider        collider    direct    NO                yes  p(Z1|x,Y)
#
# So "testing four at once" needs exactly two new pieces: `target` on make_truth
# (added in dgp.R this session) and a role-aware child factor (dag_draw.R). This
# script checks BOTH before a run is scheduled, at a small n and few reps -- it
# is a design test, NOT a validation track. Its job is to answer three questions:
#
#   Q1  Does the factorisation hold in all four cases at once, with the ORACLE
#       child model? If not, the algorithm is wrong and no amount of estimation
#       care will save it.
#   Q2  Is the LINEAR child form -- the only one that needs no new assumption --
#       enough where the arrow is linear, and insufficient where it is not?
#       That is the boundary the shipped default has to sit on.
#   Q3  Does adding the child factor HURT the fork, where Z1 is a parent and the
#       factor does not belong? A safe default must be role-aware, not always-on.
#
# NOT REGISTERED CRITERIA. Numbers here inform the design of the track; the
# track registers its own gates. PLAN §11: a gate must be able to fail for the
# reason it names.
# =============================================================================
suppressMessages({
  library(survival); library(parallel)
  for (f in c("dgp.R", "censoring.R", "dag_draw.R")) source(file.path("R", f))
})

REPS    <- as.integer(Sys.getenv("REPS",    "200"))
N_OBS   <- as.integer(Sys.getenv("N_OBS",   "2000"))
M_IMP   <- as.integer(Sys.getenv("M_IMP",   "20"))
ND      <- as.numeric(Sys.getenv("ND",      "0.4"))
WORKERS <- as.integer(Sys.getenv("WORKERS", "20"))

# REPS = 200 is chosen from the currency, not from habit. The gate is on
# bias/SE at 0.3, and the Monte-Carlo error of an R-rep mean bias/SE is
# ~1/sqrt(R) -- 0.16 at R = 40, which is half the gate and cannot resolve it;
# 0.07 at R = 200, which can. A full-wave pilot (REPS = WORKERS) measures the
# cost before the real run, per CLAUDE.md.

# The four required cases, plus the linear `pipe` as the control that isolates
# non-linearity from the factorisation itself.
CASES <- list(
  list(key = "fork",         role = "fork",     target = "direct"),
  list(key = "pipe_dir",     role = "pipe",     target = "direct"),
  list(key = "pipe_tot",     role = "pipe",     target = "total"),
  list(key = "pipenl_dir",   role = "pipe_nl",  target = "direct"),
  list(key = "pipenl_tot",   role = "pipe_nl",  target = "total"),
  list(key = "collider",     role = "collider", target = "direct")
)

# Arms, in increasing order of what they know. `noY` is the incongenial draw the
# theory predicts attenuates; `Yonly` is the congenial-but-DAG-blind draw that
# every track up to V23 used; `dag_*` add the child factor.
ARMS <- list(
  list(key = "noY",       use_y = FALSE, use_child = FALSE, form = "linear"),
  list(key = "Yonly",     use_y = TRUE,  use_child = FALSE, form = "linear"),
  list(key = "dag_lin",   use_y = TRUE,  use_child = TRUE,  form = "linear"),
  list(key = "dag_true",  use_y = TRUE,  use_child = TRUE,  form = "true")
)

# Under a fork Z1 is a PARENT, so `dag_spec()` gates the child factor off and the
# two dag_* arms are BIT-IDENTICAL to `Yonly`. Running them would spend a third
# of the fork's budget re-measuring the same number, and would then invite
# reading three matching rows as three pieces of evidence. Dropped explicitly.
arms_for <- function(cs)
  if (identical(cs$role, "fork")) ARMS[1:2] else ARMS

#' Rubin-pool the focal coefficient over m completed datasets.
pool_b1 <- function(imps, fo) {
  est <- vapply(imps, function(w) stats::coef(stats::lm(fo, data = w))[["logX1"]], 0)
  se  <- vapply(imps, function(w) {
    f <- stats::lm(fo, data = w); sqrt(diag(stats::vcov(f))[["logX1"]]) }, 0)
  m <- length(est); qbar <- mean(est); ubar <- mean(se^2); bvar <- stats::var(est)
  tot <- ubar + (1 + 1 / m) * bvar
  nu <- if (bvar > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bvar))^2 else Inf
  c(est = qbar, se = sqrt(tot), df = nu)
}

#' One replicate of one case: simulate, censor, then run every arm on the SAME
#' censored data with the SAME imputation seed, so arm differences are the draw
#' and nothing else.
one_rep <- function(rep, cs, truth, fo) {
  set.seed(20260909L + rep)              # paired across cases
  comp <- simulate_complete(N_OBS, truth)$data
  cens <- inject_left_censoring(comp, nd_frac = ND, censor_which = 1L)
  do.call(rbind, lapply(arms_for(cs), function(ar) {
    set.seed(77000L + rep)               # paired across arms within a rep
    p <- tryCatch(
      pool_b1(dag_impute_datasets(cens, truth, m = M_IMP, child_form = ar$form,
                                  use_y = ar$use_y, use_child = ar$use_child),
              fo),
      error = function(e) c(est = NA_real_, se = NA_real_, df = NA_real_))
    data.frame(case = cs$key, role = cs$role, target = cs$target, arm = ar$key,
               rep = rep, truth = truth$estimand_true,
               est = p[["est"]], se = p[["se"]], df = p[["df"]])
  }))
}

t0 <- Sys.time()
# RELOAD=1 re-aggregates the saved raw rows instead of re-simulating. The
# summary and the computed answers evolved several times while the numbers did
# not; re-running 19 minutes of simulation to re-print a table is waste, and
# hand-recomputing the summary outside the script is how a reported number stops
# matching the script that supposedly produced it.
RAW_CSV <- "results/derive_dag_four_raw.csv"
raw <- if (nzchar(Sys.getenv("RELOAD")) && file.exists(RAW_CSV)) {
  cat(sprintf("  RELOAD: re-aggregating %s (no simulation)\n", RAW_CSV))
  utils::read.csv(RAW_CSV)
} else do.call(rbind, lapply(CASES, function(cs) {
  truth <- make_truth(z_role = cs$role, target = cs$target)
  fo    <- dgp_formula(truth)
  out <- do.call(rbind, mclapply(seq_len(REPS), one_rep, cs = cs, truth = truth,
                                 fo = fo, mc.cores = WORKERS))
  cat(sprintf("  %-11s done (%d reps, %.1f min elapsed)\n", cs$key, REPS,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  utils::flush.console()
  out
}))

agg <- do.call(rbind, lapply(split(raw, list(raw$case, raw$arm), drop = TRUE),
  function(g) {
    tcrit <- stats::qt(0.975, pmax(g$df, 1))
    data.frame(case = g$case[1], arm = g$arm[1], truth = g$truth[1],
               mean_est = mean(g$est),
               bias_pct = 100 * (mean(g$est) - g$truth[1]) / g$truth[1],
               bias_se  = abs(mean(g$est) - g$truth[1]) / mean(g$se),
               mean_se  = mean(g$se),
               cover    = mean(abs(g$est - g$truth) <= tcrit * g$se),
               mcse_pct = 100 * (stats::sd(g$est) / sqrt(nrow(g))) / g$truth[1])
  }))
agg <- agg[order(match(agg$case, vapply(CASES, `[[`, "", "key")),
                 match(agg$arm,  vapply(ARMS,  `[[`, "", "key"))), ]

cat(sprintf("\nn = %d, reps = %d, m = %d, LCR = %.0f%%, %d workers, %.1f min\n",
            N_OBS, REPS, M_IMP, 100 * ND, WORKERS,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
# `se_vs_Yonly` is here because bias/SE cannot show an efficiency gain: an arm
# that removes no bias but sharpens the draw moves the SE and leaves bias/SE
# roughly alone. That is the collider's whole expected payoff.
#
# COMPUTED PAIRED, with an error bar. Every arm sees the same data and the same
# imputation seed within a rep, so an unpaired ratio of means throws that away
# and -- worse -- reports a difference with no Monte-Carlo error at all, which is
# how a 1.4% "efficiency gain" would get asserted without anything backing it.
.pair_dse <- function(cs, arm) {
  a <- raw[raw$case == cs & raw$arm == "Yonly", c("rep", "se")]
  b <- raw[raw$case == cs & raw$arm == arm,     c("rep", "se")]
  if (!nrow(b) || !nrow(a)) return(c(est = NA_real_, mcse = NA_real_))
  mg <- merge(a, b, by = "rep"); dd <- 100 * (mg$se.y / mg$se.x - 1)
  c(est = mean(dd), mcse = stats::sd(dd) / sqrt(nrow(mg)))
}
pd <- t(mapply(.pair_dse, agg$case, agg$arm))
agg$se_vs_Yonly      <- pd[, "est"]
agg$se_vs_Yonly_mcse <- pd[, "mcse"]

cat(sprintf("%-11s %-9s %7s %9s %8s %7s %15s %7s\n",
            "case", "arm", "truth", "bias%", "bias/SE", "cover",
            "dSE% (paired)", "mcse%"))
for (i in seq_len(nrow(agg))) with(agg[i, ], cat(sprintf(
  "%-11s %-9s %7.4f %+9.2f %8.2f %7.3f %+8.2f+/-%.2f %7.2f\n",
  case, arm, truth, bias_pct, bias_se, cover, se_vs_Yonly,
  se_vs_Yonly_mcse, mcse_pct)))

gv <- function(cs, ar, f) agg[[f]][agg$case == cs & agg$arm == ar]
cat("\n--- computed answers ------------------------------------------------\n")
# Q1: does the factorisation hold in all four required cases at once?
# The DAG-correct arm for the fork IS `Yonly`, because the factorisation puts no
# child factor there; for the other three it is `dag_true`.
q1 <- c(fork = abs(gv("fork", "Yonly", "bias_se")),
        vapply(c("pipenl_dir", "pipenl_tot", "collider"),
               function(c0) abs(gv(c0, "dag_true", "bias_se")), 0))
cat(sprintf("Q1 factorisation, ORACLE child: max |bias/SE| over the four = %.2f (%s)\n",
            max(q1), paste(sprintf("%s %.2f", names(q1), q1), collapse = ", ")))
cat(sprintf("   -> %s\n", if (max(q1) < 0.3)
  "HOLDS in all four simultaneously; one track can carry all four cells." else
  "does NOT hold everywhere -- fix the factorisation before designing a track."))
# Q2: is the linear child form enough where the arrow is linear, and not where not?
cat(sprintf("Q2 linear child form: pipe_dir %+.2f%% (bias/SE %.2f) | pipenl_dir %+.2f%% (%.2f)\n",
            gv("pipe_dir", "dag_lin", "bias_pct"), gv("pipe_dir", "dag_lin", "bias_se"),
            gv("pipenl_dir", "dag_lin", "bias_pct"), gv("pipenl_dir", "dag_lin", "bias_se")))
cat(sprintf("   -> %s\n", if (abs(gv("pipe_dir","dag_lin","bias_se")) < 0.3 &&
                              abs(gv("pipenl_dir","dag_lin","bias_se")) >= 0.3)
  "boundary is where expected: linear suffices iff the arrow is linear, so child_form is a real sensitivity axis." else
  "boundary is NOT where expected -- read the table before assuming child_form matters."))
# Q3: does an always-on child factor hurt the fork?
cat(sprintf("Q3 fork, DAG-correct draw (no child factor): %+.2f%% (bias/SE %.2f, cover %.3f)\n",
            gv("fork", "Yonly", "bias_pct"), gv("fork", "Yonly", "bias_se"),
            gv("fork", "Yonly", "cover")))
cat(sprintf("   -> %s\n", if (abs(gv("fork","Yonly","bias_se")) < 0.3)
  "gating the child factor on role is right: the fork needs only the Y factor." else
  "the fork is NOT clean even with the Y factor -- the parent factor is suspect."))

# Q4: is a WRONG child factor worse than no child factor? This decides whether a
# linear default can ship at all, and it is not implied by Q2 -- Q2 only says the
# linear form is insufficient, not that it is actively harmful.
q4 <- vapply(c("pipenl_dir", "pipenl_tot"), function(c0)
  c(lin = gv(c0, "dag_lin", "bias_se"), non = gv(c0, "Yonly", "bias_se")), numeric(2))
cat(sprintf("Q4 wrong child factor vs none (bias/SE): pipenl_dir %.2f vs %.2f | pipenl_tot %.2f vs %.2f\n",
            q4[1, 1], q4[2, 1], q4[1, 2], q4[2, 2]))
cat(sprintf("   -> %s\n", if (all(q4[1, ] > q4[2, ]))
  "YES -- a misspecified child factor is WORSE than omitting it, so `linear` must not be the shipped default; the form has to be declared." else
  "no -- the linear form degrades gracefully, so it is a defensible default."))
# Q5: where the child factor removes no bias, does it buy precision?
i5 <- agg$case == "collider" & agg$arm == "dag_true"
d5 <- agg$se_vs_Yonly[i5]; e5 <- agg$se_vs_Yonly_mcse[i5]
cat(sprintf("Q5 collider, child factor on: bias/SE %.2f -> %.2f, paired SE %+.2f%% +/- %.2f\n",
            gv("collider", "Yonly", "bias_se"), gv("collider", "dag_true", "bias_se"),
            d5, e5))
cat(sprintf("   -> %s\n", if (d5 + 2 * e5 < 0)
  sprintf("under a collider Z1 is congenially IGNORABLE, so the factor is an EFFICIENCY gain (%.1f%% of the SE), not a bias fix -- and %.1f%% is small against the assumption it costs.", abs(d5), abs(d5)) else
  "no precision gain distinguishable from zero -- under a collider the child factor is not worth its assumption."))

dir.create("results", showWarnings = FALSE)
utils::write.csv(raw, "results/derive_dag_four_raw.csv", row.names = FALSE)
utils::write.csv(agg, "results/derive_dag_four.csv", row.names = FALSE)
cat("\nwrote results/derive_dag_four{,_raw}.csv\n")
