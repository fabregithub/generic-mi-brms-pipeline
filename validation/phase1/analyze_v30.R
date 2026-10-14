# =============================================================================
# V30 -- does the record transfer to an ORDINAL outcome, at the derived size?
# -----------------------------------------------------------------------------
# Currency: |bias/SE| ratio, family cell over its Gaussian twin, per (cell, arm),
# the twin paired on byte-identical exposures via seed_as. Relative bias is NOT
# comparable across families (a latent-scale coefficient and a mean difference
# are different quantities), which is why V26 used this currency too.
#
# The predictions below are REGISTERED LITERALS, not recomputed -- a registered
# number that is re-derived at analysis time can drift. Reproduce them with
#     Rscript derive_ordinal_information.R
# =============================================================================
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(.here) || !nzchar(.here)) {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  .here <- if (length(f)) dirname(normalizePath(f)) else "."
}
suppressMessages({source(file.path(.here, "R", "dgp.R"))
                  source(file.path(.here, "R", "robustness.R"))})
`%||%` <- function(a, b) if (is.null(a)) b else a

TAG <- Sys.getenv("TAG", "v30")
NS  <- as.integer(strsplit(Sys.getenv("NLEVELS", "3200,12800"), ",")[[1]])

# cell -> (Gaussian twin, K, REGISTERED prediction)
CELLS <- list(
  bin_fork  = list(gau = "zr_fork", K = 2L, pred = 0.4459),
  bin_pipe  = list(gau = "zr_pipe", K = 2L, pred = 0.4488),
  ord3_fork = list(gau = "zr_fork", K = 3L, pred = 0.5091),
  ord5_fork = list(gau = "zr_fork", K = 5L, pred = 0.5449),
  ord7_fork = list(gau = "zr_fork", K = 7L, pred = 0.5559),
  ord3_pipe = list(gau = "zr_pipe", K = 3L, pred = 0.5120),
  ord5_pipe = list(gau = "zr_pipe", K = 5L, pred = 0.5477),
  ord7_pipe = list(gau = "zr_pipe", K = 7L, pred = 0.5586),
  # V30b: unbalanced categories at FIXED K = 3. These move the prediction
  # without moving K, which "ratios just rise with K" cannot do.
  ord3u_fork = list(gau = "zr_fork", K = 3L, pred = 0.4306, split = "70/20/10"),
  ord3v_fork = list(gau = "zr_fork", K = 3L, pred = 0.3358, split = "85/10/05"),
  ord3u_pipe = list(gau = "zr_pipe", K = 3L, pred = 0.4327, split = "70/20/10"),
  ord3v_pipe = list(gau = "zr_pipe", K = 3L, pred = 0.3382, split = "85/10/05"))
# The balanced K=3 cells, for the P4 ordering.
BAL3 <- c(fork = "ord3_fork", pipe = "ord3_pipe")
ARMS <- c("complete_case", "leftcens_prestep", "pipeline_bartMI",
          "pipeline_noYx", "pipeline_noYz", "pipeline_properZ",
          "pipeline_properZ_ds")
FLOOR <- 0.5          # the Gaussian twin must HAVE a defect for a ratio to mean
                      # anything. V26 registered 0.1, found it far too weak (a
                      # ratio of two near-zero arms read 23.56x) and corrected to
                      # 0.5. Registered at 0.5 here from the start.
V26_ANCHOR <- c(bin_fork = 0.4731, bin_pipe = 0.4527)

scs <- v2_scenarios(); names(scs) <- vapply(scs, `[[`, "", "name")
truth_of <- function(nm) {
  s <- scs[[nm]]
  make_truth(p = 3L, erf_form = s$erf_form, y_form = s$y_form %||% "linear",
             z_role = s$z_role %||% "precision", nl_a = s$nl_a, nl_c = s$nl_c,
             target = s$target %||% "direct", delta_xz = s$delta_xz,
             delta_zx = s$delta_zx, y_family = s$y_family %||% "gaussian",
             y_tau = s$y_tau)$estimand_true
}
load_stage <- function(n) {
  f <- file.path(.here, "results", sprintf("%s_n%d_latest.rds", TAG, n))
  if (!file.exists(f)) return(NULL)
  d <- as.data.frame(readRDS(f)$summary); d$n <- n
  d$truth   <- vapply(as.character(d$scenario), truth_of, 0)
  d$bias_se <- (d$mean_est - d$truth) / d$claimed_se
  d
}
all <- do.call(rbind, Filter(Negate(is.null), lapply(NS, load_stage)))
if (is.null(all) || !nrow(all)) stop("no V30 results for TAG=", TAG)
g <- function(d, cl, ar, f) { v <- d[[f]][d$scenario == cl & d$procedure == ar]
                              if (length(v) != 1L) NA_real_ else v }

for (n in sort(unique(all$n))) {
  d <- all[all$n == n, , drop = FALSE]
  cat(sprintf("\n=====================================================================\n"))
  cat(sprintf("  V30   n = %d   (%d reps)\n", n, max(d$n_ok, na.rm = TRUE)))
  cat(sprintf("=====================================================================\n"))

  tab <- NULL
  for (nm in names(CELLS)) {
    cf <- CELLS[[nm]]; rs <- c()
    for (ar in ARMS) {
      bb <- abs(g(d, nm, ar, "bias_se")); bg <- abs(g(d, cf$gau, ar, "bias_se"))
      if (!is.finite(bb) || !is.finite(bg) || bg < FLOOR) next
      rs <- c(rs, bb / bg)
    }
    tab <- rbind(tab, data.frame(cell = nm, K = cf$K, n_arms = length(rs),
      ratio = if (length(rs)) mean(rs) else NA_real_,
      se    = if (length(rs) > 1L) sd(rs)/sqrt(length(rs)) else NA_real_,
      med   = if (length(rs)) median(rs) else NA_real_,
      pred  = cf$pred))
  }
  tab$miss <- tab$ratio - tab$pred
  cat("\n--- ratios (arms with a real Gaussian defect, |bias/SE| > 0.5) ------\n")
  print(format(tab[, c("cell","K","n_arms","ratio","se","med","pred","miss")], digits = 3),
        row.names = FALSE)

  # ---- P1 LEVEL ------------------------------------------------------------
  cat("\n--- P1  measured ratio within 0.05 of the derived value -------------\n")
  ok <- is.finite(tab$miss) & tab$n_arms >= 3L
  if (!any(ok)) {
    cat("  VOID -- no cell has >=3 arms clearing the floor. Not evidence either way.\n")
    p1 <- NA
  } else {
    for (i in which(ok)) cat(sprintf("  %-10s %+.4f  %s\n", tab$cell[i], tab$miss[i],
                                     if (abs(tab$miss[i]) <= 0.05) "ok" else "MISS"))
    p1 <- all(abs(tab$miss[ok]) <= 0.05)
    if (any(!ok)) cat(sprintf("  (%d cell(s) excluded: fewer than 3 arms clear the floor)\n", sum(!ok)))
  }
  cat(sprintf("  P1: %s\n", if (isTRUE(p1)) "MET" else if (is.na(p1)) "VOID" else "NOT MET"))

  # ---- P2 TREND ------------------------------------------------------------
  cat("\n--- P2  the ratio RISES with K, in both roles -----------------------\n")
  p2 <- TRUE
  for (role in c("fork", "pipe")) {
    # BALANCED cells only -- the unbalanced ones share K = 3 and would both
    # duplicate a K value and drag the trend down by construction.
    s <- tab[grepl(paste0("^(bin|ord[0-9])_?", "|"), tab$cell) &
             grepl(paste0(role, "$"), tab$cell) & !grepl("[uv]_", tab$cell) &
             is.finite(tab$ratio), ]
    s <- s[order(s$K), ]
    if (nrow(s) < 3L) { cat(sprintf("  %-4s: VOID (%d usable K)\n", role, nrow(s))); p2 <- NA; next }
    inc <- all(diff(s$ratio) > 0)
    cat(sprintf("  %-4s: %s  (K=%s)  -> %s\n", role,
                paste(sprintf("%.3f", s$ratio), collapse = " -> "),
                paste(s$K, collapse = ","), if (inc) "rises" else "DOES NOT RISE"))
    if (!isTRUE(inc)) p2 <- FALSE
  }
  cat(sprintf("  P2: %s\n", if (isTRUE(p2)) "MET" else if (is.na(p2)) "VOID" else "NOT MET"))

  # ---- P3 ANCHOR -----------------------------------------------------------
  cat("\n--- P3  the K=2 cells reproduce V26 ---------------------------------\n")
  p3 <- TRUE
  for (nm in names(V26_ANCHOR)) {
    r <- tab$ratio[tab$cell == nm]
    if (!length(r) || !is.finite(r)) { p3 <- NA; next }
    d26 <- r - V26_ANCHOR[[nm]]
    cat(sprintf("  %-10s V30 %.4f  V26 %.4f  diff %+.4f  %s\n", nm, r, V26_ANCHOR[[nm]],
                d26, if (abs(d26) <= 0.06) "ok" else "DIVERGES"))
    if (abs(d26) > 0.06) p3 <- FALSE
  }
  cat(sprintf("  P3: %s\n", if (isTRUE(p3)) "MET" else if (is.na(p3)) "VOID" else "NOT MET"))
  # ---- P4  UNBALANCED categories at fixed K --------------------------------
  cat("\n--- P4  at K=3 the ratio FALLS as the split gets unbalanced ---------\n")
  p4 <- TRUE
  for (role in c("fork", "pipe")) {
    nms <- c(BAL3[[role]], paste0("ord3u_", role), paste0("ord3v_", role))
    rr  <- vapply(nms, function(z) { v <- tab$ratio[tab$cell == z]
                                     if (length(v) == 1L) v else NA_real_ }, 0)
    bin <- tab$ratio[tab$cell == paste0("bin_", role)]
    if (anyNA(rr) || !length(bin) || !is.finite(bin)) {
      cat(sprintf("  %-4s: VOID (a cell is missing or cleared no arms)\n", role)); p4 <- NA; next }
    falls <- all(diff(rr) < 0)
    below <- rr[[3]] < bin
    cat(sprintf("  %-4s balanced %.3f -> 70/20/10 %.3f -> 85/10/05 %.3f  | %s\n",
                role, rr[[1]], rr[[2]], rr[[3]], if (falls) "falls" else "DOES NOT FALL"))
    cat(sprintf("        85/10/05 (%.3f) vs its own K=2 cell (%.3f): %s\n",
                rr[[3]], bin, if (below) "BELOW, as predicted" else "NOT below"))
    if (!isTRUE(falls && below)) p4 <- FALSE
  }
  cat(sprintf("  P4: %s\n", if (isTRUE(p4)) "MET" else if (is.na(p4)) "VOID" else "NOT MET"))
}

cat("\nNOT A GATE: the fork-vs-pipe contrast at fixed K. The derivation puts\n")
cat("those within 0.003 of each other, far under any achievable MC error, so a\n")
cat("criterion on it could not fail for the reason it names. The two roles are\n")
cat("REPLICATION of the K trend, nothing more.\n")
