# =============================================================================
# Phase 1 harness -- metrics (PLAN §7.4)
# -----------------------------------------------------------------------------
# Aggregate per-replication estimates into bias / RMSE / coverage of the known
# ERF estimand, grouped by procedure × ERF form × non-detect fraction.
# =============================================================================

#' Summarise raw per-rep results against the known truth.
#'
#' @param raw A data.frame stacking `run_procedures()` output across reps, with
#'   added columns: `rep`, `erf_form`, `nd_frac`, `estimand_true`.
#' @return A tidy summary data.frame.
summarise_phase1 <- function(raw) {
  grp <- interaction(raw$erf_form, raw$nd_frac, raw$procedure, drop = TRUE)
  parts <- split(raw, grp)

  rows <- lapply(parts, function(df) {
    truth <- df$estimand_true[1]
    est <- df$estimate
    covered <- df$ci_lo <= truth & truth <= df$ci_hi
    n_ok <- sum(is.finite(est))
    data.frame(
      erf_form   = df$erf_form[1],
      nd_frac    = df$nd_frac[1],
      procedure  = df$procedure[1],
      n_rep      = nrow(df),
      n_ok       = n_ok,
      true       = truth,
      mean_est   = mean(est, na.rm = TRUE),
      bias       = mean(est - truth, na.rm = TRUE),
      rel_bias   = mean(est - truth, na.rm = TRUE) / truth,
      rmse       = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      emp_se     = stats::sd(est, na.rm = TRUE),
      mean_ci_w  = mean(df$ci_hi - df$ci_lo, na.rm = TRUE),
      coverage   = mean(covered, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  proc_order <- c("oracle", "complete_case", "leftcens_prestep",
                  "cens_mi_y", "cens_mi_y_shash", "brms_joint")
  out$procedure <- factor(out$procedure, levels = proc_order)
  out <- out[order(out$erf_form, out$nd_frac, out$procedure), ]
  rownames(out) <- NULL
  out
}

#' Pretty console print of the summary (rounded).
print_phase1 <- function(summary_df) {
  num <- vapply(summary_df, is.numeric, logical(1))
  df <- summary_df
  df[num] <- lapply(df[num], function(x) round(x, 4))
  print(df, row.names = FALSE)
  invisible(summary_df)
}
