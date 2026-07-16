source("00_config.R")
source("00_common_functions.R")

log_msg("============================================================")
log_msg("STEP 12: Export per-imputation draws for meta-analysis")
log_msg("============================================================")

export_spec <- analysis_spec$export %||% list()
cohort_id   <- export_spec$cohort_id %||% NULL

if (is.null(cohort_id) || !nzchar(trimws(cohort_id))) {
  log_msg("analysis_spec$export$cohort_id not set — skipping draw export.")
  log_msg("Set cohort_id in 00_config.R to enable export.")
  quit(save = "no", status = 0)
}

log_msg("Cohort ID:", cohort_id)

# ------------------------------------------------------------
# Determine which parameters to export
# ------------------------------------------------------------

export_scope <- export_spec$scope %||% "exposure_only"

var_dict <- readr::read_csv(
  paths$variable_dictionary,
  col_types = readr::cols(.default = "c"),
  show_col_types = FALSE
)

exposure_vars <- var_dict %>%
  dplyr::filter(role == "exposure") %>%
  dplyr::pull(var)

log_msg("Export scope:", export_scope)
if (export_scope == "exposure_only") {
  log_msg("Exposure variables:", paste(exposure_vars, collapse = ", "))
}

# ------------------------------------------------------------
# Load per-imputation draw files
# ------------------------------------------------------------

parameter_manifest <- readRDS(file.path(paths$objects, "parameter_manifest.rds"))

valid_files <- parameter_manifest %>%
  dplyr::filter(purrr::map_lgl(parameter_draw_file, rds_ok))

if (nrow(valid_files) == 0) {
  stop("No valid parameter draw files found. Run Step 6 first.")
}

log_msg("Loading", nrow(valid_files), "per-imputation draw file(s).")

# Pivot each imputation's wide draw tibble to long format, then stack
is_exposure_param <- function(param_name, exposure_vars) {
  bare  <- sub("^b_", "", param_name)
  parts <- strsplit(bare, ":", fixed = TRUE)[[1]]
  any(purrr::map_lgl(exposure_vars, function(ev) any(startsWith(parts, ev))))
}

meta_cols <- c("imputation", ".chain", ".iteration", ".draw")

export_long <- purrr::map_dfr(
  seq_len(nrow(valid_files)),
  function(i) {
    draws_i <- readRDS(valid_files$parameter_draw_file[i])

    param_cols <- setdiff(names(draws_i), meta_cols)

    if (export_scope == "exposure_only") {
      param_cols <- param_cols[
        purrr::map_lgl(param_cols, is_exposure_param, exposure_vars = exposure_vars)
      ]
    }

    if (length(param_cols) == 0) {
      warning("Imputation ", valid_files$imputation[i],
              ": no matching parameter columns after scope filter — skipped.")
      return(NULL)
    }

    draws_i %>%
      dplyr::select(dplyr::all_of(c(meta_cols, param_cols))) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(param_cols),
        names_to  = "parameter",
        values_to = "value"
      ) %>%
      dplyr::transmute(
        cohort_id  = cohort_id,
        parameter  = parameter,
        imputation = imputation,
        draw_index = .draw,
        value      = value
      )
  }
)

if (nrow(export_long) == 0) {
  stop(
    "Export produced 0 rows. Check export_scope and exposure variable roles ",
    "in 00_variable_dictionary.csv."
  )
}

n_params <- dplyr::n_distinct(export_long$parameter)
n_imp    <- dplyr::n_distinct(export_long$imputation)
n_draws  <- nrow(export_long)
log_msg(sprintf("Export: %d parameters × %d imputations = %d rows", n_params, n_imp, n_draws))

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

export_dir <- file.path(paths$results, "export")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(export_dir, "cohort_draws.rds")
csv_path <- file.path(export_dir, "cohort_draws.csv")

saveRDS(export_long, rds_path, compress = TRUE)
readr::write_csv(export_long, csv_path)

log_msg("Wrote", rds_path)
log_msg("Wrote", csv_path)

# Write a small metadata sidecar so the meta-analysis repo can sanity-check
# files without parsing the full draws
meta_sidecar <- list(
  cohort_id  = cohort_id,
  parameters = sort(unique(export_long$parameter)),
  m          = n_imp,
  n_draws    = n_draws,
  family     = analysis_spec$model$family %||% "unknown",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)
jsonlite::write_json(
  meta_sidecar,
  file.path(export_dir, "cohort_metadata.json"),
  auto_unbox = TRUE,
  pretty     = TRUE
)
log_msg("Wrote cohort_metadata.json")

log_msg("SUCCESS: STEP 12: Export per-imputation draws for meta-analysis")
