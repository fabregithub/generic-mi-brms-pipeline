# ============================================================
# 00_common_functions.R
# Generic helper functions for MI + brms pipeline
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(bayestestR)
  library(future)
  library(furrr)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ------------------------------------------------------------
# Project directories and logging
# ------------------------------------------------------------

setup_project_dirs <- function(paths) {
  dir.create("data", recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$objects, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$imputed_data, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$model_data, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$fits, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$results, recursive = TRUE, showWarnings = FALSE)
}

init_logging <- function(prefix = "pipeline") {
  assign("log_file", paste0(prefix, "_progress.log"), envir = .GlobalEnv)
  assign("heartbeat_file", paste0(prefix, "_heartbeat.txt"), envir = .GlobalEnv)
  assign("success_flag", paste0(prefix, "_success.flag"), envir = .GlobalEnv)
  assign("error_flag", paste0(prefix, "_error.flag"), envir = .GlobalEnv)
  invisible(TRUE)
}

update_heartbeat <- function(status = "running") {
  heartbeat_file <- get0("heartbeat_file", ifnotfound = "pipeline_heartbeat.txt", envir = .GlobalEnv)
  writeLines(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", status),
    heartbeat_file
  )
}

log_msg <- function(..., .sep = " ") {
  log_file <- get0("log_file", ifnotfound = "pipeline_progress.log", envir = .GlobalEnv)
  msg <- paste(..., sep = .sep)
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n", file = log_file, append = TRUE)
  cat(line, "\n")
  flush.console()
  update_heartbeat(msg)
}

log_section <- function(title) {
  log_msg("============================================================")
  log_msg(title)
  log_msg("============================================================")
}

# ------------------------------------------------------------
# Memory guard
# ------------------------------------------------------------

get_r_memory_used_gb <- function() {
  mem <- gc()
  sum(mem[, 2], na.rm = TRUE) / 1024
}

get_mac_memory_info_gb <- function() {
  if (Sys.info()[["sysname"]] != "Darwin") {
    return(list(available_gb = NA_real_, total_gb = NA_real_, used_gb = NA_real_))
  }

  vm <- tryCatch(system2("vm_stat", stdout = TRUE, stderr = TRUE), error = function(e) character(0))
  total_bytes <- tryCatch(as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE)), error = function(e) NA_real_)

  if (length(vm) == 0 || is.na(total_bytes)) {
    return(list(available_gb = NA_real_, total_gb = total_bytes / 1024^3, used_gb = NA_real_))
  }

  page_size <- suppressWarnings(as.numeric(sub(".*page size of ([0-9]+) bytes.*", "\\1", vm[1])))
  if (is.na(page_size)) page_size <- 16384

  extract_pages <- function(label) {
    line <- grep(paste0("^", label, ":"), vm, value = TRUE)
    if (length(line) == 0) return(0)
    suppressWarnings(as.numeric(gsub("[^0-9]", "", line[1])))
  }

  available_bytes <- (
    extract_pages("Pages free") +
      extract_pages("Pages inactive") +
      extract_pages("Pages speculative")
  ) * page_size

  total_gb <- total_bytes / 1024^3
  available_gb <- available_bytes / 1024^3
  list(available_gb = available_gb, total_gb = total_gb, used_gb = total_gb - available_gb)
}

log_memory <- function(memory_guard = NULL) {
  if (isTRUE(memory_guard$gc_before_check %||% TRUE)) gc()
  r_used <- get_r_memory_used_gb()
  mac <- get_mac_memory_info_gb()
  log_msg("R memory used, GB:", round(r_used, 2))
  if (!is.na(mac$available_gb)) {
    log_msg("macOS memory available, GB:", round(mac$available_gb, 2), "| used, GB:", round(mac$used_gb, 2), "| total, GB:", round(mac$total_gb, 2))
  }
}

guard_memory <- function(label = "memory check", memory_guard = NULL, min_mac_available_gb = NULL) {
  memory_guard <- memory_guard %||% list(enabled = FALSE)
  if (!isTRUE(memory_guard$enabled)) return(invisible(TRUE))

  if (isTRUE(memory_guard$gc_before_check %||% TRUE)) gc()
  r_used <- get_r_memory_used_gb()
  mac <- get_mac_memory_info_gb()

  log_msg("Memory guard at", label, "| R used GB:", round(r_used, 2))
  if (!is.na(mac$available_gb)) {
    log_msg("Memory guard at", label, "| macOS available GB:", round(mac$available_gb, 2))
  }

  max_r <- memory_guard$max_r_memory_gb %||% Inf
  min_mac <- min_mac_available_gb %||% memory_guard$min_mac_available_gb %||% -Inf

  if (is.finite(r_used) && r_used > max_r) {
    stop("Memory safeguard stopped pipeline at ", label, ": R memory used ", round(r_used, 2), " GB exceeds limit ", max_r, " GB.")
  }
  if (!is.na(mac$available_gb) && mac$available_gb < min_mac) {
    stop("Memory safeguard stopped pipeline at ", label, ": macOS available memory ", round(mac$available_gb, 2), " GB is below limit ", min_mac, " GB.")
  }
  invisible(TRUE)
}

safe_step <- function(step_name, expr, analysis_spec = NULL) {
  log_section(step_name)
  tryCatch({
    guard_memory(paste("before", step_name), analysis_spec$memory_guard %||% NULL)
    out <- force(expr)
    guard_memory(paste("after", step_name), analysis_spec$memory_guard %||% NULL)
    log_msg("SUCCESS:", step_name)
    log_memory(analysis_spec$memory_guard %||% NULL)
    out
  }, error = function(e) {
    msg <- paste("ERROR in", step_name, ":", conditionMessage(e))
    log_msg(msg)
    log_msg("Pipeline stopped.")
    error_flag <- get0("error_flag", ifnotfound = "pipeline_error.flag", envir = .GlobalEnv)
    writeLines(msg, error_flag)
    update_heartbeat(msg)
    stop(e)
  })
}

# ------------------------------------------------------------
# CmdStan / brms setup
# ------------------------------------------------------------

setup_brms_cmdstan <- function(cache_dir = file.path(path.expand("~"), ".cmdstanr-cache")) {
  cmdstan_path_main <- cmdstanr::cmdstan_path()
  cmdstanr::set_cmdstan_path(cmdstan_path_main)
  Sys.setenv(CMDSTAN = cmdstan_path_main)
  options(brms.backend = "cmdstanr")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  options(cmdstanr_write_stan_file_dir = cache_dir)
  invisible(list(cmdstan_path = cmdstan_path_main, cache_dir = cache_dir))
}

# ------------------------------------------------------------
# Basic validation and data helpers
# ------------------------------------------------------------

read_var_dict <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      across(where(is.character), stringr::str_trim),
      impute_target = as.logical(impute_target),
      use_in_model = as.logical(use_in_model),
      use_as_auxiliary = as.logical(use_as_auxiliary)
    )
}

check_required_vars <- function(df, vars, label = "variables") {
  vars <- vars[!is.na(vars) & nzchar(vars)]
  missing_vars <- setdiff(vars, names(df))
  if (length(missing_vars) > 0) {
    stop("Missing ", label, ": ", paste(missing_vars, collapse = ", "))
  }
}

prepare_raw_data <- function(raw_data, analysis_spec, var_dict) {
  out <- as_tibble(raw_data)
  row_id <- analysis_spec$data$row_id_var %||% "row_id"
  if (!row_id %in% names(out)) out[[row_id]] <- seq_len(nrow(out))

  binary_vars <- var_dict %>%
    filter(type == "binary") %>%
    pull(var) %>%
    intersect(names(out))

  categorical_vars <- var_dict %>%
    filter(type == "categorical") %>%
    pull(var) %>%
    intersect(names(out))

  ordinal_vars <- var_dict %>%
    filter(type == "ordinal") %>%
    pull(var) %>%
    intersect(names(out))

  for (v in binary_vars) {
    out[[v]] <- as.factor(out[[v]])
  }

  for (v in categorical_vars) {
    out[[v]] <- as.factor(out[[v]])
  }

  for (v in ordinal_vars) {
    out[[v]] <- ordered(out[[v]])
  }

  ref_tbl <- var_dict %>% filter(!is.na(reference), reference != "", var %in% names(out))
  if (nrow(ref_tbl) > 0) {
    for (i in seq_len(nrow(ref_tbl))) {
      v <- ref_tbl$var[i]
      ref <- as.character(ref_tbl$reference[i])
      if (is.factor(out[[v]]) && ref %in% levels(out[[v]])) {
        if (is.ordered(out[[v]])) {
          lev <- levels(out[[v]])
          out[[v]] <- ordered(out[[v]], levels = lev)
        } else {
          out[[v]] <- stats::relevel(out[[v]], ref = ref)
        }
      }
    }
  }
  out
}

make_z_stats <- function(df, scale_vars) {
  scale_vars <- intersect(scale_vars, names(df))
  if (length(scale_vars) == 0) return(tibble(variable = character(), center = numeric(), scale = numeric()))
  purrr::map_dfr(scale_vars, function(v) {
    x <- df[[v]]
    if (!is.numeric(x) && !is.integer(x)) stop("Cannot z-scale non-numeric variable: ", v)
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) stop("Cannot z-scale variable with invalid SD: ", v)
    tibble(variable = v, center = mean(x, na.rm = TRUE), scale = s)
  })
}

apply_z_stats <- function(df, z_stats) {
  out <- df
  if (is.null(z_stats) || nrow(z_stats) == 0) return(out)
  for (i in seq_len(nrow(z_stats))) {
    v <- z_stats$variable[i]
    out[[paste0(v, "_z")]] <- as.numeric((out[[v]] - z_stats$center[i]) / z_stats$scale[i])
  }
  out
}


# ------------------------------------------------------------
# Variable dictionary roles
# ------------------------------------------------------------

derive_analysis_variables <- function(var_dict, analysis_spec = NULL) {
  # 00_variable_dictionary.csv is the source of truth for variable roles,
  # types, timing, scaling, imputation targets, and model inclusion.
  #
  # analysis_spec$variables is optional. If supplied, it overrides selected
  # derived lists. This allows advanced users to make project-specific
  # adjustments without duplicating the full dictionary in 00_config.R.

  get_col <- function(df, col, default = NA_character_) {
    if (col %in% names(df)) {
      df[[col]]
    } else {
      rep(default, nrow(df))
    }
  }

  role <- get_col(var_dict, "role", "")
  type <- get_col(var_dict, "type", "")
  timing <- get_col(var_dict, "timing", "")
  scale <- get_col(var_dict, "scale", "no")
  use_as_auxiliary <- get_col(var_dict, "use_as_auxiliary", FALSE)

  use_as_auxiliary <- as.logical(use_as_auxiliary)
  use_as_auxiliary[is.na(use_as_auxiliary)] <- FALSE

  vars_from_dict <- list(
    exposure_vars = var_dict$var[role == "exposure"],
    covariate_vars = var_dict$var[role == "covariate"],
    auxiliary_vars = unique(var_dict$var[role == "auxiliary" | use_as_auxiliary]),
    continuous_vars = var_dict$var[type == "continuous"],
    categorical_vars = var_dict$var[type == "categorical"],
    ordinal_vars = var_dict$var[type == "ordinal"],
    subject_level_vars = var_dict$var[timing %in% c("single", "baseline")],
    time_varying_vars = var_dict$var[timing %in% c("repeated", "time_varying")],
    scale_vars = var_dict$var[scale == "z"]
  )

  vars_from_dict <- purrr::map(
    vars_from_dict,
    ~ unique(.x[!is.na(.x) & nzchar(.x)])
  )

  overrides <- NULL

  if (!is.null(analysis_spec) && !is.null(analysis_spec$variables)) {
    overrides <- analysis_spec$variables
  }

  if (is.null(overrides)) {
    return(vars_from_dict)
  }

  if (!is.list(overrides)) {
    stop("analysis_spec$variables must be NULL or a named list.")
  }

  unknown_override_names <- setdiff(names(overrides), names(vars_from_dict))

  if (length(unknown_override_names) > 0) {
    stop(
      "Unknown names in analysis_spec$variables: ",
      paste(unknown_override_names, collapse = ", "),
      ". Valid names are: ",
      paste(names(vars_from_dict), collapse = ", ")
    )
  }

  # For vectors, modifyList() replaces the whole vector, which is what we want.
  utils::modifyList(vars_from_dict, overrides)
}


# ------------------------------------------------------------
# Formula, family, priors
# ------------------------------------------------------------

make_brms_family <- function(outcome_spec) {
  family <- outcome_spec$family
  link <- outcome_spec$link

  # Use brmsfamily() for all supported families. This avoids namespace
  # problems such as brms::gaussian() not existing and keeps the object
  # type consistent for brms/cmdstanr.
  switch(
    family,
    gaussian = brms::brmsfamily("gaussian", link = link %||% "identity"),
    bernoulli = brms::brmsfamily("bernoulli", link = link %||% "logit"),
    poisson = brms::brmsfamily("poisson", link = link %||% "log"),
    negbinomial = brms::brmsfamily("negbinomial", link = link %||% "log"),
    beta = brms::brmsfamily("Beta", link = link %||% "logit"),
    ordinal = brms::brmsfamily("cumulative", link = link %||% "logit"),
    categorical = brms::brmsfamily("categorical", link = link %||% "logit"),
    cox = brms::brmsfamily("cox", link = link %||% "log"),
    stop("Unsupported brms family: ", family)
  )
}

make_default_priors <- function(analysis_spec, formula = NULL) {
  family <- analysis_spec$outcome$family
  priors <- switch(
    family,
    bernoulli = c(brms::prior(normal(0, 1.5), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept")),
    ordinal = c(brms::prior(normal(0, 1.5), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept")),
    categorical = c(brms::prior(normal(0, 1.5), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept")),
    gaussian = c(brms::prior(normal(0, 1), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept"), brms::prior(exponential(1), class = "sigma")),
    student = c(brms::prior(normal(0, 1), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept"), brms::prior(exponential(1), class = "sigma")),
    poisson = c(brms::prior(normal(0, 1), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept")),
    negbinomial = c(brms::prior(normal(0, 1), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept")),
    beta = c(brms::prior(normal(0, 1), class = "b"), brms::prior(student_t(3, 0, 2.5), class = "Intercept")),
    cox = c(brms::prior(normal(0, 1), class = "b")),
    stop("No default priors defined for family: ", family)
  )
  if (!is.null(formula) && grepl("\\|", paste(deparse(formula), collapse = " "))) {
    priors <- c(priors, brms::prior(exponential(1), class = "sd"))
  }
  priors
}

filter_priors_to_model <- function(priors, formula, data, family) {
  valid <- brms::get_prior(formula = formula, data = data, family = family)
  valid_classes <- unique(valid$class)
  priors_df <- as.data.frame(priors)
  keep <- priors_df$class %in% valid_classes
  dropped <- unique(priors_df$class[!keep])
  if (length(dropped) > 0) message("Dropping priors not used by this model: ", paste(dropped, collapse = ", "))
  priors[keep]
}


resolve_fixed_effects <- function(analysis_spec, var_dict, analysis_vars = NULL) {
  fixed <- analysis_spec$model$fixed_effects
  analysis_vars <- analysis_vars %||% derive_analysis_variables(var_dict, analysis_spec)

  if (length(fixed) == 1 && identical(fixed, "auto")) {
    # Default: use 00_variable_dictionary.csv.
    fixed_vars <- var_dict %>%
      dplyr::filter(
        .data$use_in_model,
        !.data$role %in% c("outcome", "binary_outcome", "id", "time")
      ) %>%
      dplyr::pull(.data$var)

    # Optional override: if users explicitly define exposure_vars and/or
    # covariate_vars in analysis_spec$variables, use those for auto fixed
    # effects. This keeps the dictionary as default while allowing advanced
    # project-specific overrides.
    if (!is.null(analysis_spec$variables) &&
        any(c("exposure_vars", "covariate_vars") %in% names(analysis_spec$variables))) {
      fixed_vars <- unique(c(
        analysis_vars$exposure_vars,
        analysis_vars$covariate_vars
      ))
    }
  } else {
    fixed_vars <- fixed
  }

  fixed_vars <- fixed_vars[!is.na(fixed_vars) & nzchar(fixed_vars)]

  scale_vars <- analysis_vars$scale_vars %||% character(0)

  fixed_vars <- ifelse(
    fixed_vars %in% scale_vars,
    paste0(fixed_vars, "_z"),
    fixed_vars
  )

  unique(fixed_vars)
}

make_brms_formula <- function(analysis_spec, var_dict) {
  if (!is.null(analysis_spec$model$custom_formula)) {
    return(analysis_spec$model$custom_formula)
  }

  y <- analysis_spec$outcome$y_var
  fixed <- resolve_fixed_effects(analysis_spec, var_dict)

  fixed_rhs <- if (length(fixed) == 0) {
    "1"
  } else {
    paste(fixed, collapse = " + ")
  }

  random_terms <- character(0)
  id_var <- analysis_spec$data$id_var

  if (isTRUE(analysis_spec$model$random_effects$subject_intercept) &&
      !is.null(id_var) &&
      !is.na(id_var) &&
      nzchar(id_var)) {
    random_terms <- c(random_terms, paste0("(1 | ", id_var, ")"))
  }

  slope_vars <- analysis_spec$model$random_effects$subject_slope_vars %||% character(0)
  slope_vars <- slope_vars[!is.na(slope_vars) & nzchar(slope_vars)]

  if (length(slope_vars) > 0 &&
      !is.null(id_var) &&
      !is.na(id_var) &&
      nzchar(id_var)) {
    random_terms <- c(
      random_terms,
      paste0("(", paste(slope_vars, collapse = " + "), " | ", id_var, ")")
    )
  }

  rhs <- paste(c(fixed_rhs, random_terms), collapse = " + ")

  stats::as.formula(paste(y, "~", rhs))
}


build_model_spec <- function(analysis_spec, var_dict, reference_data) {
  analysis_vars <- derive_analysis_variables(var_dict, analysis_spec)
  scale_vars <- analysis_vars$scale_vars %||% character(0)
  z_stats <- make_z_stats(reference_data, scale_vars)
  formula <- make_brms_formula(analysis_spec, var_dict)
  family <- make_brms_family(analysis_spec$outcome)
  prior <- if (identical(analysis_spec$model$priors, "default_weakly_regularizing")) {
    make_default_priors(analysis_spec, formula)
  } else {
    analysis_spec$model$priors
  }
  fixed_effects <- resolve_fixed_effects(analysis_spec, var_dict, analysis_vars)
  list(
    formula = formula,
    analysis_vars = analysis_vars,
    formula_vars = get_brms_formula_vars(formula),
    formula_text = get_brms_formula_text(formula),
    family = family,
    prior = prior,
    z_stats = z_stats,
    fixed_effects = fixed_effects,
    chains = analysis_spec$model$chains,
    iter = analysis_spec$model$iter,
    warmup = analysis_spec$model$warmup,
    seed = analysis_spec$model$seed,
    adapt_delta = analysis_spec$model$adapt_delta,
    max_treedepth = analysis_spec$model$max_treedepth,
    parameter_draw_regex = analysis_spec$model$parameter_draw_regex,
    allow_new_levels = analysis_spec$posterior_prediction$allow_new_levels,
    sample_new_levels = analysis_spec$posterior_prediction$sample_new_levels
  )
}


get_brms_formula_core <- function(formula) {
  if (inherits(formula, "brmsformula") && !is.null(formula$formula)) {
    return(formula$formula)
  }

  formula
}

get_brms_formula_text <- function(formula) {
  if (is.null(formula)) {
    return("not available")
  }

  paste(
    deparse(get_brms_formula_core(formula), width.cutoff = 500),
    collapse = " "
  ) %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()
}

get_brms_formula_vars <- function(formula) {
  if (is.null(formula)) {
    return(character(0))
  }

  unique(all.vars(get_brms_formula_core(formula)))
}


extract_special_term_vars <- function(formula, fun = c("s", "mo")) {
  fun <- match.arg(fun)

  if (is.null(formula)) {
    return(character(0))
  }

  formula_text <- get_brms_formula_text(formula)

  pattern <- paste0(
    "\\b",
    fun,
    "\\s*\\(\\s*`?([A-Za-z.][A-Za-z0-9._]*)`?"
  )

  matches <- gregexpr(pattern, formula_text, perl = TRUE)
  hits <- regmatches(formula_text, matches)[[1]]

  if (length(hits) == 0 || identical(hits, character(0))) {
    return(character(0))
  }

  unique(sub(pattern, "\\1", hits, perl = TRUE))
}


# ------------------------------------------------------------
# Imputation
# ------------------------------------------------------------

is_binary_like <- function(x) {
  ux <- unique(stats::na.omit(x))
  length(ux) <= 2
}

prepare_miceranger_data <- function(data) {
  # miceRanger/ranger can behave poorly with ordered factors.  For imputation,
  # treat ordered factors as ordinary factors, then restore the ordered class
  # after completeData().  This is especially important for ordinal predictors
  # that will later be used in brms::mo().
  ordered_levels <- purrr::imap(
    data,
    function(x, nm) {
      if (is.ordered(x)) {
        levels(x)
      } else {
        NULL
      }
    }
  )

  ordered_levels <- ordered_levels[!vapply(ordered_levels, is.null, logical(1))]

  out <- data

  for (v in names(ordered_levels)) {
    out[[v]] <- factor(
      as.character(out[[v]]),
      levels = ordered_levels[[v]],
      ordered = FALSE
    )
  }

  list(
    data = out,
    ordered_levels = ordered_levels
  )
}

restore_ordered_factors <- function(data, ordered_levels) {
  out <- tibble::as_tibble(data)

  if (length(ordered_levels) == 0) {
    return(out)
  }

  for (v in names(ordered_levels)) {
    if (v %in% names(out)) {
      out[[v]] <- ordered(
        as.character(out[[v]]),
        levels = ordered_levels[[v]]
      )
    }
  }

  out
}

make_row_level_imputation_spec <- function(data, analysis_spec, var_dict) {
  row_id <- analysis_spec$data$row_id_var
  y <- analysis_spec$outcome$y_var

  target_vars <- var_dict %>% filter(impute_target) %>% pull(var)
  if (isTRUE(analysis_spec$imputation$impute_y)) target_vars <- unique(c(target_vars, y))
  target_vars <- intersect(target_vars, names(data))
  target_vars <- target_vars[vapply(target_vars, function(v) anyNA(data[[v]]), logical(1))]

  analysis_vars <- derive_analysis_variables(var_dict, analysis_spec)

  base_predictors <- var_dict %>%
    filter(use_in_model | use_as_auxiliary | impute_target) %>%
    pull(var) %>%
    c(analysis_vars$auxiliary_vars %||% character(0)) %>%
    unique() %>%
    intersect(names(data))

  # Exclude variables with missingness unless they are imputation targets.
  usable_predictors <- base_predictors[vapply(base_predictors, function(v) !anyNA(data[[v]]) || v %in% target_vars, logical(1))]

  vars <- stats::setNames(
    lapply(target_vars, function(v) setdiff(usable_predictors, c(row_id, v))),
    target_vars
  )

  list(
    m = analysis_spec$imputation$m,
    maxiter = analysis_spec$imputation$maxiter,
    verbose = analysis_spec$imputation$verbose %||% FALSE,
    vars = vars,
    mean_match_k = analysis_spec$imputation$mean_match_k,
    seed = analysis_spec$imputation$seed
  )
}

run_row_level_imputation <- function(data, imputation_spec, analysis_spec) {
  vars <- imputation_spec$vars
  m <- imputation_spec$m

  if (length(vars) == 0) {
    log_msg("No imputation targets with missingness. Duplicating data m times.")
    return(rep(list(tibble::as_tibble(data)), m))
  }

  # Make the vars list robust before it reaches miceRanger.
  vars <- vars[names(vars) %in% names(data)]
  vars <- purrr::map(vars, ~ intersect(.x, names(data)))

  required_imputation_vars <- unique(c(names(vars), unlist(vars, use.names = FALSE)))
  missing_imputation_vars <- setdiff(required_imputation_vars, names(data))

  if (length(missing_imputation_vars) > 0) {
    stop(
      "The imputation specification refers to variable(s) not present in the imputation data: ",
      paste(missing_imputation_vars, collapse = ", "),
      ". Check 00_variable_dictionary.csv timing/type/impute_target settings."
    )
  }

  if (length(vars) == 0) {
    log_msg("No valid imputation targets remain after checking names. Duplicating data m times.")
    return(rep(list(tibble::as_tibble(data)), m))
  }

  mice_data_info <- prepare_miceranger_data(data)
  mice_data <- mice_data_info$data

  ordered_vars <- names(mice_data_info$ordered_levels)

  if (length(ordered_vars) > 0) {
    log_msg(
      "Temporarily treating ordered factors as unordered factors for miceRanger:",
      paste(ordered_vars, collapse = ", ")
    )
  }

  valueSelector <- vapply(names(vars), function(v) {
    x <- mice_data[[v]]
    if (is.factor(x)) {
      "value"
    } else if (is.numeric(x) && !is_binary_like(x)) {
      "meanMatch"
    } else {
      "value"
    }
  }, character(1))
  names(valueSelector) <- names(vars)

  mm_vars <- names(vars)[valueSelector == "meanMatch"]
  meanMatchCandidates <- rep(imputation_spec$mean_match_k %||% 5, length(mm_vars))
  names(meanMatchCandidates) <- mm_vars

  impute_workers <- as.integer(analysis_spec$parallel$impute_workers %||% 1L)

  num_threads_per_worker <- as.integer(
    analysis_spec$parallel$num_impute_threads_per_worker %||%
      analysis_spec$parallel$num_impute_threads %||%
      1L
  )

  use_parallel_imputation <- isTRUE(impute_workers > 1L)

  log_msg("miceRanger imputation targets:", length(names(vars)))
  log_msg("miceRanger m:", m)
  log_msg("miceRanger maxiter:", imputation_spec$maxiter)
  log_msg("miceRanger parallel:", use_parallel_imputation)
  log_msg("miceRanger impute_workers:", impute_workers)
  log_msg("miceRanger num_impute_threads_per_worker:", num_threads_per_worker)
  log_msg(
    "Approximate imputation CPU threads:",
    impute_workers * num_threads_per_worker
  )

  cl <- NULL

  if (use_parallel_imputation) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel miceRanger imputation.")
    }

    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'foreach' is required for parallel miceRanger imputation.")
    }

    cl <- parallel::makeCluster(impute_workers, type = "PSOCK")
    doParallel::registerDoParallel(cl)

    on.exit(
      {
        try(parallel::stopCluster(cl), silent = TRUE)
        try(foreach::registerDoSEQ(), silent = TRUE)
      },
      add = TRUE
    )
  }

  # miceRanger() has no seed= argument; it consumes whatever the global RNG
  # state happens to be. Seed explicitly so a given imputation_spec$seed
  # always reproduces the same batch of m datasets.
  if (!is.null(imputation_spec$seed)) {
    set.seed(imputation_spec$seed)
  }

  mice_obj <- miceRanger::miceRanger(
    data = mice_data,
    m = m,
    maxiter = imputation_spec$maxiter,
    vars = vars,
    valueSelector = valueSelector,
    meanMatchCandidates = meanMatchCandidates,
    returnModels = FALSE,
    parallel = use_parallel_imputation,
    verbose = imputation_spec$verbose %||% FALSE,
    num.threads = num_threads_per_worker
  )

  miceRanger::completeData(mice_obj, verbose = FALSE) %>%
    purrr::map(
      ~ restore_ordered_factors(
        data = .x,
        ordered_levels = mice_data_info$ordered_levels
      )
    )
}



# ------------------------------------------------------------
# Subject-wide imputation for repeated outcome data
# ------------------------------------------------------------

make_subject_wide_imputation_data <- function(data, analysis_spec, var_dict) {
  id_var <- analysis_spec$data$id_var
  time_var <- analysis_spec$data$time_var
  row_id_var <- analysis_spec$data$row_id_var
  y_var <- analysis_spec$outcome$y_var
  y_prefix <- analysis_spec$outcome$y_prefix %||% paste0(y_var, "_")

  check_required_vars(
    data,
    c(id_var, time_var, row_id_var, y_var),
    "subject-wide imputation core variables"
  )

  subject_vars <- var_dict %>%
    dplyr::filter(
      .data$timing %in% c("single", "baseline"),
      !.data$role %in% c("id", "time", "outcome", "binary_outcome")
    ) %>%
    dplyr::pull(.data$var) %>%
    unique() %>%
    intersect(names(data))

  subject_data <- data %>%
    dplyr::group_by(.data[[id_var]]) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(subject_vars), dplyr::first),
      .groups = "drop"
    )

  y_wide <- data %>%
    dplyr::select(dplyr::all_of(c(id_var, time_var, y_var))) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of(time_var),
      values_from = dplyr::all_of(y_var),
      names_prefix = y_prefix
    )

  subject_wide <- subject_data %>%
    dplyr::left_join(y_wide, by = id_var)

  long_base_vars <- var_dict %>%
    dplyr::filter(.data$timing %in% c("repeated", "time_varying")) %>%
    dplyr::pull(.data$var) %>%
    unique() %>%
    intersect(names(data))

  long_base_vars <- unique(c(row_id_var, id_var, time_var, y_var, long_base_vars))

  long_base <- data %>%
    dplyr::select(dplyr::all_of(long_base_vars))

  list(
    subject_wide = subject_wide,
    long_base = long_base,
    subject_vars = subject_vars
  )
}

make_subject_wide_imputation_spec <- function(subject_wide, analysis_spec, var_dict) {
  id_var <- analysis_spec$data$id_var
  y_wide_regex <- analysis_spec$outcome$y_wide_regex %||%
    paste0("^", analysis_spec$outcome$y_var, "_")

  y_wide_cols <- grep(y_wide_regex, names(subject_wide), value = TRUE)

  target_vars <- var_dict %>%
    dplyr::filter(.data$impute_target) %>%
    dplyr::pull(.data$var) %>%
    unique()

  if (isTRUE(analysis_spec$imputation$impute_y)) {
    target_vars <- unique(c(target_vars, y_wide_cols))
  }

  exclude_targets <- unique(c(
    id_var,
    if (!isTRUE(analysis_spec$imputation$impute_y)) y_wide_cols else character(0),
    analysis_spec$imputation$extra_exclude_targets %||% character(0)
  ))

  target_vars <- setdiff(target_vars, exclude_targets)
  target_vars <- intersect(target_vars, names(subject_wide))
  target_vars <- target_vars[
    vapply(target_vars, function(v) anyNA(subject_wide[[v]]), logical(1))
  ]

  analysis_vars <- derive_analysis_variables(var_dict, analysis_spec)

  base_predictors <- var_dict %>%
    dplyr::filter(.data$use_in_model | .data$use_as_auxiliary | .data$impute_target) %>%
    dplyr::pull(.data$var) %>%
    c(analysis_vars$auxiliary_vars %||% character(0), y_wide_cols) %>%
    unique() %>%
    intersect(names(subject_wide))

  usable_predictors <- base_predictors[
    vapply(
      base_predictors,
      function(v) !anyNA(subject_wide[[v]]) || v %in% target_vars,
      logical(1)
    )
  ]

  vars <- stats::setNames(
    lapply(target_vars, function(v) setdiff(usable_predictors, c(id_var, v))),
    target_vars
  )

  list(
    m = analysis_spec$imputation$m,
    maxiter = analysis_spec$imputation$maxiter,
    verbose = analysis_spec$imputation$verbose %||% FALSE,
    vars = vars,
    mean_match_k = analysis_spec$imputation$mean_match_k,
    seed = analysis_spec$imputation$seed
  )
}

prepare_long_imputed_from_subject_wide <- function(
    imputed_wide_list,
    long_base,
    analysis_spec
) {
  id_var <- analysis_spec$data$id_var
  y_wide_regex <- analysis_spec$outcome$y_wide_regex %||%
    paste0("^", analysis_spec$outcome$y_var, "_")

  purrr::map(
    imputed_wide_list,
    function(wide_i) {
      subject_only <- wide_i %>%
        dplyr::select(-tidyselect::matches(y_wide_regex))

      long_base %>%
        dplyr::left_join(subject_only, by = id_var)
    }
  )
}


# ------------------------------------------------------------
# Model data and fitting
# ------------------------------------------------------------

prepare_model_data_files <- function(imputation_manifest, analysis_spec, model_spec, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  y <- analysis_spec$outcome$y_var
  row_id <- analysis_spec$data$row_id_var

  manifest <- imputation_manifest %>%
    mutate(
      analysis_file = file.path(out_dir, sprintf("analysis_imp_%03d.rds", imputation)),
      pred_file = file.path(out_dir, sprintf("pred_imp_%03d.rds", imputation))
    )

  for (i in seq_len(nrow(manifest))) {
    dat_i <- readRDS(manifest$imputed_file[i])
    dat_i <- apply_z_stats(dat_i, model_spec$z_stats)

    model_vars <- model_spec$formula_vars %||% get_brms_formula_vars(model_spec$formula)
    analysis_vars <- unique(c(model_vars, row_id))
    analysis_vars <- analysis_vars[!is.na(analysis_vars) & nzchar(analysis_vars)]
    pred_vars <- unique(c(setdiff(model_vars, y), row_id))
    pred_vars <- pred_vars[!is.na(pred_vars) & nzchar(pred_vars)]

    check_required_vars(dat_i, analysis_vars, "analysis variables")
    check_required_vars(dat_i, pred_vars, "prediction variables")

    analysis_i <- dat_i %>% filter(!is.na(.data[[y]])) %>% select(all_of(analysis_vars))
    pred_i <- dat_i %>% filter(is.na(.data[[y]])) %>% select(all_of(pred_vars))

    saveRDS(analysis_i, manifest$analysis_file[i], compress = FALSE)
    saveRDS(pred_i, manifest$pred_file[i], compress = FALSE)
    rm(dat_i, analysis_i, pred_i); gc()
  }
  manifest
}

rds_ok <- function(file) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({ obj <- readRDS(file); rm(obj); gc(); TRUE }, error = function(e) FALSE)
}

fit_one_brm_file <- function(ii, fit_manifest, model_spec, analysis_spec, setup_info, log_file = NULL) {
  imp_i <- fit_manifest$imputation[ii]
  fit_file_i <- fit_manifest$fit_file[ii]
  analysis_file_i <- fit_manifest$analysis_file[ii]

  worker_log <- function(...) {
    msg <- paste(..., sep = " ")
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
    if (!is.null(log_file)) cat(line, "\n", file = log_file, append = TRUE)
    cat(line, "\n")
    flush.console()
    invisible(line)
  }

  worker_log("Worker started for imputation", imp_i)
  worker_log("Analysis file:", analysis_file_i)
  worker_log("Fit file:", fit_file_i)

  if (rds_ok(fit_file_i)) {
    worker_log("Existing valid fit found; skipping", imp_i)
    return(tibble(imputation = imp_i, status = "skipped_existing_valid_fit", fit_file = fit_file_i))
  }

  if (file.exists(fit_file_i) && !rds_ok(fit_file_i)) {
    worker_log("Existing fit file is invalid; removing", fit_file_i)
    file.remove(fit_file_i)
  }

  worker_log("Setting CmdStan path")
  cmdstanr::set_cmdstan_path(setup_info$cmdstan_path)
  options(brms.backend = "cmdstanr")
  options(cmdstanr_write_stan_file_dir = setup_info$cache_dir)

  worker_log("Reading analysis data")
  dat_i <- readRDS(analysis_file_i)
  worker_log("Analysis rows:", nrow(dat_i), "columns:", ncol(dat_i))
  worker_log("Formula:", paste(deparse(model_spec$formula), collapse = " "))
  worker_log("Family:", paste(capture.output(print(model_spec$family)), collapse = " | "))

  if (anyNA(dat_i)) {
    na_counts <- vapply(dat_i, function(x) sum(is.na(x)), integer(1))
    na_counts <- na_counts[na_counts > 0]
    worker_log("WARNING: analysis data contains missing values:", paste(names(na_counts), na_counts, sep = "=", collapse = ", "))
  }

  worker_log("Filtering priors to actual model")
  prior_i <- filter_priors_to_model(model_spec$prior, model_spec$formula, dat_i, model_spec$family)
  worker_log("Prior classes used:", paste(unique(as.data.frame(prior_i)$class), collapse = ", "))

  worker_log("Starting brm()")
  fit_i <- tryCatch(
    brms::brm(
      formula = model_spec$formula,
      data = dat_i,
      family = model_spec$family,
      prior = prior_i,
      chains = model_spec$chains,
      iter = model_spec$iter,
      warmup = model_spec$warmup,
      cores = analysis_spec$parallel$cores_per_fit %||% 1,
      seed = model_spec$seed + imp_i,
      init = analysis_spec$model$init %||% "random",
      refresh = analysis_spec$model$refresh %||% 10,
      silent = analysis_spec$model$silent %||% 0,
      backend = "cmdstanr",
      control = list(adapt_delta = model_spec$adapt_delta, max_treedepth = model_spec$max_treedepth)
    ),
    error = function(e) {
      worker_log("ERROR during brm() for imputation", imp_i, ":", conditionMessage(e))
      stop(e)
    }
  )

  worker_log("brm() finished; saving fit")
  saveRDS(fit_i, fit_file_i, compress = FALSE)
  rm(dat_i, fit_i)
  gc()
  worker_log("Worker completed for imputation", imp_i)

  tibble(imputation = imp_i, status = "completed", fit_file = fit_file_i)
}

# ------------------------------------------------------------
# Diagnostics / posterior summaries / predictions
# ------------------------------------------------------------

diagnose_one_fit <- function(fit_file, imputation, max_treedepth) {
  fit_i <- readRDS(fit_file)
  np <- brms::nuts_params(fit_i)
  out <- tibble(
    imputation = imputation,
    divergent = sum(np$Parameter == "divergent__" & np$Value == 1, na.rm = TRUE),
    treedepth_hits = sum(np$Parameter == "treedepth__" & np$Value >= max_treedepth, na.rm = TRUE)
  )
  rm(fit_i, np); gc()
  out
}

summarise_parameter_draws_table <- function(draws_df, summary_spec) {
  meta_cols <- c("imputation", ".chain", ".iteration", ".draw")
  param_cols <- setdiff(names(draws_df), meta_cols)
  param_cols <- grep("^(b_|sd_|sigma)", param_cols, value = TRUE)
  if (identical(summary_spec$effects, "fixed")) param_cols <- grep("^b_", param_cols, value = TRUE)
  if (length(param_cols) == 0) stop("No parameter columns found for posterior summary.")

  args <- list(
    x = draws_df[, param_cols, drop = FALSE],
    centrality = summary_spec$centrality,
    ci = summary_spec$ci,
    ci_method = summary_spec$ci_method,
    verbose = FALSE
  )
  if (!is.null(summary_spec$test) && length(summary_spec$test) > 0) args$test <- summary_spec$test
  if (!is.null(summary_spec$rope$fixed_range)) args$rope_range <- summary_spec$rope$fixed_range

  suppressWarnings(do.call(bayestestR::describe_posterior, args)) %>% as_tibble()
}

predict_missing_y_draws_one <- function(fit_i, pred_i, analysis_spec, model_spec, ndraws = 1000) {
  row_id <- analysis_spec$data$row_id_var
  if (nrow(pred_i) == 0) return(tibble())
  new_i <- pred_i %>% select(-all_of(row_id))
  if (anyNA(new_i)) stop("Prediction rows still contain missing predictors.")

  draw_mat <- brms::posterior_predict(
    fit_i,
    newdata = new_i,
    ndraws = ndraws,
    allow_new_levels = model_spec$allow_new_levels %||% TRUE,
    sample_new_levels = model_spec$sample_new_levels %||% "gaussian"
  ) %>% as.matrix()

  as_tibble(draw_mat, .name_repair = "minimal") %>%
    stats::setNames(paste0("row_", pred_i[[row_id]])) %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(cols = -.draw, names_to = row_id, values_to = "y_draw") %>%
    mutate(!!row_id := as.integer(sub("^row_", "", .data[[row_id]])))
}

summarise_missing_y_draws <- function(pred_draws_long, analysis_spec, summary_spec) {
  row_id <- analysis_spec$data$row_id_var
  if (nrow(pred_draws_long) == 0) return(tibble())
  pred_draws_long %>%
    group_by(.data[[row_id]]) %>%
    group_split() %>%
    map_dfr(function(df_i) {
      suppressWarnings(bayestestR::describe_posterior(
        df_i$y_draw,
        centrality = summary_spec$centrality,
        ci = summary_spec$ci,
        ci_method = summary_spec$ci_method,
        verbose = FALSE
      )) %>% as_tibble() %>% mutate(!!row_id := unique(df_i[[row_id]]), .before = 1)
    })
}

# ------------------------------------------------------------
# Imputation-count stability: lightweight two-batch check
# ------------------------------------------------------------
#
# A minimal, dependency-light comparison restricted to exactly the two batch
# sizes supplied, used by run_all.R's automatic m-increment loop to decide
# whether to stop fitting more imputations. The full multi-batch diagnostic
# report (11_check_imputation_stability.R) still runs once, at whatever m
# the loop settles on, and is the source of truth for the published report;
# this function only gates the loop and is not itself a publication output.
evaluate_mi_stability_batches <- function(paths, analysis_spec, m_previous, m_final) {
  cfg <- analysis_spec$mi_stability %||% list()

  parameter_regex <- cfg$parameter_regex %||% "^b_"
  primary_parameters <- cfg$primary_parameters %||% NULL
  exclude_intercept <- cfg$exclude_intercept %||% TRUE

  ci <- analysis_spec$summary$ci %||% 0.95
  alpha <- (1 - ci) / 2

  estimate_tolerance <- cfg$estimate_tolerance %||% 0.05
  ci_endpoint_tolerance <- cfg$ci_endpoint_tolerance %||% 0.05
  pd_tolerance <- cfg$pd_tolerance %||% 0.02

  manifest_file <- file.path(paths$objects, "parameter_manifest.rds")

  if (!rds_ok(manifest_file)) {
    stop("parameter_manifest.rds not found. Run Step 6 first.")
  }

  pm <- readRDS(manifest_file)
  file_col <- intersect(
    c("parameter_draw_file", "draw_file", "file", "parameter_file"),
    names(pm)
  )[1]

  draw_manifest <- pm %>%
    dplyr::transmute(
      imputation = as.integer(.data$imputation),
      parameter_draw_file = as.character(.data[[file_col]])
    ) %>%
    dplyr::filter(!is.na(.data$imputation)) %>%
    dplyr::filter(purrr::map_lgl(.data$parameter_draw_file, rds_ok)) %>%
    dplyr::arrange(.data$imputation)

  if (nrow(draw_manifest) < m_final) {
    stop(
      "Requested batch size ", m_final,
      " exceeds the number of valid per-imputation draw files (",
      nrow(draw_manifest), ")."
    )
  }

  first_draws <- readRDS(draw_manifest$parameter_draw_file[1])
  meta_cols <- c("imputation", ".chain", ".iteration", ".draw")
  parameter_cols <- setdiff(names(first_draws), meta_cols)
  parameter_cols <- parameter_cols[vapply(first_draws[parameter_cols], is.numeric, logical(1))]

  if (!is.null(primary_parameters) && length(primary_parameters) > 0) {
    parameter_cols <- intersect(parameter_cols, primary_parameters)
  } else if (!is.null(parameter_regex) && nzchar(parameter_regex)) {
    parameter_cols <- parameter_cols[stringr::str_detect(parameter_cols, parameter_regex)]
  }

  if (isTRUE(exclude_intercept)) {
    parameter_cols <- parameter_cols[
      !parameter_cols %in% c("b_Intercept", "Intercept", "(Intercept)")
    ]
  }

  rm(first_draws)

  if (length(parameter_cols) == 0) {
    stop("No parameter columns selected for the imputation-stability check.")
  }

  summarise_one <- function(x) {
    x <- x[is.finite(x)]

    if (length(x) == 0) {
      return(c(Median = NA_real_, CI_low = NA_real_, CI_high = NA_real_, pd = NA_real_))
    }

    c(
      Median = stats::median(x),
      CI_low = as.numeric(stats::quantile(x, alpha, names = FALSE)),
      CI_high = as.numeric(stats::quantile(x, 1 - alpha, names = FALSE)),
      pd = max(mean(x > 0), mean(x < 0))
    )
  }

  summarise_batch <- function(batch_n) {
    batch_files <- draw_manifest$parameter_draw_file[seq_len(batch_n)]

    batch_draws <- purrr::map_dfr(batch_files, function(f) {
      d <- readRDS(f)
      d[, intersect(names(d), parameter_cols), drop = FALSE]
    })

    purrr::map_dfr(parameter_cols, function(p) {
      s <- summarise_one(batch_draws[[p]])
      tibble::tibble(
        Parameter = p,
        Median = s[["Median"]],
        CI_low = s[["CI_low"]],
        CI_high = s[["CI_high"]],
        pd = s[["pd"]]
      )
    })
  }

  summary_previous <- summarise_batch(m_previous)
  summary_final <- summarise_batch(m_final)

  comparison <- summary_final %>%
    dplyr::rename_with(~ paste0(.x, "_final"), -Parameter) %>%
    dplyr::left_join(
      summary_previous %>% dplyr::rename_with(~ paste0(.x, "_previous"), -Parameter),
      by = "Parameter"
    ) %>%
    dplyr::mutate(
      abs_Median_change = abs(.data$Median_final - .data$Median_previous),
      max_abs_CI_endpoint_change = pmax(
        abs(.data$CI_low_final - .data$CI_low_previous),
        abs(.data$CI_high_final - .data$CI_high_previous),
        na.rm = TRUE
      ),
      abs_pd_change = abs(.data$pd_final - .data$pd_previous),
      stable_by_thresholds =
        .data$abs_Median_change <= estimate_tolerance &
        .data$max_abs_CI_endpoint_change <= ci_endpoint_tolerance &
        .data$abs_pd_change <= pd_tolerance
    )

  list(
    m_previous = m_previous,
    m_final = m_final,
    n_parameters = nrow(comparison),
    n_stable = sum(comparison$stable_by_thresholds, na.rm = TRUE),
    all_stable = nrow(comparison) > 0 && all(comparison$stable_by_thresholds, na.rm = TRUE),
    detail = comparison
  )
}

# ------------------------------------------------------------
# Optional runtime overrides
# ------------------------------------------------------------
#
# Every numbered pipeline script re-sources 00_config.R as its first step,
# and 00_config.R itself is routinely replaced wholesale per example/study
# (see test/test_example_common.sh's prepare_*_example() helpers, which copy
# an example's own 00_config_*.R over the project's 00_config.R). So an
# override hook cannot live in 00_config.R itself without being silently
# lost the moment a different config file replaces it.
#
# This file, 00_common_functions.R, is shared and never replaced per-example,
# and is always sourced immediately after 00_config.R -- so analysis_spec
# and paths already exist in the calling environment by the time this runs.
# That makes it the right place for an orchestrator (currently only
# run_all.R's automatic m-increment loop) to persist a small override to
# disk and have it reapplied after every fresh source() of 00_config.R.
# Absent the override file, behaviour is unchanged from a plain config read.
#
# The override file holds a nested list keyed by top-level analysis_spec
# section, e.g. list(imputation = list(m = 8), model = list(run_smoke_fit
# = FALSE)), so an orchestrator can patch fields in more than one section
# (e.g. disabling the repeated smoke fit across loop batches, in addition
# to incrementing m).
mi_runtime_override_file <- file.path(paths$objects, "mi_runtime_override.rds")

if (file.exists(mi_runtime_override_file)) {
  mi_runtime_override <- readRDS(mi_runtime_override_file)

  for (mi_override_section in names(mi_runtime_override)) {
    section_overrides <- mi_runtime_override[[mi_override_section]]

    for (mi_override_name in names(section_overrides)) {
      analysis_spec[[mi_override_section]][[mi_override_name]] <- section_overrides[[mi_override_name]]
    }
  }

  rm(mi_runtime_override)
}

rm(mi_runtime_override_file)
