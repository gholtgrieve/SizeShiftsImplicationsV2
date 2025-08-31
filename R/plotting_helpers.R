#' Internal plotting helpers (not exported)
#' @keywords internal

#' Helper: Standardize scenario labels and factor orders (internal)
#' Ensures consistent printed labels across all figures and summaries.
#' - **factorMSY**: `0.75` → "liberal", `1` → "MSY", `1.5` → "precautionary"
#' - **trends**:
#'   - "no trends" (and "No trends") → "no trends"
#'   - "age-sex-length trends" → "ASL trends stabilized"
#'   - "age-length trends" → "AL trends stabilized"
#'   - "continuing trends" → "ASL trends continued"
#' - **mgmt**: map common long-form to short codes: "smsy_goal" → "TRM", "s_eq_goal" → "YPR", "smsy_dlm_goal" → "DLM"
#' - **selectivity**: standardize to "small-mesh", "unselective", "large-mesh"
#' If a column is absent, it is ignored. Factors are (re)leveled for stable facet/legend order:
#' - `factorMSY`: `c("liberal","MSY","precautionary")`
#' - `trends`:    `c("no trends","ASL trends stabilized","AL trends stabilized","ASL trends continued")`
#' - `mgmt`:      `c("TRM","YPR","DLM")` (plus any other levels appended at end, if present)
#' - `selectivity`: `c("small-mesh","unselective","large-mesh")`
#' @param df A data.frame/tibble with any of the columns above.
#' @return The same data.frame with standardized character values and ordered factors.

.standardize_scenario_labels <- function(df) {
  # factorMSY --------------------------------------------------------------
  if ("factorMSY" %in% names(df)) {
    fmsy_chr <- trimws(as.character(df$factorMSY))
    fmsy_chr <- tolower(fmsy_chr)
    fmsy_chr <- dplyr::case_when(
      fmsy_chr %in% c("0.75", "0.750", "liberal") ~ "liberal",
      fmsy_chr %in% c("1", "1.0", "msy") ~ "MSY",
      fmsy_chr %in% c("1.5", "1.50", "precautionary") ~ "precautionary",
      TRUE ~ fmsy_chr
    )
    df$factorMSY <- factor(fmsy_chr, levels = c("liberal", "MSY", "precautionary"))
  }

  # trends -----------------------------------------------------------------
  if ("trends" %in% names(df)) {
    tr_chr <- trimws(as.character(df$trends))
    tr_chr <- tolower(tr_chr)
    tr_chr <- dplyr::case_when(
      tr_chr %in% c("no trends") ~ "no trends",
      tr_chr %in% c("age-sex-length trends") ~ "ASL trends stabilized",
      tr_chr %in% c("age-length trends") ~ "AL trends stabilized",
      tr_chr %in% c("continuing trends") ~ "ASL trends continued",
      TRUE ~ tr_chr
    )
    df$trends <- factor(
      tr_chr,
      levels = c("no trends", "ASL trends stabilized", "AL trends stabilized", "ASL trends continued")
    )
  }

  # mgmt -------------------------------------------------------------------
  if ("mgmt" %in% names(df)) {
    mg_chr <- trimws(as.character(df$mgmt))
    mg_chr <- tolower(mg_chr)
    mg_chr <- dplyr::case_when(
      mg_chr %in% c("smsy_goal", "trm") ~ "TRM",
      mg_chr %in% c("s_eq_goal", "ypr") ~ "YPR",
      mg_chr %in% c("smsy_dlm_goal", "dlm") ~ "DLM",
      TRUE ~ mg_chr
    )
    base_levels <- c("TRM", "YPR", "DLM")
    mg_levels <- unique(c(base_levels, setdiff(mg_chr, base_levels)))
    df$mgmt <- factor(mg_chr, levels = mg_levels)
  }

  # selectivity ------------------------------------------------------------
  if ("selectivity" %in% names(df)) {
    sel_chr <- trimws(as.character(df$selectivity))
    sel_chr <- tolower(sel_chr)
    sel_chr <- dplyr::case_when(
      sel_chr %in% c("6 inch gillnet", "small-mesh", "small mesh") ~ "small-mesh",
      sel_chr %in% c("unselective") ~ "unselective",
      sel_chr %in% c("8.5 inch gillnet", "large-mesh", "large mesh") ~ "large-mesh",
      TRUE ~ sel_chr
    )
    df$selectivity <- factor(sel_chr, levels = c("small-mesh", "unselective", "large-mesh"))
  }

  df
}


# Standardize scenario labels and add IDs that match results list names
.scenarios_from_ssirun <- function(data) {
  # 'data' must be an ssi_run object, so 'data$scenarios' is expected to be present.
  df <- .standardize_scenario_labels(tibble::as_tibble(data$scenarios))

  # row-order id (useful generally)
  df$scen <- paste0("scenario_", seq_len(nrow(df)))

  # ensure scen_num and scen_key ("scen_<num>" used in run$results$S_msy & sr_sim)
  if (!"scen_num" %in% names(df)) df$scen_num <- seq_len(nrow(df))
  df$scen_key <- paste0("scen_", df$scen_num)

  df
}

# Create ./<file_basename>/ under output_dir and return that folder path
.ensure_outdir <- function(output_dir, file_basename) {
  out <- file.path(output_dir, file_basename)
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  out
}

# Safe pull from parameters with default
.get_param <- function(data, name, default = NULL) {
  val <- tryCatch(data$parameters[[name]], error = function(e) NULL)
  if (is.null(val)) default else val
}

.get_nyh <- function(data, default = 50L) as.integer(.get_param(data, "nyh", default))
.get_nyi <- function(data, default = 50L) as.integer(.get_param(data, "nyi", default))
.get_ny  <- function(data, default = NA_integer_) as.integer(.get_param(data, "ny", default))
.get_goalfreq <- function(data, default = NA_integer_) as.integer(.get_param(data, "goalfreq", default))
.get_firstrev <- function(data, default = 20L) as.integer(.get_param(data, "firstrev", default)) # run_model uses 20

# Nested obs list: [[scenario]][[iteration]] -> data.frame with obsEsc/obsHarv/obsRet
.get_obs_list <- function(data) {
  obs_list <- data$results$obs
  if (is.null(obs_list) || !length(obs_list)) stop("No observation lists found at data$results$obs.")
  obs_list
}

# SR params list (alpha, beta, ...) per scenario/iteration (data.frame per iter)
.get_sr_params_list <- function(data) {
  sr <- data$results$sr_sim
  if (is.null(sr)) stop("SR parameter list not found at data$results$sr_sim.")
  sr
}

# S_MSY series per scenario/iteration (numeric vector over review years)
.get_smsy_series_list <- function(data) {
  sm <- data$results$S_msy
  if (is.null(sm)) stop("S_MSY series not found at data$results$S_msy.")
  sm
}

# Review years as in run_model(): goalrev <- seq(nyi + firstrev, ny, by = goalfreq)
.get_review_years <- function(data, series_length = NULL) {
  nyi      <- .get_nyi(data, 50L)
  ny       <- .get_ny(data, NA_integer_)
  goalfreq <- .get_goalfreq(data, NA_integer_)
  firstrev <- .get_firstrev(data, 20L)

  if (is.finite(nyi) && is.finite(ny) && is.finite(goalfreq) && goalfreq > 0) {
    return(seq(from = nyi + firstrev, to = ny, by = goalfreq))
  }

  # Fallback: 1:N if we can’t construct calendar years
  if (is.null(series_length)) stop("Cannot infer review_years: provide series_length.")
  seq_len(series_length)
}

# Index of the first review *at or after* end of history.
# In run_model(), calendar review years are absolute (e.g., nyi+20, ...).
# The “post-historical review” target is at calendar year (nyi + nyh).
.get_review_index_post_history <- function(data, review_years) {
  nyi <- .get_nyi(data, 50L)
  nyh <- .get_nyh(data, 50L)
  target <- nyi + nyh
  idx <- which(review_years == target)
  if (length(idx) == 0L) {
    idx <- suppressWarnings(min(which(review_years >= target)))
  }
  if (!is.finite(idx)) idx <- length(review_years) # fallback to last
  idx
}

# Per-iteration mean over a window (utility for other figs)
.iter_window_means <- function(obs_list, year_index, col) {
  nscen <- length(obs_list)
  niter <- length(obs_list[[1L]])
  mat <- matrix(NA_real_, nrow = nscen, ncol = niter)
  for (j in seq_len(nscen)) {
    for (k in seq_len(niter)) {
      obs <- obs_list[[j]][[k]]
      if (!is.data.frame(obs) || !nrow(obs) || !(col %in% names(obs))) next
      v <- obs[[col]][year_index]
      mat[j, k] <- mean(v, na.rm = TRUE)
    }
  }
  mat
}



#' Summarize 50-year average metrics from ssi_run model output
#'
#' Accepts an `ssi_run` object returned by `run_scenarios()`.
#'
#' For each scenario and iteration, computes the mean across the last
#' `years` rows (default 50) for `obsEsc`, `obsHarv`, and `obsRet`.
#' Then returns iteration-level summaries (means, SDs, quantiles)
#' grouped by scenario.
#'
#' @param data `ssi_run` object from `run_scenarios()`.
#' @param years Integer. Number of most recent years to summarize
#'   (default 50).
#'
#' @return A named list with elements `escapement`, `harvest`, and `return`.
#'   Each is a list (one element per scenario) containing:
#'   - `means`: numeric vector (length = niter) of iteration means
#'   - `sd`: numeric vector (length = niter) of iteration SDs
#'   - `quantiles`: matrix niter x 9 with iteration-wise quantiles
#'
#' @keywords internal

.summarize_50_year_avg <- function(data, years = 50L) {
  # Validate input
  if (!inherits(data, "ssi_run")) {
    stop("`data` must be an 'ssi_run' object from run_scenarios().")
  }

  obs_list <- data$results$obs
  if (is.null(obs_list) || !length(obs_list)) {
    stop("No observation lists found in `ssi_run` (data$results$obs).")
  }

  # Iterations
  nscen <- length(obs_list)
  niter <- length(obs_list[[1L]])
  iter_names <- paste0("iter_", seq_len(niter))

  # Years: determine common window
  first_obs <- obs_list[[1L]][[1L]]
  if (!is.data.frame(first_obs) || !nrow(first_obs)) {
    stop("Empty observation data in first scenario/iteration.")
  }
  tail_n <- min(nrow(first_obs), as.integer(years))
  year_index <- seq.int(nrow(first_obs) - tail_n + 1L, nrow(first_obs))

  # Helper: per-iteration stats across selected years
  compute_stats <- function(data_matrix) {
    means <- round(rowMeans(data_matrix, na.rm = TRUE), 0)
    names(means) <- iter_names

    sd_val <- round(apply(data_matrix, 1L, stats::sd, na.rm = TRUE), 0)
    names(sd_val) <- iter_names

    qs <- apply(
      data_matrix, 1L, stats::quantile,
      probs = c(0.025, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975),
      na.rm = TRUE
    )
    qs <- round(t(qs), 0)
    rownames(qs) <- iter_names
    colnames(qs) <- c("q2.5","q5","q10","q25","q50","q75","q90","q95","q97.5")

    list(means = means, sd = sd_val, quantiles = qs)
  }

  # Build summaries for each scenario
  escapement  <- vector("list", nscen)
  harvest     <- vector("list", nscen)
  return_list <- vector("list", nscen)

  for (j in seq_len(nscen)) {
    esc_mat  <- matrix(NA_real_, niter, length(year_index))
    harv_mat <- matrix(NA_real_, niter, length(year_index))
    ret_mat  <- matrix(NA_real_, niter, length(year_index))

    for (k in seq_len(niter)) {
      obs <- obs_list[[j]][[k]]
      if (!is.data.frame(obs) || !nrow(obs)) next

      esc_jk <- obs$obsEsc[year_index]
      harv_jk <- obs$obsHarv[year_index]
      ret_jk  <- obs$obsRet[year_index]

      if (!is.null(esc_jk))  esc_mat[k, ]  <- esc_jk
      if (!is.null(harv_jk)) harv_mat[k, ] <- harv_jk
      if (!is.null(ret_jk))  ret_mat[k, ]  <- ret_jk
    }

    escapement[[j]]   <- compute_stats(esc_mat)
    harvest[[j]]      <- compute_stats(harv_mat)
    return_list[[j]]  <- compute_stats(ret_mat)
  }

  names(escapement) <- names(harvest) <- names(return_list) <- paste0("scenario_", seq_len(nscen))

  list(
    escapement = escapement,
    harvest    = harvest,
    return     = return_list
  )
}



#' Summarize observations by year across scenarios and iterations
#'
#' Accepts an `ssi_run` object from `run_scenarios()`.
#'
#' For each scenario and each year position, computes summary statistics
#' (Mean, SD, and quantiles) across iterations for `obsEsc`, `obsHarv`,
#' and `obsRet`.
#'
#' @param data `ssi_run` object.
#' @param years Integer or `NULL`. The function uses the most recent `years`
#'   (default 50) rows from each iteration. If an iteration has fewer than
#'   `years` rows, the largest common window across all iterations is used.
#'
#' @return A nested list `MeanByYear` with components:
#'   - `escapement`, `harvest`, `return`: each is a list of scenarios;
#'     each scenario contains a list of length `ny` (years), where each element
#'     is a named numeric vector with stats:
#'     `Mean`, `SD`, `2.5%`, `5%`, `10%`, `25%`, `50%`, `75%`, `90%`, `95%`, `97.5%`.
#' @keywords internal

.summarize_by_year <- function(data, years = 50L) {
  if (!inherits(data, "ssi_run")) {
    stop("`data` must be an 'ssi_run' object.")
  }

  obs_list <- data$results$obs
  if (is.null(obs_list) || !length(obs_list)) {
    stop("No observation lists found in `ssi_run` (data$results$obs).")
  }

  nscen <- length(obs_list)
  niter <- length(obs_list[[1L]])

  if (is.data.frame(data$scenarios) && "scen_num" %in% names(data$scenarios)) {
    scen_names <- paste0("scenario_", data$scenarios$scen_num[seq_len(nscen)])
  } else {
    scen_names <- paste0("scenario_", seq_len(nscen))
  }
  iter_names <- paste0("iter_", seq_len(niter))

  lens <- unlist(lapply(obs_list, function(iters) {
    vapply(iters, function(df) if (is.data.frame(df)) nrow(df) else 0L, integer(1))
  }), use.names = FALSE)
  if (!length(lens) || all(lens == 0L)) stop("Observation data frames are empty.")

  years <- as.integer(years)
  if (is.na(years) || years <= 0L) years <- 50L
  nyh <- min(min(lens[lens > 0L]), years)
  year_names <- paste0("t", seq_len(nyh))

  stat_names <- c("Mean","SD","2.5%","5%","10%","25%","50%","75%","90%","95%","97.5%")
  summarize_vec <- function(x) {
    stats::setNames(round(c(
      mean(x, na.rm = TRUE),
      stats::sd(x, na.rm = TRUE),
      stats::quantile(x, probs = c(0.025,0.05,0.10,0.25,0.5,0.75,0.90,0.95,0.975), na.rm = TRUE)
    ), 0), stat_names)
  }

  make_nested_list <- function() {
    setNames(
      replicate(nscen,
                replicate(nyh, rep(NA_real_, length(stat_names)), simplify = FALSE),
                simplify = FALSE
      ),
      scen_names
    )
  }

  MeanByYear <- list(
    escapement = make_nested_list(),
    harvest    = make_nested_list(),
    return     = make_nested_list()
  )

  for (j in seq_len(nscen)) {
    sims <- obs_list[[j]]
    esc_mat  <- matrix(NA_real_, nrow = nyh, ncol = niter)
    harv_mat <- matrix(NA_real_, nrow = nyh, ncol = niter)
    ret_mat  <- matrix(NA_real_, nrow = nyh, ncol = niter)

    for (k in seq_len(niter)) {
      sim <- sims[[k]]
      if (!is.null(sim) && is.data.frame(sim) && nrow(sim) > 0L) {
        rows <- seq.int(nrow(sim) - nyh + 1L, nrow(sim))
        esc <- sim$obsEsc[rows]
        harv <- sim$obsHarv[rows]
        ret  <- sim$obsRet[rows]
        if (!is.null(esc))  esc_mat[,  k] <- esc
        if (!is.null(harv)) harv_mat[, k] <- harv
        if (!is.null(ret))  ret_mat[,  k] <- ret
      }
    }

    for (i in seq_len(nyh)) {
      MeanByYear$escapement[[j]][[i]] <- summarize_vec(esc_mat[i, ])
      MeanByYear$harvest[[j]][[i]]    <- summarize_vec(harv_mat[i, ])
      MeanByYear$return[[j]][[i]]     <- summarize_vec(ret_mat[i, ])
    }
  }

  for (var in names(MeanByYear)) {
    for (scen in scen_names) {
      names(MeanByYear[[var]][[scen]]) <- year_names
    }
  }

  return(MeanByYear)
}

