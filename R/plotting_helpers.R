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
#' @keywords internal

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


#' Standardize scenario labels and add IDs that match results list names
#'
#' Used by Ohlberger figures. Adds scen_key matching the names in results lists.
#'
#' @param data An `ssi_run` object
#' @return Tibble with standardized scenario labels and keys
#' @keywords internal
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


#' Create output directory for figures
#'
#' Used by Ohlberger figures.
#'
#' @param output_dir Base output directory
#' @param file_basename Subfolder name to create
#' @return Full path to created directory
#' @keywords internal
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
.get_firstrev <- function(data, default = 20L) as.integer(.get_param(data, "firstrev", default))

# Nested obs list: [[scenario]][[iteration]] -> data.frame with obsEsc/obsHarv/obsRet
.get_obs_list <- function(data) {
  obs_list <- data$results$obs
  if (is.null(obs_list) || !length(obs_list)) {
    stop("No observation lists found at data$results$obs.")
  }
  obs_list
}

# SR params list (alpha, beta, ...) per scenario/iteration (data.frame per iter)
.get_sr_params_list <- function(data) {
  sr <- data$results$sr_sim
  if (is.null(sr)) {
    stop("SR parameter list not found at data$results$sr_sim.")
  }
  sr
}

# S_MSY series per scenario/iteration (numeric vector over review years)
.get_smsy_series_list <- function(data) {
  sm <- data$results$S_msy
  if (is.null(sm)) {
    stop("S_MSY series not found at data$results$S_msy.")
  }
  sm
}


#' Access review years from an ssi_run
#'
#' Returns the sequence of review years used during the simulation.
#' This is saved into `ssi_run$parameters$review_years` by [run_scenarios()].
#'
#' @param data An `ssi_run` object.
#' @return Integer vector of review years.
#' @keywords internal
.get_review_years <- function(data) {
  ry <- tryCatch(data$parameters$review_years, error = function(e) NULL)
  if (is.null(ry)) {
    stop("`review_years` not found in `ssi_run$parameters`. ",
         "Re-run scenarios with a recent version of the package.")
  }
  ry
}


#' Index of the first post-historical review
#'
#' Uses `ssi_run$parameters$hist_end_year` and `ssi_run$parameters$review_years`
#' to locate the review at or just after the end of the historical period.
#'
#' @param data An `ssi_run` object.
#' @param review_years (optional) supply review years explicitly; by default
#'   uses `.get_review_years(data)`.
#' @return Integer scalar index into `review_years`.
#' @keywords internal
.get_review_index_post_history <- function(data, review_years = .get_review_years(data)) {
  hist_end <- tryCatch(data$parameters$hist_end_year, error = function(e) NULL)
  if (is.null(hist_end)) {
    stop("`hist_end_year` not found in `ssi_run$parameters`. ",
         "Re-run scenarios with a recent version of the package or upgrade legacy runs.")
  }

  idx <- which(review_years == hist_end)
  if (length(idx) == 0L) {
    idx <- suppressWarnings(min(which(review_years >= hist_end)))
  }

  if (!is.finite(idx) || idx < 1L) {
    stop("Could not locate first post-historical review in `review_years`.")
  }
  idx
}


#' Summarize 50-year average metrics from ssi_run model output
#'
#' Accepts an `ssi_run` object returned by `run_scenarios()`.
#'
#' For each scenario and iteration, computes the mean across the last 50 years
#' (determined by `nyh` parameter stored in the ssi_run object) for `obsEsc`,
#' `obsHarv`, and `obsRet`. Then returns iteration-level summaries (means, SDs,
#' quantiles) grouped by scenario.
#'
#' @param data `ssi_run` object from `run_scenarios()`.
#'
#' @return A named list with elements `escapement`, `harvest`, and `return`.
#'   Each is a list (one element per scenario) containing:
#'   - `means`: numeric vector (length = niter) of iteration means
#'   - `sd`: numeric vector (length = niter) of iteration SDs
#'   - `quantiles`: matrix niter x 9 with iteration-wise quantiles
#'   - `across_iter`: list with median and quantile bands across iterations
#'
#' @keywords internal

.summarize_50_year_avg <- function(data) {
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
  nyh      <- .get_nyh(data, 50L)
  ny_obs   <- nrow(as.data.frame(obs_list[[1]][[1]]))
  year_index <- (nyh + 1L):ny_obs  # fixed window used for *all* scen/iter

  # Helper: per-iteration stats across selected years
  compute_stats <- function(data_matrix) {
    means <- rowMeans(data_matrix, na.rm = TRUE)
    names(means) <- iter_names

    sd_val <- apply(data_matrix, 1L, stats::sd, na.rm = TRUE)
    names(sd_val) <- iter_names

    # (Not used by paper figs, but keep for completeness: within-iter, across time)
    qs <- t(apply(
      data_matrix, 1L, stats::quantile,
      probs = c(0.025, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975),
      na.rm = TRUE
    ))

    # Add across-iteration summaries of iteration means: median, 50% & 80% bands
    across_iter <- list(
      median = stats::median(means, na.rm = TRUE),
      q50_lo = stats::quantile(means, 0.25, na.rm = TRUE),
      q50_hi = stats::quantile(means, 0.75, na.rm = TRUE),
      q80_lo = stats::quantile(means, 0.10, na.rm = TRUE),
      q80_hi = stats::quantile(means, 0.90, na.rm = TRUE)
    )

    rownames(qs) <- iter_names
    colnames(qs) <- c("q2.5","q5","q10","q25","q50","q75","q90","q95","q97.5")

    list(means = means, across_iter = across_iter, sd = sd_val, quantiles = qs)
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
  if (!length(lens) || all(lens == 0L)) {
    stop("Observation data frames are empty.")
  }

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


#' Helper: Compute closure metrics (probability and cumulative)
#'
#' Calculates both probability of zero harvest AND cumulative closures by year.
#' This allows easy switching between metrics without recomputing.
#'
#' @param obs_list Simulation results from `run_scenarios()`.
#' @param scen_names Scenario labels.
#' @param iter_names Iteration labels.
#' @param year_rows_fun Function to extract last 50 years from a simulation data.frame.
#' @param scen_df Scenario metadata.
#' @return Tidy data.frame with `scen`, `year`, `prob_zero_harvest`,
#'   `mean_cumulative_closures`, `median_cumulative_closures`, etc.
#' @keywords internal
.compute_closure_metrics <- function(obs_list, scen_names, iter_names, year_rows_fun, scen_df) {
  dt_list <- lapply(seq_along(obs_list), function(i) {
    sims <- obs_list[[i]]
    data.table::rbindlist(lapply(seq_along(sims), function(j) {
      sim <- sims[[j]]
      if (!is.data.frame(sim) || is.null(sim$obsHarv)) return(NULL)
      rows <- year_rows_fun(sim)

      # Identify zero harvest years (closures)
      harvest_values <- sim$obsHarv[rows]
      zero_harvest_indicator <- harvest_values == 0

      # Calculate cumulative count of closures
      cumulative_closures <- cumsum(zero_harvest_indicator)

      data.table::data.table(
        scen   = scen_names[i],
        iter   = iter_names[j],
        year   = seq_len(length(rows)),
        obsHarv = harvest_values,
        cumulative_closures = cumulative_closures
      )
    }), fill = TRUE, use.names = TRUE)
  })

  dt_all <- data.table::rbindlist(dt_list, fill = TRUE, use.names = TRUE)
  dt_all$zeroHarvest <- ifelse(is.na(dt_all$obsHarv), NA, dt_all$obsHarv == 0)

  # Aggregate: compute BOTH probability AND cumulative metrics
  tibble::as_tibble(dt_all) |>
    dplyr::group_by(.data$scen, .data$year) |>
    dplyr::summarise(
      n_zero_harvest = sum(.data$zeroHarvest == TRUE, na.rm = TRUE),
      n_obs          = sum(!is.na(.data$zeroHarvest)),
      prob_zero_harvest = n_zero_harvest / n_obs,
      mean_cumulative_closures = mean(.data$cumulative_closures, na.rm = TRUE),
      median_cumulative_closures = median(.data$cumulative_closures, na.rm = TRUE),
      q25_cumulative_closures = stats::quantile(.data$cumulative_closures, 0.25, na.rm = TRUE),
      q75_cumulative_closures = stats::quantile(.data$cumulative_closures, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(scen_df, by = "scen")
}


#' Helper: Build tidy 50-yr summaries with quantiles for specified metrics
#'
#' Generalized function to extract point estimates and quantiles for any
#' metrics from `.summarize_50_year_avg()` output.
#'
#' @param summary_list Output from `.summarize_50_year_avg()`
#' @param scen_df Scenario metadata tibble with `scen` column
#' @param summary_fn "mean" or "median" for point estimates
#' @param metrics Character vector of metric names to extract (e.g., c("escapement", "harvest"))
#' @param prefix Optional prefix for column names (e.g., "esc_" produces "esc_point", "esc_q10", etc.)
#'   If NULL, uses metric name as prefix (e.g., "escapement_point")
#' @return Tidy dataframe (one row per scenario) with columns for each metric's
#'   point estimate and quantiles (point, q10, q25, q50, q75, q90)
#' @keywords internal
.make_50yr_metrics_tidy <- function(summary_list, scen_df, summary_fn = "mean",
                                    metrics = c("escapement", "harvest", "return"),
                                    prefix = NULL) {

  # Validate summary_fn
  if (!summary_fn %in% c("mean", "median")) {
    stop("`summary_fn` must be 'mean' or 'median'.")
  }

  # Internal extractor returning a named vector of point + quantiles
  .extract_point_and_intervals <- function(metric_list, summary_fn) {
    extract_one <- function(scn) {
      # central point
      point <- if (identical(summary_fn, "mean")) {
        base::mean(scn$means, na.rm = TRUE)
      } else {
        stats::median(scn$quantiles[, "q50"], na.rm = TRUE)
      }

      qs <- c("q10", "q25", "q50", "q75", "q90")
      # average each quantile across the 50-year window
      qvals <- vapply(qs, function(q) base::mean(scn$quantiles[, q], na.rm = TRUE), numeric(1))

      c(point = point, qvals)
    }

    out <- t(vapply(metric_list, extract_one, numeric(6)))
    colnames(out) <- c("point", "q10", "q25", "q50", "q75", "q90")
    out
  }

  # Start with scenario metadata
  df <- scen_df

  # Extract stats for each requested metric
  for (metric in metrics) {
    if (!metric %in% names(summary_list)) {
      stop("Metric '", metric, "' not found in summary_list.")
    }

    stats <- .extract_point_and_intervals(summary_list[[metric]], summary_fn)

    # Determine column prefix
    col_prefix <- if (!is.null(prefix) && length(prefix) == length(metrics)) {
      prefix[match(metric, metrics)]
    } else if (!is.null(prefix) && length(prefix) == 1) {
      prefix
    } else {
      paste0(metric, "_")
    }

    # Add columns with appropriate names
    df[[paste0(col_prefix, "point")]] <- round(stats[, "point"], 0)
    df[[paste0(col_prefix, "q10")]]   <- round(stats[, "q10"], 0)
    df[[paste0(col_prefix, "q25")]]   <- round(stats[, "q25"], 0)
    df[[paste0(col_prefix, "q50")]]   <- round(stats[, "q50"], 0)
    df[[paste0(col_prefix, "q75")]]   <- round(stats[, "q75"], 0)
    df[[paste0(col_prefix, "q90")]]   <- round(stats[, "q90"], 0)
  }

  df
}
