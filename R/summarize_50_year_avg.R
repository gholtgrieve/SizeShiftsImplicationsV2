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
#' @export
summarize_50_year_avg <- function(data, years = 50L) {
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
