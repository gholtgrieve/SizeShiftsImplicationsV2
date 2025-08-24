#' Summarize 50-year average metrics from model data
#'
#' Accepts either:
#' - an `ssi_run` object returned by `run_scenarios()`,
#' - a character path to an `.rds` saved by `run_scenarios()`, or
#' - a legacy list with elements: `obs_list`, `iter_names`, `year_index`.
#'
#' The function computes, for each scenario and iteration, the mean across the
#' last 50 years (or supplied indices) for `obsEsc`, `obsHarv`, and `obsRet`,
#' then returns iteration-level summaries (means, SDs, quantiles) grouped by
#' scenario, preserving the original return structure.
#'
#' @param data `ssi_run` object, file path to a saved run `.rds`, or legacy list.
#' @param years Integer. If using an `ssi_run`, number of most recent years to
#'   summarize (default 50). Ignored if legacy `year_index` is supplied.
#'
#' @return A named list with elements `escapement`, `harvest`, and `return`.
#'   Each is a list (one element per scenario) containing:
#'   - `means`: numeric vector (length = niter) of iteration means,
#'   - `sd`: numeric vector (length = niter) of iteration SDs,
#'   - `quantiles`: matrix niter x 9 with iteration-wise quantiles.
#' @export
summarize_50_year_avg <- function(data, years = 50L) {
  # --------------------------- normalize inputs ---------------------------
  # Allow passing a path to the saved run
  if (is.character(data) && length(data) == 1L && file.exists(data)) {
    data <- readRDS(data)
  }

  obs_list   <- NULL
  iter_names <- NULL
  year_index <- NULL

  if (inherits(data, "ssi_run")) {
    # New run object shape
    obs_list <- data$results$obs
    if (is.null(obs_list) || !length(obs_list)) {
      stop("No observation lists found in `ssi_run` (run$results$obs).")
    }
    # Derive iteration count and names from first scenario
    niter <- length(obs_list[[1L]])
    iter_names <- paste0("iter_", seq_len(niter))

    # Derive year_index = last `years` rows of obs for the first iteration
    first_obs <- obs_list[[1L]][[1L]]
    if (!is.data.frame(first_obs) || !nrow(first_obs)) {
      stop("Empty observation data in first scenario/iteration.")
    }
    tail_n <- min(nrow(first_obs), as.integer(years))
    year_index <- seq.int(nrow(first_obs) - tail_n + 1L, nrow(first_obs))
  } else if (is.list(data) && !is.null(data$obs_list)) {
    # Legacy bespoke structure
    obs_list   <- data$obs_list
    iter_names <- if (!is.null(data$iter_names)) data$iter_names else NULL
    year_index <- if (!is.null(data$year_index)) data$year_index else NULL
    if (is.null(year_index)) stop("Legacy input must include `year_index`.")
    if (is.null(iter_names)) {
      niter <- length(obs_list[[1L]])
      iter_names <- paste0("iter_", seq_len(niter))
    }
  } else {
    stop("`data` must be an 'ssi_run', a path to its .rds, or a legacy list with `obs_list`.")
  }

  # ------------------------------- validate -------------------------------
  nscen <- length(obs_list)
  if (nscen == 0L) stop("No scenarios found in `obs_list`.")
  niter <- length(obs_list[[1L]])
  if (niter == 0L) stop("No iterations found in first scenario of `obs_list`.")
  if (length(iter_names) != niter) {
    iter_names <- paste0("iter_", seq_len(niter))  # fallback
  }

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

  # ------------------------- build scenario summaries -------------------------
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

      # Column names expected from run_model()$obs (dataObs):
      # obsEsc, obsHarv, obsRet
      esc_jk <- obs$obsEsc[year_index]
      harv_jk <- obs$obsHarv[year_index]
      ret_jk <- obs$obsRet[year_index]

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
