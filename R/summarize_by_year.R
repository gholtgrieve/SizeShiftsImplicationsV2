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
#' @export
summarize_by_year <- function(data, years = 50L) {
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
