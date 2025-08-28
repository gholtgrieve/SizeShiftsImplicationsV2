#' Internal plotting helpers (not exported)
#' @keywords internal
#' @noRd

# Standardize scenario labels and add IDs that match results list names
.scenarios_from_ssirun <- function(data) {
  df <- standardize_scenario_labels(tibble::as_tibble(data$scenarios))

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
