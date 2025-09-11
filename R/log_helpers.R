# ---- log_helpers.R -----------------------------------------------------------
# Minimal-deps logging helpers for run_model() review sub-blocks 7a–7h
# Each fun_log_7*() appends exactly one row to a dedicated CSV per review.

# ============== tiny utilities ===============================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

fun_log_dir <- function(config = NULL) {
  # Priority: config$log_dir -> option -> ./logs
  d <- config$log_dir %||% getOption("ssi.log_dir")
  if (is.null(d) || is.na(d)) d <- file.path(getwd(), "logs")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

fun_log_paths <- function(config = NULL) {
  d <- fun_log_dir(config)
  list(
    a = file.path(d, "log_7a_yearloop.csv"),
    b = file.path(d, "log_7b_reconstruct_window.csv"),
    c = file.path(d, "log_7c_sim_gls.csv"),
    d = file.path(d, "log_7d_reconstruct_obs.csv"),
    e = file.path(d, "log_7e_obs_gls.csv"),
    f = file.path(d, "log_7f_dlm.csv"),
    g = file.path(d, "log_7g_ypr.csv"),
    h = file.path(d, "log_7h_goal_update.csv")
  )
}

fun_reason <- function(...) {
  x <- Filter(function(z) is.character(z) && nzchar(z), unlist(list(...), use.names = FALSE))
  if (!length(x)) NA_character_ else paste(unique(x), collapse = ";")
}


# ======== parallel-safe writer: shard per worker (PID) ========================

fun_shard_path <- function(base_csv_path) {
  d    <- dirname(base_csv_path)
  stem <- tools::file_path_sans_ext(basename(base_csv_path))  # "log_7a_yearloop"
  file.path(d, sprintf("%s.pid%s.csv", stem, Sys.getpid()))
}

fun_write_row <- function(path, row_list) {
  df <- as.data.frame(row_list, check.names = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  shard <- fun_shard_path(path)
  has_file <- file.exists(shard)
  utils::write.table(
    df, file = shard, sep = ",",
    row.names = FALSE, col.names = !has_file,
    append   = has_file, qmethod = "double"
  )
}

# ======== merge shards -> final CSV (run once after all workers finish) =======

fun_merge_one_log <- function(base_csv_path, remove_shards = TRUE) {
  d    <- dirname(base_csv_path)
  stem <- tools::file_path_sans_ext(basename(base_csv_path))
  pat  <- paste0("^", stem, "\\.pid[0-9]+\\.csv$")
  shards <- list.files(d, pattern = pat, full.names = TRUE)
  if (!length(shards)) return(invisible(NULL))

  # Read all shards, harmonize columns (union), then rbind
  dfs <- lapply(shards, function(f) utils::read.csv(f, check.names = FALSE))
  all_cols <- Reduce(union, lapply(dfs, names))
  dfs2 <- lapply(dfs, function(x) {
    missing <- setdiff(all_cols, names(x))
    if (length(missing)) x[missing] <- NA
    # reorder columns to union order
    x[all_cols]
  })

  out <- do.call(rbind, dfs2)
  utils::write.table(out, file = base_csv_path, sep = ",",
                     row.names = FALSE, col.names = TRUE, qmethod = "double")

  if (isTRUE(remove_shards)) unlink(shards)
  invisible(base_csv_path)
}

fun_merge_all_logs <- function(log_dir = getOption("ssi.log_dir", file.path(getwd(), "logs")),
                               remove_shards = TRUE) {
  stems <- c("log_7a_yearloop",
             "log_7b_reconstruct_window",
             "log_7c_sim_gls",
             "log_7d_reconstruct_obs",
             "log_7e_obs_gls",
             "log_7f_dlm",
             "log_7g_ypr",
             "log_7h_goal_update")
  for (s in stems) {
    fun_merge_one_log(file.path(log_dir, paste0(s, ".csv")), remove_shards = remove_shards)
  }
  invisible(log_dir)
}

fun_id_cols <- function(config, scen_num, iteration, review_i, y_rev, seed) {
  # Try to use a readable scenario_id if present
  scenario_id <- config$scenario_id %||% config$scen_id %||%
    config$scenario %||% (if (!is.null(scen_num)) paste0("scen_", scen_num) else NA_character_)
  list(
    timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    scenario_id  = scenario_id,
    scen_num     = scen_num %||% NA_integer_,
    iteration    = iteration %||% NA_integer_,
    review_i     = review_i %||% NA_integer_,
    y_rev        = y_rev %||% NA_integer_,
    seed         = seed %||% NA_integer_
  )
}

fun_rng_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(list(min = NA_real_, max = NA_real_))
  list(min = min(x), max = max(x))
}

# Binary "is this YPR result usable?" test
# Default rule: all finite, positive H+S, interior optimum (not near bounds).
fun_ypr_converged <- function(alpha_rep_out, beta_rep_out,
                              zPR0, R0_val,
                              F_eq, fit_value, H_eq, S_eq, U_eq,
                              fit_convergence,
                              lower = -10, upper = 2,
                              tol = getOption("ssi.ypr_boundary_tol", 1e-4),
                              allow_boundary = getOption("ssi.ypr_allow_boundary", FALSE)) {
  finite_fit <- all(is.finite(c(alpha_rep_out, beta_rep_out, zPR0, R0_val,
                                F_eq, fit_value, H_eq, S_eq, U_eq)))
  if (!finite_fit) return(FALSE)

  den_HS <- H_eq + S_eq
  if (!is.finite(den_HS) || den_HS <= 0) return(FALSE)

  near_lower <- F_eq <= (lower + tol)
  near_upper <- F_eq >= (upper - tol)
  at_boundary <- near_lower || near_upper

  # Brent returns 0 even on boundary; keep the check for completeness
  ok_alg <- isTRUE(fit_convergence == 0)

  if (!ok_alg) return(FALSE)
  if (at_boundary && !allow_boundary) return(FALSE)

  TRUE
}


# ============== 7a: year loop aggregate ======================================

fun_log_7a_yearloop <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    PopDat, HarvRates, selectivities_by_age, yindex,
    sel_fallback_count = 0L, sample_fallback_count = 0L
) {
  pths <- fun_log_paths(config)
  n_years <- length(yindex)

  # Regime shifts across the window
  reg_y <- PopDat$reg[yindex]
  reg_y <- reg_y[is.finite(reg_y)]
  n_regime_shifts <- if (length(reg_y) > 1) sum(diff(reg_y) != 0) else 0

  # Ranges
  rEsc  <- fun_rng_range(PopDat$Esc[yindex])
  rHarv <- fun_rng_range(PopDat$Harv[yindex])
  rRet  <- fun_rng_range(PopDat$Ret[yindex])
  rU    <- fun_rng_range(HarvRates[yindex])

  # Selectivity sanity per year
  if (is.null(dim(selectivities_by_age))) {
    n_sel_zero_max <- NA_integer_
  } else {
    sel_sub <- selectivities_by_age[yindex, , drop = FALSE]
    row_max <- suppressWarnings(apply(sel_sub, 1L, function(v) {
      m <- suppressWarnings(max(v, na.rm = TRUE))
      if (!is.finite(m)) NA_real_ else m
    }))
    n_sel_zero_max <- sum(!is.finite(row_max) | row_max <= 0, na.rm = TRUE)
  }

  # SR non-finite counts (Rec as proxy)
  n_sr_nonfinite <- sum(!is.finite(PopDat$Rec[yindex]), na.rm = TRUE)

  row <- c(
    fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
    list(
      n_years = n_years,
      n_regime_shifts = n_regime_shifts,
      harv_rate_min = rU$min, harv_rate_max = rU$max,
      Esc_min = rEsc$min,     Esc_max = rEsc$max,
      Harv_min = rHarv$min,   Harv_max = rHarv$max,
      Ret_min  = rRet$min,    Ret_max  = rRet$max,
      n_years_selectivity_zeroMax = n_sel_zero_max,
      n_years_sr_nonfinite       = n_sr_nonfinite,
      sel_fallback_count = as.integer(sel_fallback_count),
      sample_fallback_count = as.integer(sample_fallback_count)
    )
  )
  fun_write_row(pths$a, row)
}

# ============== 7b: reconstructed window info =================================

fun_log_7b_reconstruct_window <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    nyrec, year_index_start, year_index_end
) {
  pths <- fun_log_paths(config)
  row <- c(
    fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
    list(
      nyrec = nyrec,
      year_index_start = year_index_start,
      year_index_end   = year_index_end
    )
  )
  fun_write_row(pths$b, row)
}

# ============== 7c: GLS on simulated data =====================================

fun_log_7c_sim_gls <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    sig_sim, phi_sim, log_a_sim, b_sim, alpha_sim, beta_sim,
    S_msy_sim, U_msy_sim,
    sim_gls_guard = 0L, sim_gls_reason = NA_character_
  ) {
    pths <- fun_log_paths(config)
    finite_fit <- all(is.finite(c(sig_sim, phi_sim, log_a_sim, b_sim, alpha_sim, beta_sim,
                                  S_msy_sim, U_msy_sim)))
    row <- c(
      fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
      list(
        sig_sim = sig_sim, phi_sim = phi_sim, log_a_sim = log_a_sim, b_sim = b_sim,
        alpha_sim = alpha_sim, beta_sim = beta_sim, S_msy_sim = S_msy_sim, U_msy_sim = U_msy_sim,
        finite_sim_fit = as.integer(all(is.finite(c(sig_sim, phi_sim, log_a_sim, b_sim, alpha_sim, beta_sim, S_msy_sim, U_msy_sim)))),
        sim_gls_guard = as.integer(sim_gls_guard),
        sim_gls_reason = sim_gls_reason
      )
    )
    fun_write_row(pths$c, row)
}

# ============== 7d: reconstruction from observations ===========================

fun_log_7d_reconstruct_obs <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    dataObs
) {
  pths <- fun_log_paths(config)
  n_rows <- tryCatch(nrow(dataObs), error = function(e) NA_integer_)
  min_obsEsc <- suppressWarnings(min(dataObs$obsEsc[is.finite(dataObs$obsEsc)], na.rm = TRUE))
  min_recRec <- suppressWarnings(min(dataObs$recRec[is.finite(dataObs$recRec)], na.rm = TRUE))
  n_bad_ln   <- sum(!is.finite(dataObs$obslnRS), na.rm = TRUE)
  too_few    <- is.finite(n_rows) && n_rows < 3L

  row <- c(
    fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
    list(
      n_rows_obs_kept = n_rows,
      min_obsEsc = if (is.infinite(min_obsEsc)) NA_real_ else min_obsEsc,
      min_recRec = if (is.infinite(min_recRec)) NA_real_ else min_recRec,
      n_nonfinite_obslnRS = n_bad_ln,
      too_few_rows_flag = as.integer(too_few)
    )
  )
  fun_write_row(pths$d, row)
}

# ============== 7e: GLS on reconstructed data (MSY goals) ======================

fun_log_7e_obs_gls <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    sig_obs, phi_obs, log_a_obs, b_obs, alpha_obs, beta_obs,
    S_msy_obs, U_msy_obs,
    obs_gls_guard = 0L, obs_gls_reason = NA_character_
) {
  pths <- fun_log_paths(config)
  finite_fit <- all(is.finite(c(sig_obs, phi_obs, log_a_obs, b_obs, alpha_obs, beta_obs,
                                S_msy_obs, U_msy_obs)))
  row <- c(
    fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
    list(
      sig_obs = sig_obs, phi_obs = phi_obs, log_a_obs = log_a_obs, b_obs = b_obs,
      alpha_obs = alpha_obs, beta_obs = beta_obs, S_msy_obs = S_msy_obs, U_msy_obs = U_msy_obs,
      finite_obs_fit = as.integer(all(is.finite(c(sig_obs, phi_obs, log_a_obs, b_obs, alpha_obs, beta_obs, S_msy_obs, U_msy_obs)))),
      obs_gls_guard = as.integer(obs_gls_guard),
      obs_gls_reason = obs_gls_reason
    )
  )
  fun_write_row(pths$e, row)
}

# ============== 7f: DLM branch =================================================

fun_log_7f_dlm <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    var_mode, nrd, sig_dlm, alpha_dlm_recent, beta_dlm_recent,
    Smsy_dlm, Umsy_dlm,
    dlm_guard = 0L, dlm_guard_reason = NA_character_
) {
  pths <- fun_log_paths(config)
  finite_fit <- all(is.finite(c(sig_dlm, alpha_dlm_recent, beta_dlm_recent, Smsy_dlm, Umsy_dlm)))
  row <- c(
    fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
    list(
      var_mode = as.character(var_mode),
      n_rows_dlm = nrd %||% NA_integer_,
      sig_dlm = sig_dlm,
      alpha_dlm_recent = alpha_dlm_recent,
      beta_dlm_recent  = beta_dlm_recent,
      Smsy_dlm = Smsy_dlm, Umsy_dlm = Umsy_dlm,
      finite_dlm_fit = as.integer(all(is.finite(c(sig_dlm, alpha_dlm_recent, beta_dlm_recent, Smsy_dlm, Umsy_dlm)))),
      dlm_guard = as.integer(dlm_guard),
      dlm_guard_reason = dlm_guard_reason
    )
  )
  fun_write_row(pths$f, row)
}

# ============== 7g: YPR branch (post-optim) ===================================

fun_log_7g_ypr <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    sig_rep_out, phi_rep_out, log_a_rep_out, alpha_rep_out, beta_rep_out,
    zPR0, R0_val, F_eq, fit_value, fit_convergence, H_eq, S_eq, U_eq,
    lower = -10, upper = 2,
    tol = getOption("ssi.ypr_boundary_tol", 1e-4),
    convergence = NULL,
    ypr_guard = NULL, ypr_guard_reason = NA_character_
) {
  pths <- fun_log_paths(config)
  near_lower <- is.finite(F_eq) && (F_eq <= (lower + tol))
  near_upper <- is.finite(F_eq) && (F_eq >= (upper - tol))
  boundary_side <- if (near_lower) "lower" else if (near_upper) "upper" else "interior"

  if (is.null(convergence)) {
    convergence <- fun_ypr_converged(alpha_rep_out, beta_rep_out,
                                     zPR0, R0_val, F_eq, fit_value, H_eq, S_eq, U_eq,
                                     fit_convergence, lower, upper, tol)
  }
  if (is.null(ypr_guard)) ypr_guard <- as.integer(!isTRUE(convergence))

  finite_fit <- all(is.finite(c(alpha_rep_out, beta_rep_out, zPR0, R0_val, F_eq, fit_value, H_eq, S_eq, U_eq)))

  row <- c(
    fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
    list(
      sig_rep_out = sig_rep_out, phi_rep_out = phi_rep_out, log_a_rep_out = log_a_rep_out,
      alpha_rep_out = alpha_rep_out, beta_rep_out = beta_rep_out,
      zPR0 = zPR0, R0 = R0_val,
      F_eq = F_eq, fit_value = fit_value, fit_convergence = fit_convergence,
      H_eq = H_eq, S_eq = S_eq, U_eq = U_eq,
      boundary_side = boundary_side,
      finite_ypr = as.integer(finite_fit),
      convergence = as.integer(isTRUE(convergence)),
      ypr_guard = as.integer(ypr_guard),
      ypr_guard_reason = ypr_guard_reason
    )
  )
  fun_write_row(pths$g, row)
}


# ============== 7h: goal update ===============================================

fun_log_7h_goal_update <- function(
    config, scen_num, iteration, review_i, y_rev, seed,
    harvmgmt, msygoal_before, msygoal_after, factorMSY,
    S_msy_obs = NA_real_, U_msy_obs = NA_real_,
    S_eq = NA_real_, U_eq = NA_real_,
    Smsy_dlm = NA_real_, Umsy_dlm = NA_real_
) {
  pths <- fun_log_paths(config)
  proposal_type <- switch(as.character(harvmgmt),
                          "smsy_goal"="S_msy_obs", "umsy_goal"="U_msy_obs",
                          "s_eq_goal"="S_eq", "u_eq_goal"="U_eq",
                          "smsy_dlm_goal"="Smsy_dlm", "umsy_dlm_goal"="Umsy_dlm",
                          "fixed_rate")
  proposal_value <- switch(proposal_type,
                           "S_msy_obs"=S_msy_obs, "U_msy_obs"=U_msy_obs, "S_eq"=S_eq, "U_eq"=U_eq,
                           "Smsy_dlm"=Smsy_dlm, "Umsy_dlm"=Umsy_dlm, NA_real_)
  accepted <- is.finite(msygoal_before) && is.finite(msygoal_after) && !isTRUE(all.equal(msygoal_before, msygoal_after))

  # Try to infer why it wasn’t accepted
  reason <- NA_character_
  if (!accepted) {
    if (!is.finite(proposal_value)) {
      reason <- "nonfinite_proposal"
    } else if (proposal_type %in% c("S_msy_obs","S_eq","Smsy_dlm")) {
      if (!(proposal_value > msygoal_before*0.5 && proposal_value < msygoal_before*2)) reason <- "out_of_range"
    } else if (proposal_type %in% c("U_msy_obs","U_eq","Umsy_dlm")) {
      if (!(proposal_value > 0 && proposal_value < 0.85)) reason <- "out_of_range"
    }
    if (is.na(reason) && review_i == 1L) reason <- "fallback_to_initial"
  }

  row <- c(
    fun_id_cols(config, scen_num, iteration, review_i, y_rev, seed),
    list(
      harvmgmt = as.character(harvmgmt),
      proposal_type  = proposal_type,
      proposal_value = proposal_value,
      msygoal_before = msygoal_before,
      msygoal_after  = msygoal_after,
      factorMSY      = factorMSY,
      accepted       = as.integer(accepted),
      goal_skip_reason = reason
    )
  )
  fun_write_row(pths$h, row)
}


# ---- end log_helpers.R -------------------------------------------------------
