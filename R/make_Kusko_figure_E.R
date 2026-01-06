
#' Generate Figure E: P(Zero harvest) vs Harvest with uncertainty (50-year summaries)
#'
#' Same layout and styling as Figure D. Builds a facet grid of scatterplots with
#' x = Harvest (in 1000s) and y = Probability of zero harvest (50-year average).
#' Points represent factorMSY (liberal, MSY, precautionary) and error bars show
#' 80% (thin; q10–q90) and 50% (thick; q25–q75) intervals across Monte Carlo
#' iterations of the 50-year means. Facets: rows = trends; columns = mgmt
#' (TRM, YPR, DLM). Labels and theme follow Figure A/D.
#'
#' @param data `ssi_run` from `run_scenarios()`
#' @param output_dir Directory to save outputs (default `"."`)
#' @return (invisible) list with `data` and `plot`
#' @noRd
#' @keywords internal
.make_Kusko_figure_E <- function(data, output_dir = ".") {
  statistic <- "mean"
  selectivity_filter <- "unselective"
  file_basename <- "FigureE"

  # Colors by mgmt (TRM, YPR, DLM)
  colors_mgmt <- c("TRM" = "darkgray",
                   "YPR" = "deepskyblue3",
                   "DLM" = "orange")

  # 50-year summaries for harvest (and we’ll compute P0 via helper)
  summary_list <- .summarize_50_year_avg(data)

  # Scenario metadata
  scen_df <- .standardize_scenario_labels(tibble::as_tibble(data$scenarios))
  scen_df$scen <- paste0("scenario_", seq_len(nrow(scen_df)))

  # Build tidy DF (Harvest + ProbZero) and apply filters / factor orders
  plot_df <- .make_50yr_hp_tidy(data, summary_list, scen_df, statistic) |>
    dplyr::filter(
      .data$selectivity == !!selectivity_filter,
      .data$factorMSY == "precautionary"                 # (1) keep only precautionary
    ) |>
    dplyr::mutate(
      trends = forcats::fct_relevel(                      # (3) facet columns order
        .data$trends,
        "no trends", "ASL trends stabilized", "ASL trends continued"
      ),
      mgmt = forcats::fct_relevel(.data$mgmt, "TRM", "YPR", "DLM")
    )

  # Output folder
  output_dir <- file.path(output_dir, file_basename)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Plot
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$harv_mean,               # x = Harvest (1000s)
      y = .data$.data$p0_mean,            # y = Average probability of zero harvest in any given year
      color = .data$mgmt                 # (4) color by mgmt
    )
  ) +
    # 80% intervals (thin)
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$p0_q10, ymax = .data$p0_q90),
      linewidth = 0.4, alpha = 0.9, width = 0
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$harv_q10, xmax = .data$harv_q90),
      width = 0, linewidth = 0.4, alpha = 0.9
    ) +
    # 50% intervals (thick)
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$p0_q25, ymax = .data$p0_q75),
      linewidth = 0.9, width = 0
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$harv_q25, xmax = .data$harv_q75),
      width = 0, linewidth = 0.9
    ) +
    # Points
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colors_mgmt, name = "Method") +  # legend = mgmt
    ggplot2::scale_x_continuous(
      expand = c(0.07, 0.07),
      limits = c(0, NA),
      labels = function(x) x / 1000
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0.0, 0.07),
      limits = c(0, NA),
      labels = scales::percent_format(accuracy = 1)     # show probability as %
    ) +
    # (2) remove factorMSY from facets; (3) columns = trends, single row
    ggplot2::facet_grid(
      rows = NULL,
      cols = ggplot2::vars(.data$trends),
      scales = "fixed",
      switch = "y"
    ) +
    ggplot2::labs(
      x = "Harvest (1000s)",
      y = "P(Zero harvest, 50-yr median) (%)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10),
      axis.line = ggplot2::element_line(linewidth = 0.1),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      strip.placement = "outside"
    )

  # Save
  plot_path <- file.path(output_dir, paste0(file_basename, ".pdf"))
  data_path <- file.path(output_dir, paste0(file_basename, ".csv"))
  ggplot2::ggsave(plot_path, p, width = 9, height = 4, units = "in")  # shorter since 1 row
  readr::write_csv(plot_df, data_path)

  message("Figure E saved to: ", plot_path)
  message("Data saved to: ", data_path)

  invisible(list(data = plot_df, plot = p))
}


#' Helper: Build tidy HP (Harvest / ProbZeroHarvest) 50-yr summaries with quantiles
#'
#' Harvest stats are read from `summary_list$harvest` (output of
#' `.summarize_50_year_avg()`). Probability of zero harvest (P0) is computed
#' directly from `data$results$obs` by:
#' \enumerate{
#'   \item Identifying the common last-50-year window across iterations (same rule
#'         as `.summarize_by_year()`): the minimum available length across all
#'         iterations and scenarios, capped at 50.
#'   \item For each scenario and iteration, computing the 50-yr mean of the
#'         zero-harvest indicator: `as.numeric(obsHarv <= 0)`.
#'   \item Across iterations, summarizing to mean and quantiles (10, 25, 50, 75, 90\%).
#' }
#'
#' @param data `ssi_run` object (from `run_scenarios()`).
#' @param summary_list Output from `.summarize_50_year_avg(data)`.
#' @param scen_df Scenario metadata tibble with standardized labels and column `scen`
#'        matching the order of scenarios (e.g., "scenario_1", ...).
#' @param statistic Point estimate for harvest axis: `"mean"` or `"median"`.
#'
#' @return Tidy data.frame (one row per scenario) with columns:
#'   `mgmt, selectivity, trends, factorMSY, scen,
#'    harv_mean, harv_q10, harv_q25, harv_q50, harv_q75, harv_q90,
#'    p0_mean,  p0_q10,  p0_q25,  p0_q50,  p0_q75,  p0_q90`
#' @keywords internal
.make_50yr_hp_tidy <- function(data, summary_list, scen_df, statistic = "mean") {
  if (!inherits(data, "ssi_run")) {
    stop("`data` must be an 'ssi_run' object.")
  }
  obs_list <- data$results$obs
  if (is.null(obs_list) || !length(obs_list)) {
    stop("No observation lists found in `ssi_run` (data$results$obs).")
  }

  # ----- harvest stats from 50-yr summaries (matches Figures A/D) ----------------
  .extract_point_and_intervals <- function(metric_list, statistic) {
    extract_one <- function(scn) {
      point <- if (identical(statistic, "mean")) {
        base::mean(scn$means, na.rm = TRUE)
      } else if (identical(statistic, "median")) {
        stats::median(scn$quantiles[, "q50"], na.rm = TRUE)
      } else {
        stop("`statistic` must be 'mean' or 'median'.")
      }
      qs <- c("q10","q25","q50","q75","q90")
      qvals <- vapply(qs, function(q) base::mean(scn$quantiles[, q], na.rm = TRUE), numeric(1))
      c(point = point, qvals)
    }
    out <- t(vapply(metric_list, extract_one, numeric(6)))
    colnames(out) <- c("mean","q10","q25","q50","q75","q90")
    out
  }
  har_stats <- .extract_point_and_intervals(summary_list$harvest, statistic)

  # ----- establish common 50-yr window across iterations (like .summarize_by_year) -----
  nscen <- length(obs_list)
  niter <- length(obs_list[[1L]])
  # lengths of each iteration DF across all scenarios
  lens <- unlist(lapply(obs_list, function(iters) {
    vapply(iters, function(df) if (is.data.frame(df)) nrow(df) else 0L, integer(1))
  }), use.names = FALSE)
  if (!length(lens) || all(lens == 0L)) stop("Observation data frames are empty.")
  nyh <- min(min(lens[lens > 0L]), 50L)  # common window length, capped at 50
  if (!is.finite(nyh) || nyh <= 0L) stop("Could not determine a positive common window length.")

  scen_names <- scen_df$scen
  if (length(scen_names) != nscen) {
    # Fallback to scenario_1.. if sizes mismatch
    scen_names <- paste0("scenario_", seq_len(nscen))
  }
  iter_names <- paste0("iter_", seq_len(niter))

  # helper to return row indices for the last nyh rows of a simulation DF
  year_rows_fun <- function(sim_df) {
    n <- nrow(sim_df)
    if (is.null(n) || n <= 0L) integer(0)
    else seq.int(n - nyh + 1L, n)
  }

  # ----- build per-iteration, per-scenario 0/1 zero-harvest and reduce to 50-yr mean -----
  # produce a data.frame with columns: scen, iter, p0_iter_mean
  per_iter_list <- vector("list", nscen)
  for (i in seq_len(nscen)) {
    sims <- obs_list[[i]]
    rows_list <- lapply(seq_len(niter), function(j) {
      sim <- sims[[j]]
      if (!is.data.frame(sim) || is.null(sim$obsHarv)) return(NULL)
      rows <- year_rows_fun(sim)
      if (length(rows) != nyh) return(NULL)
      # 50-yr mean probability of zero harvest for this iteration
      p0_iter_mean <- mean(as.numeric(sim$obsHarv[rows] <= 0), na.rm = TRUE)
      data.frame(scen = scen_names[i], iter = iter_names[j], p0_iter_mean = p0_iter_mean)
    })
    per_iter_list[[i]] <- dplyr::bind_rows(rows_list)
  }
  per_iter <- dplyr::bind_rows(per_iter_list)
  if (!nrow(per_iter)) {
    stop("Failed to compute per-iteration zero-harvest means; check that obsHarv exists and has rows.")
  }

  # ----- summarise across iterations to mean & quantiles for each scenario -------
  p0_stats_df <- per_iter |>
    dplyr::group_by(.data$scen) |>
    dplyr::summarise(
      mean = mean(.data$p0_iter_mean, na.rm = TRUE),
      q10  = as.numeric(stats::quantile(.data$p0_iter_mean, 0.10, na.rm = TRUE)),
      q25  = as.numeric(stats::quantile(.data$p0_iter_mean, 0.25, na.rm = TRUE)),
      q50  = as.numeric(stats::quantile(.data$p0_iter_mean, 0.50, na.rm = TRUE)),
      q75  = as.numeric(stats::quantile(.data$p0_iter_mean, 0.75, na.rm = TRUE)),
      q90  = as.numeric(stats::quantile(.data$p0_iter_mean, 0.90, na.rm = TRUE)),
      .groups = "drop"
    )

  # align to scenario order of scen_df and coerce to matrix-like access
  p0_stats <- as.matrix(
    p0_stats_df[match(scen_df$scen, p0_stats_df$scen),
                c("mean","q10","q25","q50","q75","q90")]
  )

  # ----- assemble tidy dataframe -------------------------------------------------
  df <- scen_df

  # Harvest (x-axis) — integer counts
  df$harv_mean <- round(har_stats[, "mean"], 0)
  df$harv_q10  <- round(har_stats[, "q10"],  0)
  df$harv_q25  <- round(har_stats[, "q25"],  0)
  df$harv_q50  <- round(har_stats[, "q50"],  0)
  df$harv_q75  <- round(har_stats[, "q75"],  0)
  df$harv_q90  <- round(har_stats[, "q90"],  0)

  # Probability (y-axis) — keep on 0–1 scale (plot formats as %)
  df$p0_mean <- as.numeric(p0_stats[, "mean"])
  df$p0_q10  <- as.numeric(p0_stats[, "q10"])
  df$p0_q25  <- as.numeric(p0_stats[, "q25"])
  df$p0_q50  <- as.numeric(p0_stats[, "q50"])
  df$p0_q75  <- as.numeric(p0_stats[, "q75"])
  df$p0_q90  <- as.numeric(p0_stats[, "q90"])

  df
}


