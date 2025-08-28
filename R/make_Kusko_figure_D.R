#' Generate Figure D: Harvest vs Escapement with uncertainty (50-year summaries)
#'
#' Builds a 3×3 facet grid of scatterplots of Harvest (y, in 1000s) versus
#' Escapement (x, in 1000s). Each panel shows TRM (gray), YPR (blue), and DLM
#' (orange). Points are the averages across the last 50 simulated years and all
#' Monte Carlo iterations; error bars depict 50% (thick; q25–q75) and 80%
#' (thin; q10–q90) intervals for both axes. Facets: rows = trends
#' ("no trends", "ASL trends stabilized", "ASL trends continued"); columns =
#' factorMSY ("liberal", "MSY", "precautionary"). Labels and styling follow
#' Figure A helpers and theme.
#'
#' @param data `ssi_run` from `run_scenarios()`
#' @param output_dir Directory to save outputs (default `"."`)
#' @return (invisible) list with `data` and `plot`
#' @noRd
#' @keywords internal

.make_Kusko_figure_D <- function(data, output_dir = ".") {
  # Fixed settings
  statistic <- "mean"
  selectivity_filter <- "unselective"

  # Colors now map to factorMSY (liberal, MSY, precautionary)
  colors_fMSY <- c("liberal" = "darkgray",
                   "MSY" = "deepskyblue3",
                   "precautionary" = "orange")

  file_basename <- "FigureD"

  # 50-year summaries (same as Figure A)
  summary_list <- .summarize_50_year_avg(data)

  # Scenario metadata + standardized labels
  scen_df <- standardize_scenario_labels(tibble::as_tibble(data$scenarios))
  scen_df$scen <- paste0("scenario_", seq_len(nrow(scen_df)))

  # Build tidy dataframe with means + quantiles for Harvest and Escapement
  tidy_df <- .make_50yr_xy_tidy(summary_list, scen_df, statistic)

  # Select, order factors for facets and legend
  plot_df <-
    tidy_df |>
    dplyr::filter(.data$selectivity == !!selectivity_filter) |>
    dplyr::mutate(
      # rows (unchanged): trends
      trends = forcats::fct_relevel(
        .data$trends,
        "no trends", "ASL trends stabilized", "ASL trends continued"
      ),
      # columns (NEW): mgmt panels in this order
      mgmt = forcats::fct_relevel(.data$mgmt, "TRM", "YPR", "DLM"),
      # points (NEW): factorMSY mapped to color, in this order
      factorMSY = forcats::fct_relevel(.data$factorMSY, "liberal", "MSY", "precautionary")
    )

  # --- Make x & y share the same global limits (>= 0) ---
  max_val <- max(
    plot_df$esc_mean, plot_df$harv_mean,
    plot_df$esc_q90,  plot_df$harv_q90,
    na.rm = TRUE
  )
  lims_equal <- c(0, max_val)

  # Output folder
  output_dir <- file.path(output_dir, file_basename)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Plot
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$esc_mean,
      y = .data$harv_mean,
      color = .data$factorMSY
    )
  ) +
    # 80% intervals (thin): q10–q90
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$harv_q10, ymax = .data$harv_q90),
      linewidth = 0.4, alpha = 0.9
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$esc_q10, xmax = .data$esc_q90),
      height = 0, linewidth = 0.4, alpha = 0.9
    ) +
    # 50% intervals (thick): q25–q75
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$harv_q25, ymax = .data$harv_q75),
      linewidth = 0.9
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$esc_q25, xmax = .data$esc_q75),
      height = 0, linewidth = 0.9
    ) +
    # 1:1 diagonal line (drawn within each panel)
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", linewidth = 0.4, alpha = 0.7, color = "slategray") +
    # Points (central)
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colors_fMSY, guide = "none") +
    ggplot2::scale_y_continuous(
      limits = lims_equal, expand = c(0.07, 0.07),
      labels = function(y) y / 1000
    ) +
    ggplot2::scale_x_continuous(
      limits = lims_equal, expand = c(0.07, 0.07),
      labels = function(x) x / 1000
    ) +
    # Keep axes 1:1 visually too
    ggplot2::coord_equal() +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$trends),      # rows: trends (unchanged)
      cols = ggplot2::vars(.data$mgmt),        # cols: mgmt (TRM, YPR, DLM)
      scales = "fixed",
      switch = "y"
    ) +
    ggplot2::labs(
      x = "Escapement (1000s)",
      y = "Harvest (1000s)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10),
      strip.text.y = ggplot2::element_text(size = 10),
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
  ggplot2::ggsave(plot_path, p, width = 9, height = 9, units = "in")
  readr::write_csv(plot_df, data_path)

  message("Figure D saved to: ", plot_path)
  message("Data saved to: ", data_path)

  return(invisible(list(data = plot_df, plot = p)))
}


#' Helper: Build tidy XY (Harvest/Escapement) 50-yr summaries with quantiles
#'
#' Converts the output of `.summarize_50_year_avg()` into a tidy data frame that
#' includes central points (means or medians) and quantile ranges (q10, q25, q50,
#' q75, q90) for Harvest and Escapement, joined with standardized scenario
#' metadata.
#'
#' @param summary_list Output from `.summarize_50_year_avg()`.
#' @param scen_df Scenario metadata tibble (with `scen` column).
#' @param statistic Which statistic for the point estimate: "mean" or "median".
#' @return Tidy dataframe (one row per scenario) with columns:
#'   mgmt, selectivity, trends, factorMSY, scen, and
#'   esc_mean, esc_q10, esc_q25, esc_q50, esc_q75, esc_q90,
#'   harv_mean, harv_q10, harv_q25, harv_q50, harv_q75, harv_q90.
#' @keywords internal
.make_50yr_xy_tidy <- function(summary_list, scen_df, statistic = "mean") {

  # Internal extractor returning a named vector of point + quantiles
  .extract_point_and_intervals <- function(metric_list, statistic) {
    # For each scenario, compute the point (mean or median of q50) and
    # quantile summaries averaged across the 50-yr window (over iterations).
    extract_one <- function(scn) {
      # central point
      point <- if (identical(statistic, "mean")) {
        base::mean(scn$means, na.rm = TRUE)
      } else if (identical(statistic, "median")) {
        stats::median(scn$quantiles[, "q50"], na.rm = TRUE)
      } else {
        stop("`statistic` must be 'mean' or 'median'.")
      }

      qs <- c("q10", "q25", "q50", "q75", "q90")
      # average each quantile across the 50-year window
      qvals <- vapply(qs, function(q) base::mean(scn$quantiles[, q], na.rm = TRUE), numeric(1))

      c(point = point, qvals)
    }

    out <- t(vapply(metric_list, extract_one, numeric(6)))
    colnames(out) <- c("mean", "q10", "q25", "q50", "q75", "q90")
    out
  }

  esc_stats <- .extract_point_and_intervals(summary_list$escapement, statistic)
  har_stats <- .extract_point_and_intervals(summary_list$harvest,    statistic)

  df <- scen_df
  # Bind stats with informative names
  df$esc_mean <- round(esc_stats[, "mean"], 0)
  df$esc_q10  <- round(esc_stats[, "q10"],  0)
  df$esc_q25  <- round(esc_stats[, "q25"],  0)
  df$esc_q50  <- round(esc_stats[, "q50"],  0)
  df$esc_q75  <- round(esc_stats[, "q75"],  0)
  df$esc_q90  <- round(esc_stats[, "q90"],  0)

  df$harv_mean <- round(har_stats[, "mean"], 0)
  df$harv_q10  <- round(har_stats[, "q10"],  0)
  df$harv_q25  <- round(har_stats[, "q25"],  0)
  df$harv_q50  <- round(har_stats[, "q50"],  0)
  df$harv_q75  <- round(har_stats[, "q75"],  0)
  df$harv_q90  <- round(har_stats[, "q90"],  0)

  df
}
