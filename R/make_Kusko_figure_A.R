#' Generate Figure A: Summary plot of 50-year metrics
#'
#' Builds a facet grid of median/mean 50-year metrics by trends (x), colored by
#' management method, with columns for factorMSY. Returns a list containing plot
#' object and tidy data used in plotting invisibly.
#'
#' @param data `ssi_run` object from `run_scenarios()`
#' @param output_dir Directory to save plots and data (default: current working directory).
#' @param summary_fn Character: "mean" (default) or "median" for point estimates
#'   across Monte Carlo iterations.
#'
#' @return A list with:
#'   - `data`: tidy data.frame used for plotting
#'   - `plot`: ggplot object
#' @noRd
#' @keywords internal

.make_Kusko_figure_A <- function(data, output_dir = ".", summary_fn = "mean") {
  # Validate summary_fn parameter
  summary_fn <- match.arg(summary_fn, choices = c("mean", "median"))

  # Fixed settings
  selectivity_filter <- "unselective"
  colors <- c("darkgray", "deepskyblue3", "orange")
  file_basename <- "FigureA"

  # Summarize 50-year averages
  summary_list <- .summarize_50_year_avg(data)

  # Scenario metadata
  scen_df <- .standardize_scenario_labels(tibble::as_tibble(data$scenarios))
  scen_df$scen <- paste0("scenario_", seq_len(nrow(scen_df)))

  # Build tidy dataframe using generalized helper
  summary_df <- .make_50yr_metrics_tidy(
    summary_list,
    scen_df,
    summary_fn,
    metrics = c("harvest", "return", "escapement"),
    prefix = NULL  # Will produce "harvest_point", "return_point", etc.
  )

  # Extract just the point estimates and rename for Figure A
  summary_df$harvest    <- summary_df$harvest_point
  summary_df$return     <- summary_df$return_point
  summary_df$escapement <- summary_df$escapement_point

  # Prepare final plotting dataframe
  plot_df <- .prepare_plot_A_df(summary_df, selectivity_filter)

  # Output folder: FigureA inside output_dir
  output_dir <- file.path(output_dir, file_basename)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Plot
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = forcats::fct_inorder(.data$trends),
                 y = .data$value,
                 color = forcats::fct_inorder(.data$mgmt),
                 group = .data$mgmt)
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_y_continuous(
      expand = c(0.07, 0.07),
      limits = c(0, NA),
      labels = function(y) y / 1000
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL, color = "Method") +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$metric_label),
      cols = ggplot2::vars(.data$factorMSY),
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10),
      strip.text.y = ggplot2::element_text(size = 10),
      axis.line = ggplot2::element_line(linewidth = 0.1),
      axis.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      strip.placement = "outside"
    )

  # Save outputs
  plot_path <- file.path(output_dir, paste0(file_basename, ".pdf"))
  data_path <- file.path(output_dir, paste0(file_basename, ".csv"))
  ggplot2::ggsave(plot_path, p, width = 9, height = 9, units = "in")
  readr::write_csv(plot_df, data_path)

  message("Figure A saved to: ", plot_path)
  message("Data saved to: ", data_path)
  message("Summary function used: ", summary_fn)

  return(invisible(list(data = plot_df, plot = p)))
}

#' Helper: Reshape tidy dataframe for plotting Figure A
#'
#' Filters to desired selectivity and reshapes to long format
#' with metric labels for facetting.
#'
#' @param df Tidy dataframe with metrics and metadata.
#' @param selectivity_filter Which selectivity to keep.
#' @return Long-format dataframe ready for plotting.
#' @keywords internal
.prepare_plot_A_df <- function(df, selectivity_filter) {
  df |>
    dplyr::filter(.data$selectivity == !!selectivity_filter) |>
    tidyr::pivot_longer(
      cols = c("harvest", "return", "escapement"),
      names_to = "metric", values_to = "value"
    ) |>
    dplyr::mutate(
      metric_label = dplyr::case_when(
        .data$metric == "harvest"    ~ "Harvest (1000s)",
        .data$metric == "return"     ~ "Return (1000s)",
        .data$metric == "escapement" ~ "Escapement (1000s)"
      ),
      metric_label = factor(
        .data$metric_label,
        levels = c("Harvest (1000s)", "Escapement (1000s)", "Return (1000s)")
      )
    )
}
