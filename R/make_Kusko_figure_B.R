#' Generate Figure B: Ellipses of Harvest vs Escapement (1000 iters)
#'
#' Builds 80% ellipses of (escapement, harvest) using per-iteration
#' 50-year means for each scenario. Facets by trends (rows) × factorMSY (cols),
#' colors by management method. Saves a single PDF and CSV, returns results invisibly.
#'
#' @param data `ssi_run` object from `run_scenarios()`.
#' @param output_dir Directory to save outputs (default: current working directory).
#'
#' @return invisible(list(plot = ggplot, data = tidy_plotting_data))
#' @noRd
#' @keywords internal

.make_Kusko_figure_B <- function(data, output_dir = ".") {
  # Fixed settings (match Figures A/C)
  selectivity_filter <- "unselective"
  colors <- c("darkgray", "deepskyblue3", "orange")
  file_basename <- "FigureB"

  # 1) Summarize per-iteration 50y means for each scenario
  summary_list <- .summarize_50_year_avg(data)  # uses last 50y by default

  # 2) Scenario metadata (standardize + stable IDs)
  scen_df <- standardize_scenario_labels(tibble::as_tibble(data$scenarios))
  scen_df$scen <- paste0("scenario_", seq_len(nrow(scen_df)))

  # 3) Build tidy per-iteration dataset for ellipses
  ellipse_df <- .make_figure_B_ellipse_data(summary_list, scen_df, selectivity_filter)

  # 4) Output folder: ./<output_dir>/FigureB
  output_dir <- file.path(output_dir, file_basename)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 5) Plot (style consistent with A/C)
  p <- ggplot2::ggplot(
    ellipse_df,
    ggplot2::aes(
      x = .data$escapement,
      y = .data$harvest,
      color = forcats::fct_inorder(.data$mgmt),
      group = .data$mgmt
    )
  ) +
    ggplot2::stat_ellipse(level = 0.80) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_y_continuous(
      name   = "Harvest (1000s)",
      limits = c(0, 210000),
      breaks = seq(0, 250000, 50000),
      labels = function(y) y / 1000,
      expand = c(0.05, 0.05)
    ) +
    ggplot2::scale_x_continuous(
      name   = "Escapement (1000s)",
      limits = c(0, 210000),
      breaks = seq(0, 250000, 50000),
      labels = function(x) x / 1000
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$trends),
      cols = ggplot2::vars(.data$factorMSY),
      switch = "y"
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = "Method") +
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

  # 6) Save outputs
  plot_path <- file.path(output_dir, paste0(file_basename, ".pdf"))
  data_path <- file.path(output_dir, paste0(file_basename, ".csv"))
  ggplot2::ggsave(plot_path, p, width = 9, height = 9, units = "in")
  readr::write_csv(ellipse_df, data_path)

  message("Figure B saved to: ", plot_path)
  message("Data saved to: ", data_path)

  return(invisible(list(plot = p, data = ellipse_df)))
}

#' Helper: Build per-iteration ellipse data for Figure B
#'
#' Uses per-iteration 50y means from `.summarize_50_year_avg()` for harvest/escapement.
#' Keeps only scenarios matching `selectivity_filter`, and joins scenario metadata.
#'
#' @param summary_list Output from `.summarize_50_year_avg()`
#' @param scen_df Scenario metadata tibble with standardized labels and `scen`
#' @param selectivity_filter Character selectivity label to include
#' @return Tidy dataframe with columns: scen, iter, escapement, harvest, trends, mgmt, factorMSY, selectivity
#' @keywords internal
.make_figure_B_ellipse_data <- function(summary_list, scen_df, selectivity_filter = "unselective") {
  # filter scenarios after standardization
  scen_keep <- dplyr::filter(scen_df, .data$selectivity == !!selectivity_filter)

  # safety: names of summary_list elements should be "scenario_#"
  # create a look-up from scen -> index
  idx_lookup <- setNames(seq_len(nrow(scen_df)), scen_df$scen)

  df_list <- purrr::map(scen_keep$scen, function(sid) {
    i <- idx_lookup[[sid]]
    if (is.null(i)) stop("Scenario id '", sid, "' not found in scenario metadata.")

    # per-iteration 50y means
    esc_means <- summary_list$escapement[[sid]]$means
    har_means <- summary_list$harvest[[sid]]$means

    if (!is.numeric(esc_means) || !is.numeric(har_means)) {
      stop("Non-numeric per-iteration means for scenario ", sid, ".")
    }
    tibble::tibble(
      scen       = sid,
      iter       = seq_along(esc_means),
      escapement = esc_means,
      harvest    = har_means
    ) |>
      dplyr::left_join(scen_keep, by = c("scen" = "scen"))
  })

  df <- dplyr::bind_rows(df_list)
  # idempotent standardization post-join (keeps factor levels stable)
  standardize_scenario_labels(df)
}
