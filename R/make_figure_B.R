#' Generate Figure B: Tradeoff plot of 50-year average escapement vs harvest by scenario and harvest strategy.
#'
#' @param data A list returned by `summarize_50_year_avg()`.
#' @param scenarios_all A tibble describing the scenarios (typically from `load_model_outputs()`).
#' @param selectivity_filter A string for filtering scenarios by selectivity type (default = "unselective").
#' @param colors A character vector of colors for plotting (default = c("darkgray", "deepskyblue3", "orange")).
#' @param output_dir Directory where plot and CSV are saved (default = current directory).
#' @param file_basename Base filename (without extension) for saved files (default = "FigureB").
#'
#' @return A list with:
#'   - `plot`: the ggplot object
#'   - `data`: the tidy data.frame used for plotting
#' @export
make_figure_B <- function(data,
                          scenarios_all,
                          selectivity_filter = "unselective",
                          colors = c("darkgray", "deepskyblue3", "orange"),
                          output_dir = ".",
                          file_basename = "FigureB") {
  
  # Filter scenarios_all
  filtered_meta <- scenarios_all %>%
    dplyr::filter(selectivity == selectivity_filter) %>%
    dplyr::mutate(scenario_id = paste0("scenario_", dplyr::row_number()))
  
  # Build long data for plotting
  df_list <- purrr::map(filtered_meta$scenario_id, function(sid) {
    i <- which(names(data$escapement) == sid)
    tibble::tibble(
      scenario = sid,
      escapement = data$escapement[[i]]$mean,
      harvest = data$harvest[[i]]$mean,
      iter = seq_along(data$escapement[[i]]$mean),
      trends = filtered_meta$trends[i],
      mgmt = filtered_meta$mgmt[i],
      factorMSY = filtered_meta$factorMSY[i]
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  
  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = escapement, y = harvest,
                                        color = forcats::fct_inorder(mgmt), group = mgmt)) +
    ggplot2::stat_ellipse(level = 0.80) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_y_continuous(expand=c(0.05,0.05),limits=c(0,210000), name = "Harvest (1000s)", 
                                breaks = seq(0,250000,50000), labels = function(y) y / 1000) +
    ggplot2::scale_x_continuous(name = "Escapement (1000s)",limits=c(0,210000), 
                                breaks = seq(0,250000,50000), labels = function(x) x / 1000) +
    ggplot2::facet_grid(rows = ggplot2::vars(trends), cols = ggplot2::vars(factorMSY),
                        switch = "y") +
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
  
  # Save outputs
  plot_path <- file.path(output_dir, paste0(file_basename, ".pdf"))
  data_path <- file.path(output_dir, paste0(file_basename, "_data.csv"))
  ggplot2::ggsave(plot_path, p, width = 9, height = 9, units = "in")
  readr::write_csv(df, data_path)
  
  return(list(
    plot = p,
    data = df
  ))
}
