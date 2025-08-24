#' Generate Figure A: Median summary plot of 50-year metrics
#'
#' @param data A list returned by `summarize_50_year_avg()`.
#' @param scenarios_all A tibble describing the scenarios (typically from `load_model_outputs()`).
#' @param use_quant Either "median" (default) or "mean", determines summary statistic used.
#' @param selectivity_filter A string for filtering scenarios by selectivity type (default = "unselective").
#' @param colors A character vector of colors for plotting (default = c("darkgray", "deepskyblue3", "orange")).
#' @param output_dir Directory where plot and CSV are saved (default = current directory).
#' @param file_basename Base filename (without extension) for saved files (default = "FigureA").
#'
#' @return A list with:
#'   - `plot`: the ggplot object
#'   - `data`: the tidy data.frame used for plotting
#' @export
make_figure_A <- function(data,
                          scenarios_all,
                          use_quant = "median",
                          selectivity_filter = "unselective",
                          colors = c("darkgray", "deepskyblue3", "orange"),
                          output_dir = ".",
                          file_basename = "FigureA") {

  # Extract and prepare scenario-level data
  df <- scenarios_all %>%
    dplyr::select(trends, factorMSY, mgmt, selectivity)

  # Helper to extract values based on use_quant
  extract_values <- function(metric_list, use_quant) {
    sapply(metric_list, function(scenario) {
      if (use_quant == "mean") {
        mean(scenario$mean, na.rm = TRUE)
      } else if (use_quant == "median") {
        if (!is.null(scenario$quantiles) && "q50" %in% colnames(scenario$quantiles)) {
          median(scenario$quantiles[, "q50"], na.rm = TRUE)
        } else {
          stop("Missing 'q50' in scenario$quantiles for one or more scenarios.")
        }
      } else {
        stop("use_quant must be either 'mean' or 'median")
      }
    })
  }

  # Build metric columns
  df$mean_50yr_harvest <- round(extract_values(data$harvest, use_quant), 0)
  df$mean_50yr_return <- round(extract_values(data$return, use_quant), 0)
  df$mean_50yr_escapement <- round(extract_values(data$escapement, use_quant), 0)

  # Filter selectivity
  df <- df %>% dplyr::filter(selectivity == selectivity_filter)

  # Tidy format for plotting
  df_tidy <- df %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("mean_50yr"),
                        names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric_label = dplyr::case_when(
        metric == "mean_50yr_harvest" ~ "Harvest (1000s)",
        metric == "mean_50yr_return" ~ "Return (1000s)",
        metric == "mean_50yr_escapement" ~ "Escapement (1000s)"
      ),
      metric_label = factor(
        metric_label,
        levels = c("Harvest (1000s)", "Escapement (1000s)", "Return (1000s)")
      ),
      factorMSY = as.character(.data$factorMSY),
      factorMSY = dplyr::case_when(
        factorMSY %in% c("0.75", "liberal")       ~ "liberal",
        factorMSY %in% c("1", "MSY")              ~ "MSY",
        factorMSY %in% c("1.5", "precautionary")  ~ "precautionary",
        TRUE                                       ~ factorMSY
      ),
      factorMSY = factor(factorMSY, levels = c("liberal", "MSY", "precautionary"))
    )

  # Build plot
  p <- ggplot2::ggplot(df_tidy, ggplot2::aes(x = forcats::fct_inorder(trends), y = value,
                                             color = forcats::fct_inorder(mgmt), group = mgmt)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_y_continuous(expand = c(0.07, 0.07), limits = c(0, NA),
                                labels = function(y) y / 1000) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL, color = "Method") +
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
      legend.text = ggplot2::element_text(size = 8)
    )

  g <- p + ggplot2::facet_grid(rows = ggplot2::vars(metric_label), cols = ggplot2::vars(factorMSY),
                               scales = "free_y", switch = "y") +
    ggplot2::theme(strip.placement = "outside")

  # Save plot and data
  plot_path <- file.path(output_dir, paste0(file_basename, ".pdf"))
  data_path <- file.path(output_dir, paste0(file_basename, "_data.csv"))
  ggplot2::ggsave(plot_path, g, width = 9, height = 9, units = "in")
  readr::write_csv(df_tidy, data_path)

  return(list(
    plot = g,
    data = df_tidy
  ))
}
