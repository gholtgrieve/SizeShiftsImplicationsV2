#' Generate Figure C: Time series of escapement, harvest, and closure metrics
#'   by trend, harvest strategy, and management
#'
#' Accepts an `ssi_run` object from `run_scenarios()`.
#'
#' Internally calls `.summarize_by_year()` (which defaults to the last 50 years
#' for `ssi_run` inputs). Produces ribbons for central intervals and overlays
#' either probability of zero harvest or cumulative closures on the secondary axis.
#'
#' @param data `ssi_run` object.
#' @param output_dir Directory to save plots and data (default: current working directory).
#' @param closure_metric Character: "probability" (default, original behavior with 5-year smoothing)
#'   or "cumulative" (cumulative count of closures, no smoothing).
#'
#' @return A list with:
#'   - `data`: tidy data.frame used for plotting
#'   - `plots`: named list of ggplot objects by management strategy
#' @noRd
#' @keywords internal

.make_Kusko_figure_C <- function(data, output_dir = ".", closure_metric = "cumulative", statistic = "median") {

  # Validate parameters
  closure_metric <- match.arg(closure_metric, choices = c("probability", "cumulative"))
  statistic <- match.arg(statistic, choices = c("mean", "median"))

  summary_stats <- c(mean = "Mean", ymin = "25%", ymax = "75%")
  selectivity_filter <- "unselective"
  colors <- c("darkgray", "deepskyblue3", "orange")
  file_basename <- "FigureC"

  # Make and set output directory
  output_dir <- file.path(output_dir, file_basename)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Summarize and prep metadata
  summary_list <- .summarize_by_year(data)
  obs_list <- data$results$obs
  scen_df_all <- .standardize_scenario_labels(tibble::as_tibble(data$scenarios))

  if (!"scen_num" %in% names(scen_df_all)) {
    scen_df_all$scen_num <- seq_len(length(obs_list))
  }
  scen_df_all$scen <- paste0("scenario_", scen_df_all$scen_num)

  scen_names <- scen_df_all$scen
  niter <- length(obs_list[[1]])
  iter_names <- paste0("iter_", seq_len(niter))

  known_order <- c("TRM", "YPR", "DLM")
  mgmt_names <- intersect(known_order, unique(scen_df_all$mgmt))
  if (length(mgmt_names) == 0) mgmt_names <- unique(scen_df_all$mgmt)

  year_names <- names(summary_list$escapement[[1]])
  nyh <- length(year_names)
  year_rows_fun <- function(sim) {
    rows <- seq.int(nrow(sim) - nyh + 1L, nrow(sim))
    rows[rows >= 1L]
  }

  # Create summary dataframes
  summary_df <- .make_tidy_summary(summary_list, scen_df_all, scen_names, summary_stats)
  closure_df <- .compute_closure_metrics(obs_list, scen_names, iter_names, year_rows_fun, scen_df_all)
  plot_df <- .prepare_plot_df(summary_df, closure_df, selectivity_filter)

  # Determine scaling and smoothing based on closure_metric
  if (closure_metric == "probability") {
    #  5-year rolling mean on probability
    metric_df <- plot_df |>
      dplyr::filter(.data$label == "Harvest (1000s)") |>
      dplyr::arrange(.data$scen, .data$mgmt, .data$trends, .data$factorMSY, .data$year) |>
      dplyr::group_by(.data$scen, .data$mgmt, .data$trends, .data$factorMSY) |>
      dplyr::mutate(plot_value = zoo::rollmean(.data$prob_zero_harvest, k = 5, fill = NA, align = "center")) |>
      dplyr::ungroup()

    half_primary <- 220000 / 2
    scaling_factor <- half_primary  # Scale 0-1 probability to 0-110000
    secondary_breaks <- seq(0, 1, 0.5)
    secondary_labels <- scales::percent_format(accuracy = 1)
    line_label <- "Probability of zero harvest\n(right axis)"

  } else {  # cumulative
    # New behavior: cumulative closures, no smoothing
    if (statistic == "mean") {
      metric_df <- plot_df |>
        dplyr::filter(.data$label == "Harvest (1000s)") |>
        dplyr::arrange(.data$scen, .data$mgmt, .data$trends, .data$factorMSY, .data$year) |>
        dplyr::mutate(plot_value = .data$mean_cumulative_closures)
    } else {
      metric_df <- plot_df |>
        dplyr::filter(.data$label == "Harvest (1000s)") |>
        dplyr::arrange(.data$scen, .data$mgmt, .data$trends, .data$factorMSY, .data$year) |>
        dplyr::mutate(plot_value = .data$median_cumulative_closures)
    }

    max_closures <- max(metric_df$plot_value, na.rm = TRUE)
    secondary_max <- if(max_closures <= 25) 25 else 50
    scaling_factor <- 220000 / secondary_max
    secondary_breaks <- seq(0, secondary_max, by = if(secondary_max == 25) 5 else 10)
    secondary_labels <- function(x) as.character(round(x))
    line_label <- "Cumulative fishery closures\n(right axis)"
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plot_component <- function(df, mgmt_name, metric_data) {
    # Build base plot
    p <- ggplot2::ggplot(df |> dplyr::filter(.data$mgmt == mgmt_name),
                         ggplot2::aes(x = .data$year, y = .data$mean, ymin = .data$ymin, ymax = .data$ymax,
                                      color = forcats::fct_inorder(.data$label),
                                      group = .data$label, fill = .data$label)) +
      ggplot2::geom_ribbon(alpha = 0.3, color = NA, show.legend = FALSE) +
      ggplot2::geom_line(linewidth = 0.8)

    # Add the metric line
    p <- p + ggplot2::geom_line(
      data = metric_data |> dplyr::filter(.data$mgmt == mgmt_name),
      ggplot2::aes(x = .data$year,
                   y = .data$plot_value * scaling_factor,
                   color = line_label),
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
      ggplot2::scale_colour_manual(
        values = stats::setNames(
          c(colors[3], colors[2], colors[1]),
          c(line_label, "Escapement (1000s)", "Harvest (1000s)")
        ),
        breaks = c("Harvest (1000s)", "Escapement (1000s)", line_label)
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Escapement (1000s)" = colors[2],
          "Harvest (1000s)"    = colors[1],
          line_label           = colors[3]
        )
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0.05, 0.05),
        limits = c(0, 220000),
        breaks = seq(0, 220000, 50000),
        labels = ~ .x / 1000,
        sec.axis = ggplot2::sec_axis(
          ~ . / scaling_factor,
          name = NULL,
          breaks = secondary_breaks,
          labels = secondary_labels
        )
      ) +
      ggplot2::scale_x_continuous(name = "Year", limits = c(1, length(year_names)), breaks = seq(0, length(year_names), 10)) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$trends), cols = ggplot2::vars(.data$factorMSY),
                          scales = "free_y", switch = "y") +
      ggplot2::theme_classic(base_family = "sans") +
      ggplot2::labs(x = "", y = "", color = "") +
      ggplot2::theme(
        strip.background   = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "dashed", linewidth = 0.3),
        strip.text.x       = ggplot2::element_text(size = 10),
        strip.text.y       = ggplot2::element_text(size = 10),
        axis.line          = ggplot2::element_line(linewidth = 0.3),
        axis.text          = ggplot2::element_text(size = 10),
        axis.text.x        = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title         = ggplot2::element_text(size = 12),
        panel.border       = ggplot2::element_rect(fill = NA, linewidth = 0.5),
        legend.key.size    = grid::unit(0.5, 'cm'),
        legend.title       = ggplot2::element_text(size = 10),
        legend.text        = ggplot2::element_text(size = 8),
        strip.placement    = "outside"
      ) -> p

    pdf_path <- file.path(output_dir, paste0(file_basename, "_", mgmt_name, ".pdf"))
    ggplot2::ggsave(pdf_path, p, width = 9, height = 9, units = "in")
    p
  }

  plots <- purrr::map(mgmt_names, ~ plot_component(plot_df, .x, metric_df))
  names(plots) <- mgmt_names

  csv_path <- file.path(output_dir, paste0(file_basename, ".csv"))
  readr::write_csv(dplyr::arrange(plot_df, .data$mgmt, .data$trends, .data$factorMSY, .data$scen, .data$year), csv_path)

  # Report all saved plot paths
  pdf_paths <- file.path(output_dir, paste0(file_basename, "_", mgmt_names, ".pdf"))
  for (path in pdf_paths) {
    message("Figure C plot saved to: ", path)
  }

  # Report saved CSV path and metrics reported
  message("Data saved to: ", csv_path)
  message("Closure metric used: ", closure_metric)
  if (closure_metric == "cumulative") {
    message("Statistic used for cumulative closures: ", statistic)
  }

  return(invisible(list(data = plot_df, plots = plots, closure_metric = closure_metric)))
}




#' Helper: Tidy simulation summaries for plotting
#'
#' Converts list-of-lists summary output from `.summarize_by_year()` into a single
#' tidy data frame for escapement and harvest, with ribbons and means.
#'
#' @param summary_list Output from `.summarize_by_year()`.
#' @param scen_df Full scenario metadata.
#' @param scen_names Vector of scenario labels.
#' @param summary_stats Named vector mapping new names to original stats (e.g., mean = "Mean").
#' @return Tidy data.frame with `scen`, `year`, `mean`, `ymin`, `ymax`, `label`, and scenario metadata.
#' @keywords internal
.make_tidy_summary <- function(summary_list, scen_df, scen_names, summary_stats) {
  year_names <- names(summary_list$escapement[[1]])
  summary_tbl <- tibble::tibble(
    stat_in_data = unname(summary_stats),
    new_col_name = names(summary_stats)
  )
  data_pluck <- function(component, stat) {
    purrr::map_dfr(scen_names, function(scen) {
      purrr::map_dfr(year_names, function(yr) {
        tibble::tibble(
          scen = scen,
          year = which(year_names == yr),
          value = summary_list[[component]][[scen]][[yr]][[stat]]
        )
      })
    })
  }
  tidy_list <- purrr::map(c("escapement", "harvest"), function(component) {
    dfs <- purrr::pmap(summary_tbl, function(stat_in_data, new_col_name) {
      df <- data_pluck(component, stat_in_data)
      df <- dplyr::rename(df, !!new_col_name := value)
      df$label <- ifelse(component == "escapement", "Escapement (1000s)", "Harvest (1000s)")
      df
    })
    Reduce(function(left, right) dplyr::full_join(left, right, by = c("scen", "year", "label")), dfs)
  })
  dplyr::bind_rows(tidy_list) |>
    dplyr::left_join(scen_df, by = "scen") |>
    dplyr::mutate(label = factor(.data$label, levels = c("Harvest (1000s)", "Escapement (1000s)")))
}



#' Helper: Filter and finalize tidy plot data
#'
#' Joins harvest summaries with closure metrics and filters by selectivity.
#'
#' @param summary_df Output of `.make_tidy_summary()`.
#' @param closure_df Output of `.compute_closure_metrics()`.
#' @param selectivity_filter Which selectivity level to keep (e.g., "unselective").
#' @return Final tidy data.frame for plotting.
#' @keywords internal
.prepare_plot_df <- function(summary_df, closure_df, selectivity_filter) {
  dplyr::left_join(
    summary_df,
    closure_df |> dplyr::select("scen", "year", "n_obs", "n_zero_harvest",
                                "prob_zero_harvest", "mean_cumulative_closures",
                                "median_cumulative_closures", "q25_cumulative_closures",
                                "q75_cumulative_closures"),
    by = c("scen", "year")
  ) |> dplyr::filter(.data$selectivity == selectivity_filter)
}
