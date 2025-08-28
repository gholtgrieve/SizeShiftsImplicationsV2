#' Generate Ohlberger Figure 2: S_MSY at first post-historical review
#'
#' Reproduces the published Figure 2: 50% and 80% intervals for \eqn{S_{MSY}}
#' at the first review year after the historical period. Columns = selectivity,
#' x = trends, color = management method (TRM/YPR/DLM). Uses per-iteration
#' \eqn{S_{MSY}} values at that review year across scenarios.
#'
#' @section Caption (from the publication):
#' The plot shows the distribution of \eqn{S_{MSY}} estimates at the first
#' post-historical review year across simulation iterations. Points indicate
#' medians, with vertical bars showing central 50% and 80% intervals. Panels are
#' split by gear selectivity, the x-axis groups demographic trend scenarios, and
#' colors denote management method.
#'
#' @section Data & assumptions:
#' - Input is an **`ssi_run`** from [run_scenarios()].
#' - `data$results$S_msy` is a nested list named by original scenario ids
#'   (`"scen_<scen_num>"`), each containing a list of numeric vectors (one value
#'   per **review year**) for each iteration.
#' - Review years follow the simulation logic:
#'   `seq(nyi + 20L, ny, by = goalfreq)`. The figure uses the **first review at
#'   or after** calendar year `nyi + nyh` (i.e., the first post-historical review).
#' - Scenario labels and factor orders are standardized via
#'   `standardize_scenario_labels()`.
#'
#' @section Computation:
#' 1. Identify the review year index at/after `nyi + nyh`.
#' 2. For each scenario/iteration, extract that single \eqn{S_{MSY}} value.
#' 3. For each (trend × selectivity × mgmt) cell, compute the **median** and the
#'    central **50%** and **80%** intervals across iterations.
#' 4. Plot medians with 50% and 80% intervals; facet by selectivity; x = trends,
#'    color = mgmt (TRM/YPR/DLM).
#'
#' @section Files written:
#' - `file.path(output_dir, "Figure2", "Figure2.pdf")`
#' - `file.path(output_dir, "Figure2", "Figure2.csv")`
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for
#' salmon body size declines in fishery management can reduce conservation
#' risks. *Fish and Fisheries*, 26, 113–130. https://doi.org/10.1111/faf.12869
#'
#' @param data `ssi_run` object from [run_scenarios()].
#' @param output_dir Directory to save outputs (default: current working directory).
#'
#' @return Invisible list with `plot` (ggplot) and `data` (tidy plotting data).
#' @noRd
#' @keywords internal

.make_Ohlberger_figure_2 <- function(data, output_dir = ".") {
  colors <- c("darkgray", "deepskyblue3", "orange")
  file_basename <- "Figure2"

  # Summary functions (match publication: 80% and 50% intervals around median)
  summary_CI80 <- function(x) data.frame(
    y = stats::median(x, na.rm = TRUE),
    ymin = stats::quantile(x, 0.10, na.rm = TRUE),
    ymax = stats::quantile(x, 0.90, na.rm = TRUE)
  )
  summary_CI50 <- function(x) data.frame(
    y = stats::median(x, na.rm = TRUE),
    ymin = stats::quantile(x, 0.25, na.rm = TRUE),
    ymax = stats::quantile(x, 0.75, na.rm = TRUE)
  )

  # Scenarios (standardized); original script excludes "continuing trends"
  scen_df_all <- .scenarios_from_ssirun(data)
  scen_use <- dplyr::filter(scen_df_all, .data$trends != "ASL trends continued")

  # Locate S_MSY series and align to the first post-historical review
  smsy_series <- .get_smsy_series_list(data)              # data$results$S_msy
  series_len  <- length(smsy_series[[ scen_use$scen_key[1] ]][[1]])
  review_years <- .get_review_years(data, series_length = series_len)
  idx <- .get_review_index_post_history(data, review_years)

  # Build long per-iteration vector at target review index, using scen_key
  df_list <- lapply(seq_len(nrow(scen_use)), function(r) {
    key <- scen_use$scen_key[r]          # e.g., "scen_13"
    iters <- smsy_series[[key]]
    if (is.null(iters)) {
      stop("S_MSY series not found for ", key, ". Available: ",
           paste0(names(smsy_series), collapse = ", "))
    }
    vals <- vapply(iters, function(v) if (length(v) >= idx) v[[idx]] else NA_real_, numeric(1))
    tibble::tibble(
      scen  = scen_use$scen[r],
      value = vals
    ) |>
      dplyr::left_join(scen_use[r, , drop = FALSE], by = "scen")
  })
  plot_df <- dplyr::bind_rows(df_list)
  plot_df <- plot_df[stats::complete.cases(plot_df$value), , drop = FALSE]

  # Output folder
  outdir <- .ensure_outdir(output_dir, file_basename)

  # Plot (preserve published styling)
  jig <- ggplot2::position_dodge(width = 0.6)
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = forcats::fct_inorder(.data$trends),
                 y = .data$value,
                 color = forcats::fct_inorder(.data$mgmt))
  ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::stat_summary(fun.data = summary_CI80, position = jig, linewidth = 0.2) +
    ggplot2::stat_summary(fun.data = summary_CI50, position = jig, linewidth = 0.6) +
    ggplot2::scale_y_continuous(breaks = seq(0, 20000, 4000), expand = c(0.02, 0.02)) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = "Method", shape = "Method") +
    ggplot2::labs(x = "", y = expression("" * S[MSY] * "")) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$selectivity), scales = "free_y") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10),
      axis.line = ggplot2::element_line(linewidth = 0.1),
      axis.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8)
    )

  # Save
  plot_path <- file.path(outdir, paste0(file_basename, ".pdf"))
  data_path <- file.path(outdir, paste0(file_basename, ".csv"))
  ggplot2::ggsave(plot_path, p, width = 6, height = 4, units = "in")
  readr::write_csv(plot_df, data_path)

  message("Figure 2 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  return(invisible(list(data = plot_df, plot = p)))
}
