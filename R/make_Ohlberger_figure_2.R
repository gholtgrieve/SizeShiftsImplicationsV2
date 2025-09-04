#' Generate Ohlberger Figure 2: \eqn{S_{MSY}} at first post-historical review
#'
#' Reproduces the published Figure 2. For each scenario, we extract the
#' \eqn{S_{MSY}} value at the **first review at or after** the end of the
#' historical period (\eqn{nyi + nyh}), then summarize **across iterations**
#' within each (trends × selectivity × mgmt) cell. The plot shows medians
#' (points) with central **50%** and **80%** intervals (vertical bars).
#' Columns = gear selectivity; x = demographic trend scenarios; color =
#' management method (TRM/YPR/DLM).
#'
#' @section Data source:
#' - Input is an **`ssi_run`** from [run_scenarios()].
#' - Required: ssi_run$parameters must include nyi, nyh, ny, goalfreq,
#'   firstrev, review_years. The function validates these and errors if missing.
#' - Uses `data$results$S_msy`: a nested list keyed by scenario (e.g.,
#'   `"scen_<num>"`), where each element is a list of numeric vectors
#'   (one value per **review year**) for each iteration.
#'
#' @section Review year selection:
#' Review years follow the simulation setup:
#' \code{seq(nyi + firstrev, ny, by = goalfreq)}. The figure uses the **first
#' review at or after** calendar year \code{nyi + nyh} (i.e., the first
#' post-historical review). Internally this is resolved via
#' \code{.get_review_years()} and \code{.get_review_index_post_history()}.
#'
#' @section Scenario filtering & labeling:
#' - Scenario labels and factor orders are standardized via
#'   \code{.standardize_scenario_labels()}.
#' - Keep only \code{factorMSY == "MSY"}.
#' - Exclude \code{trends == "ASL trends continued"}.
#' - Facet by \code{selectivity} with \code{scales = "free_y"}.
#'
#' @section Statistics (across iterations):
#' For each (trends × selectivity × mgmt) cell, compute the **median** and the
#' central **50%** (\eqn{0.25, 0.75}) and **80%** (\eqn{0.10, 0.90}) intervals
#' of the per-iteration \eqn{S_{MSY}} values at the selected review.
#'
#' @section Plot design:
#' - x: \code{trends} (in order); y: \eqn{S_{MSY}} (native units; no scaling).
#' - color: \code{mgmt} with palette \code{c("darkgray","deepskyblue3","orange")}.
#' - Two \code{stat_summary()} layers: 80% band (thin) and 50% band (thicker),
#'   both centered on the median.
#' - y label: \eqn{S_{MSY}} rendered via \code{expression(""*S[MSY]*"")}.
#'
#' @section Files written:
#' - \code{file.path(output_dir, "Figure2", "Figure2.pdf")}
#' - \code{file.path(output_dir, "Figure2", "Figure2.csv")}
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for
#' salmon body size declines in fishery management can reduce conservation
#' risks. *Fish and Fisheries*, 26, 113–130. https://doi.org/10.1111/faf.12869
#'
#' @param data `ssi_run` object from [run_scenarios()].
#' @param output_dir Directory to save the figure (default: temp directory).
#' @param file_basename Base name for the output folder/file (default: "Figure2").
#' @param width_in,height_in Dimensions (inches) for the saved PDF.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{plot}: the ggplot object
#'   \item \code{data}: tidy plotting data (one value per iteration)
#'   \item \code{outdir}: output directory used
#' }
#' @keywords internal


#' Make Ohlberger Figure 2 (internal)
#' @keywords internal
.make_Ohlberger_figure_2 <- function(data, output_dir = tempdir(),
                                     file_basename = "Figure2",
                                     width_in = 6, height_in = 4) {
  stopifnot(inherits(data, "ssi_run"))

  # Hard requirement: params must exist and validate
  .validate_run_params(data$parameters)

  # Helpers and inputs
  review_years <- .get_review_years(data)
  idx_hist     <- .get_review_index_post_history(data, review_years)
  scenarios_df <- .scenarios_from_ssirun(data)  # adds scen_key, standardized labels
  smsy_list    <- .get_smsy_series_list(data)   # data$results$S_msy
  if (!is.finite(idx_hist) || idx_hist < 1L || idx_hist > length(review_years)) {
    stop("Could not locate first post-historical review in `review_years`.")
  }

  # Keep only MSY factor and drop "continued" trends (per paper)
  scen_use <- scenarios_df |>
    dplyr::filter(.data$factorMSY == "MSY") |>
    dplyr::filter(.data$trends != "ASL trends continued")

  # Build long table: one S_MSY per (scenario, iteration) at review index
  pull_smsy_at <- function(j) {
    iters <- smsy_list[[j]]
    vapply(iters, function(v) {
      if (is.null(v) || !length(v)) return(NA_real_)
      x <- as.numeric(v)
      if (!is.finite(idx_hist) || idx_hist < 1L || idx_hist > length(x)) return(NA_real_)
      x[idx_hist]
    }, numeric(1))
  }

  # Map standardized scenarios back to positions in results list
  # Assume scen_key == names(smsy_list); if not, fall back to row order
  smsy_names <- names(smsy_list)
  if (!is.null(smsy_names) && "scen_key" %in% names(scen_use) &&
      all(scen_use$scen_key %in% smsy_names)) {
    j_idx <- match(scen_use$scen_key, smsy_names)
  } else {
    # fallback to row order
    j_idx <- scen_use$scen_num
  }

  vals <- lapply(j_idx, pull_smsy_at)
  names(vals) <- scen_use$scen

  plot_df <- scen_use |>
    dplyr::mutate(values = vals) |>
    tidyr::unnest_wider(values, names_sep = "_iter_") |>
    tidyr::pivot_longer(dplyr::starts_with("values_iter_"),
                        names_to = "iteration",
                        values_to = "value",
                        values_drop_na = TRUE) |>
    dplyr::select(.data$trends, .data$selectivity, .data$mgmt, .data$value)

  # Colors per the paper
  colors <- c("darkgray", "deepskyblue3", "orange")

  # Summary functions (median + central intervals)
  summary_CI80 <- function(x) {
    data.frame(
      y    = stats::median(x, na.rm = TRUE),
      ymin = stats::quantile(x, probs = 0.10, na.rm = TRUE, names = FALSE),
      ymax = stats::quantile(x, probs = 0.90, na.rm = TRUE, names = FALSE)
    )
  }
  summary_CI50 <- function(x) {
    data.frame(
      y    = stats::median(x, na.rm = TRUE),
      ymin = stats::quantile(x, probs = 0.25, na.rm = TRUE, names = FALSE),
      ymax = stats::quantile(x, probs = 0.75, na.rm = TRUE, names = FALSE)
    )
  }

  jig <- ggplot2::position_dodge(width = 0.6)

  p <- plot_df |>
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_inorder(.data$trends),
                                 y = .data$value,
                                 color = forcats::fct_inorder(.data$mgmt))) +
    ggplot2::scale_color_manual(
      values = c("TRM"="darkgray","YPR"="deepskyblue3","DLM"="orange"),
      breaks = c("TRM","YPR","DLM"),
      drop   = FALSE) +
    ggplot2::stat_summary(fun.data = summary_CI80, position = jig, linewidth = 0.2, na.rm = TRUE) +
    ggplot2::stat_summary(fun.data = summary_CI50, position = jig, linewidth = 0.6, na.rm = TRUE) +
    ggplot2::scale_y_continuous(breaks = seq(0, 20000, 4000), expand = c(0.02, 0.02)) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = "Method", x = "", y = expression(""*S[MSY]*"")) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$selectivity), scales = "free_y") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x     = ggplot2::element_text(size = 10),
      axis.line        = ggplot2::element_line(linewidth = 0.1),
      axis.text        = ggplot2::element_text(size = 10),
      axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title       = ggplot2::element_text(size = 10),
      panel.border     = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.key.size  = grid::unit(0.5, "cm"),
      legend.title     = ggplot2::element_text(size = 10),
      legend.text      = ggplot2::element_text(size = 8)
    )

  # Ensure output directory
  outdir <- .ensure_outdir(output_dir, file_basename)
  # Save .pdf
  ggplot2::ggsave(
    file.path(outdir, paste0(file_basename, ".pdf")),
    p, width = width_in, height = height_in, units = "in"
  )
    # Save .csv (tidy per-iteration values for reproducibility / inspection)
  readr::write_csv(
    plot_df,
    file.path(outdir, paste0(file_basename, ".csv"))
  )

  invisible(list(plot = p, data = plot_df, outdir = outdir))
}
