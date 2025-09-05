#' Generate Ohlberger Figure 4: % difference in run size vs TRM
#'
#' Reproduces the published Figure 4: median **percent difference** in mean
#' run size relative to the traditional Ricker method (TRM), by trends (x),
#' colored by method (YPR/DLM only), with columns for selectivity.
#'
#' @section Caption (from the publication):
#' Median % difference in average total run size from the time-invariant model.
#' Shown are differences for the two alternative estimation methods, the
#' yield-per-recruit analysis (YPR, blue) and the Dynamic Linear Model (DLM,
#' yellow), compared to the time-invariant Ricker model for three different
#' fishery selectivities (columns) given the simulated demographic trends
#' (x-axis). This figure is based on numbers presented in the second row of
#' Figure 3.
#'
#' @section Data & assumptions:
#' - Input is an **`ssi_run`** from [run_scenarios()].
#' - Uses `data$results$obs` for time series and `data$scenarios` for labels.
#' - Window is the **last `nyh` years** per iteration (default 50; uses up to
#'   50 from the tail of available rows).
#' - For each scenario we identify the **TRM reference** (same trends,
#'   selectivity, and factorMSY) and compute iteration-wise % differences
#'   for YPR/DLM vs TRM.
#' - Scenario labels are standardized via `.standardize_scenario_labels()`.
#'
#' @section Computation:
#' 1. For each scenario/iteration, compute mean run size over the last `nyh` years:
#'    \code{mean(obsRet)}.
#' 2. For each non-TRM scenario, pair with its TRM reference and compute
#'    \code{100 * (mean_ret - mean_ret_TRM) / mean_ret_TRM} per iteration.
#' 3. Take the **median across iterations** by scenario.
#' 4. Plot medians by trends; color = mgmt (YPR, DLM); facet by selectivity.
#'
#' @section Files written:
#' - \code{file.path(output_dir, "Figure4", "Figure4.pdf")}
#' - \code{file.path(output_dir, "Figure4", "Figure4.csv")}
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for
#' salmon body size declines in fishery management can reduce conservation
#' risks. *Fish and Fisheries*, 26, 113–130. https://doi.org/10.1111/faf.12869
#'
#' @param data `ssi_run` object from [run_scenarios()].
#' @param output_dir Directory for outputs (default: current working directory).
#'
#' @return Invisible list with `plot` (ggplot) and `data` (tidy plotting data).
#' @noRd
#' @keywords internal

#' @noRd
#' @keywords internal
.make_Ohlberger_figure_4 <- function(data,
                                     output_dir = tempdir(),
                                     file_basename = "Figure4",
                                     width_in = 6,
                                     height_in = 4) {
  stopifnot(inherits(data, "ssi_run"))
  .validate_run_params(data$parameters)

  # --- access recorded params & lists -----------------------------------------
  nyh      <- .get_nyh(data)
  obs_list <- .get_obs_list(data)

  # --- scenarios, standardized, MSY only; drop AL trend to match paper --------
  scen_df <- .scenarios_from_ssirun(data) |>
    dplyr::filter(.data$factorMSY == "MSY") |>
    dplyr::mutate(trends = droplevels(.data$trends)) |>
    dplyr::filter(.data$trends != "AL trends stabilized") |>
    droplevels()

  # prefer name-based mapping to results
  obs_names <- names(obs_list)
  j_idx <- if (!is.null(obs_names) && all(scen_df$scen_key %in% obs_names)) {
    match(scen_df$scen_key, obs_names)
  } else {
    scen_df$scen_num
  }
  if (anyNA(j_idx)) stop("Scenario key mapping to obs list failed.")

  # --- fixed post-history window (same for all scen/iter) ----------------------
  first_obs <- as.data.frame(obs_list[[ j_idx[1] ]][[1]])
  ny_obs    <- nrow(first_obs)
  year_index <- seq.int(nyh + 1L, ny_obs)
  if (length(year_index) <= 0) stop("Empty post-history window; check nyh/obs length.")

  # --- per-scenario, per-iteration mean return over the window -----------------
  mean_ret_last <- function(iter_obs_list) {
    vapply(iter_obs_list, function(df) {
      if (!is.data.frame(df) || !nrow(df) || is.null(df$obsRet)) return(NA_real_)
      mean(df$obsRet[year_index], na.rm = TRUE)
    }, numeric(1))
  }

  # compute vector of means for each scenario (keyed by scen_key)
  scen_keys <- scen_df$scen_key
  mean_ret_by_scen <- setNames(vector("list", length(scen_keys)), scen_keys)
  for (i in seq_along(scen_keys)) {
    key <- scen_keys[i]
    mean_ret_by_scen[[key]] <- mean_ret_last(obs_list[[ key ]])
  }

  # --- pair each YPR/DLM to TRM within the same cell ---------------------------
  # cell definition matches the paper's grouping
  cell_vars <- c("trends", "selectivity", intersect("factorMSY", names(scen_df)))
  scen_df$cell <- do.call(paste, c(scenv = scen_df[cell_vars], sep = "||"))

  trm_map <- scen_df |>
    dplyr::filter(.data$mgmt == "TRM") |>
    dplyr::select(.data$cell, trm_key = .data$scen_key)

  cand_df <- scen_df |>
    dplyr::filter(.data$mgmt %in% c("YPR", "DLM")) |>
    dplyr::left_join(trm_map, by = "cell")

  # iteration-wise % diff versus TRM reference
  diff_rows <- lapply(seq_len(nrow(cand_df)), function(i) {
    row  <- cand_df[i, ]
    keym <- row$scen_key
    keyr <- row$trm_key
    if (is.na(keyr)) return(NULL)

    vm <- mean_ret_by_scen[[keym]]
    vr <- mean_ret_by_scen[[keyr]]
    if (is.null(vm) || is.null(vr) || length(vm) != length(vr)) return(NULL)

    pct <- 100 * (vm - vr) / vr
    tibble::tibble(
      trends      = row$trends,
      selectivity = row$selectivity,
      mgmt        = row$mgmt,
      scen_key    = keym,
      iter        = seq_along(pct),
      pct_diff_run = pct
    )
  })
  diff_long <- dplyr::bind_rows(diff_rows)
  if (!nrow(diff_long)) stop("No paired TRM references found to compute differences.")

  # scenario-level median across iterations (paper)
  df_plot <- diff_long |>
    dplyr::group_by(.data$trends, .data$selectivity, .data$mgmt, .data$scen_key) |>
    dplyr::summarise(median = stats::median(.data$pct_diff_run, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(.data$mgmt %in% c("YPR", "DLM")) |>
    droplevels()

  # --- plot aesthetics: match paper -------------------------------------------
  colors <- c(TRM = "darkgray", YPR = "deepskyblue3", DLM = "orange")

  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = forcats::fct_inorder(trends),
      y = median,
      color = forcats::fct_inorder(mgmt),
      group = mgmt
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.1) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_point(shape = 1, size = 3, fill = NA, color = "black") +
    ggplot2::scale_colour_manual(values = colors[-1], breaks = c("YPR","DLM"), drop = FALSE) +
    ggplot2::scale_y_continuous(expand = c(0.1, 0.1), breaks = seq(-20, 30, 5)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "", y = "Median difference in run size (%)", color = "Method") +
    ggplot2::facet_grid(. ~ selectivity) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10),
      axis.line = ggplot2::element_line(linewidth = 0.1),
      axis.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 8),
      legend.background = ggplot2::element_blank()
    )

  # --- save & announce ---------------------------------------------------------
  outdir    <- .ensure_outdir(output_dir, file_basename)
  plot_path <- file.path(outdir, paste0(file_basename, ".pdf"))
  data_path <- file.path(outdir, paste0(file_basename, ".csv"))

  ggplot2::ggsave(plot_path, p, width = width_in, height = height_in, units = "in")
  readr::write_csv(df_plot, data_path)

  message("Figure 4 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  invisible(list(plot = p, data = df_plot, outdir = outdir))

}
