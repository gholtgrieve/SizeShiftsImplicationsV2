#' Generate Ohlberger Figure 3: Performance metrics over last 50 years
#'
#' Reproduces the published Figure 3: for each scenario, compute iteration-wise
#' metrics over the last 50 years of observations, summarize across iterations
#' by the **median**, and plot by *trends* (x), colored by management method,
#' with columns for *selectivity* and rows for four metrics.
#'
#' @section Caption (from the publication):
#' Evaluation of fishery management performance under different scenarios. Shown
#' are median values of four performance metrics (rows) under three different
#' selectivity regimes (columns) for each estimation method (colours) as a
#' function of the demographic trends. Evaluated performance metrics were mean
#' harvest, mean run size, probability that recruitment was above 50% of the
#' maximum recruitment (Rmax), and the probability that the spawner escapement
#' was above 50% of the equilibrium spawner abundance (S0). Simulated
#' selectivity regimes were a small-mesh (6-in.), unselective or large-mesh
#' (8.5-in.) fishery. The three estimation methods were the traditional Ricker
#' model (TRM, grey), the yield-per-recruit analysis (YPR, blue) and the Dynamic
#' Linear Model (DLM, yellow). Age-sex-length (ASL) trends were either not
#' included, stabilised after 50 years or were assumed to continue into the
#' future (x-axis). Simulations were run for 100 years, and fishery performance
#' was evaluated for the last 50 years.
#'
#' @section Data & assumptions:
#' - Input is an **`ssi_run`** from [run_scenarios()].
#' - Uses `data$results$obs` (observed time series per iteration) and
#'   `data$results$sr_sim` (per-iteration Ricker parameter fits).
#' - Window is the **last `nyh` rows** per iteration (typically 50).
#' - Scenario labels and factor orders are standardized via
#'   `.standardize_scenario_labels()`.
#'
#' @section Computation:
#' For each scenario \eqn{j} and iteration \eqn{k}, over the last `nyh` years:
#' \itemize{
#'   \item \strong{Mean harvest}: \code{mean(obsHarv)}
#'   \item \strong{Mean run size}: \code{mean(obsRet)}
#'   \item \strong{Prob. > 0.5 × Rmax}: \code{mean(recRec > 0.5 * Rmax)}
#'   \item \strong{Prob. > 0.5 × S0}:   \code{mean(obsEsc > 0.5 * S0)}
#' }
#' with \eqn{R_{\max} = (\alpha/\beta)\exp(-1)} and \eqn{S_{0} = \log(\alpha)/\beta},
#' taking \eqn{\alpha,\beta} from `sr_sim`. Then take the **median across iterations**
#' per scenario and facet rows by metric, columns by selectivity; x = trends;
#' color = mgmt (TRM/YPR/DLM).
#'
#' @section Files written:
#' - `file.path(output_dir, "Figure3", "Figure3.pdf")`
#' - `file.path(output_dir, "Figure3", "Figure3.csv")`
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for
#' salmon body size declines in fishery management can reduce conservation
#' risks. *Fish and Fisheries*, 26, 113–130. https://doi.org/10.1111/faf.12869
#'
#' @param data `ssi_run` object from [run_scenarios()].
#' @param output_dir Directory where outputs are written (default: current working directory).
#'
#' @return Invisible list with `plot` (ggplot) and `data` (tidy plotting data).
#' @noRd
#' @keywords internal

#' @noRd
#' @keywords internal
.make_Ohlberger_figure_3 <- function(data,
                                     output_dir = tempdir(),
                                     file_basename = "Figure3",
                                     width_in = 5.5,
                                     height_in = 7) {
  stopifnot(inherits(data, "ssi_run"))
  .validate_run_params(data$parameters)

  # Accessors (recorded in run)
  nyh      <- .get_nyh(data)
  obs_list <- .get_obs_list(data)
  sr_list  <- .get_sr_params_list(data)

  # Scenario frame with standardized labels/keys
  scen_df <- .scenarios_from_ssirun(data) |>
    dplyr::filter(.data$factorMSY == "MSY") |>
    dplyr::mutate(trends = droplevels(.data$trends)) |>
    # After standardization, exclude "AL trends stabilized"
    dplyr::filter(.data$trends != "AL trends stabilized") |>
    droplevels()

  # Map scenarios to results (prefer names)
  obs_names <- names(obs_list)
  j_idx <- if (!is.null(obs_names) && all(scen_df$scen_key %in% obs_names)) {
    match(scen_df$scen_key, obs_names)
  } else {
    scen_df$scen_num
  }
  if (anyNA(j_idx)) stop("Scenario key mapping to obs/sr_sim failed.")

  # Determine fixed post-history window from a non-empty obs df
  first_obs <- as.data.frame(obs_list[[ j_idx[1] ]][[1]])
  ny_obs    <- nrow(first_obs)
  year_index <- seq.int(nyh + 1L, ny_obs)  # post-history window
  if (length(year_index) <= 0) stop("Empty post-history window; check nyh/obs length.")

  # Preallocate per-metric matrices (rows = scenarios, cols = iterations)
  nscen <- nrow(scen_df)
  niter <- length(obs_list[[ j_idx[1] ]])
  m_av_harv <- m_av_ret <- m_pRmax50 <- m_pSeq50 <- matrix(NA_real_, nscen, niter)

  # Core loop (paper-faithful)
  for (ii in seq_len(nscen)) {
    j <- j_idx[ii]
    for (k in seq_len(niter)) {
      ob <- as.data.frame(obs_list[[j]][[k]])
      if (!is.data.frame(ob) || !nrow(ob)) next

      # Harvest / Return time-means over window (observed)
      if ("obsHarv" %in% names(ob)) m_av_harv[ii, k] <- mean(ob$obsHarv[year_index], na.rm = TRUE)
      if ("obsRet"  %in% names(ob)) m_av_ret [ii, k] <- mean(ob$obsRet [year_index], na.rm = TRUE)

      # Probabilities relative to Rmax/Seq thresholds (NA-safe denominators)
      sr <- as.data.frame(sr_list[[j]][[k]])
      if (is.data.frame(sr) && all(c("alpha","beta") %in% names(sr))) {
        R_max <- (sr$alpha / sr$beta) * exp(-1)
        S_eq  <- log(sr$alpha) / sr$beta

        rec <- if ("recRec" %in% names(ob)) ob$recRec[year_index] else NULL
        esc <- if ("obsEsc" %in% names(ob)) ob$obsEsc[year_index] else NULL

        if (!is.null(rec)) {
          den <- sum(!is.na(rec))
          if (den > 0) m_pRmax50[ii, k] <- sum(rec > 0.5 * R_max, na.rm = TRUE) / den
        }
        if (!is.null(esc)) {
          den <- sum(!is.na(esc))
          if (den > 0) m_pSeq50[ii, k] <- sum(esc > 0.5 * S_eq, na.rm = TRUE) / den
        }
      }
    }
  }

  # Across-iteration medians (paper)
  safemed <- function(v) stats::median(v, na.rm = TRUE)
  df <- scen_df |>
    dplyr::mutate(
      av_harv       = apply(m_av_harv, 1, safemed) / 1e3,  # thousands
      av_ret        = apply(m_av_ret , 1, safemed) / 1e3,  # thousands
      p_over_Rmax50 = apply(m_pRmax50, 1, safemed),
      p_over_Seq50  = apply(m_pSeq50 , 1, safemed)
    )

  # Pivot only metric columns (avoid mixing types)
  df_metrics <- df |>
    dplyr::select(trends, selectivity, mgmt,
                  av_harv, av_ret, p_over_Rmax50, p_over_Seq50)

  dfp <- df_metrics |>
    tidyr::pivot_longer(
      c(av_harv, av_ret, p_over_Rmax50, p_over_Seq50),
      names_to  = "metric",
      values_to = "median"
    ) |>
    dplyr::mutate(
      metric_label = dplyr::case_when(
        metric == "av_harv"       ~ "Mean harvest\n(thousands)",
        metric == "av_ret"        ~ "Mean run size\n(thousands)",
        metric == "p_over_Rmax50" ~ "Probability\nabove 50% Rmax",
        metric == "p_over_Seq50"  ~ "Probability\nabove 50% S0",
        TRUE ~ metric
      ),
      # lock facet row order to published layout
      metric_label = factor(
        metric_label,
        levels = c("Mean harvest\n(thousands)",
                   "Mean run size\n(thousands)",
                   "Probability\nabove 50% Rmax",
                   "Probability\nabove 50% S0")
      )
    ) |>
    droplevels()

  # Colors and plot (match paper)
  colors <- c(TRM = "darkgray", YPR = "deepskyblue3", DLM = "orange")

  p <- dfp |>
    ggplot2::ggplot(ggplot2::aes(
      x = forcats::fct_inorder(trends),
      y = median,
      color = forcats::fct_inorder(mgmt),
      group = mgmt
    )) +
    ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = colors,
                                 breaks = c("TRM","YPR","DLM"),
                                 drop = FALSE) +
    ggplot2::scale_y_continuous(expand = c(0.07, 0.07), limits = c(0, NA)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "", y = "", color = "Method") +
    # use formula syntax (portable) and place y-strips outside panels
    ggplot2::facet_grid(
      metric_label ~ selectivity,
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10),
      strip.text.y = ggplot2::element_text(size = 10),
      strip.placement = "outside",
      axis.line   = ggplot2::element_line(linewidth = 0.1),
      axis.text   = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title  = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 8)
    )

  outdir <- .ensure_outdir(output_dir, file_basename)

  ggplot2::ggsave(
    file.path(outdir, paste0(file_basename, ".pdf")),
    p, width = width_in, height = height_in, units = "in"
  )

  # Write summarized medians (what the figure shows)
  readr::write_csv(dfp, file.path(outdir, paste0(file_basename, ".csv")))

  invisible(list(plot = p, data = dfp, outdir = outdir))
}
