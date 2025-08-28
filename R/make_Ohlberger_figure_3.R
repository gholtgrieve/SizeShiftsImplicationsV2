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
#'   `standardize_scenario_labels()`.
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

.make_Ohlberger_figure_3 <- function(data, output_dir = ".") {
  stopifnot(inherits(data, "ssi_run"))

  colors <- c("darkgray", "deepskyblue3", "orange")  # TRM, YPR, DLM
  file_basename <- "Figure3"

  # Scenario metadata (standardized) and nested results
  scen_df <- .scenarios_from_ssirun(data)           # adds scen, scen_num, scen_key
  obs_list <- .get_obs_list(data)                   # data$results$obs
  sr_list  <- .get_sr_params_list(data)             # data$results$sr_sim

  # Keep only the three trends used in the paper AND factorMSY == 1 ("MSY"),
  # then lock factor orders for stable facets/legend.
  scen_df <- scen_df |>
    dplyr::filter(
      .data$trends %in% c("no trends", "ASL trends stabilized", "ASL trends continued")
    ) |>
    # robust filter for factorMSY == 1 regardless of encoding ("MSY", 1, "1")
    dplyr::mutate(factorMSY_chr = as.character(.data$factorMSY)) |>
    dplyr::filter(.data$factorMSY_chr %in% c("MSY", "1", "1.0")) |>
    dplyr::select(-.data$factorMSY_chr) |>
    dplyr::mutate(
      trends = forcats::fct_relevel(
        .data$trends,
        "no trends", "ASL trends stabilized", "ASL trends continued"
      ),
      selectivity = forcats::fct_relevel(
        .data$selectivity,
        "small-mesh", "unselective", "large-mesh"
      ),
      mgmt = factor(.data$mgmt, levels = c("TRM", "YPR", "DLM"))
    )



  # Ensure we have matching keys in results
  missing_keys <- setdiff(scen_df$scen_key, names(obs_list))
  if (length(missing_keys)) {
    stop("Missing obs_list entries for scenario keys: ", paste(missing_keys, collapse = ", "))
  }
  missing_keys_sr <- setdiff(scen_df$scen_key, names(sr_list))
  if (length(missing_keys_sr)) {
    stop("Missing sr_sim entries for scenario keys: ", paste(missing_keys_sr, collapse = ", "))
  }

  # Window = last nyh rows per iteration
  nyh <- .get_nyh(data, 50L)

  scen_keys <- scen_df$scen_key
  nscen <- length(scen_keys)
  niter <- length(obs_list[[scen_keys[1L]]])

  av_harv_mat <- matrix(NA_real_, nrow = nscen, ncol = niter)
  av_ret_mat  <- matrix(NA_real_, nrow = nscen, ncol = niter)
  p_over_Rmax <- matrix(NA_real_, nrow = nscen, ncol = niter)
  p_over_S0   <- matrix(NA_real_, nrow = nscen, ncol = niter)

  for (j in seq_len(nscen)) {
    key <- scen_keys[j]
    iter_obs <- obs_list[[key]]
    iter_sr  <- sr_list[[key]]

    for (k in seq_len(niter)) {
      obs <- iter_obs[[k]]
      srk <- iter_sr[[k]]

      if (!is.data.frame(obs) || !nrow(obs) || !is.data.frame(srk) || !nrow(srk)) next
      if (!all(c("alpha", "beta") %in% names(srk))) next

      nrows <- nrow(obs)
      win_n <- min(nyh, nrows)
      rows  <- seq.int(nrows - win_n + 1L, nrows)

      esc <- obs$obsEsc[rows]
      har <- obs$obsHarv[rows]
      ret <- obs$obsRet[rows]
      rec <- obs$recRec[rows]

      alpha <- as.numeric(srk$alpha[1L])
      beta  <- as.numeric(srk$beta[1L])
      if (is.na(alpha) || is.na(beta) || beta <= 0) next

      Rmax <- (alpha / beta) * exp(-1)
      S0   <- log(alpha) / beta

      av_harv_mat[j, k] <- mean(har, na.rm = TRUE)
      av_ret_mat[j, k]  <- mean(ret, na.rm = TRUE)
      p_over_Rmax[j, k] <- mean(rec > 0.5 * Rmax, na.rm = TRUE)
      p_over_S0[j, k]   <- mean(esc > 0.5 * S0,   na.rm = TRUE)
    }
  }

  # Scenario-level median across iterations, then tidy
  df <- scen_df |>
    dplyr::select(.data$trends, .data$selectivity, .data$mgmt)

  df$av_harv_thousands <- apply(av_harv_mat, 1L, function(x) stats::median(x, na.rm = TRUE) / 1e3)
  df$av_ret_thousands  <- apply(av_ret_mat,  1L, function(x) stats::median(x, na.rm = TRUE) / 1e3)
  df$p_over_Rmax50     <- apply(p_over_Rmax, 1L, function(x) stats::median(x, na.rm = TRUE))
  df$p_over_S0_50      <- apply(p_over_S0,   1L, function(x) stats::median(x, na.rm = TRUE))

  df_long <- df |>
    tidyr::pivot_longer(
      cols = c(.data$av_harv_thousands, .data$av_ret_thousands,
               .data$p_over_Rmax50, .data$p_over_S0_50),
      names_to = "metric",
      values_to = "median"
    ) |>
    dplyr::mutate(
      metric_label = dplyr::case_when(
        .data$metric == "av_harv_thousands" ~ "Mean harvest\n(thousands)",
        .data$metric == "av_ret_thousands"  ~ "Mean run size\n(thousands)",
        .data$metric == "p_over_Rmax50"     ~ "Probability\nabove 50% Rmax",
        .data$metric == "p_over_S0_50"      ~ "Probability\nabove 50% S0"
      ),
      metric_label = factor(
        .data$metric_label,
        levels = c("Mean harvest\n(thousands)",
                   "Mean run size\n(thousands)",
                   "Probability\nabove 50% Rmax",
                   "Probability\nabove 50% S0")
      )
    )

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = forcats::fct_inorder(.data$trends),
                 y = .data$median,
                 color = forcats::fct_inorder(.data$mgmt),
                 group = .data$mgmt)
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_y_continuous(expand = c(0.07, 0.07), limits = c(0, NA)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "", y = "", color = "Method") +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$metric_label),
      cols = ggplot2::vars(.data$selectivity),
      scales = "free_y",
      switch = "y",
      drop = FALSE
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

  outdir <- .ensure_outdir(output_dir, file_basename)
  plot_path <- file.path(outdir, paste0(file_basename, ".pdf"))
  data_path <- file.path(outdir, paste0(file_basename, ".csv"))
  ggplot2::ggsave(plot_path, p, width = 5.5, height = 7, units = "in")
  readr::write_csv(df_long, data_path)

  message("Figure 3 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  return(invisible(list(plot = p, data = df_long)))
}
