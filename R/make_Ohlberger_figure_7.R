#' Generate Ohlberger Figure 7: Performance for liberal vs precautionary goals
#'
#' Reproduces the published Figure 7: median values of four performance metrics
#' (rows) for each estimation method (colors) using management targets of
#' either 0.75 × S_MSY (“liberal”) or 1.5 × S_MSY (“precautionary”),
#' assuming unselective fishing (columns = factorMSY).
#'
#' @section Caption (from the publication):
#' Performance metrics for liberal and precautionary management strategies.
#' Shown are median values of four performance metrics (rows) for each
#' estimation method (colours) using management targets of either 0.75 S_MSY
#' (left, liberal) or 1.5 S_MSY (right, precautionary), assuming unselective
#' fishing. Details on the performance metrics, the estimation methods, and
#' time period are provided in the caption to Figure 3.
#'
#' @section Data & assumptions:
#' - Input is an **`ssi_run`** from [run_scenarios()].
#' - Uses `data$results$obs` (time series) and `data$results$sr_sim` (α, β)
#'   and the scenario table `data$scenarios`.
#' - Window = **post-historical last 50 years** per iteration (clamped to
#'   available rows).
#' - Metrics per iteration:
#'   * mean harvest, mean return (individuals);
#'   * P(recruitment > 0.5·Rmax), with Rmax = (α/β)·exp(-1);
#'   * P(escapement > 0.5·S0), with S0 = log(α)/β.
#' - Summary across iterations: **median** (modifiable inside function).
#'
#' @section Files written:
#' - `file.path(output_dir, "Figure7", "Figure7.pdf")`
#' - `file.path(output_dir, "Figure7", "Figure7.csv")`
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for
#' salmon body size declines in fishery management can reduce conservation
#' risks. *Fish and Fisheries*, 26, 113–130. <doi:10.1111/faf.12869>.
#'
#' @param data `ssi_run` object.
#' @param output_dir Directory to save outputs (default: current working directory).
#'
#' @return Invisible list with `plot` (ggplot) and `data` (tidy data.frame).
#' @noRd
#' @keywords internal

#' @noRd
#' @keywords internal
.make_Ohlberger_figure_7 <- function(data,
                                     output_dir = tempdir(),
                                     file_basename = "Figure7",
                                     width_in = 4.5,
                                     height_in = 7) {
  stopifnot(inherits(data, "ssi_run"))
  .validate_run_params(data$parameters)

  # ---- access recorded params / results ---------------------------------------
  nyh       <- .get_nyh(data)
  obs_list  <- .get_obs_list(data)
  sr_list   <- .get_sr_params_list(data)
  colors    <- c(TRM = "darkgray", YPR = "deepskyblue3", DLM = "orange")

  # ---- scenarios: standardize, filter per paper -------------------------------
  scen <- .scenarios_from_ssirun(data) |>
    dplyr::filter(.data$selectivity == "unselective") |>
    dplyr::filter(.data$factorMSY %in% c("liberal","precautionary")) |>
    # drop the AL-only trend (matches paper’s `filter(trends!="age-length trends")`)
    dplyr::filter(.data$trends %in% c("no trends","ASL trends stabilized","ASL trends continued")) |>
    dplyr::mutate(
      factorMSY = factor(.data$factorMSY, levels = c("liberal","precautionary")),
      mgmt      = factor(.data$mgmt, levels = c("TRM","YPR","DLM")),
      trends    = factor(.data$trends,
                         levels = c("no trends","ASL trends stabilized","ASL trends continued"))
    ) |>
    droplevels()

  if (!nrow(scen)) stop("No scenarios after filtering to unselective & (liberal/precautionary).")

  # prefer name-based mapping (already provided by .scenarios_from_ssirun)
  obs_names <- names(obs_list)
  if (!is.null(obs_names) && !all(scen$scen_key %in% obs_names)) {
    stop("Scenario key mapping to obs list failed.")
  }

  # ---- fixed post-history window: (nyh+1):ny_obs ------------------------------
  # derive ny_obs from a representative obs df
  first_obs <- as.data.frame(obs_list[[ scen$scen_key[1] ]][[1]])
  ny_obs    <- nrow(first_obs)
  from_idx  <- nyh + 1L
  if (from_idx > ny_obs) stop("Empty post-history window; check nyh/obs length.")
  # window = (nyh+1):ny_obs; in Ohlberger config that’s exactly 50 years.

  # helpers that use the post-history window explicitly
  mean_post_hist <- function(iter_list, col) {
    vapply(iter_list, function(df) {
      if (!is.data.frame(df) || !nrow(df) || is.null(df[[col]])) return(NA_real_)
      n <- nrow(df); end <- ny_obs
      if (from_idx > n || end < from_idx) return(NA_real_)
      mean(df[[col]][from_idx:end], na.rm = TRUE)
    }, numeric(1))
  }
  prob_post_hist <- function(iter_list, pred_fun) {
    vapply(iter_list, function(df) {
      if (!is.data.frame(df) || !nrow(df)) return(NA_real_)
      n <- nrow(df); end <- ny_obs
      if (from_idx > n || end < from_idx) return(NA_real_)
      v <- df[from_idx:end, , drop = FALSE]
      x <- pred_fun(v)
      if (!length(x)) return(NA_real_)
      mean(x, na.rm = TRUE)
    }, numeric(1))
  }

  # ---- build per-scenario summaries (median across iterations) ----------------
  rows <- lapply(seq_len(nrow(scen)), function(i) {
    r   <- scen[i, ]
    key <- r$scen_key

    it_obs <- obs_list[[ key ]]
    it_fit <- sr_list [[ key ]]
    if (is.null(it_obs) || is.null(it_fit)) return(NULL)

    # α, β per iteration from sr_sim
    get_param <- function(it_fit, name) {
      vapply(it_fit, function(x) {
        if (!is.data.frame(x) || is.null(x[[name]])) return(NA_real_)
        as.numeric(x[[name]][1])
      }, numeric(1))
    }
    alpha <- get_param(it_fit, "alpha")
    beta  <- get_param(it_fit, "beta")

    # reference points per iteration
    Rmax <- (alpha / beta) * exp(-1)
    S0   <- log(alpha) / beta

    # time-window means and probabilities per iteration (post-history)
    mean_h <- mean_post_hist(it_obs, "obsHarv")
    mean_r <- mean_post_hist(it_obs, "obsRet")
    prob_R <- prob_post_hist(it_obs, function(df) df$recRec > 0.5 * Rmax)
    prob_S <- prob_post_hist(it_obs, function(df) df$obsEsc > 0.5 * S0)

    tibble::tibble(
      trends      = r$trends,
      factorMSY   = r$factorMSY,
      mgmt        = r$mgmt,
      selectivity = r$selectivity,
      mean_h      = stats::median(mean_h, na.rm = TRUE),
      mean_r      = stats::median(mean_r, na.rm = TRUE),
      p_over_Rmax50 = stats::median(prob_R, na.rm = TRUE),
      p_over_Seq50  = stats::median(prob_S, na.rm = TRUE)
    )
  })

  df <- dplyr::bind_rows(rows)
  if (!nrow(df)) stop("No summaries computed for Figure 7.")

  # ---- tidy for plotting (thousands scaling for harvest/return) ---------------
  df_tidy <- df |>
    tidyr::pivot_longer(
      cols = c(.data$mean_h, .data$mean_r, .data$p_over_Rmax50, .data$p_over_Seq50),
      names_to = "metric", values_to = "value"
    ) |>
    dplyr::mutate(
      metric_label = dplyr::case_when(
        .data$metric == "mean_h"        ~ "Mean harvest\n(thousands)",
        .data$metric == "mean_r"        ~ "Mean run size\n(thousands)",
        .data$metric == "p_over_Rmax50" ~ "Probability\nabove 50% Rmax",
        .data$metric == "p_over_Seq50"  ~ "Probability\nabove 50% S0"
      ),
      value_plot = dplyr::case_when(
        .data$metric %in% c("mean_h","mean_r") ~ .data$value / 1000,
        TRUE ~ .data$value
      ),
      metric_label = factor(
        .data$metric_label,
        levels = c("Mean harvest\n(thousands)",
                   "Mean run size\n(thousands)",
                   "Probability\nabove 50% Rmax",
                   "Probability\nabove 50% S0")
      )
    ) |>
    droplevels()

  # ---- plot -------------------------------------------------------------------
  p <- ggplot2::ggplot(
    df_tidy,
    ggplot2::aes(x = forcats::fct_inorder(.data$trends),
                 y = .data$value_plot,
                 color = forcats::fct_inorder(.data$mgmt),
                 group = .data$mgmt)
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colors,
                                 breaks = c("TRM","YPR","DLM"),
                                 drop   = FALSE) +
    ggplot2::scale_y_continuous(expand = c(0.07, 0.07), limits = c(0, NA)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL, color = "Method") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 10),
      strip.text.y = ggplot2::element_text(size = 10),
      axis.line     = ggplot2::element_line(linewidth = 0.1),
      axis.text     = ggplot2::element_text(size = 10),
      axis.text.x   = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title    = ggplot2::element_text(size = 12),
      panel.border  = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title    = ggplot2::element_text(size = 10),
      legend.text     = ggplot2::element_text(size = 8)
    )

  g <- p +
    ggplot2::facet_grid(
      metric_label ~ factorMSY,
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::theme(strip.placement = "outside")

  # ---- save & announce ---------------------------------------------------------
  outdir    <- .ensure_outdir(output_dir, file_basename)
  plot_path <- file.path(outdir, paste0(file_basename, ".pdf"))
  data_path <- file.path(outdir, paste0(file_basename, ".csv"))

  ggplot2::ggsave(plot_path, g, width = width_in, height = height_in, units = "in")
  readr::write_csv(df_tidy, data_path)

  message("Figure 7 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  invisible(list(plot = g, data = df_tidy, outdir = outdir))
}
