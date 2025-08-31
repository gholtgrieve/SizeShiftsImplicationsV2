#' Generate Ohlberger Figure 7: Performance for liberal vs precautionary goals
#'
#' Reproduces the published Figure 7: median values of four performance metrics
#' (rows) for each estimation method (colours) using management targets of
#' either **0.75 × S[MSY]** (“liberal”) or **1.5 × S[MSY]** (“precautionary”),
#' assuming **unselective** fishing (columns = factorMSY).
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

.make_Ohlberger_figure_7 <- function(data, output_dir = ".") {
  stopifnot(inherits(data, "ssi_run"))

  # --- tunables / style (edit here if needed) ----------------------------------
  summary_fun <- stats::median   # statistic across iterations
  nyh_default <- 50L             # last-50-year window
  colors <- c(TRM = "darkgray", YPR = "deepskyblue3", DLM = "orange")
  file_basename <- "Figure7"

  # --- scenario frame ----------------------------------------------------------
  scen <- tibble::as_tibble(data$scenarios)
  scen <- .standardize_scenario_labels(scen)

  # stable key to reach nested lists
  if ("scen_num" %in% names(scen)) {
    scen$scen_key <- paste0("scen_", scen$scen_num)
  } else {
    nm <- names(data$results$obs)
    scen$scen_key <- if (!is.null(nm) && length(nm) == nrow(scen)) nm else paste0("scenario_", seq_len(nrow(scen)))
  }

  # filter to unselective, factorMSY in {liberal, precautionary}, and drop AL-only trend
  scen <- scen |>
    dplyr::filter(.data$selectivity == "unselective",
                  .data$factorMSY %in% c("liberal","precautionary"),
                  .data$trends %in% c("no trends", "ASL trends stabilized", "ASL trends continued"))

   if (!nrow(scen)) stop("No scenarios after filtering to unselective & (liberal/precautionary).")

  # ensure factor order for facets and legend stability
  scen$factorMSY <- factor(scen$factorMSY, levels = c("liberal","precautionary"))
  scen$mgmt      <- factor(scen$mgmt, levels = c("TRM","YPR","DLM"))
  scen$trends <- factor(scen$trends,
                        levels = c("no trends", "ASL trends stabilized", "ASL trends continued"))

  # --- access results ----------------------------------------------------------
  obs_list   <- data$results$obs
  sr_simlist <- data$results$sr_sim
  if (is.null(obs_list) || !length(obs_list))   stop("Missing run$results$obs.")
  if (is.null(sr_simlist) || !length(sr_simlist)) stop("Missing run$results$sr_sim.")

  # helper: last-50-year mean of a column per iteration
  mean_last <- function(iter_list, col, nyh = nyh_default) {
    vapply(iter_list, function(df) {
      if (!is.data.frame(df) || !nrow(df) || is.null(df[[col]])) return(NA_real_)
      n <- nrow(df); start <- max(1L, n - nyh + 1L)
      mean(df[[col]][start:n], na.rm = TRUE)
    }, numeric(1))
  }

  # helper: last-50-year probability of condition per iteration
  prob_last <- function(iter_list, predicate, nyh = nyh_default) {
    vapply(iter_list, function(df) {
      if (!is.data.frame(df) || !nrow(df)) return(NA_real_)
      n <- nrow(df); start <- max(1L, n - nyh + 1L)
      v <- df[start:n, , drop = FALSE]
      x <- predicate(v)
      if (!length(x)) return(NA_real_)
      mean(x, na.rm = TRUE)
    }, numeric(1))
  }

  # build per-scenario summaries (median across iterations)
  rows <- lapply(seq_len(nrow(scen)), function(i) {
    r <- scen[i, ]
    key <- r$scen_key

    it_obs <- obs_list[[key]]
    it_fit <- sr_simlist[[key]]
    if (is.null(it_obs) || is.null(it_fit)) return(NULL)

    # α, β per iteration (from sr_sim fit)
    get_param <- function(it_fit, name) {
      vapply(it_fit, function(x) {
        if (is.null(x) || !is.data.frame(x) || is.null(x[[name]])) return(NA_real_)
        as.numeric(x[[name]][1])
      }, numeric(1))
    }
    alpha <- get_param(it_fit, "alpha")
    beta  <- get_param(it_fit, "beta")

    # reference points per iteration
    Rmax <- (alpha / beta) * exp(-1)
    S0   <- log(alpha) / beta

    # time-window means per iteration
    mean_h <- mean_last(it_obs, "obsHarv", nyh_default)
    mean_r <- mean_last(it_obs, "obsRet",  nyh_default)

    # probabilities per iteration
    prob_R <- prob_last(it_obs, function(df) df$recRec > 0.5 * Rmax, nyh_default)
    prob_S <- prob_last(it_obs, function(df) df$obsEsc > 0.5 * S0,   nyh_default)

    tibble::tibble(
      trends      = r$trends,
      factorMSY   = r$factorMSY,
      mgmt        = r$mgmt,
      selectivity = r$selectivity,
      mean_50yr_harvest    = summary_fun(mean_h, na.rm = TRUE),
      mean_50yr_return     = summary_fun(mean_r, na.rm = TRUE),
      p_over_Rmax50        = summary_fun(prob_R, na.rm = TRUE),
      p_over_Seq50         = summary_fun(prob_S, na.rm = TRUE)
    )
  })

  df <- dplyr::bind_rows(rows)
  if (!nrow(df)) stop("No summaries computed for Figure 7.")

  # tidy for plotting
  df_tidy <- df |>
    tidyr::pivot_longer(
      cols = c(.data$mean_50yr_harvest, .data$mean_50yr_return,
               .data$p_over_Rmax50, .data$p_over_Seq50),
      names_to = "metric", values_to = "value"
    ) |>
    dplyr::mutate(
      metric_label = dplyr::case_when(
        .data$metric == "mean_50yr_harvest" ~ "Mean harvest\n(thousands)",
        .data$metric == "mean_50yr_return"  ~ "Mean run size\n(thousands)",
        .data$metric == "p_over_Rmax50"     ~ "Probability\nabove 50% Rmax",
        .data$metric == "p_over_Seq50"      ~ "Probability\nabove 50% S0"
      ),
      # thousands for the two abundance metrics
      value_plot = dplyr::case_when(
        .data$metric %in% c("mean_50yr_harvest","mean_50yr_return") ~ .data$value / 1000,
        TRUE ~ .data$value
      ),
      metric_label = factor(
        .data$metric_label,
        levels = c("Mean harvest\n(thousands)",
                   "Mean run size\n(thousands)",
                   "Probability\nabove 50% Rmax",
                   "Probability\nabove 50% S0")
      )
    )

  # plot
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
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_y_continuous(expand = c(0.07, 0.07), limits = c(0, NA)) +
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

  g <- p +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$metric_label),
      cols = ggplot2::vars(.data$factorMSY),
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::theme(strip.placement = "outside")

  # save
  outdir <- file.path(output_dir, "Figure7")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plot_path <- file.path(outdir, paste0(file_basename, ".pdf"))
  data_path <- file.path(outdir, paste0(file_basename, ".csv"))
  ggplot2::ggsave(plot_path, g, width = 4.5, height = 7, units = "in")
  readr::write_csv(df_tidy, data_path)

  message("Figure 7 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  return(invisible(list(plot = g, data = df_tidy)))
}
