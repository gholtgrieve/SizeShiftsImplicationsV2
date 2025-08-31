#' Generate Ohlberger Figure 5: Decadal % difference in run size vs TRM
#'
#' Reproduces the published Figure 5: **median percent difference** in average
#' run size for each decade after the historical period, comparing YPR and DLM
#' to the time-invariant Ricker model (TRM). The x-axis shows the decade
#' endpoints (10, 20, 30, 40, 50), facets are selectivity regimes, and we
#' condition on **ASL trends continued**.
#'
#' @section Caption (from the publication):
#' Median % difference in run size from the time-invariant model across years.
#' Shown is the median difference for the yield-per-recruit analysis (YPR, blue)
#' and the Dynamic Linear Model (DLM, yellow) compared to the time-invariant
#' Ricker model in each decade starting after the historical period (e.g., year
#' 10 is the difference for future years 1–10) assuming continued ASL trends.
#'
#' @section Data & assumptions:
#' - Input is an **`ssi_run`** from `run_scenarios()`.
#' - Uses `data$results$obs` for time series and `data$scenarios` for labels.
#' - We restrict to **trends = "ASL trends continued"** and **factorMSY = "MSY"**
#'   (base-case), then pair each YPR/DLM scenario with its TRM reference within
#'   the same (selectivity, factorMSY) cell.
#' - For each iteration and each decade endpoint `t ∈ {10,20,30,40,50}`, we
#'   average `obsRet` over years `(nyh+1) : (nyh+t)` (clamped to available rows),
#'   compute `% diff` vs TRM, and then take the **median across iterations**.
#'
#' @section Files written:
#' - `file.path(output_dir, "Figure5", "Figure5.pdf")`
#' - `file.path(output_dir, "Figure5", "Figure5.csv")`
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for
#' salmon body size declines in fishery management can reduce conservation
#' risks. *Fish and Fisheries*, 26, 113–130. https://doi.org/10.1111/faf.12869
#'
#' @param data `ssi_run` object from `run_scenarios()`.
#' @param output_dir Directory to save outputs (default: current working directory).
#'
#' @return Invisible list with `plot` (ggplot) and `data` (tidy plotting data).
#' @noRd
#' @keywords internal

.make_Ohlberger_figure_5 <- function(data, output_dir = ".") {
  stopifnot(inherits(data, "ssi_run"))

  # ---- constants / style (match other figs) -----------------------------------
  colors <- c(TRM = "darkgray", YPR = "deepskyblue3", DLM = "orange")
  decades <- c(10L, 20L, 30L, 40L, 50L)
  nyh_default <- 50L
  file_basename <- "Figure5"

  # ---- scenarios (standardize, filter to continued trends & base-case MSY) ----
  scen_df <- tibble::as_tibble(data$scenarios)
  scen_df <- .standardize_scenario_labels(scen_df)

  # attach stable key matching results lists
  if ("scen_num" %in% names(scen_df)) {
    scen_df$scen_key <- paste0("scen_", scen_df$scen_num)
  } else {
    obs_names <- names(data$results$obs)
    scen_df$scen_key <- if (!is.null(obs_names) && length(obs_names) == nrow(scen_df)) obs_names else paste0("scenario_", seq_len(nrow(scen_df)))
  }

  # filter: ASL trends continued, base-case MSY
  scen_df <- scen_df |>
    dplyr::filter(.data$trends == "ASL trends continued") |>
    dplyr::filter(.data$factorMSY == "MSY")

  if (!nrow(scen_df)) stop("No scenarios after filtering to 'ASL trends continued' and factorMSY == 'MSY'.")

  # ---- pull obs and helper to compute future means by decade -------------------
  obs_list <- data$results$obs
  if (is.null(obs_list) || !length(obs_list)) {
    stop("No observation lists found in `ssi_run` (run$results$obs).")
  }

  compute_future_means <- function(iter_obs, nyh = nyh_default, decades = decades) {
    niter <- length(iter_obs)
    res <- matrix(NA_real_, nrow = niter, ncol = length(decades),
                  dimnames = list(paste0("iter_", seq_len(niter)), as.character(decades)))
    for (k in seq_len(niter)) {
      df <- iter_obs[[k]]
      if (!is.data.frame(df) || !nrow(df) || is.null(df$obsRet)) next
      n <- nrow(df)
      start <- min(nyh + 1L, n)
      if (start > n) next
      for (j in seq_along(decades)) {
        end <- min(nyh + decades[j], n)
        if (end < start) { res[k, j] <- NA_real_; next }
        res[k, j] <- mean(df$obsRet[start:end], na.rm = TRUE)
      }
    }
    res
  }

  # compute per-scenario (niter x ndecades) means
  scen_keys <- scen_df$scen_key
  mean_mats <- setNames(vector("list", length(scen_keys)), scen_keys)
  for (i in seq_along(scen_keys)) {
    key <- scen_keys[i]
    mean_mats[[key]] <- compute_future_means(obs_list[[key]], nyh_default, decades)
  }

  # ---- pair YPR/DLM with TRM within each (selectivity, factorMSY) cell ---------
  scen_df$cell <- do.call(paste, c(scenv = scen_df[c("selectivity", "factorMSY")], sep = "||"))

  trm_map <- scen_df |>
    dplyr::filter(.data$mgmt == "TRM") |>
    dplyr::select(.data$cell, trm_key = .data$scen_key)

  cand_df <- scen_df |>
    dplyr::filter(.data$mgmt %in% c("YPR", "DLM")) |>
    dplyr::left_join(trm_map, by = "cell")

  if (!nrow(cand_df)) stop("No YPR/DLM scenarios to compare against TRM for Figure 5.")

  # build tidy output of median % diffs per decade
  rows <- lapply(seq_len(nrow(cand_df)), function(i) {
    row <- cand_df[i, ]
    key_m <- row$scen_key
    key_r <- row$trm_key
    if (is.na(key_r)) return(NULL)

    Mm <- mean_mats[[key_m]]
    Mr <- mean_mats[[key_r]]
    if (is.null(Mm) || is.null(Mr)) return(NULL)
    if (!all(colnames(Mm) == colnames(Mr))) return(NULL)

    pct <- 100 * (Mm - Mr) / Mr   # matrix [niter x ndecades]
    med <- apply(pct, 2L, stats::median, na.rm = TRUE)

    tibble::tibble(
      selectivity = row$selectivity,
      mgmt        = row$mgmt,
      period      = as.integer(names(med)),
      median      = as.numeric(med)
    )
  })

  df_plot <- dplyr::bind_rows(rows)
  if (!nrow(df_plot)) stop("No paired differences computed for Figure 5.")

  df_plot <- df_plot |>
    .standardize_scenario_labels() |>
    dplyr::mutate(
      # keep only the two alternatives for this figure
      mgmt = factor(.data$mgmt, levels = c("YPR", "DLM"))
    ) |>
    droplevels()

  # ---- plot --------------------------------------------------------------------
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$period, y = .data$median,
                 color = forcats::fct_inorder(.data$mgmt), group = .data$mgmt)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.1) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = c(YPR = "deepskyblue3", DLM = "orange")) +
    ggplot2::scale_x_continuous(breaks = decades, limits = range(decades), expand = c(0.1, 0.1)) +
    ggplot2::scale_y_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Year", y = "Median difference in run size (%)", color = "Method") +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$selectivity)) +
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
      legend.text = ggplot2::element_text(size = 8),
      legend.background = ggplot2::element_blank()
    )

  # ---- save & return -----------------------------------------------------------
  outdir <- file.path(output_dir, "Figure5")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plot_path <- file.path(outdir, paste0(file_basename, ".pdf"))
  data_path <- file.path(outdir, paste0(file_basename, ".csv"))

  ggplot2::ggsave(plot_path, p, width = 6.5, height = 3, units = "in")
  readr::write_csv(df_plot, data_path)

  message("Figure 5 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  return(invisible(list(plot = p, data = df_plot)))
}
