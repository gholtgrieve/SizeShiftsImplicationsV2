#' Generate Ohlberger Figure 6: Δ harvest vs Δ escapement (per iteration)
#'
#' Reproduces the published Figure 6: percent differences (YPR/DLM vs TRM)
#' in **mean harvest** (x) and **mean escapement** (y) for each stochastic
#' iteration under a **large-mesh** fishery with **ASL trends continued**,
#' shown as a scatter with marginal density plots.
#'
#' @section Caption (from the publication):
#' Difference in mean harvest and mean escapement from time-invariant model.
#' Shown is the % difference in harvest and escapement in each stochastic run
#' (circles) for the yield-per-recruit analysis (YPR, blue) and the Dynamic
#' Linear Model (DLM, yellow) compared to the time-invariant Ricker model for a
#' large-mesh fishery assuming continued ASL trends. Marginal plots along the
#' secondary axes show the respective density distributions.
#'
#' @section Data & assumptions:
#' - Input is an **`ssi_run`** produced by [run_scenarios()].
#' - Uses `data$results$obs` for time series and `data$scenarios` for labels.
#' - Window = **post-historical years** (rows `nyh+1 : n`, with `nyh = 50`
#'   by default, clamped to available rows).
#' - We restrict to `selectivity == "large-mesh"`, `trends == "ASL trends continued"`,
#'   and `factorMSY == "MSY"` (base case).
#' - For each YPR/DLM scenario we compute **per-iteration** % differences
#'   against the TRM scenario in the same (selectivity, trends, factorMSY) cell,
#'   pairing the **same iteration index**.
#'
#' @section Files written:
#' - `file.path(output_dir, "Figure6", "Figure6.pdf")`
#' - `file.path(output_dir, "Figure6", "Figure6.csv")`
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for
#' salmon body size declines in fishery management can reduce conservation
#' risks. *Fish and Fisheries*, 26, 113–130. https://doi.org/10.1111/faf.12869
#'
#' @param data `ssi_run` object from [run_scenarios()].
#' @param output_dir Directory to save outputs (default: current working directory).
#'
#' @return Invisible list with `plot` (ggplot/ggExtra object) and `data`
#'   (tidy data frame of per-iteration percent differences).
#' @noRd
#' @keywords internal

.make_Ohlberger_figure_6 <- function(data,
                                     output_dir = tempdir(),
                                     file_basename = "Figure6",
                                     width_in = 4,
                                     height_in = 4) {
  stopifnot(inherits(data, "ssi_run"))
  .validate_run_params(data$parameters)

  # ---- constants / style ------------------------------------------------------
  nyh <- .get_nyh(data)           # historical length (e.g., 50)
  obs_list <- .get_obs_list(data)
  colors <- c(YPR = "deepskyblue3", DLM = "orange")  # TRM not plotted

  # ---- scenarios (standardize + keys) ----------------------------------------
  scen <- .scenarios_from_ssirun(data)

  # attach stable key matching nested results (already present from helper)
  # filter to the captioned setting (paper): large-mesh, continued trends, MSY
  scen <- scen |>
    dplyr::filter(.data$selectivity == "large-mesh",
                  .data$trends      == "ASL trends continued",
                  .data$factorMSY   == "MSY") |>
    droplevels()

  if (!nrow(scen)) {
    stop("No scenarios after filtering to large-mesh, ASL trends continued, factorMSY == 'MSY'.")
  }

  # build cell id for robust TRM pairing (matches paper's intent)
  scen$cell <- do.call(paste, c(scen[c("trends","selectivity","factorMSY")], sep = "||"))

  trm_map <- scen |>
    dplyr::filter(.data$mgmt == "TRM") |>
    dplyr::select(.data$cell, trm_key = .data$scen_key)

  cand <- scen |>
    dplyr::filter(.data$mgmt %in% c("YPR","DLM")) |>
    dplyr::left_join(trm_map, by = "cell")

  if (!nrow(cand)) stop("No YPR/DLM scenarios to compare against TRM in the filtered set.")

  # ---- helper: per-iteration mean over fixed post-history window -------------
  mean_post_hist <- function(iter_obs, nyh, col) {
    vapply(iter_obs, function(df) {
      if (!is.data.frame(df) || !nrow(df) || is.null(df[[col]])) return(NA_real_)
      n    <- nrow(df)
      from <- nyh + 1L
      to   <- nyh + 50L
      end  <- min(to, n)
      if (from > n || end < from) return(NA_real_)
      mean(df[[col]][from:end], na.rm = TRUE)
    }, numeric(1))
  }

  # ---- compute per-iteration % diffs vs TRM for selected cells ----------------
  rows <- lapply(seq_len(nrow(cand)), function(i) {
    r <- cand[i, ]
    key_m <- r$scen_key
    key_r <- r$trm_key
    if (is.na(key_r)) return(NULL)

    it_m <- obs_list[[key_m]]
    it_r <- obs_list[[key_r]]
    if (is.null(it_m) || is.null(it_r)) return(NULL)

    # per-iteration means over post-history window
    mh_m <- mean_post_hist(it_m, nyh, "obsHarv")
    me_m <- mean_post_hist(it_m, nyh, "obsEsc")
    mh_r <- mean_post_hist(it_r, nyh, "obsHarv")
    me_r <- mean_post_hist(it_r, nyh, "obsEsc")

    # align lengths defensively
    n <- min(length(mh_m), length(mh_r), length(me_m), length(me_r))
    if (n == 0) return(NULL)

    d_h <- 100 * (mh_m[seq_len(n)] - mh_r[seq_len(n)]) / mh_r[seq_len(n)]
    d_e <- 100 * (me_m[seq_len(n)] - me_r[seq_len(n)]) / me_r[seq_len(n)]

    tibble::tibble(
      mgmt        = r$mgmt,
      selectivity = r$selectivity,
      trends      = r$trends,
      factorMSY   = r$factorMSY,
      iter        = paste0("iter_", seq_len(n)),
      diff_harvest_pct    = as.numeric(d_h),
      diff_escapement_pct = as.numeric(d_e)
    )
  })

  df <- dplyr::bind_rows(rows)
  if (!nrow(df)) stop("No paired iteration-wise differences computed for Figure 6.")

  # Keep only the two alternatives for this figure and lock legend order
  df <- df |>
    dplyr::filter(.data$mgmt %in% c("YPR","DLM")) |>
    dplyr::mutate(mgmt = factor(.data$mgmt, levels = c("YPR","DLM"))) |>
    droplevels()

  # ---- plot (scatter with marginal densities) ---------------------------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$diff_harvest_pct,
                 y = .data$diff_escapement_pct,
                 color = .data$mgmt)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.1, color = "black") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.1, color = "black") +
    ggplot2::geom_point(size = 0.75, shape = 16) +
    ggplot2::scale_colour_manual(values = colors, breaks = c("YPR","DLM"), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(-80, 190), breaks = seq(-50, 150, 50)) +
    ggplot2::scale_y_continuous(limits = c(-80, 190), breaks = seq(-50, 150, 50)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Difference in mean harvest (%)",
                  y = "Difference in mean escapement (%)",
                  color = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.position = c(0.1, 0.95),
      legend.key.size = grid::unit(0.25, "cm"),
      legend.text = ggplot2::element_text(size = 8)
    )

  # Add marginals (ggExtra)
  g <- ggExtra::ggMarginal(p, groupColour = TRUE, groupFill = TRUE,
                           margins = "both", size = 5, type = "density")

  # ---- save & return ----------------------------------------------------------
  outdir <- .ensure_outdir(output_dir, file_basename)
  plot_path <- file.path(outdir, paste0(file_basename, ".pdf"))
  data_path <- file.path(outdir, paste0(file_basename, ".csv"))

  ggplot2::ggsave(plot_path, g, width = width_in, height = height_in, units = "in")
  readr::write_csv(df, data_path)

  message("Figure 6 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  invisible(list(plot = g, data = df, outdir = outdir))
}

