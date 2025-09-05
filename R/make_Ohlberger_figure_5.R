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

#' @noRd
#' @keywords internal
.make_Ohlberger_figure_5 <- function(data,
                                     output_dir = tempdir(),
                                     file_basename = "Figure5",
                                     width_in = 6.5,
                                     height_in = 3) {
  stopifnot(inherits(data, "ssi_run"))
  .validate_run_params(data$parameters)

  # ---- constants / style ------------------------------------------------------
  colors  <- c(YPR = "deepskyblue3", DLM = "orange")  # (TRM not plotted)
  decades <- c(10L, 20L, 30L, 40L, 50L)
  nyh     <- .get_nyh(data)            # post-history length recorded with run
  obs_list <- .get_obs_list(data)

  # ---- scenarios (standardize) ------------------------------------------------
  scen_df <- .scenarios_from_ssirun(data) |>
    # published Fig 5 operates on base cases (MSY)
    dplyr::filter(.data$factorMSY == "MSY") |>
    droplevels()

  if (!nrow(scen_df)) stop("No scenarios after filtering to factorMSY == 'MSY'.")

  # prefer name-based mapping to results
  obs_names <- names(obs_list)
  j_idx <- if (!is.null(obs_names) && all(scen_df$scen_key %in% obs_names)) {
    match(scen_df$scen_key, obs_names)
  } else {
    scen_df$scen_num
  }
  if (anyNA(j_idx)) stop("Scenario key mapping to obs list failed.")

  # ---- fixed post-history window (same for all scen/iter) ---------------------
  first_obs <- as.data.frame(obs_list[[ j_idx[1] ]][[1]])
  ny_obs    <- nrow(first_obs)
  from_idx  <- nyh + 1L
  if (from_idx > ny_obs) stop("Empty post-history window; check nyh/obs length.")
  # For decade i, window is (nyh+1) : (nyh + i*10)

  # ---- helper: per-iteration future means by decade ---------------------------
  compute_future_means <- function(iter_obs, nyh, decades) {
    niter <- length(iter_obs)
    res <- matrix(NA_real_, nrow = niter, ncol = length(decades),
                  dimnames = list(paste0("iter_", seq_len(niter)), as.character(decades)))
    for (k in seq_len(niter)) {
      df <- iter_obs[[k]]
      if (!is.data.frame(df) || !nrow(df) || is.null(df$obsRet)) next
      n  <- nrow(df)
      from <- nyh + 1L
      if (from > n) next
      for (j in seq_along(decades)) {
        to   <- nyh + decades[j]
        end  <- min(to, n)
        if (end < from) { res[k, j] <- NA_real_; next }
        res[k, j] <- mean(df$obsRet[from:end], na.rm = TRUE)
      }
    }
    res
  }

  # compute per-scenario (niter x ndecades) mean returns
  scen_keys <- scen_df$scen_key
  mean_mats <- setNames(vector("list", length(scen_keys)), scen_keys)
  for (i in seq_along(scen_keys)) {
    key <- scen_keys[i]
    mean_mats[[key]] <- compute_future_means(
      iter_obs = obs_list[[ key ]],
      nyh      = nyh,
      decades  = decades
    )
  }

  # ---- pair YPR/DLM with TRM within (trends, selectivity, factorMSY) ----------
  scen_df$cell <- do.call(paste, c(scenv = scen_df[c("trends","selectivity","factorMSY")], sep = "||"))

  trm_map <- scen_df |>
    dplyr::filter(.data$mgmt == "TRM") |>
    dplyr::select(.data$cell, trm_key = .data$scen_key)

  cand_df <- scen_df |>
    dplyr::filter(.data$mgmt %in% c("YPR", "DLM")) |>
    dplyr::left_join(trm_map, by = "cell")

  if (!nrow(cand_df)) stop("No YPR/DLM scenarios to compare against TRM for Figure 5.")

  # ---- median % differences per decade (vs TRM) --------------------------------
  rows <- lapply(seq_len(nrow(cand_df)), function(i) {
    row  <- cand_df[i, ]
    keym <- row$scen_key
    keyr <- row$trm_key
    if (is.na(keyr)) return(NULL)

    Mm <- mean_mats[[keym]]
    Mr <- mean_mats[[keyr]]
    if (is.null(Mm) || is.null(Mr)) return(NULL)
    if (!identical(colnames(Mm), colnames(Mr))) return(NULL)

    pct <- 100 * (Mm - Mr) / Mr   # matrix [niter x ndecades]
    med <- apply(pct, 2L, stats::median, na.rm = TRUE)

    tibble::tibble(
      trends      = row$trends,
      selectivity = row$selectivity,
      mgmt        = row$mgmt,
      period      = as.integer(names(med)),
      median      = as.numeric(med)
    )
  })

  df_plot <- dplyr::bind_rows(rows)
  if (!nrow(df_plot)) stop("No paired differences computed for Figure 5.")

  # Published filter timing: keep only ASL trends stabilized, drop TRM before plotting
  df_plot <- df_plot |>
    dplyr::filter(.data$trends == "ASL trends stabilized") |>
    dplyr::filter(.data$mgmt %in% c("YPR","DLM")) |>
    dplyr::select(-.data$trends) |>
    dplyr::mutate(mgmt = factor(.data$mgmt, levels = c("YPR","DLM"))) |>
    droplevels()

  # ---- plot -------------------------------------------------------------------
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$period, y = .data$median,
                 color = forcats::fct_inorder(.data$mgmt), group = .data$mgmt)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.1) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(shape = 1, size = 2.5, fill = NA, color = "black") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colors, breaks = c("YPR","DLM"), drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = decades, limits = range(decades), expand = c(0.1, 0.1)) +
    ggplot2::scale_y_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Year", y = "Median difference in run size (%)", color = "Method") +
    # Portable faceting (avoids vars())
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

  message("Figure 5 saved to: ", plot_path)
  message("Data saved to: ", data_path)

  invisible(list(plot = p, data = df_plot, outdir = outdir))

}

