#' Generate Figure C: Time series of escapement, harvest, and probability of zero harvest
#'   by trend, harvest strategy, and management
#'
#' Accepts:
#' - an `ssi_run` object from `run_scenarios()` (recommended), or
#' - legacy flat outputs compatible with older code.
#'
#' Internally calls `summarize_by_year()` (which defaults to the last 50 years
#' for `ssi_run` inputs). Produces ribbons for central intervals and overlays a
#' smoothed (5-year rolling mean) probability of zero harvest on the secondary axis.
#'
#' @param outputs `ssi_run` (preferred) or legacy outputs list.
#' @param selectivity_filter Character; filter for selectivity type (default: "unselective").
#' @param summary_stats Named vector mapping plot columns to summary statistic names
#'   (default: c(mean = "Mean", ymin = "25%", ymax = "75%")).
#' @param colors Character vector of colors (default: c("darkgray", "deepskyblue3", "orange")).
#' @param output_dir Directory to save plots and data (default: current working directory).
#' @param file_basename Base filename for saved files (default: "FigureC").
#'
#' @return A list with:
#'   - `data`: tidy data.frame used for plotting
#'   - `plots`: named list of ggplot objects by management strategy
#' @export
make_figure_C <- function(
    outputs,
    selectivity_filter = "unselective",
    summary_stats = c(mean = "Mean", ymin = "25%", ymax = "75%"),
    colors = c("darkgray", "deepskyblue3", "orange"),
    output_dir = ".",
    file_basename = "FigureC"
) {
  # 1) Summaries (uses last 50y by default for ssi_run)
  summary_list <- summarize_by_year(outputs)

  # Detect ssi_run vs legacy and extract needed bits
  if (inherits(outputs, "ssi_run")) {
    # ssi_run path
    obs_list <- outputs$results$obs
    stopifnot(is.list(obs_list), length(obs_list) > 0)

    scen_df <- tibble::as_tibble(outputs$scenarios)
    if (!"scen_num" %in% names(scen_df)) {
      scen_df$scen_num <- seq_len(length(obs_list))
    }
    scen_names <- paste0("scenario_", scen_df$scen_num)
    scen_df$scen <- scen_names

    niter <- length(obs_list[[1]])
    iter_names <- paste0("iter_", seq_len(niter))

    # management names (stable ordering if present)
    known_order <- c("smsy_goal", "s_eq_goal", "smsy_dlm_goal")
    mgmt_names <- intersect(known_order, unique(scen_df$mgmt))
    if (length(mgmt_names) == 0) mgmt_names <- unique(scen_df$mgmt)

    year_names <- names(summary_list$escapement[[1]])
    # From summarize_by_year.ssi_run these are t1..tN; use 1..N for plotting
    year_nums <- seq_along(year_names)
    nyh <- length(year_names)
    year_rows_fun <- function(sim) {
      rows <- seq.int(nrow(sim) - nyh + 1L, nrow(sim))
      rows[rows >= 1L]
    }

  } else {
    # legacy path (kept for backward compatibility)
    obs_list <- outputs$obs_list
    scen_names <- names(summary_list$escapement)
    year_names <- names(summary_list$escapement[[1]])
    # If legacy used "year=####", strip; else just coerce to positions
    if (grepl("^year=", year_names[1])) {
      year_nums <- as.numeric(gsub("^year=", "", year_names))
    } else if (grepl("^t\\d+$", year_names[1])) {
      year_nums <- seq_along(year_names)
    } else {
      year_nums <- seq_along(year_names)
    }
    mgmt_names <- outputs$mgmt_names %||% unique(outputs$scenarios_all$mgmt)
    if (is.null(mgmt_names)) mgmt_names <- c("smsy_goal", "s_eq_goal", "smsy_dlm_goal")

    iter_names <- outputs$iter_names %||% paste0("iter_", seq_len(length(obs_list[[1]])))

    scen_df <- tibble::as_tibble(outputs$scenarios_all) |>
      dplyr::mutate(scen = scen_names)

    year_index <- outputs$year_index
    year_rows_fun <- function(sim) year_index
    nyh <- length(year_names)
  }

  # 2) Tidy the summaries for escapement & harvest
  summary_tbl <- tibble::tibble(
    stat_in_data = unname(summary_stats),
    new_col_name = names(summary_stats)
  )

  data_pluck <- function(component, stat) {
    purrr::map_dfr(scen_names, function(scen) {
      purrr::map_dfr(year_names, function(yr) {
        tibble::tibble(
          scen = scen,
          year = which(year_names == yr), # 1..nyh for plotting
          value = summary_list[[component]][[scen]][[yr]][[stat]]
        )
      })
    })
  }

  tidy_list <- purrr::map(c("escapement", "harvest"), function(component) {
    dfs <- purrr::pmap(summary_tbl, function(stat_in_data, new_col_name) {
      df <- data_pluck(component, stat_in_data)
      df <- dplyr::rename(df, !!new_col_name := value)
      df$label <- ifelse(component == "escapement", "Escapement (1000s)", "Harvest (1000s)")
      df
    })
    Reduce(function(left, right) dplyr::full_join(left, right, by = c("scen", "year", "label")), dfs)
  })

  scenarios_all_df <- scen_df

  combined_df <- dplyr::bind_rows(tidy_list) |>
    dplyr::left_join(scenarios_all_df, by = "scen") |>
    dplyr::mutate(
      mgmt = factor(.data$mgmt, levels = mgmt_names),
      trends = dplyr::case_when(
        .data$trends == "age-sex-length trends" ~ "ASL trends stabilized",
        .data$trends == "continuing trends"     ~ "ASL trends continued",
        TRUE ~ as.character(.data$trends)
      ),
      trends = factor(.data$trends, levels = c("no trends", "ASL trends stabilized", "ASL trends continued")),
      label  = factor(.data$label, levels = c("Harvest (1000s)", "Escapement (1000s)"))
    ) |>
    dplyr::mutate(
      factorMSY = as.character(.data$factorMSY),
      factorMSY = dplyr::case_when(
        factorMSY %in% c("0.75", "liberal")       ~ "liberal",
        factorMSY %in% c("1", "MSY")              ~ "MSY",
        factorMSY %in% c("1.5", "precautionary")  ~ "precautionary",
        TRUE                                       ~ factorMSY
      ),
      factorMSY = factor(factorMSY, levels = c("liberal", "MSY", "precautionary"))
    ) |>
    dplyr::arrange(.data$label) |>
    dplyr::select(-dplyr::any_of(c("ageT", "sizeT", "sexT", "futureT", "sdsel", "maxsel", "scenario_name")))

  # 3) Compute probability of zero harvest on the same 50-year window
  dt_list <- lapply(seq_along(obs_list), function(i) {
    sims <- obs_list[[i]]
    data.table::rbindlist(lapply(seq_along(sims), function(j) {
      sim <- sims[[j]]
      if (!is.data.frame(sim) || is.null(sim$obsHarv)) return(NULL)
      rows <- year_rows_fun(sim)
      data.table::data.table(
        scen   = scen_names[i],
        iter   = iter_names[j],
        year   = seq_len(length(rows)),     # 1..nyh aligned to summaries
        obsHarv = sim$obsHarv[rows]
      )
    }), fill = TRUE, use.names = TRUE)
  })

  dt_all <- data.table::rbindlist(dt_list, fill = TRUE, use.names = TRUE)
  dt_all$zeroHarvest <- ifelse(is.na(dt_all$obsHarv), NA, dt_all$obsHarv == 0)

  zero_harvest_df <- tibble::as_tibble(dt_all) |>
    dplyr::group_by(.data$scen, .data$year) |>
    dplyr::summarise(
      n_zero_harvest = sum(.data$zeroHarvest == TRUE, na.rm = TRUE),
      n_obs          = sum(!is.na(.data$zeroHarvest)),
      prob_zero_harvest = n_zero_harvest / n_obs,
      .groups = "drop"
    ) |>
    dplyr::left_join(scenarios_all_df, by = "scen") |>
    dplyr::select(-dplyr::any_of(c("ageT", "sizeT", "sexT", "futureT", "sdsel", "maxsel", "scenario_name")))

  plot_df <- combined_df |>
    dplyr::left_join(
      zero_harvest_df |> dplyr::select(.data$scen, .data$year, .data$n_obs, .data$n_zero_harvest, .data$prob_zero_harvest),
      by = c("scen", "year")
    ) |>
    dplyr::filter(.data$selectivity == !!selectivity_filter)

  # 4) Pre-smooth probability (5-year centered rolling mean)
  half_primary <- 220000 / 2
  prob_smooth_df <- plot_df |>
    dplyr::filter(.data$label == "Harvest (1000s)") |>
    dplyr::arrange(.data$scen, .data$mgmt, .data$trends, .data$factorMSY, .data$year) |>
    dplyr::group_by(.data$scen, .data$mgmt, .data$trends, .data$factorMSY) |>
    dplyr::mutate(prob_smooth = zoo::rollmean(.data$prob_zero_harvest, k = 5, fill = NA, align = "center")) |>
    dplyr::ungroup()

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plot_component <- function(df, mgmt_name) {
    ggplot2::ggplot(df |> dplyr::filter(.data$mgmt == mgmt_name),
                    ggplot2::aes(x = .data$year, y = .data$mean, ymin = .data$ymin, ymax = .data$ymax,
                                 color = forcats::fct_inorder(.data$label),
                                 group = .data$label, fill = .data$label)) +
      ggplot2::geom_ribbon(alpha = 0.3, color = NA, show.legend = FALSE) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_line(
        data = prob_smooth_df |> dplyr::filter(.data$mgmt == mgmt_name),
        ggplot2::aes(x = .data$year, y = .data$prob_smooth * half_primary, color = "Probability of zero harvest"),
        linewidth = 0.7,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          "Probability of zero harvest" = colors[3],
          "Escapement (1000s)"         = colors[2],
          "Harvest (1000s)"            = colors[1]
        ),
        breaks = c("Harvest (1000s)", "Escapement (1000s)", "Probability of zero harvest")
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Escapement (1000s)" = colors[2],
          "Harvest (1000s)"    = colors[1]
        )
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0.05, 0.05),
        limits = c(0, 220000),
        breaks = seq(0, 220000, 50000),
        labels = ~ .x / 1000,
        sec.axis = ggplot2::sec_axis(
          ~ . / half_primary,
          name = NULL,
          breaks = seq(0, 1, 0.5),
          labels = scales::percent_format(accuracy = 1)
        )
      ) +
      ggplot2::scale_x_continuous(name = "Year", limits = c(1, length(year_names)), breaks = seq(0, length(year_names), 10)) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$trends), cols = ggplot2::vars(.data$factorMSY),
                          scales = "free_y", switch = "y") +
      ggplot2::theme_classic(base_family = "sans") +
      ggplot2::labs(x = "", y = "", color = "") +
      ggplot2::theme(
        strip.background   = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "dashed", linewidth = 0.3),
        strip.text.x       = ggplot2::element_text(size = 10),
        strip.text.y       = ggplot2::element_text(size = 10),
        axis.line          = ggplot2::element_line(linewidth = 0.3),
        axis.text          = ggplot2::element_text(size = 10),
        axis.text.x        = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title         = ggplot2::element_text(size = 12),
        panel.border       = ggplot2::element_rect(fill = NA, linewidth = 0.5),
        legend.key.size    = grid::unit(0.5, 'cm'),
        legend.title       = ggplot2::element_text(size = 10),
        legend.text        = ggplot2::element_text(size = 8),
        strip.placement    = "outside"
      ) -> p

    pdf_path <- file.path(output_dir, paste0(file_basename, "_", mgmt_name, ".pdf"))
    ggplot2::ggsave(pdf_path, p, width = 9, height = 9, units = "in")
    p
  }

  plots <- purrr::map(mgmt_names, ~ plot_component(plot_df, .x))
  names(plots) <- mgmt_names

  csv_path <- file.path(output_dir, paste0(file_basename, ".csv"))
  readr::write_csv(dplyr::arrange(combined_df, .data$mgmt, .data$trends, .data$factorMSY, .data$scen, .data$year), csv_path)

  message("Plots saved to: ", output_dir)
  message("Data saved to: ", csv_path)

  list(data = combined_df, plots = plots)
}
