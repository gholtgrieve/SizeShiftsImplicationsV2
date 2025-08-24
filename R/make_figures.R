#' Wrapper to Generate Figures A, B, and C from Simulation Outputs
#'
#' This function loads model outputs and generates figure panels A, B, and/or C.
#' You can pass either the in-memory run object returned by `run_scenarios()`,
#' a path to a saved `.rds` from `run_scenarios()`, or (for back-compat) a
#' `timestamp` + `data_dir` pair to locate the saved `.rds`.
#'
#' @param figures Character vector of figures to generate. Options: "A","B","C". Default: all.
#' @param run Either an `ssi_run` object (returned by `run_scenarios()`) or a path to its `.rds`. Optional.
#' @param timestamp (Back-compat) Timestamp string used to find `run_<timestamp>.rds` in `data_dir`.
#'   Accepts either `"YYYY-MM-DD_HH-MM-SS"` or `"YYYY-MM-DD HH-MM-SS"` (space will be normalized to underscore).
#' @param selectivity_filter Character filter for selectivity type (default: "unselective").
#' @param data_dir Directory containing saved `.rds` runs (default: `here()/outputs`).
#' @param figure_dir Directory to save figures (default: `here()/figures`).
#' @param colors Optional vector of colors for plotting.
#'
#' @return A named list containing results for each figure generated.
#' @export
make_figures <- function(
    figures = c("A", "B", "C"),
    run = NULL,
    timestamp = NULL,
    selectivity_filter = "unselective",
    data_dir = paste0(here::here(), "/outputs"),
    figure_dir = paste0(here::here(), "/figures"),
    colors = c("darkgray", "deepskyblue3", "orange")
) {
  # Validate figures
  figures <- toupper(figures)
  valid_figs <- c("A","B","C")
  if (!all(figures %in% valid_figs)) {
    stop("Invalid `figures`. Allowed: 'A','B','C'.")
  }

  # Ensure figure dir
  if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

  # ---- Load/normalize inputs into a common `outputs` object (must be ssi_run) ----
  if (!is.null(run) && inherits(run, "ssi_run")) {
    outputs <- run                                    # keep full ssi_run
  } else if (is.character(run) && length(run) == 1L && file.exists(run)) {
    run_obj <- readRDS(run)
    if (!inherits(run_obj, "ssi_run")) stop("`run` path did not contain an 'ssi_run' object.")
    outputs <- run_obj
  } else {
    if (is.null(timestamp) || !is.character(timestamp) || nchar(timestamp) == 0) {
      stop("Provide either `run` (ssi_run or .rds path) OR a valid `timestamp` to locate the saved run.")
    }
    ts_norm <- gsub(" ", "_", timestamp)
    rds_path <- file.path(data_dir, paste0("run_", ts_norm, ".rds"))
    if (!file.exists(rds_path)) stop("Could not find saved run at: ", rds_path)
    run_obj <- readRDS(rds_path)
    if (!inherits(run_obj, "ssi_run")) stop("Saved file is not an 'ssi_run' object: ", rds_path)
    outputs <- run_obj
  }

  # Alias scenarios for helpers that expect `scenarios_all`
  scenarios_all <- outputs$scenarios

  results <- list(outputs = outputs)

  # If A or B, compute summary (pass the ssi_run directly)
  if (any(c("A","B") %in% figures)) {
    MeanByScenario <- SizeShiftsImplicationsV2::summarize_50_year_avg(data = outputs)
    results$MeanByScenario <- MeanByScenario
  }

  # Figure A
  if ("A" %in% figures) {
    results$figureA <- SizeShiftsImplicationsV2::make_figure_A(
      data = results$MeanByScenario,
      scenarios_all = scenarios_all,
      use_quant = "mean",
      selectivity_filter = selectivity_filter,
      colors = colors,
      output_dir = figure_dir,
      file_basename = "FigureA"
    )
  }

  # Figure B
  if ("B" %in% figures) {
    results$figureB <- SizeShiftsImplicationsV2::make_figure_B(
      data = results$MeanByScenario,
      scenarios_all = scenarios_all,
      selectivity_filter = selectivity_filter,
      colors = colors,
      output_dir = figure_dir,
      file_basename = "FigureB"
    )
  }

  # Figure C
  if ("C" %in% figures) {
    results$figureC <- SizeShiftsImplicationsV2::make_figure_C(
      outputs = outputs,
      selectivity_filter = selectivity_filter,
      summary_stats = c(mean = "Mean", ymin = "25%", ymax = "75%"),
      colors = colors,
      output_dir = figure_dir,
      file_basename = "FigureC"
    )
  }

  return(results)
}
