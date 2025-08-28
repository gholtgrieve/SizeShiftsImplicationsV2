#' Build Selected Figures from an `ssi_run`
#'
#' Convenience wrapper to generate either the **Ohlberger** figures (2–7),
#' the **Kuskokwim** figures (A–C), or a user-specified subset.
#'
#' - `"Ohlberger"` builds Figures **2, 3, 4, 5, 6, 7**.
#' - `"Kuskokwim"` builds Figures **A, B, C**.
#' - You can also pass explicit figure codes, e.g. `c("A","C")` or `c("2","5")`.
#'   (Case-insensitive; `"Figure5"` also works.)
#'
#' @param figures Character vector indicating what to build.
#'   Examples: `"Ohlberger"`, `"Kuskokwim"`, `c("A","C")`, `c("2","5")`.
#'   Defaults to `"Kuskokwim"`.
#' @param data Either an `ssi_run` object (from [run_scenarios()]) **or**
#'   a file path to a saved `.rds` containing an `ssi_run`.
#' @param figure_dir Directory to save figures (default: `here()/figures`).
#'
#' @return Invisible named list of results for each figure generated.
#' @export

make_figures <- function(
    figures    = "Kuskokwim",
    data,
    figure_dir = file.path(here::here(), "figures")
) {
  # ---- normalize figure selection ---------------------------------------------
  normalize_tokens <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("^figure\\s*", "", x, ignore.case = TRUE)  # "Figure5" -> "5"
    tokens <- toupper(x)

    expand_one <- function(tok) {
      if (grepl("^OHL", tok))  return(c("2","3","4","5","6","7"))      # Ohlberger group
      if (grepl("^KUSK", tok)) return(c("A","B","C","D","E"))          # Kuskokwim group
      if (tok %in% c("A","B","C","D","E")) return(tok)
      if (tok %in% c("2","3","4","5","6","7")) return(tok)
      character(0)
    }

    out <- unlist(lapply(tokens, expand_one), use.names = FALSE)
    out[!duplicated(out)]
  }

  wanted <- normalize_tokens(figures)
  if (length(wanted) == 0L) {
    stop("No valid figures selected. Use 'Ohlberger', 'Kuskokwim', or specific codes like 'A','C','E','2','5'.")
  }

  # ---- load/normalize the run into `run_obj` -----------------------------------
  if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  if (inherits(data, "ssi_run")) {
    run_obj <- data
  } else if (is.character(data) && length(data) == 1L && file.exists(data)) {
    run_obj <- readRDS(data)
    if (!inherits(run_obj, "ssi_run")) stop("`data` path did not contain an 'ssi_run' object.")
  } else {
    stop("`data` must be an 'ssi_run' object or a path to a saved '.rds' run.")
  }

  # ---- dispatch map ------------------------------------------------------------
  # All functions should accept (data = run_obj, output_dir = figure_dir)
  builders <- list(
    A   = function() .make_Kusko_figure_A(data = run_obj, output_dir = figure_dir),
    B   = function() .make_Kusko_figure_B(data = run_obj, output_dir = figure_dir),
    C   = function() .make_Kusko_figure_C(data = run_obj, output_dir = figure_dir),
    D   = function() .make_Kusko_figure_D(data = run_obj, output_dir = figure_dir),
    E   = function() .make_Kusko_figure_E(data = run_obj, output_dir = figure_dir),
    `2` = function() .make_Ohlberger_figure_2(data = run_obj, output_dir = figure_dir),
    `3` = function() .make_Ohlberger_figure_3(data = run_obj, output_dir = figure_dir),
    `4` = function() .make_Ohlberger_figure_4(data = run_obj, output_dir = figure_dir),
    `5` = function() .make_Ohlberger_figure_5(data = run_obj, output_dir = figure_dir),
    `6` = function() .make_Ohlberger_figure_6(data = run_obj, output_dir = figure_dir),
    `7` = function() .make_Ohlberger_figure_7(data = run_obj, output_dir = figure_dir)
  )

  missing <- setdiff(wanted, names(builders))
  if (length(missing)) stop("No builder available for: ", paste(missing, collapse = ", "))

  # ---- build requested figures -------------------------------------------------
  out <- list(data = run_obj)
  for (code in wanted) {
    res <- builders[[code]]()
    out[[paste0("Figure", code)]] <- res
  }

  return(invisible(out))
}
