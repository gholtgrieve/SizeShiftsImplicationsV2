#' Select scenarios to run
#'
#' Filter the package's built-in scenarios table using a single selector argument.
#' The selector can be:
#' \itemize{
#'   \item `"all"` (or `"Ohlberger"`) — return all rows (subject to \code{enforce_constraints}).
#'   \item `"Kuskokwim"` — a convenience alias that keeps only
#'         \code{selectivity == "unselective"} and \code{trends != "age-length trends"}
#'         (still subject to \code{enforce_constraints}).
#'   \item A single dplyr-style filter \strong{string}, parsed with \code{rlang::parse_expr()},
#'         e.g. \code{"selectivity == 'unselective' & factorMSY \%in\% c(0.75, 1.50)"}.
#'   \item A \strong{numeric} vector — treated as \code{scen_num} values if that
#'         column exists; otherwise as row indices.
#' }
#'
#' By default, the function enforces your design constraint: for gillnets
#' (\code{"8.5 inch gillnet"} or \code{"6 inch gillnet"}) only \code{factorMSY == 1.00}
#' is allowed; for \code{"unselective"} nets, \code{factorMSY} may be
#' \code{1.00}, \code{0.75}, or \code{1.50}.
#'
#' @param selector See details above. Case-insensitive for the special aliases
#'   \code{"all"}, \code{"Ohlberger"}, and \code{"Kuskokwim"}.
#' @param enforce_constraints Logical, default \code{TRUE}. If \code{TRUE}, keep only
#'   the allowed combinations: for `"unselective"`, \code{factorMSY} in
#'   \code{c(1.00, 0.75, 1.50)}; for `"8.5 inch gillnet"` or `"6 inch gillnet"`,
#'   \code{factorMSY == 1.00}.
#'
#' @return A tibble with the selected scenarios (zero rows if nothing matches).
#'
#' @details
#' \itemize{
#'   \item Uses the package's built-in \emph{scenarios} dataset; no argument is required.
#'   \item The function coerces \code{factorMSY} to numeric internally to make filtering robust.
#'   \item The filter string is parsed with \code{rlang::parse_expr()} and evaluated in the
#'         context of the scenarios table, so column names should be referenced directly
#'         (e.g., \code{mgmt == 'smsy_goal'}).
#'   \item When \code{selector} is numeric and \code{scen_num} exists, selection is by ID and
#'         does not depend on row order.
#' }
#'
#' @examples
#' # 1) Return all allowed scenarios (built-in dataset)
#' head(select_scenarios("all"))
#'
#' # 2) Convenience aliases
#' kusk <- select_scenarios("Kuskokwim")   # unselective & trends != "age-length trends"
#' ohl  <- select_scenarios("Ohlberger")   # same as "all"
#'
#' # 3) Filter by expression: unselective nets, no DLM mgmt, conservative or aggressive MSY
#' sel <- "selectivity == 'unselective' & mgmt %in% c('smsy_goal','s_eq_goal') & factorMSY %in% c(0.75, 1.50)"
#' out <- select_scenarios(sel)
#' nrow(out)
#'
#' # 4) Select by scenario numbers (IDs)
#' select_scenarios(c(1, 12, 37, 60))
#'
#' # 5) Inspect everything in the table (including intentionally omitted combos)
#' #    by disabling constraints (useful for QA)
#' all_raw <- select_scenarios("all", enforce_constraints = FALSE)
#' nrow(all_raw)
#'
#' @keywords internal
#' @export
select_scenarios <- function(selector = "all", enforce_constraints = TRUE) {
  scens <- scenarios  # built-in dataset

  # Coerce factorMSY robustly
  if ("factorMSY" %in% names(scens)) {
    scens$factorMSY <- as.numeric(scens$factorMSY)
  }

  # Enforce default constraints (allowed combinations)
  if (isTRUE(enforce_constraints)) {
    scens <- dplyr::filter(
      scens,
      (selectivity == "unselective" &
         factorMSY %in% c(1.00, 0.75, 1.50)) |
        (selectivity %in% c("8.5 inch gillnet", "6 inch gillnet") &
           factorMSY == 1.00)
    )
  }

  # Handle character selectors (aliases or filter string)
  if (is.character(selector) && length(selector) == 1L) {
    key <- tolower(selector)

    # Aliases
    if (key %in% c("all", "ohlberger")) {
      return(tibble::as_tibble(scens))
    }
    if (key == "kuskokwim") {
      out <- dplyr::filter(
        scens,
        selectivity == "unselective",
        trends != "age-length trends"
      )
      return(tibble::as_tibble(out))
    }

    # Otherwise: treat as a filter expression string
    expr <- try(rlang::parse_expr(selector), silent = TRUE)
    if (inherits(expr, "try-error")) {
      cond <- attr(expr, "condition")
      stop("Could not parse selector: ", conditionMessage(cond))
    }
    out <- dplyr::filter(scens, !!expr)
    return(tibble::as_tibble(out))
  }

  # Numeric IDs (prefer scen_num if present, else row indices)
  if (is.numeric(selector)) {
    idx <- unique(as.integer(selector))
    if ("scen_num" %in% names(scens)) {
      out <- dplyr::filter(scens, scen_num %in% idx)
    } else {
      idx <- idx[idx >= 1 & idx <= nrow(scens)]
      out <- dplyr::slice(scens, idx)
    }
    return(tibble::as_tibble(out))
  }

  stop("`selector` must be 'all'/'Ohlberger', 'Kuskokwim', a single filter string, or numeric scen_num/row indices.")
}
