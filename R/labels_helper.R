#' Standardize scenario labels and factor orders (internal)
#'
#' Ensures consistent printed labels across all figures and summaries.
#'
#' - **factorMSY**: `0.75` → `"liberal"`, `1` → `"MSY"`, `1.5` → `"precautionary"`
#' - **trends**:
#'   - `"no trends"` (and `"No trends"`) → `"no trends"`
#'   - `"age-sex-length trends"` → `"ASL trends stabilized"`
#'   - `"age-length trends"` → `"AL trends stabilized"`
#'   - `"continuing trends"` → `"ASL trends continued"`
#' - **mgmt**: map common long-form to short codes: `"smsy_goal"` → `"TRM"`, `"s_eq_goal"` → `"YPR"`, `"smsy_dlm_goal"` → `"DLM"`
#' - **selectivity**: standardize to `"small-mesh"`, `"unselective"`, `"large-mesh"`
#'
#' If a column is absent, it is ignored. Factors are (re)leveled for stable facet/legend order:
#' - `factorMSY`: `c("liberal","MSY","precautionary")`
#' - `trends`:    `c("no trends","ASL trends stabilized","AL trends stabilized","ASL trends continued")`
#' - `mgmt`:      `c("TRM","YPR","DLM")` (plus any other levels appended at end, if present)
#' - `selectivity`: `c("small-mesh","unselective","large-mesh")`
#'
#' @param df A data.frame/tibble with any of the columns above.
#' @return The same data.frame with standardized character values and ordered factors.
#' @keywords internal
#' @export
standardize_scenario_labels <- function(df) {
  # factorMSY --------------------------------------------------------------
  if ("factorMSY" %in% names(df)) {
    fmsy_chr <- as.character(df$factorMSY)
    fmsy_chr <- dplyr::case_when(
      fmsy_chr %in% c("0.75", "0.750", "0.75 ", "liberal") ~ "liberal",
      fmsy_chr %in% c("1", "1.0", "MSY")                   ~ "MSY",
      fmsy_chr %in% c("1.5", "1.50", "1.5 ", "precautionary") ~ "precautionary",
      TRUE ~ fmsy_chr
    )
    df$factorMSY <- factor(fmsy_chr, levels = c("liberal", "MSY", "precautionary"))
  }

  # trends -----------------------------------------------------------------
  if ("trends" %in% names(df)) {
    tr_chr <- as.character(df$trends)
    tr_chr <- dplyr::case_when(
      tr_chr %in% c("no trends", "No trends")                          ~ "no trends",
      tr_chr %in% c("age-sex-length trends", "ASL trends stabilized")  ~ "ASL trends stabilized",
      tr_chr %in% c("age-length trends", "AL trends stabilized")       ~ "AL trends stabilized",
      tr_chr %in% c("continuing trends", "ASL trends continued")       ~ "ASL trends continued",
      TRUE ~ tr_chr
    )
    df$trends <- factor(
      tr_chr,
      levels = c("no trends", "ASL trends stabilized", "AL trends stabilized", "ASL trends continued")
    )
  }

  # mgmt -------------------------------------------------------------------
  if ("mgmt" %in% names(df)) {
    mg_chr <- as.character(df$mgmt)
    mg_chr <- dplyr::case_when(
      mg_chr %in% c("smsy_goal", "TRM")     ~ "TRM",
      mg_chr %in% c("s_eq_goal", "YPR")     ~ "YPR",
      mg_chr %in% c("smsy_dlm_goal", "DLM") ~ "DLM",
      TRUE ~ mg_chr
    )
    base_levels <- c("TRM", "YPR", "DLM")
    # preserve base order; append any extras (rare)
    mg_levels <- unique(c(base_levels, setdiff(mg_chr, base_levels)))
    df$mgmt <- factor(mg_chr, levels = mg_levels)
  }

  # selectivity ------------------------------------------------------------
  if ("selectivity" %in% names(df)) {
    sel_chr <- as.character(df$selectivity)
    sel_chr <- dplyr::case_when(
      sel_chr %in% c("6 inch gillnet", "small-mesh", "small mesh") ~ "small-mesh",
      sel_chr %in% c("unselective", "Unselective")                 ~ "unselective",
      sel_chr %in% c("8.5 inch gillnet", "large-mesh", "large mesh") ~ "large-mesh",
      TRUE ~ sel_chr
    )
    df$selectivity <- factor(sel_chr, levels = c("small-mesh", "unselective", "large-mesh"))
  }

  return(df)
}
