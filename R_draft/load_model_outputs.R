#' Load model outputs and scenario metadata
#'
#' Loads simulation output (`.Rdata`), parameter objects (`_parms.Rdata`),
#' and scenario metadata (`_scen.xlsx`) from a specified directory.
#' Applies consistent labeling and prepares indexing variables to support
#' tidy downstream analysis and plotting.
#'
#' @param main_dir Character. Path to the directory containing the simulation files.
#' @param timestamp Character. Timestamp suffix used in file names (e.g., "2025-08-04").
#'
#' @return A named list containing:
#' \describe{
#'   \item{obs_list}{List of simulation outputs loaded from `.Rdata`.}
#'   \item{parameters}{List of simulation parameters such as \code{nscen}, \code{niter}, and \code{nyh}.}
#'   \item{scenarios_all}{Cleaned data frame of scenario metadata with labeled factor columns.}
#'   \item{scenario_names}{Character vector combining trend, mgmt, and MSY labels.}
#'   \item{trend_names}{Character vector of unique trend types.}
#'   \item{mgmt_names}{Character vector of unique management types.}
#'   \item{factorMSY_names}{Character vector of unique MSY factor levels.}
#'   \item{scen_names}{Character vector of scenario index labels (e.g., \code{"scen=1"}).}
#'   \item{iter_names}{Character vector of iteration index labels (e.g., \code{"iter=1"}).}
#'   \item{year_names}{Character vector of hindcast year labels (e.g., \code{"year=1"}).}
#'   \item{year_index}{Integer vector indexing the most recent 50 years.}
#' }
#'
#' @export
load_model_outputs <- function(main_dir, timestamp) {
  # Define paths and files names
  path <- fs::path(main_dir)
  rdata_path <- fs::path(path, paste0("run_", timestamp, ".Rdata"))
  parameter_path <- fs::path(path, paste0("run_", timestamp, "_parms.Rdata"))
  scenarios_path <- fs::path(path, paste0("run_", timestamp, "_scen.xlsx"))
  
  # Load model output (obs_list)
  load(rdata_path, envir = rlang::current_env())
  
  # Load and unpack parameter list
  parameters_raw <- readRDS(parameter_path)
  parameters <- parameters_raw |>
    purrr::map(~ unlist(.x, recursive = TRUE, use.names = TRUE)) |>
    purrr::list_flatten() |>
    as.list()
  
  # Inject global variables into parameters if present
  if (exists("nscen", inherits = FALSE)) parameters$nscen <- nscen
  if (exists("niter", inherits = FALSE)) parameters$niter <- niter
  if (exists("nyh", inherits = FALSE))   parameters$nyh   <- nyh
  
  # Load and clean scenario metadata
  scenarios_all <- readxl::read_excel(scenarios_path) |>
    dplyr::select(-1) |>
    dplyr::mutate(
      selectivity = factor(
        selectivity,
        levels = c("6 inch gillnet", "unselective", "8.5 inch gillnet"),
        labels = c("small-mesh", "unselective", "large-mesh")
      ),
      mgmt = dplyr::case_when(
        mgmt == "smsy_goal" ~ "TRM",
        mgmt == "s_eq_goal" ~ "YPR",
        mgmt == "smsy_dlm_goal" ~ "DLM",
        TRUE ~ mgmt
      ),
      factorMSY = factor(
        factorMSY,
        levels = c("0.75", "1", "1.5"),
        labels = c("liberal", "MSY", "precautionary")
      ),
      trends = dplyr::case_when(
        trends == "age-sex-length trends" ~ "ASL trends stabilized",
        trends == "continuing trends" ~ "ASL trends continued",
        TRUE ~ trends
      ),
      scenario_name = paste(trends, mgmt, factorMSY)
    )
  
  # Validate scenario count
  stopifnot(nrow(scenarios_all) == parameters$nscen)
  
  # Extract factor levels
  trend_names <- scenarios_all |> dplyr::distinct(trends) |> dplyr::pull(trends)
  mgmt_names <- scenarios_all |> dplyr::distinct(mgmt) |> dplyr::pull(mgmt)
  factorMSY_names <- scenarios_all |> dplyr::distinct(factorMSY) |> dplyr::pull(factorMSY)
  
  # Create index labels
  scen_names <- paste0("scen=", seq_len(parameters$nscen))
  iter_names <- paste0("iter=", seq_len(parameters$niter))
  year_names <- paste0("year=", seq_len(parameters$nyh))
  
  # Define index for most recent 50-year period
  ny_obs <- nrow(obs.list[[1]][[1]])
  stopifnot(ny_obs >= parameters$nyh)
  year_index <- (parameters$nyh + 1):ny_obs
  
  return(list(
    obs_list = obs.list,
    parameters = parameters,
    scenarios_all = scenarios_all,
    scenario_names = scenarios_all$scenario_name,
    trend_names = trend_names,
    mgmt_names = mgmt_names,
    factorMSY_names = factorMSY_names,
    scen_names = scen_names,
    iter_names = iter_names,
    year_names = year_names,
    year_index = year_index
  ))
}
