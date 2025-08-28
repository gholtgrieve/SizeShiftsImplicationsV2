#' Scenarios list
#'
#' A tibble of scenarios from Ohlberger et al. 2025.
#'
#' @references
#' Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025).
#' *Accounting for salmon body size declines in fishery management can reduce conservation risks*.
#' Fish and Fisheries, 26, 113–130. \doi{10.1111/faf.12869}
#'
#' @format A tibble with rows × 11 columns:
#' \describe{
#'   \item{ageT}{integer}
#'   \item{sizeT}{integer}
#'   \item{sexT}{integer}
#'   \item{futureT}{character}
#'   \item{trends}{character}
#'   \item{sdsel}{double}
#'   \item{maxsel}{double}
#'   \item{selectivity}{factor}
#'   \item{mgmt}{character}
#'   \item{factorMSY}{factor}
#'   \item{scenario_name}{character}
#' }
#' @source Generated from `data-raw/scenarios.R`.
#' @docType data
#' @keywords datasets
#' @name scenarios
#' @usage data(scenarios)
NULL
