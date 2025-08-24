#' Generate recruits from spawners using a Ricker stock-recruitment model
#'
#' Internal function to simulate recruitment from a given number of spawners
#' using the Ricker model, with log-normal environmental variability and
#' temporal autocorrelation in residuals.
#'
#' @param spawn Numeric. Number of spawners in the current year.
#' @param sigma Numeric. Standard deviation of the random error term.
#' @param alpha Numeric. Productivity parameter of the Ricker model.
#' @param beta Numeric. Density-dependence parameter of the Ricker model.
#' @param rho Numeric. Temporal autocorrelation coefficient (between 0 and 1).
#' @param last.eps Numeric. Recruitment residual (in log-space) from the
#'   previous year.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{rec}}{Numeric. Recruits in the current year.}
#'   \item{\code{eps}}{Numeric. Residual (in log-space) for the current year.}
#' }
#'
#' @examples
#' ricker(spawn = 100, sigma = 0.3, alpha = 5, beta = 0.01,
#'         rho = 0.5, last.eps = 0.1)
#'
#' @keywords internal
#' @export


ricker <- function(spawn, sigma, alpha, beta, rho, last.eps) {

  ## normal random error with bias correction for the mean
  delta <- rnorm(1, mean = -(sigma^2)/2, sd = sigma)

  ## residual (year y) in log-space with temporal autocorrelation
  eps <- rho * last.eps + sqrt(1 - rho^2) * delta

  ## ln(recruits/spawner) with autocorrelated residuals
  lnrs <- log(alpha) - beta * spawn + eps

  ## recruits
  rec <- exp(lnrs) * spawn

  return(list(
    rec = rec,
    eps = eps
  ))
}
