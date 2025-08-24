#' Harvest selectivity based on fish size and mesh size
#'
#' Internal function to calculate the proportion of fish retained by a gillnet
#' of a given mesh size, using the selectivity model described by Bromaghin (2005).
#'
#' @param size Numeric. Fish size in millimeters.
#' @param meshsize Numeric. Stretch mesh size in inches.
#' @param s Numeric. Standard deviation of selectivity, typically \code{0.204}.
#'
#' @return Numeric. Proportion of fish of the given size retained by the net.
#'
#' @references
#' Bromaghin, J. F. (2005). The statistical evaluation of relative gillnet selectivity
#' estimates and their application to the management of Alaskan sockeye salmon.
#' \emph{Alaska Fishery Research Bulletin}, 11(2), 258–266.

#'
#' @examples
#' selectivity(size = 550, meshsize = 5.0, s = 0.204)
#'
#' @keywords internal
#' @export

selectivity <- function(size, meshsize, s) {

  t <- 1.920
  h <- 0.622
  l <- -0.547
  perimeter <- 25.4 * meshsize * 2
  x <- size / perimeter

  select <- (1 + l^2 / (4 * h^2))^h *
    (1 + (x - s * l / 2 * h - t)^2 / (s^2))^-h *
    exp(-l * (atan((x - s * l / 2 * h - t) / s) + atan(l / (2 * h))))

  return(select)
}
