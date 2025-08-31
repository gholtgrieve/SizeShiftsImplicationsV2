#' Calculate fecundity or egg mass based on female size
#'
#' Internal function to compute fecundity and egg mass from female length using
#' log-log allometric relationships.
#'
#' @param size Numeric. Female length in millimeters.
#' @param allometry Named numeric vector of allometric parameters:
#'   \describe{
#'     \item{\code{size_ref}}{Reference size.}
#'     \item{\code{fec_ref}}{Fecundity at reference size.}
#'     \item{\code{b_fec}}{Slope of log-log fecundity relationship.}
#'     \item{\code{egg_ref}}{Egg mass at reference size.}
#'     \item{\code{b_eggs}}{Slope of log-log egg mass relationship.}
#'   }
#'
#' @return Named list with:
#' \describe{
#'   \item{\code{fecundity}}{Numeric vector of fecundity (log-scale) for each size.}
#'   \item{\code{eggmass}}{Numeric vector of egg mass (log-scale) for each size.}
#' }
#'
#' @examples
#' \dontrun{
#' allom <- c(size_ref = 500, fec_ref = 1000, b_fec = 3.0,
#'            egg_ref = 0.2, b_eggs = 1.0)
#' reprod_output(size = 550, allom)
#' }
#' @keywords internal
#' @export

.calc_reprod_output <- function(size, allometry) {

  size_ref <- allometry[names(allometry) == "size_ref"]

  ## allometry for fecundity
  fec_ref <- allometry[names(allometry) == "fec_ref"]
  b_fec <- allometry[names(allometry) == "b_fec"]
  a_fec <- exp(log(fec_ref) - b_fec * log(size_ref))
  fecundity <- log(a_fec) + b_fec * log(size)

  ## allometry for egg mass
  egg_ref <- allometry[names(allometry) == "egg_ref"]
  b_eggs <- allometry[names(allometry) == "b_eggs"]
  a_eggs <- exp(log(egg_ref) - b_eggs * log(size_ref))
  eggmass <- log(a_eggs) + b_eggs * log(size)

  return(list(
    fecundity = fecundity,
    eggmass = eggmass
  ))
}
