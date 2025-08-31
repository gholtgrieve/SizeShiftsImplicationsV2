#' Generate age composition of brood year recruits
#'
#' Internal function to simulate the age distribution of recruits from a given brood year.
#'
#' @param ages Numeric vector. All simulated age classes.
#' @param recruits Numeric. Number of recruits from that brood year.
#' @param meanage Numeric. Mean age of recruits.
#' @param sdage Numeric. Standard deviation of the age distribution.
#'
#' @return Numeric vector. Number of recruits in each age class.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' ages <- 2:6
#' agecomp(ages, recruits = 100, meanage = 4, sdage = 0.8)
#' }
#' @export


.calc_agecomp <- function(ages, recruits, meanage, sdage) {

  ## probabilities by age given mean age
  probs_a <- dnorm(ages, meanage, sdage)
  probs_a <- probs_a / sum(probs_a)

  ## draw Dirichlet based on probabilities by age
  d <- 50
  probs <- as.vector(gtools::rdirichlet(1, probs_a * d))

  ## age counts
  agecounts <- round(probs * recruits)

  return(agecounts)
}
