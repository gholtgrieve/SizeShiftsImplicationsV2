#' Fit dynamic linear models (DLMs) with time-varying intercept and/or slope
#'
#' Internal function to fit a dynamic linear model to escapement and recruitment data,
#' allowing for time-varying intercept (\eqn{\alpha}) and/or slope (\eqn{\beta})
#' parameters in a linearized Ricker model.
#'
#' @param data A data frame containing:
#'   \describe{
#'     \item{\code{Esc}}{Numeric. Escapement values.}
#'     \item{\code{Rec}}{Numeric. Recruitment values.}
#'   }
#' @param var_alpha Logical. If \code{TRUE}, \eqn{\alpha} is estimated as a time-varying parameter.
#' @param var_beta Logical. If \code{TRUE}, \eqn{\beta} is estimated as a time-varying parameter.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{results}}{Data frame of original data with columns \code{alpha_y} and \code{beta_y} (smoothed parameters).}
#'   \item{\code{AICc}}{Numeric. Corrected Akaike Information Criterion.}
#'   \item{\code{sigma}}{Numeric. Standard deviation of observation error.}
#' }
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(Esc = c(100, 120, 150, 130),
#'                   Rec = c(500, 550, 600, 580))
#' DLMfit(dat, var_alpha = TRUE, var_beta = FALSE)
#' }
#'
#' @keywords internal
#' @export

DLMfit <- function(data, var_alpha, var_beta) {

  lnRS <- log(data$Rec / data$Esc)

  alpha_y <- beta_y <- NULL
  mod <- dlm::dlmModReg(data$Esc)
  npara <- 3 + sum(c(var_alpha, var_beta))

  build_mod <- function(parm) {
    mod$V <- exp(parm[1])
    mod$W[1,1] <- mod$W[2,2] <- 0

    if (var_alpha) mod$W[1,1] <- exp(parm[2])
    if (var_beta)  mod$W[2,2] <- exp(parm[2])
    if (var_alpha & var_beta) {
      mod$W[1,1] <- exp(parm[2])
      mod$W[2,2] <- exp(parm[3])
    }
    return(mod)
  }

  dlm_out <- dlm::dlmMLE(y = lnRS, build = build_mod, parm = c(-.1,-.1,-.1), method = "L-BFGS-B")
  lls <- dlm_out$value
  dlmMod <- build_mod(dlm_out$par)
  sigma <- sqrt(dlmMod$V)

  outsFilter <- dlm::dlmFilter(y = lnRS, mod = dlmMod)
  outsSmooth <- dlm::dlmSmooth(outsFilter)

  alpha_y <- signif(outsSmooth$s[-1,1], 4)
  beta_y <- signif(outsSmooth$s[-1,2], 5)

  results_df <- data
  results_df$alpha_y <- alpha_y
  results_df$beta_y <- beta_y

  AICc <- 2 * lls + 2 * npara + (2 * npara * (npara + 1) / (length(data$Rec) - npara - 1))

  output <- list(results = results_df, AICc = AICc, sigma = sigma)
  return(output)
}
