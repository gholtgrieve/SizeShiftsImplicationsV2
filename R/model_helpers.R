# Internal helpers for run-model()


#' Generate age composition of brood year recruits
#'
#' Internal function to simulate the age distribution of recruits from a given brood year.
#'
#' @param ages Numeric vector. All simulated age classes.
#' @param recruits Numeric. Number of recruits from that brood year.
#' @param meanage Numeric. Mean age of recruits.
#' @param sdage Numeric. Standard deviation of the age distribution.
#' @return Numeric vector. Number of recruits in each age class.
#' @keywords internal
#' @examples
#' ages <- 2:6
#' .calc_agecomp(ages, recruits = 100, meanage = 4, sdage = 0.8)
.calc_agecomp <- function(ages, recruits, meanage, sdage) {

  # Probabilities by age given mean age
  probs_a <- dnorm(ages, meanage, sdage)
  probs_a <- probs_a / sum(probs_a)

  # Draw Dirichlet based on probabilities by age
  d <- 50
  probs <- as.vector(gtools::rdirichlet(1, probs_a * d))

  # Age counts
  agecounts <- round(probs * recruits)

  agecounts
}


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
#'   \describe{
#'     \item{\code{results}}{Data frame of original data with columns \code{alpha_y} and \code{beta_y} (smoothed parameters).}
#'     \item{\code{AICc}}{Numeric. Corrected Akaike Information Criterion.}
#'     \item{\code{sigma}}{Numeric. Standard deviation of observation error.}
#'   }
#' @examples
#' dat <- data.frame(Esc = c(100, 120, 150, 130), Rec = c(500, 550, 600, 580))
#' .calc_DLMfit(dat, var_alpha = TRUE, var_beta = FALSE)
#' @keywords internal
.calc_DLMfit <- function(data, var_alpha, var_beta) {

  lnRS <- log(data$Rec / data$Esc)
  alpha_y <- beta_y <- NULL
  mod <- dlm::dlmModReg(data$Esc)
  npara <- 3 + sum(c(var_alpha, var_beta))

  build_mod <- function(parm) {
    mod$V <- exp(parm[1])
    mod$W[1, 1] <- mod$W[2, 2] <- 0

    if (var_alpha) mod$W[1, 1] <- exp(parm[2])
    if (var_beta) mod$W[2, 2] <- exp(parm[2])
    if (var_alpha & var_beta) {
      mod$W[1, 1] <- exp(parm[2])
      mod$W[2, 2] <- exp(parm[3])
    }
    mod
  }

  dlm_out <- dlm::dlmMLE(y = lnRS, build = build_mod, parm = c(-0.1, -0.1, -0.1), method = "L-BFGS-B")
  lls <- dlm_out$value
  dlmMod <- build_mod(dlm_out$par)
  sigma <- sqrt(dlmMod$V)

  outsFilter <- dlm::dlmFilter(y = lnRS, mod = dlmMod)
  outsSmooth <- dlm::dlmSmooth(outsFilter)

  alpha_y <- signif(outsSmooth$s[-1, 1], 4)
  beta_y <- signif(outsSmooth$s[-1, 2], 5)

  results_df <- data
  results_df$alpha_y <- alpha_y
  results_df$beta_y <- beta_y

  AICc <- 2 * lls + 2 * npara + (2 * npara * (npara + 1) / (length(data$Rec) - npara - 1))

  list(results = results_df, AICc = AICc, sigma = sigma)
}


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
#' @return Named list with:
#'   \describe{
#'     \item{\code{fecundity}}{Numeric vector of fecundity (log-scale) for each size.}
#'     \item{\code{eggmass}}{Numeric vector of egg mass (log-scale) for each size.}
#'   }
#' @examples
#' allom <- c(size_ref = 500, fec_ref = 1000, b_fec = 3.0, egg_ref = 0.2, b_eggs = 1.0)
#' .calc_reprod_output(size = 550, allometry = allom)
#' @keywords internal
.calc_reprod_output <- function(size, allometry) {

  size_ref <- allometry["size_ref"]

  # Allometry for fecundity
  fec_ref <- allometry["fec_ref"]
  b_fec <- allometry["b_fec"]
  a_fec <- exp(log(fec_ref) - b_fec * log(size_ref))
  fecundity <- log(a_fec) + b_fec * log(size)

  # Allometry for egg mass
  egg_ref <- allometry["egg_ref"]
  b_eggs <- allometry["b_eggs"]
  a_eggs <- exp(log(egg_ref) - b_eggs * log(size_ref))
  eggmass <- log(a_eggs) + b_eggs * log(size)

  list(fecundity = fecundity, eggmass = eggmass)
}


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
#' @param last.eps Numeric. Recruitment residual (in log-space) from the previous year.
#' @return A named list with:
#'   \describe{
#'     \item{\code{rec}}{Numeric. Recruits in the current year.}
#'     \item{\code{eps}}{Numeric. Residual (in log-space) for the current year.}
#'   }
#' @examples
#' .calc_ricker(spawn = 100, sigma = 0.3, alpha = 5, beta = 0.01,
#'              rho = 0.5, last.eps = 0.1)
#' @keywords internal
.calc_ricker <- function(spawn, sigma, alpha, beta, rho, last.eps) {

  # Normal random error with bias correction for the mean
  delta <- rnorm(1, mean = -(sigma^2) / 2, sd = sigma)

  # Residual (year y) in log-space with temporal autocorrelation
  eps <- rho * last.eps + sqrt(1 - rho^2) * delta

  # ln(recruits/spawner) with autocorrelated residuals
  lnrs <- log(alpha) - beta * spawn + eps

  # Recruits
  rec <- exp(lnrs) * spawn

  list(rec = rec, eps = eps)
}


#' Harvest selectivity based on fish size and mesh size
#'
#' Internal function to calculate the proportion of fish retained by a gillnet
#' of a given mesh size, using the selectivity model described by Bromaghin (2005).
#'
#' @param size Numeric. Fish size in millimeters.
#' @param meshsize Numeric. Stretch mesh size in inches.
#' @param s Numeric. Standard deviation of selectivity, typically \code{0.204}.
#' @return Numeric. Proportion of fish of the given size retained by the net.
#' @references
#' Bromaghin, J. F. (2005). The statistical evaluation of relative gillnet selectivity
#' estimates and their application to the management of Alaskan sockeye salmon.
#' \emph{Alaska Fishery Research Bulletin}, 11(2), 258–266.
#' @examples
#' .calc_selectivity(size = 550, meshsize = 5.0, s = 0.204)
#' @keywords internal
.calc_selectivity <- function(size, meshsize, s) {

  t <- 1.920  # Location parameter from Bromaghin (2005)
  h <- 0.622
  l <- -0.547
  perimeter <- 25.4 * meshsize * 2
  x <- size / perimeter

  select <- (1 + l^2 / (4 * h^2))^h *
    (1 + (x - s * l / (2 * h) - t)^2 / (s^2))^-h *
    exp(-l * (atan((x - s * l / (2 * h) - t) / s) + atan(l / (2 * h))))

  select
}


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
#'         e.g. \code{"selectivity == 'unselective' & factorMSY %in% c(0.75, 1.50)"}.
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
#'   the allowed combinations: for "unselective", \code{factorMSY} in
#'   \code{c(1.00, 0.75, 1.50)}; for "8.5 inch gillnet" or "6 inch gillnet",
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
.select_scenarios <- function(selector = "all", enforce_constraints = TRUE) {
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
    key <- tolower(trimws(selector))

    if (nzchar(key) == FALSE) {
      stop("Empty selector string is not valid.")
    }

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

  stop("Invalid `selector`: must be one of 'all', 'Ohlberger', 'Kuskokwim', a filter expression string, or a numeric vector of scenario IDs or row indices.")
}
