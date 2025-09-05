# Internal helpers for run-model()


#' Generate age composition of brood year recruits
#'
#' Internal function to simulate the age distribution of recruits from a given brood year
#' using a normal-distribution-based Dirichlet draw. Accepts an optional RNG seed for reproducibility.
#'
#' @param ages Numeric vector. All simulated age classes.
#' @param recruits Numeric. Number of recruits from that brood year.
#' @param meanage Numeric. Mean age of recruits.
#' @param sdage Numeric. Standard deviation of the age distribution.
#' @param seed Optional integer. If provided, sets the random seed for reproducible draws.
#'
#' @return Numeric vector. Number of recruits in each age class (same length as `ages`).
#'
#' @keywords internal
#'
#' @examples
#' ages <- 2:6
#' .calc_agecomp(ages, recruits = 100, meanage = 4, sdage = 0.8, seed = 123)
.calc_agecomp <- function(ages, recruits, meanage, sdage, seed = NULL) {
  if (!is.null(seed)) {
    withr::with_seed(seed, {
      probs_a <- dnorm(ages, meanage, sdage)
      probs_a <- probs_a / sum(probs_a)
      d <- 50
      probs <- as.vector(gtools::rdirichlet(1, probs_a * d))
      round(probs * recruits)
    })
  } else {
    probs_a <- dnorm(ages, meanage, sdage)
    probs_a <- probs_a / sum(probs_a)
    d <- 50
    probs <- as.vector(gtools::rdirichlet(1, probs_a * d))
    round(probs * recruits)
  }
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
#' temporal autocorrelation in residuals. Accepts an optional RNG seed for reproducibility.
#'
#' @param spawn Numeric. Number of spawners in the current year.
#' @param sigma Numeric. Standard deviation of the random error term.
#' @param alpha Numeric. Productivity parameter of the Ricker model.
#' @param beta Numeric. Density-dependence parameter of the Ricker model.
#' @param rho Numeric. Temporal autocorrelation coefficient (between 0 and 1).
#' @param last.eps Numeric. Recruitment residual (in log-space) from the previous year.
#' @param seed Optional integer. If provided, sets the random seed for reproducible draws.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{rec}}{Numeric. Recruits in the current year.}
#'   \item{\code{eps}}{Numeric. Residual (in log-space) for the current year.}
#' }
#'
#' @keywords internal
#'
#' @examples
#' .calc_ricker(spawn = 100, sigma = 0.3, alpha = 5, beta = 0.01,
#'              rho = 0.5, last.eps = 0.1, seed = 456)
.calc_ricker <- function(spawn, sigma, alpha, beta, rho, last.eps, seed = NULL) {
  if (!is.null(seed)) {
    withr::with_seed(seed, {
      delta <- rnorm(1, mean = -(sigma^2) / 2, sd = sigma)
      eps <- rho * last.eps + sqrt(1 - rho^2) * delta
      lnrs <- log(alpha) - beta * spawn + eps
      rec <- exp(lnrs) * spawn
      list(rec = rec, eps = eps)
    })
  } else {
    delta <- rnorm(1, mean = -(sigma^2) / 2, sd = sigma)
    eps <- rho * last.eps + sqrt(1 - rho^2) * delta
    lnrs <- log(alpha) - beta * spawn + eps
    rec <- exp(lnrs) * spawn
    list(rec = rec, eps = eps)
  }
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
#' Filter the package's built-in scenarios table using a single selector argument.
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


# Internal: required fields for a valid ssi_run
.required_run_params <- c(
  # timeline / reviews
  "nyi", "nyh", "ny", "goalfreq", "firstrev", "review_years",
  # reproducibility / provenance
  "rng_kind", "blas_threads",
  "pkg_version", "timestamp",
  # (optional but recommended) scenario knobs used to build scenarios
  "scenario_grid_hash"
)


# Internal: validate parameters are present and consistent
.validate_run_params <- function(params) {
  miss <- setdiff(.required_run_params, names(params %||% list()))
  if (length(miss)) {
    stop("Missing required run parameters: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  # Basic consistency checks
  if (!is.numeric(params$nyi) || !is.numeric(params$nyh) || !is.numeric(params$ny)) {
    stop("nyi/nyh/ny must be numeric.", call. = FALSE)
  }
  if (!is.numeric(params$goalfreq) || params$goalfreq <= 0) {
    stop("goalfreq must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(params$firstrev)) stop("firstrev must be numeric.", call. = FALSE)

  expected_reviews <- seq(params$nyi + params$firstrev, params$ny, by = params$goalfreq)
  if (!isTRUE(all.equal(unname(expected_reviews), unname(params$review_years)))) {
    stop("`review_years` does not match seq(nyi + firstrev, ny, by = goalfreq).", call. = FALSE)
  }
  invisible(TRUE)
}


#' Upgrade a legacy ssi_run by injecting required parameters
#'
#' Use this when an older `ssi_run` is missing timeline/review parameters.
#' It records nyi, nyh, ny, goalfreq, firstrev, computed review_years,
#' and hist_end_year, plus light provenance for reproducibility.
#'
#' @param run       legacy `ssi_run` object
#' @param nyi       integer(1) initial years (e.g., 10)
#' @param nyh       integer(1) historical length in years (e.g., 50)
#' @param ny        integer(1) total simulated years (e.g., 110)
#' @param goalfreq  integer(1) review cadence in years (e.g., 10)
#' @param firstrev  integer(1) first review offset after nyi (e.g., 20)
#' @param seednum   optional integer seed recorded for provenance
#' @return          upgraded `ssi_run` (invisible)
#' @keywords internal
.upgrade_ssi_run_params <- function(run, nyi, nyh, ny, goalfreq, firstrev, seednum = NULL) {
  stopifnot(inherits(run, "ssi_run"))

  # ---- validate inputs ------------------------------------------------------
  req <- list(nyi = nyi, nyh = nyh, ny = ny, goalfreq = goalfreq, firstrev = firstrev)
  bad <- vapply(req, function(x) !is.numeric(x) || length(x) != 1L || !is.finite(x), logical(1))
  if (any(bad)) {
    stop("All of nyi, nyh, ny, goalfreq, firstrev must be numeric scalars.", call. = FALSE)
  }
  if (goalfreq <= 0) stop("`goalfreq` must be > 0.", call. = FALSE)
  if (nyi < 0 || nyh <= 0 || ny <= 0) stop("`nyi`>=0, `nyh`>0, `ny`>0 required.", call. = FALSE)

  # ---- compute derived values ----------------------------------------------
  review_years  <- seq(from = nyi + firstrev, to = ny, by = goalfreq)
  if (!length(review_years)) {
    stop("Computed `review_years` is empty. Check nyi/ny/goalfreq/firstrev.", call. = FALSE)
  }
  hist_end_year <- nyi + nyh

  # ---- light provenance (best-effort; no hard deps) -------------------------
  rng_kind    <- try(paste(utils::capture.output(RNGkind()), collapse = " "), silent = TRUE)
  rng_kind    <- if (inherits(rng_kind, "try-error")) NA_character_ else rng_kind
  blas_threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)))
  pkg_version <- try(as.character(utils::packageVersion("SizeShiftsImplicationsV2")), silent = TRUE)
  pkg_version <- if (inherits(pkg_version, "try-error")) NA_character_ else pkg_version
  scen_hash   <- try({
    if (requireNamespace("digest", quietly = TRUE)) digest::digest(run$scenarios) else NA_character_
  }, silent = TRUE)
  scen_hash   <- if (inherits(scen_hash, "try-error")) NA_character_ else scen_hash

  # ---- merge into run$parameters (preserve existing where present) ----------
  if (is.null(run$parameters)) run$parameters <- list()
  add_if_null <- function(x, name, value) {
    if (is.null(x[[name]])) x[[name]] <- value
    x
  }

  run$parameters$nyi         <- nyi
  run$parameters$nyh         <- nyh
  run$parameters$ny          <- ny
  run$parameters$goalfreq    <- goalfreq
  run$parameters$firstrev    <- firstrev
  run$parameters$review_years <- review_years
  run$parameters$hist_end_year <- hist_end_year

  run$parameters <- add_if_null(run$parameters, "seednum", seednum)
  run$parameters <- add_if_null(run$parameters, "rng_kind", rng_kind)
  run$parameters <- add_if_null(run$parameters, "blas_threads", blas_threads)
  run$parameters <- add_if_null(run$parameters, "pkg_version", pkg_version)
  run$parameters <- add_if_null(run$parameters, "timestamp", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  run$parameters <- add_if_null(run$parameters, "scenario_grid_hash", scen_hash)

  # ---- internal consistency checks -----------------------------------------
  expected_reviews <- seq(nyi + firstrev, ny, by = goalfreq)
  if (!isTRUE(all.equal(unname(expected_reviews), unname(run$parameters$review_years)))) {
    stop("`review_years` mismatch with seq(nyi + firstrev, ny, by = goalfreq).", call. = FALSE)
  }
  if (!is.numeric(run$parameters$hist_end_year) || length(run$parameters$hist_end_year) != 1L) {
    stop("`hist_end_year` must be a numeric scalar.", call. = FALSE)
  }

  # Reassert class tag
  class(run) <- unique(c("ssi_run", class(run)))

  invisible(run)
}


# Pick a sensible default worker count: (available cores - 1), capped by nscen, min 1.
.auto_workers <- function(nscen, requested = NULL) {
  # If the user explicitly provided workers, sanitize & cap by scenarios
  if (!is.null(requested) && is.finite(requested)) {
    return(as.integer(max(1L, min(nscen, requested))))
  }

  # Auto-detect cores robustly in local & scheduler environments
  cores <- 1L
  cores <- tryCatch(
    parallelly::availableCores(omit = c("system", "fallback")),
    error = function(e) 1L
  )
  # Leave one core for the main/session/OS; keep at least 1
  as.integer(max(1L, min(nscen, cores - 1L)))
}


