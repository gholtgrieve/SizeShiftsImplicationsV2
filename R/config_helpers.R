# R/config_helpers.R

#' Minimal defaults layered under profiles
#' @keywords internal
default_config <- function() {
  list(
    # filled per run:
    j = NA_integer_, k = NA_integer_, seednum = NA_integer_,
    # optional diagnostics:
    print = FALSE, plot = FALSE
  )
}

#' Validate that all required fields exist for run_model()
#' @keywords internal

validate_config <- function(cfg) {
  required <- c(
    "agediff","agetrend","ages","allometry","alpha_mean","alt_sr_param",
    "beta_mean","factorMSY","futureT","goalfreq","harverr","harvmgmt",
    "harvrate","hobserr","maxsel","meanageini","ny","nyh","nyi","obserr",
    "ocean0s","procerr","propF","propFtrend","reglength","regstr","rho",
    "ricker_type","sdSaA","sdage","sdsel","seednum","sim_recruits",
    "sizetrends","sr_corr","sr_parms_sd","var","vonB_k","vonB_Linf"
  )
  miss <- setdiff(required, names(cfg))
  if (length(miss)) {
    stop("Parameter set is missing required field(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  cfg
}

#' Apply scenario-driven overrides to a profile config
#' @param cfg A config list (typically default_config() %<-% profile)
#' @param scen_row One-row data.frame/tibble from scenarios
#' @keywords internal
apply_scenario_overrides <- function(cfg, scen_row) {
  if (!is.data.frame(scen_row) || nrow(scen_row) != 1L) {
    stop("`scen_row` must be a one-row data.frame/tibble.", call. = FALSE)
  }
  # helpers to coerce safely
  chr1 <- function(x) as.character(x[[1]])
  num1 <- function(x) as.numeric(x[[1]])

  # direct fields from scenarios
  if ("mgmt"      %in% names(scen_row)) cfg$harvmgmt  <- chr1(scen_row$mgmt)
  if ("factorMSY" %in% names(scen_row)) cfg$factorMSY <- num1(scen_row$factorMSY)
  if ("futureT"   %in% names(scen_row)) cfg$futureT   <- chr1(scen_row$futureT)
  if ("sdsel"     %in% names(scen_row)) cfg$sdsel     <- num1(scen_row$sdsel)
  if ("maxsel"    %in% names(scen_row)) cfg$maxsel    <- num1(scen_row$maxsel)

  # switches multiply profile baselines (0/1 flags typically)
  if ("ageT"  %in% names(scen_row)) cfg$agetrend   <- cfg$agetrend   * num1(scen_row$ageT)
  if ("sizeT" %in% names(scen_row)) cfg$sizetrends <- cfg$sizetrends * num1(scen_row$sizeT)
  if ("sexT"  %in% names(scen_row)) cfg$propFtrend <- cfg$propFtrend * num1(scen_row$sexT)

  cfg
}

#' Resolve a parameter profile (by name in param_configs or as a list)
#' @param params Character scalar (name in `param_configs`) or a list
#' @keywords internal
resolve_profile <- function(params) {
  if (is.character(params) && length(params) == 1L) {
    avail <- names(param_configs)
    if (!(params %in% avail)) {
      stop("`params` must be one of ", paste(shQuote(avail), collapse = ", "),
           " or a list.", call. = FALSE)
    }
    param_configs[[params]]
  } else if (is.list(params)) {
    params
  } else {
    stop("`params` must be a profile name (in `param_configs`) or a list.",
         call. = FALSE)
  }
}

#' Build a ready-to-run config for run_model()
#' @param params Profile name or list (see resolve_profile())
#' @param scen_row One-row scenarios slice to apply
#' @param j Scenario index (integer)
#' @param k Iteration index (integer)
#' @param seednum Integer seed for the run
#' @return A validated config list
#' @keywords internal
build_config <- function(params, scen_row, j, k, seednum, log_dir) {
  profile <- resolve_profile(params)
  cfg <- utils::modifyList(default_config(), profile)
  cfg <- apply_scenario_overrides(cfg, scen_row)
  cfg$j <- as.integer(j)
  cfg$k <- as.integer(k)
  cfg$seednum <- as.integer(seednum)
  cfg$log_dir <- log_dir    # ensure run_model() logging goes to outputs/logs
  validate_config(cfg)
}
