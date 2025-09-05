#' Convert a legacy Ohlberger .RData environment to a package-compliant `ssi_run` object.
#'
#' This function loads a saved legacy `.RData` file (the whole run environment)
#' and converts it into a strict `ssi_run` object with the correct top-level
#' keys: **meta**, **scenarios**, **results**, and **parameters**.
#'
#' The sublists under `$results` are direct copies of the legacy objects,
#' with only the scenario level renamed to `scen_1, scen_2, …`.
#'
#' @param rdata_file Character. Path to the legacy `.RData` file containing
#'   the saved environment from a run (must include objects like `data.list`,
#'   `obs.list`, `egg.list`, `fec.list`, `sr_sim.list`, `parms.list`, etc.).
#' @param env_yaml Optional. Path to a YAML file with environment metadata
#'   (e.g. RNG kind, BLAS threads, pkg version) to enrich `meta` and `parameters`.
#'
#' @return A list of class `ssi_run` with the following top-level components:
#'   - `meta`       : list with run-level metadata (timestamp, niter, nscen, etc.)
#'   - `scenarios`  : tibble of scenario metadata (with `scen_num`, no `No`)
#'   - `results`    : list with sublists
#'       (`para.list`, `sr_sim.list`, `fec.list`, `egg.list`, `S_msy.list`,
#'        `data.list`, `obs.list`, `ret_age.list`, `meanSaA.list`,
#'        `propfemale.list`, `select_age.list`, `MSY_Goals`, `impl_errors.list`)
#'   - `parameters` : list of run parameters (ny, nyh, nyi, goalfreq, firstrev,
#'       review_years, seednum, …)
#'
#' @export
convert_legacy_to_ssi <- function(rdata_file, env_yaml = NULL) {
  if (!file.exists(rdata_file)) stop("File not found: ", rdata_file, call. = FALSE)

  # ---- load into isolated environment
  legacy_env <- new.env(parent = emptyenv())
  load(rdata_file, envir = legacy_env)
  g <- as.list(legacy_env)

  `%||%` <- function(x,y) if (is.null(x)) y else x

  # core objects
  data_list <- g[["data.list"]] %||% g[["data"]]
  parms     <- g[["parms.list"]] %||% g[["parms"]]
  scen_df   <- g[["scen"]] %||% g[["scenarios"]]
  if (is.null(data_list) || !length(data_list)) stop("legacy$data.list missing/empty.", call.=FALSE)
  if (is.null(parms)) stop("legacy$parms.list is required in the .RData file.", call.=FALSE)

  nscen <- length(data_list)
  scen_names <- paste0("scen_", seq_len(nscen))

  # helper: rename scenario level
  rename_scenarios <- function(x) {
    if (is.null(x)) return(setNames(vector("list", nscen), scen_names))
    if (!is.list(x)) return(x)
    names(x) <- scen_names
    x
  }

  # ---- results: copy + rename only
  search_map <- list(
    "para.list"       = c("para.list","para"),
    "sr_sim.list"     = c("sr_sim.list","sr_sim","sr.list"),
    "fec.list"        = c("fec.list","fec"),
    "egg.list"        = c("egg.list","egg"),
    "S_msy.list"      = c("S_msy.list","S_msy"),
    "data.list"       = c("data.list","data"),
    "obs.list"        = c("obs.list","obs"),
    "ret_age.list"    = c("ret_age.list","ret_by_age"),
    "meanSaA.list"    = c("meanSaA.list","meanSaA"),
    "propfemale.list" = c("propfemale.list","propfemale"),
    "select_age.list" = c("select_age.list","selectivity.list","selectivity"),
    "MSY_Goals"       = c("MSY_Goals","MSY_Goals.list"),
    "impl_errors.list"= c("impl_errors.list","impl_errors")
  )

  results <- list()
  for (target_name in names(search_map)) {
    cand_keys <- search_map[[target_name]]
    obj <- NULL
    for (k in cand_keys) { if (!is.null(g[[k]])) { obj <- g[[k]]; break } }
    results[[target_name]] <- rename_scenarios(obj)
  }

  # enforce order and fill missing
  required_order <- names(search_map)
  for (nm in required_order) if (is.null(results[[nm]])) {
    results[[nm]] <- setNames(vector("list", nscen), scen_names)
  }
  results <- results[required_order]

  # ---- scenarios: copy, drop No
  scenarios <- scen_df %||% data.frame(scen_num = seq_len(nscen))
  if (is.data.frame(scenarios)) {
    if ("No" %in% names(scenarios)) scenarios$No <- NULL
    if (!"scen_num" %in% names(scenarios)) {
      scenarios$scen_num <- seq_len(NROW(scenarios))
    }
  }
  scenarios <- tibble::as_tibble(scenarios, .name_repair = "minimal")

  # ---- meta
  niter <- max(vapply(data_list, function(s) if (is.list(s)) length(s) else 0L, integer(1)))
  meta <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
    niter     = as.integer(niter),
    nscen     = as.integer(nscen),
    seed_mode = if (!is.null(parms$seednum)) "reproducible" else "unknown"
  )
  if (!is.null(env_yaml) && requireNamespace("yaml", quietly = TRUE)) {
    env <- try(yaml::read_yaml(env_yaml), silent = TRUE)
    if (is.list(env)) {
      if (!is.null(env$rng_kind))     meta$rng_kind     <- env$rng_kind
      if (!is.null(env$blas_threads)) meta$blas_threads <- env$blas_threads
      if (!is.null(env$pkg_version))  meta$pkg_version  <- env$pkg_version
    }
  }

  # ---- parameters
  ny       <- as.integer(parms$ny   %||% NA)
  nyh      <- as.integer(parms$nyh  %||% ny)
  nyi      <- as.integer(parms$nyi  %||% 0L)
  goalfreq <- as.integer(parms$goalfreq %||% NA)
  firstrev <- as.integer(parms$firstrev %||% 1L)
  seednum  <- as.integer(parms$seednum %||% NA)
  review_years <- if (!is.na(ny) && !is.na(goalfreq) && !is.na(firstrev)) {
    seq.int(nyi + firstrev, ny, by = goalfreq)
  } else NA

  parameters <- list(
    ny           = ny,
    nyh          = nyh,
    nyi          = nyi,
    goalfreq     = goalfreq,
    firstrev     = firstrev,
    review_years = review_years,
    seednum      = seednum
  )

  # ---- assemble
  out <- list(
    meta       = meta,
    scenarios  = scenarios,
    results    = results,
    parameters = parameters
  )
  class(out) <- c("ssi_run", class(out))
  out
}
