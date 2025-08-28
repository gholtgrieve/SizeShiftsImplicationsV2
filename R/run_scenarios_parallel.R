
#' Run a batch of Management Strategy Evaluation (MSE) scenarios
#'
#' High-level interface to run the simulation across a set of scenarios and
#' Monte Carlo iterations. This function selects which scenarios to run,
#' seeds each run for reproducibility or randomness, executes the simulation,
#' and saves outputs to disk.
#'
#' @param scenarios Selector passed to [select_scenarios()]. One of:
#'   \itemize{
#'     \item `"all"` — run all allowed scenarios;
#'     \item a single dplyr-style filter **string**, e.g.
#'       `"selectivity == 'unselective' & factorMSY %in% c(0.75, 1.50)"`;
#'     \item an integer vector of scenario IDs, e.g. `c(1, 12, 37)`.
#'   }
#'   See [select_scenarios()] for details.
#' @param niter Integer. Number of Monte Carlo iterations **per scenario**. Default `1000`.
#' @param seed Character. `"reproducible"` (iteration seeds = `1:niter`, shared across scenarios)
#'   or `"random"` (iteration seeds sampled from `1:1e6`, shared across scenarios). Any other value errors.
#' @param params Either a **named profile** (`"Ohlberger"` or `"Kuskokwim"`) for [resolve_profile()],
#'   or a fully-specified **list** passed to [build_config()]. Default `"Ohlberger"`.
#' @param output_dir String. Directory where outputs are written. Default `file.path(here::here(), "outputs")`.
#' @param parallel Logical. Use parallel workers for iterations? Default `TRUE`.
#' @param workers Integer or `NULL`. Number of parallel workers if `parallel = TRUE`.
#'   Default: `min(parallelly::availableCores(), niter)`.
#'
#' @return An **`ssi_run`** object; also saved as `run_<timestamp>.rds` under `output_dir`.
#' @export
run_scenarios <- function(scenarios,
                          niter = 1000,
                          seed = "reproducible",
                          params = "Ohlberger",
                          output_dir = file.path(here::here(), "outputs"),
                          parallel = TRUE,
                          workers = NULL) {

  ## ---- scenarios to run -------------------------------------------------------
  scen  <- select_scenarios(selector = scenarios, enforce_constraints = TRUE)
  nscen <- nrow(scen)

  ## ---- seeds: common random numbers across scenarios -------------------------
  if (seed == "reproducible") {
    seeds_k <- seq_len(niter)                # iteration seeds: 1,2,...,niter
  } else if (seed == "random") {
    seeds_k <- sample.int(1e6, niter, replace = TRUE)
  } else {
    stop("Please specify 'reproducible' or 'random' for seed.")
  }

  ## ---- resolve profile once ---------------------------------------------------
  profile <- resolve_profile(params)

  ## ---- preallocate result lists ----------------------------------------------
  para.list        <- vector("list", nscen)
  sr_sim.list      <- vector("list", nscen)
  fec.list         <- vector("list", nscen)
  egg.list         <- vector("list", nscen)
  S_msy.list       <- vector("list", nscen)
  data.list        <- vector("list", nscen)
  obs.list         <- vector("list", nscen)
  ret_age.list     <- vector("list", nscen)
  meanSaA.list     <- vector("list", nscen)
  propfemale.list  <- vector("list", nscen)
  select_age.list  <- vector("list", nscen)
  MSY_Goals.list   <- vector("list", nscen)
  impl_errors.list <- vector("list", nscen)

  if ("scen_num" %in% names(scen)) {
    nm <- paste0("scen_", scen$scen_num)
    names(para.list)        <- nm
    names(sr_sim.list)      <- nm
    names(fec.list)         <- nm
    names(egg.list)         <- nm
    names(S_msy.list)       <- nm
    names(data.list)        <- nm
    names(obs.list)         <- nm
    names(ret_age.list)     <- nm
    names(meanSaA.list)     <- nm
    names(propfemale.list)  <- nm
    names(select_age.list)  <- nm
    names(MSY_Goals.list)   <- nm
    names(impl_errors.list) <- nm
  }

  ## ---- inner runner (one iteration for one scenario row) ----------------------
  .run_one_iter <- function(k, scen_row, profile) {
    seednum <- seeds_k[k]  # CRN: same seed for iteration k across all scenarios
    cfg <- build_config(
      params   = profile,
      scen_row = scen_row,
      j        = NA_integer_,
      k        = k,
      seednum  = seednum
    )
    # run_model() itself uses withr::with_seed(seednum, { ... }) internally
    try(run_model(config = cfg), silent = TRUE)
  }

  ## ---- parallel plan (if requested and available) ----------------------------
  use_parallel <- isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)
  if (use_parallel) {
    if (is.null(workers)) workers <- min(parallelly::availableCores(), niter)
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
  }

  ## ---- progress bar (progressr) ----------------------------------------------
  # Tip: run once per session to enable a nice bar: progressr::handlers(global = TRUE)
  total_steps <- nscen * niter
  t_start <- Sys.time()

  progressr::with_progress({
    p <- progressr::progressor(steps = total_steps)

    for (j in seq_len(nscen)) {
      scen_row <- scen[j, , drop = FALSE]

      # 2nd-level preallocation per scenario
      para.list[[j]]        <- vector("list", niter)
      sr_sim.list[[j]]      <- vector("list", niter)
      fec.list[[j]]         <- vector("list", niter)
      egg.list[[j]]         <- vector("list", niter)
      S_msy.list[[j]]       <- vector("list", niter)
      data.list[[j]]        <- vector("list", niter)
      obs.list[[j]]         <- vector("list", niter)
      ret_age.list[[j]]     <- vector("list", niter)
      meanSaA.list[[j]]     <- vector("list", niter)
      propfemale.list[[j]]  <- vector("list", niter)
      select_age.list[[j]]  <- vector("list", niter)
      MSY_Goals.list[[j]]   <- vector("list", niter)
      impl_errors.list[[j]] <- vector("list", niter)

      iter_indices <- seq_len(niter)

      if (use_parallel) {
        res_iter <- future.apply::future_lapply(
          iter_indices,
          FUN = function(k) {
            out <- .run_one_iter(k, scen_row, profile)
            p()  # one tick per completed iteration
            out
          },
          future.seed = TRUE
        )
      } else {
        res_iter <- lapply(iter_indices, function(k) {
          out <- .run_one_iter(k, scen_row, profile)
          p()
          out
        })
      }

      # collect results
      for (k in iter_indices) {
        mod.out <- res_iter[[k]]
        if (!inherits(mod.out, "try-error") && !is.null(mod.out)) {
          para.list[[j]][[k]]        <- mod.out$para
          sr_sim.list[[j]][[k]]      <- mod.out$sr_sim
          fec.list[[j]][[k]]         <- mod.out$fec
          egg.list[[j]][[k]]         <- mod.out$egg
          S_msy.list[[j]][[k]]       <- mod.out$S_msy
          data.list[[j]][[k]]        <- mod.out$data
          obs.list[[j]][[k]]         <- mod.out$obs
          ret_age.list[[j]][[k]]     <- mod.out$ret_by_age
          meanSaA.list[[j]][[k]]     <- mod.out$meanSaA
          propfemale.list[[j]][[k]]  <- mod.out$propfemale
          select_age.list[[j]][[k]]  <- mod.out$selectivities_by_age
          MSY_Goals.list[[j]][[k]]   <- mod.out$MSY_Goals
          impl_errors.list[[j]][[k]] <- mod.out$impl_errors
        }
      }
    } # end scenarios loop
  }) # end with_progress

  t_end <- Sys.time()
  message("Elapsed time: ", round(t_end - t_start, 2))

  ## ---- package run object -----------------------------------------------------
  time.save <- format(t_end, "%Y-%m-%d_%H-%M-%S")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  run_out <- list(
    meta = list(
      timestamp   = time.save,
      niter       = niter,
      nscen       = nscen,
      seed_mode   = seed,
      seed_list   = seeds_k,
      params_spec = params,
      output_dir  = output_dir
    ),
    scenarios = scen,
    results = list(
      para         = para.list,
      sr_sim       = sr_sim.list,
      fec          = fec.list,
      egg          = egg.list,
      S_msy        = S_msy.list,
      data         = data.list,
      obs          = obs.list,
      ret_by_age   = ret_age.list,
      meanSaA      = meanSaA.list,
      propfemale   = propfemale.list,
      selectivity  = select_age.list,
      MSY_Goals    = MSY_Goals.list,
      impl_errors  = impl_errors.list
    )
  )
  class(run_out) <- c("ssi_run", "list")

  rds_path <- file.path(output_dir, paste0("run_", time.save, ".rds"))
  saveRDS(run_out, rds_path)
  openxlsx::write.xlsx(scen, file = file.path(output_dir, paste0("run_", time.save, "_scen.xlsx")))

  return(invisible(run_out))
}
