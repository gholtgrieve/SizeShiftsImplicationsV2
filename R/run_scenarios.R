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
#' @param seed Character. Seeding strategy across the scenario–iteration grid.
#'   Must be exactly `"reproducible"` (seeds are `1:(niter*nscen)`) or `"random"`
#'   (seeds are sampled from `1:1000000`). Any other value triggers an error.
#' @param params Either a **named profile** (`"Ohlberger"` or `"Kuskokwim"`) that will be
#'   resolved via [resolve_profile()], or a fully-specified **list** of parameter
#'   values to pass through to [build_config()]. Default `"Ohlberger"`.
#' @param output_dir String. Directory where the `.rds` run artifact and other
#'   outputs will be written. Default `file.path(here::here(), "outputs")`.
#'
#' @return An object of class **`ssi_run`** with components:
#'   \itemize{
#'     \item `meta`: run metadata (timestamp, seeds, niter, nscen, dirs, profile);
#'     \item `scenarios`: the filtered scenario tibble that was run;
#'     \item `results`: nested per-scenario/per-iteration outputs (lists of lists).
#'   }
#'   A copy is also saved as an `.rds` under `output_dir`.
#'
#' @section Files written:
#' Files are written under `output_dir` (created if needed):
#' \itemize{
#'   \item `run_<timestamp>.rds` — serialized `ssi_run` object;
#'   \item `run_<timestamp>_scen.xlsx` — scenario table (via **openxlsx**).
#' }
#'
#' @section Dependencies:
#' Requires packages **here**, **progress**, and **openxlsx** for file paths,
#' progress bar, and Excel export (all imported at call sites with `pkg::fun()`).
#'
#' @seealso [select_scenarios()] for scenario selection,
#'   [run_model()] for the single-simulation routine, [build_config()] and
#'   [resolve_profile()] for configuration.
#'
#' @examples
#' \dontrun{
#' # All allowed scenarios, reproducible seeding
#' run_scenarios("all", niter = 1000, seed = "reproducible")
#'
#' # Subset by filter
#' sel <- "selectivity == 'unselective' & mgmt %in% c('smsy_goal','s_eq_goal') & factorMSY %in% c(0.75, 1.50)"
#' run_scenarios(sel, niter = 1000, seed = "reproducible", params = "Ohlberger")
#'
#' # Specific scenario IDs
#' run_scenarios(c(1, 12, 37), niter = 500, seed = "random", params = "Ohlberger")
#' }
#' @export


run_scenarios <- function(scenarios,
                          niter = 1000,
                          seed = "reproducible",
                          params = "Ohlberger",
                          output_dir = paste0(here::here(), "/outputs")) {

  ## scenarios to run
  scen  <- select_scenarios(selector = scenarios, enforce_constraints = TRUE)
  nscen <- nrow(scen)

  ## seeds
  if (seed == "reproducible") {
    seed.list <- seq_len(niter * nscen)
  } else if (seed == "random") {
    seed.list <- sample.int(1e6, niter * nscen, replace = TRUE)
  } else {
    stop("Please specify 'reproducible' or 'random' for seed.")
  }
  stopifnot(length(seed.list) == niter * nscen)

  message("time estimate: ", round((niter * nscen / 20) / 60, 2), " hours")
  pb <- progress::progress_bar$new(total = niter * nscen)

  ## resolve profile ONCE (either a name in param_configs or a list)
  profile <- resolve_profile(params)

  ##=====================================================## output lists (top-level preallocation)
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

  ## Name by scenario id
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

  start.time <- Sys.time()

  for (j in seq_len(nscen)) {
    scen_row <- scen[j, , drop = FALSE] # single-row data frame for this scenario
    ## second-level pre-allocation: one slot per iteration
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

    for (k in seq_len(niter)) {
      idx     <- (j - 1L) * niter + k
      seednum <- seed.list[idx]

      # Build a complete, validated config for this run:
      cfg <- build_config(params  = profile,   # you can pass the profile object here
                          scen_row = scen_row,
                          j        = j,
                          k        = k,
                          seednum  = seednum)

      mod.out <- try(run_model(config = cfg), silent = TRUE)
      if (!inherits(mod.out, "try-error")) {
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
      pb$tick()
    }
  }

  ##====================================================## true run time
  end.time <- Sys.time()
  print(round(end.time - start.time, 2))
  time.save <- format(end.time, "%Y-%m-%d_%H-%M-%S")  # safe timestamp

  # Create output dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Package the run into a single structured object ----
  run_out <- list(
    meta = list(
      timestamp   = time.save,
      niter       = niter,
      nscen       = nscen,
      seed_mode   = seed,
      seed_list   = seed.list,
      params_spec = params,          # name or list the user passed
      output_dir  = output_dir
    ),
    scenarios = scen,                # the filtered scenario table
    results = list(                  # all your per-scenario/per-iter outputs
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

  # Save a single compact file for later use
  rds_path <- file.path(output_dir, paste0("run_", time.save, ".rds"))
  saveRDS(run_out, rds_path)

  # Keep your Excel summary if you like
  openxlsx::write.xlsx(scen, file = file.path(output_dir, paste0("run_", time.save, "_scen.xlsx")))

  # Return the run object so you can pass it directly to make_figures()
  return(run_out)

} #end function

