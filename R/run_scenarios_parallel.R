#' Run a batch of Management Strategy Evaluation (MSE) scenarios
#'
#' High-level interface to run the simulation across a set of scenarios and
#' Monte Carlo iterations. Selects scenarios, seeds iterations (CRN across
#' scenarios), runs the model (optionally in parallel), and saves outputs.
#'
#' @param scenarios See [.select_scenarios()] (e.g., "all", filter string, or IDs).
#' @param niter Integer. Iterations per scenario. Default 1000.
#' @param seed Character. `"reproducible"` seeds each (scenario, iteration) pair
#'   deterministically (base seeds 1, 2, ..., niter); two runs with the **same
#'   backend** (both parallel or both sequential) will produce bitwise-identical
#'   results. `"random"` draws random base seeds so results differ between runs.
#'   **Parallel and sequential runs with the same `seed` value will differ
#'   numerically** due to different seednum offsets; see Details.
#' @param params Profile name ("Ohlberger","Kuskokwim") or full list for [build_config()].
#' @param output_dir Directory to write outputs (default: here()/outputs).
#' @param parallel Logical. Run scenarios in parallel. Default = TRUE.
#'   When enabled, each worker handles one scenario.
#' @param workers Integer or NULL. Number of parallel workers. If NULL
#'   (the default), the function uses \code{parallelly::availableCores(omit=c("system","fallback")) - 1},
#'   capped at the number of scenarios and never below 1. This leaves one core free
#'   for the main R session / operating system.
#'
#' @details
#' ## Reproducibility
#'
#' All model random draws are wrapped in `withr::with_seed(seednum, ...)` where
#' `seednum` is a deterministic function of the scenario index `j` and iteration
#' index `k`:
#' * **Sequential** (`parallel = FALSE`): `seednum = seeds_k[k]`
#' * **Parallel** (`parallel = TRUE`): `seednum = seeds_k[k] + 10000 * j`
#'
#' Using `seed = "reproducible"` guarantees that two runs with the **same
#' backend** produce bitwise-identical results. Parallel and sequential runs
#' with the same `seed` value **will differ numerically** because the seednum
#' offsets differ. This is intentional: the `+10000 * j` offset in parallel
#' mode ensures no two workers share an RNG subsequence even when base seeds
#' overlap across scenarios.
#'
#' `future.seed = TRUE` is retained as a safety net for any future code changes
#' that introduce un-seeded draws outside `withr::with_seed`; it does not
#' currently contribute to model reproducibility.
#'
#' The seed strategy for each run is recorded in `meta$seed_strategy`
#' (`"CRN"` for sequential, `"per_scenario_iter"` for parallel) and the base
#' seeds are stored in `meta$seed_list`.
#' @export
run_scenarios <- function(scenarios,
                          niter      = 1000,
                          seed       = "reproducible",
                          params     = "Ohlberger",
                          output_dir = file.path(here::here(), "outputs"),
                          parallel   = TRUE,
                          workers    = getOption("ssi.workers", NULL)) {
    ## --- Select scenarios
  scen  <- .select_scenarios(selector = scenarios, enforce_constraints = TRUE)
  nscen <- nrow(scen)
  if (nscen == 0L) stop("No scenarios selected.")

  ## --- Base seeds per iteration (used directly for CRN or offset for parallel)
  seeds_k <- switch(seed,
                    reproducible = seq_len(niter),
                    random       = sample.int(1e6, niter, replace = TRUE),
                    stop("Please specify 'reproducible' or 'random' for `seed`.")
  )

  ## --- Parallel plan (by scenario)
  use_future <- isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)

  if (use_future) {
    # Save current thread settings so they can be restored when the function exits.
    old_omp  <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
    old_mkl  <- Sys.getenv("MKL_NUM_THREADS", unset = NA_character_)
    old_blas <- if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_get_num_procs()
    } else {
      NULL
    }
    on.exit({
      if (is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS") else
        Sys.setenv(OMP_NUM_THREADS = old_omp)
      if (is.na(old_mkl)) Sys.unsetenv("MKL_NUM_THREADS") else
        Sys.setenv(MKL_NUM_THREADS = old_mkl)
      if (!is.null(old_blas) &&
          requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(old_blas)
      }
    }, add = TRUE)

    # Avoid BLAS over-subscription when forking/spawning parallel R sessions
    # (This is safe to do repeatedly; packages can still override if they want)
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    }

    # Compute default: (cores - 1), capped by nscen, min 1
    workers <- .auto_workers(nscen = nrow(scen), requested = workers)

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    # multisession works across OSes (RStudio safe). If you ever want multicore on
    # non-RStudio Unix, you can switch based on interactive/session checks.
    future::plan(future::multisession, workers = workers)
    message(sprintf("Parallel ON with %d worker(s).", future::nbrOfWorkers()))
  } else {
    # Ensure sequential fallback is clear and reproducible
    message("Parallel OFF (running sequentially).")
  }

  # --- Logging directory (shared across workers) --------------------------------
  log_dir <- file.path(output_dir, "logs")
  options(ssi.log_dir = log_dir)  # used by fun_log_dir()
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  # ------------------------------------------------------------------------------


  # Extract config parameters for recording later
  profile0 <- resolve_profile(params)  # or however you turn "Ohlberger" into a config list

  # Normalize possible field names coming from build_config()/param_configs
  nyi      <- profile0$nyi
  nyh      <- profile0$nyh %||% profile0$nyh_len
  ny       <- profile0$ny
  goalfreq <- profile0$goalfreq %||% profile0$goal_freq
  firstrev <- profile0$firstrev %||% 20L
  review_years  <- seq(from = nyi + firstrev, to = ny, by = goalfreq) # derived
  hist_end_year <- nyi + nyh # derived

  missing <- c(
    if (is.null(nyi)) "nyi",
    if (is.null(nyh)) "nyh/nyh_len",
    if (is.null(ny)) "ny",
    if (is.null(goalfreq)) "goalfreq/goal_freq",
    if (is.null(firstrev)) "firstrev/first_rev"
  )
  if (length(missing)) {
    stop("Param profile is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  review_years  <- seq(from = nyi + firstrev, to = ny, by = goalfreq)
  if (!length(review_years)) {
    stop("Computed review_years is empty. Check nyi/ny/goalfreq/firstrev.", call. = FALSE)
  }
  hist_end_year <- nyi + nyh


  ## --- Progress setup
  total_steps <- nscen * niter
  # global handler conflicts with testthat's handler stack; skip silently if needed
  try(progressr::handlers(global = TRUE), silent = TRUE)
  # If you want a simple console bar, uncomment:
  # progressr::handlers("txtprogressbar")

  start.time <- Sys.time()

  ## --- Run scenarios (each worker owns one whole scenario)
  res_all <- progressr::with_progress({
    p <- progressr::progressor(steps = total_steps)

    run_one_scenario <- function(scen_row, j, niter, seeds_k, params, p) {
      # Resolve profile inside worker (keeps exported globals tiny)
      profile <- resolve_profile(params)
      out <- vector("list", niter)
      for (k in seq_len(niter)) {
        # Abandon CRN when parallel (avoid cross-worker randomness collisions)
        seednum <- if (isTRUE(parallel)) {
          seeds_k[k] + 10000 * j
        } else {
          seeds_k[k]
        }

        cfg <- build_config(
          params   = profile,
          scen_row = scen_row,
          j        = j,
          k        = k,
          seednum  = seednum,
          log_dir = log_dir
        )

        # run_model() should use withr::with_seed(seednum, { ... }) internally
        out[[k]] <- try(run_model(cfg), silent = FALSE)
        p() # tick progress once per iteration
      }
      out
    }

    if (use_future && future::nbrOfWorkers() > 1L) {
      future.apply::future_lapply(
        X = seq_len(nscen),
        FUN = function(j) {
          scen_row <- scen[j, , drop = FALSE]
          run_one_scenario(scen_row, j, niter, seeds_k, params, p)
        },
        future.seed     = TRUE,
        future.packages = "SizeShiftsImplicationsV2"
      )
    } else {
      # Fallback to sequential if no parallel workers available
      lapply(seq_len(nscen), function(j) {
        scen_row <- scen[j, , drop = FALSE]
        run_one_scenario(scen_row, j, niter, seeds_k, params, p)
      })
    }
  })

  ## --- Collect results
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
    names(para.list)        <- nm; names(sr_sim.list)   <- nm
    names(fec.list)         <- nm; names(egg.list)      <- nm
    names(S_msy.list)       <- nm; names(data.list)     <- nm
    names(obs.list)         <- nm; names(ret_age.list)  <- nm
    names(meanSaA.list)     <- nm; names(propfemale.list) <- nm
    names(select_age.list)  <- nm; names(MSY_Goals.list)  <- nm
    names(impl_errors.list) <- nm
  }

  for (j in seq_len(nscen)) {
    res_iter <- res_all[[j]]
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

    n_fail_j <- sum(vapply(res_iter, inherits, logical(1L), "try-error"))
    if (n_fail_j > 0L) {
      warning(sprintf(
        "Scenario %d: %d of %d iteration(s) failed and were excluded from results.",
        j, n_fail_j, length(res_iter)
      ), call. = FALSE)
    }
  }

  ## --- Save
  end.time  <- Sys.time()
  time.save <- format(end.time, "%Y-%m-%d_%H-%M-%S")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  run_out <- list(
    meta = list(
      timestamp      = time.save,
      niter          = niter,
      nscen          = nscen,
      seed_mode      = seed,              # "reproducible" or "random"
      seed_strategy  = if (isTRUE(parallel)) "per_scenario_iter" else "CRN",
      seed_list      = seeds_k,           # base seeds used per iteration
      parallel       = isTRUE(parallel),
      workers        = workers,
      plan           = if (isTRUE(parallel)) as.character(attr(future::plan(), "call")) else "sequential",
      params_spec    = params,
      output_dir     = output_dir
    ),
    parameters = list(
      nyi          = nyi,
      nyh          = nyh,
      ny           = ny,
      goalfreq     = goalfreq,
      firstrev     = firstrev,     # hard coded; recorded for transparency
      review_years = review_years,
      hist_end_year = hist_end_year,
      # provenance
      rng_kind     = paste(utils::capture.output(RNGkind()), collapse = " "),
      blas_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)),
      pkg_version  = as.character(utils::packageVersion("SizeShiftsImplicationsV2")),
      timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      scenario_grid_hash = if (requireNamespace("digest", quietly = TRUE))
        digest::digest(scen) else NA_character_
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
  .validate_run_params(run_out$parameters)
  class(run_out) <- c("ssi_run", "list")

  rds_path <- file.path(output_dir, paste0("run_", time.save, ".rds"))
  saveRDS(run_out, rds_path)
  openxlsx::write.xlsx(scen, file = file.path(output_dir, paste0("run_", time.save, "_scen.xlsx")))

  # --- Merge per-worker log shards into final CSVs ------------------------------
  if (exists("fun_merge_all_logs", mode = "function")) {
    try(fun_merge_all_logs(log_dir, remove_shards = TRUE), silent = TRUE)
  }
  # ------------------------------------------------------------------------------

  message(sprintf("Completed %d scenarios \u00d7 %d iterations in %s",
                  nscen, niter, format(end.time - start.time)))
  invisible(run_out)
}
