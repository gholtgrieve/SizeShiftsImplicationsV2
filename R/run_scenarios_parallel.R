#' Run a batch of Management Strategy Evaluation (MSE) scenarios
#'
#' High-level interface to run the simulation across a set of scenarios and
#' Monte Carlo iterations. Selects scenarios, seeds iterations (CRN across
#' scenarios), runs the model (optionally in parallel), and saves outputs.
#'
#' @param scenarios See [.select_scenarios()] (e.g., "all", filter string, or IDs).
#' @param niter Integer. Iterations per scenario. Default 1000.
#' @param seed "reproducible" (1:niter) or "random" (sampled seeds).
#' @param params Profile name ("Ohlberger","Kuskokwim") or full list for [build_config()].
#' @param output_dir Directory to write outputs (default: here()/outputs).
#' @param parallel Logical. Run scenarios in parallel (one worker per scenario).
#' @param workers Integer or NULL. Number of workers; default auto-detect.
#' @export
run_scenarios <- function(scenarios,
                          niter      = 1000,
                          seed       = "reproducible",
                          params     = "Ohlberger",
                          output_dir = file.path(here::here(), "outputs"),
                          parallel   = TRUE,
                          workers    = NULL) {

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
    # pick a safe default if user didn't supply workers
    if (is.null(workers) || !is.numeric(workers) || length(workers) != 1L || workers < 1L) {
      workers <- tryCatch(parallelly::availableCores(omit = c("system","fallback")),
                          error = function(e) 1L)
      workers <- max(1L, min(workers - 1L, nscen))  # leave one core for the main session
      if (!is.finite(workers) || workers < 1L) workers <- 1L
    } else {
      workers <- as.integer(min(max(workers, 1L), nscen))
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    # multisession works on macOS/Windows/Linux
    future::plan(future::multisession, workers = workers)
    message(sprintf("Parallel ON with %d worker(s).", future::nbrOfWorkers()))
  } else {
    message("Parallel OFF (running sequentially).")
  }

  ## --- Progress setup
  total_steps <- nscen * niter
  progressr::handlers(global = TRUE)
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
          seednum  = seednum
        )

        # run_model() should use withr::with_seed(seednum, { ... }) internally
        out[[k]] <- try(run_model(cfg), silent = TRUE)
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
        future.seed = TRUE
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
    run$parameters <- list(
      # timeline / reviews
      nyi = nyi,
      nyh = nyh,
      ny = ny,
      goalfreq = goalfreq,
      firstrev = firstrev,
      review_years = seq(from = nyi + firstrev, to = ny, by = goalfreq),
      # reproducibility / provenance
      seednum = seednum,
      rng_kind = paste(utils::capture.output(RNGkind()), collapse = " "),
      blas_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)),
      pkg_version = as.character(utils::packageVersion("SizeShiftsImplicationsV2")),
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      scenario_grid_hash = digest::digest(run$scenarios)
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
  .validate_run_params(run$parameters)
  class(run_out) <- c("ssi_run", "list")

  rds_path <- file.path(output_dir, paste0("run_", time.save, ".rds"))
  saveRDS(run_out, rds_path)
  openxlsx::write.xlsx(scen, file = file.path(output_dir, paste0("run_", time.save, "_scen.xlsx")))

  message(sprintf("Completed %d scenarios × %d iterations in %s",
                  nscen, niter, format(end.time - start.time)))
  invisible(run_out)
}
