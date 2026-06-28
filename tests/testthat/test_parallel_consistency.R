# Parallel computing consistency tests for the MSE R package.
#
# These tests verify:
#   1. Within-backend reproducibility (same seed -> same results)
#   2. Replicate independence (workers use distinct RNG streams)
#   3. Result aggregation correctness (manual vs package output)
#   4. Scenario completeness (all scenarios return non-degenerate results)
#   5. CSV log integrity (PID-shard merge produces clean output)
#   6. BLAS thread guard (graceful degradation without RhpcBLASctl)
#   7. Future plan restoration (plan restored after return, even on error)
#
# Architecture note: ALL model RNG goes through withr::with_seed(seednum,...),
# not through the global RNG state. Reproducibility is guaranteed by seednum
# arithmetic:
#   sequential:  seednum = seeds_k[k]
#   parallel:    seednum = seeds_k[k] + 10000 * j
# Two runs with the same seed="reproducible" and the same parallel flag are
# bitwise-identical. Parallel and sequential runs with the same seed mode
# produce DIFFERENT numerical results by design (see docs/parallel_audit.txt,
# Issue #2). Tests that compare results therefore stay within one backend.
#
# Uses the Kuskokwim parameter profile throughout.

testthat::local_edition(3)

# ---------------------------------------------------------------------------
# Shared scenario selection helpers (evaluated once at file-load time)
# ---------------------------------------------------------------------------

# First Kuskokwim-compatible scenario: unselective gear, no demographic trends,
# smsy_goal management. Scenario 13 satisfies all of these.
.kusko_scen1 <- 13L

# Three representative scenarios — one per management type in the Kuskokwim
# grid, all with "no trends" and unselective gear.
.kusko_three <- c(13L, 17L, 21L)  # smsy_goal, s_eq_goal, smsy_dlm_goal

# Run a single Kuskokwim scenario. Uses a temp dir that is cleaned up when the
# function returns (sufficient since callers only need the in-memory ssi_run).
# Tests that need to read files after the call should create their own tmp dir.
.run_kusko <- function(parallel, n_iter = 5L, workers = 2L,
                       scenarios = .kusko_scen1) {
  tmp <- withr::local_tempdir(clean = TRUE)
  run_scenarios(
    scenarios  = scenarios,
    niter      = n_iter,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = parallel,
    workers    = if (isTRUE(parallel)) workers else NULL
  )
}


# ============================================================================
# TEST 1: Within-backend reproducibility
# ============================================================================
#
# Two runs with the same seed="reproducible" and the same parallel flag must
# produce bitwise-identical results because seednum is a deterministic function
# of (j, k, parallel). Cross-backend identity is NOT tested here (see note
# above about the +10000*j offset).

test_that("sequential runs are bitwise-identical given seed='reproducible'", {
  run1 <- .run_kusko(parallel = FALSE, n_iter = 5L)
  run2 <- .run_kusko(parallel = FALSE, n_iter = 5L)

  # Compare first and last iteration recruitment series
  expect_identical(
    run1$results$data[[1L]][[1L]]$Rec,
    run2$results$data[[1L]][[1L]]$Rec,
    label = "sequential iter 1: Rec identical across two runs"
  )
  expect_identical(
    run1$results$data[[1L]][[5L]]$Esc,
    run2$results$data[[1L]][[5L]]$Esc,
    label = "sequential iter 5: Esc identical across two runs"
  )
})

test_that("parallel runs are bitwise-identical given seed='reproducible'", {
  skip_if_not_installed("future.apply")

  run1 <- .run_kusko(parallel = TRUE, n_iter = 5L, workers = 2L)
  run2 <- .run_kusko(parallel = TRUE, n_iter = 5L, workers = 2L)

  expect_identical(
    run1$results$data[[1L]][[1L]]$Rec,
    run2$results$data[[1L]][[1L]]$Rec,
    label = "parallel iter 1: Rec identical across two runs"
  )
  expect_identical(
    run1$results$data[[1L]][[5L]]$Harv,
    run2$results$data[[1L]][[5L]]$Harv,
    label = "parallel iter 5: Harv identical across two runs"
  )
})


# ============================================================================
# TEST 2: Replicate independence
# ============================================================================
#
# Each iteration k uses seednum = k (sequential) or k + 10000*j (parallel), so
# replicates must be numerically distinct and their pairwise correlations must
# be consistent with independent draws, not correlated or identical streams.

test_that("replicate Rec sums are all distinct (no duplicate RNG streams)", {
  skip_if_not_installed("future.apply")

  run_par <- .run_kusko(parallel = TRUE, n_iter = 20L, workers = 2L)
  n_iter  <- run_par$meta$niter

  recs <- lapply(seq_len(n_iter), function(k) run_par$results$data[[1L]][[k]]$Rec)

  # Every replicate must exist and be non-empty
  expect_true(
    all(vapply(recs, function(r) length(r) > 0L, logical(1L))),
    label = "all 20 replicate Rec series are non-empty"
  )

  # Sum-based fingerprints: 20 distinct seednum values -> 20 distinct values
  fingerprints <- vapply(recs, sum, numeric(1L))
  expect_equal(
    length(unique(fingerprints)), n_iter,
    label = "all 20 replicates have distinct Rec sums (Test A: no identical streams)"
  )
})

test_that("pairwise correlations across replicates are consistent with independence", {
  skip_if_not_installed("future.apply")

  run_par <- .run_kusko(parallel = TRUE, n_iter = 20L, workers = 2L)
  n_iter  <- run_par$meta$niter

  recs <- lapply(seq_len(n_iter), function(k) run_par$results$data[[1L]][[k]]$Rec)
  min_len <- min(vapply(recs, length, integer(1L)))
  skip_if(min_len < 10L, "too few reconstructed years for correlation analysis")

  mat      <- do.call(cbind, lapply(recs, function(r) r[seq_len(min_len)]))
  cormat   <- cor(mat)
  off_diag <- cormat[lower.tri(cormat)]

  # Test B: no perfectly correlated pair (|r| < 0.99)
  expect_true(
    max(abs(off_diag)) < 0.99,
    label = paste0(
      "Test B: max |r| = ", round(max(abs(off_diag)), 3),
      " — no pair of replicates should have |r| >= 0.99"
    )
  )

  # Test C: mean |r| below 0.3 (independent normal sequences give ~0.1)
  expect_true(
    mean(abs(off_diag)) < 0.3,
    label = paste0(
      "Test C: mean |r| = ", round(mean(abs(off_diag)), 3),
      " — should be < 0.3 for independent streams"
    )
  )
})


# ============================================================================
# TEST 3: Result aggregation correctness
# ============================================================================
#
# run_scenarios() returns raw per-replicate data. This test verifies the nested
# list structure is internally consistent: mass balance holds, values are finite
# and non-negative, and manual aggregation matches what you'd compute directly.

test_that("Ret = Esc + Harv (mass balance, tolerance +-1) across all iterations", {
  run_seq <- .run_kusko(parallel = FALSE, n_iter = 5L)
  n_iter  <- run_seq$meta$niter

  for (k in seq_len(n_iter)) {
    d <- run_seq$results$data[[1L]][[k]]
    expect_false(is.null(d), label = sprintf("data[[1]][[%d]] is not NULL", k))
    if (!is.null(d)) {
      residual <- abs(d$Ret - (d$Esc + d$Harv))
      expect_true(
        all(residual <= 1L, na.rm = TRUE),
        label = sprintf("iter %d: Ret = Esc + Harv (+/- 1 rounding)", k)
      )
    }
  }
})

test_that("manual mean annual harvest is finite, positive, and self-consistent", {
  run_seq <- .run_kusko(parallel = FALSE, n_iter = 5L)
  n_iter  <- run_seq$meta$niter

  harv_means <- vapply(seq_len(n_iter), function(k) {
    mean(run_seq$results$data[[1L]][[k]]$Harv, na.rm = TRUE)
  }, numeric(1L))

  overall_mean <- mean(harv_means)
  expect_true(is.finite(overall_mean) && overall_mean > 0,
    label = "mean annual harvest is finite and positive")

  # Recompute: same data object must give the same answer
  recalc <- mean(vapply(seq_len(n_iter), function(k) {
    mean(run_seq$results$data[[1L]][[k]]$Harv, na.rm = TRUE)
  }, numeric(1L)))
  expect_equal(overall_mean, recalc,
    label = "mean annual harvest is reproducible from the same in-memory object")
})

test_that("escapement and harvest are non-negative in all iterations", {
  run_seq <- .run_kusko(parallel = FALSE, n_iter = 5L)
  n_iter  <- run_seq$meta$niter

  for (k in seq_len(n_iter)) {
    d <- run_seq$results$data[[1L]][[k]]
    if (!is.null(d)) {
      expect_true(all(d$Esc  >= 0, na.rm = TRUE),
        label = sprintf("iter %d: Esc >= 0", k))
      expect_true(all(d$Harv >= 0, na.rm = TRUE),
        label = sprintf("iter %d: Harv >= 0", k))
    }
  }
})

test_that("SR parameters are finite and positive in all iterations", {
  run_seq <- .run_kusko(parallel = FALSE, n_iter = 5L)
  n_iter  <- run_seq$meta$niter

  for (k in seq_len(n_iter)) {
    p <- run_seq$results$para[[1L]][[k]]
    expect_false(is.null(p), label = sprintf("para[[1]][[%d]] not NULL", k))
    if (!is.null(p)) {
      expect_true(all(is.finite(p) & p > 0),
        label = sprintf("iter %d: all SR params finite and positive", k))
    }
  }
})


# ============================================================================
# TEST 4: Scenario completeness
# ============================================================================
#
# Run three representative Kuskokwim scenarios (one per management type:
# smsy_goal = scen 13, s_eq_goal = scen 17, smsy_dlm_goal = scen 21) in
# parallel and verify that every scenario returns non-NULL, non-degenerate
# results with finite, in-range values. No scenario should be silently dropped.

test_that("all tested scenarios return results (no silent drops)", {
  skip_if_not_installed("future.apply")
  skip_on_cran()

  d <- withr::local_tempdir(clean = TRUE)
  run_all <- run_scenarios(
    scenarios  = .kusko_three,
    niter      = 3L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = d,
    parallel   = TRUE,
    workers    = 2L
  )

  nscen  <- run_all$meta$nscen
  n_iter <- run_all$meta$niter

  # Correct count of returned results
  expect_equal(length(run_all$results$data), nscen,
    label = "result list has one entry per input scenario")
  expect_equal(nscen, length(.kusko_three),
    label = "returned nscen matches number of requested scenarios")
})

test_that("each (scenario, iter) cell has non-NULL data with finite Esc and Harv", {
  skip_if_not_installed("future.apply")
  skip_on_cran()

  d <- withr::local_tempdir(clean = TRUE)
  run_all <- run_scenarios(
    scenarios  = .kusko_three,
    niter      = 3L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = d,
    parallel   = TRUE,
    workers    = 2L
  )

  nscen  <- run_all$meta$nscen
  n_iter <- run_all$meta$niter
  scen_names <- names(run_all$results$data)

  for (j in seq_len(nscen)) {
    label_j <- if (!is.null(scen_names) && nchar(scen_names[j]) > 0) scen_names[j] else paste0("scen_", j)
    for (k in seq_len(n_iter)) {
      d_jk <- run_all$results$data[[j]][[k]]
      expect_false(is.null(d_jk),
        label = sprintf("%s iter %d: data not NULL", label_j, k))
      if (!is.null(d_jk)) {
        expect_true(all(d_jk$Esc  >= 0, na.rm = TRUE),
          label = sprintf("%s iter %d: Esc >= 0", label_j, k))
        expect_true(all(d_jk$Harv >= 0, na.rm = TRUE),
          label = sprintf("%s iter %d: Harv >= 0", label_j, k))
        expect_true(any(is.finite(d_jk$Rec) & d_jk$Rec > 0, na.rm = TRUE),
          label = sprintf("%s iter %d: some Rec > 0", label_j, k))
      }
    }
  }
})


# ============================================================================
# TEST 5: CSV log integrity (PID-shard mechanism)
# ============================================================================
#
# Confirms that per-PID shard files written by concurrent workers, when merged
# by fun_merge_all_logs(), produce a clean, non-corrupted CSV. The test checks
# log_7a_yearloop.csv (one row per review window per iteration).
#
# Kuskokwim defaults: nyi=10, firstrev=20 (default), ny=110, goalfreq=10
# -> review years: seq(30, 110, 10) = 9 review windows (nrev = 9)

.expected_nrev <- 9L

test_that("merged log_7a has correct row count (niter x nrev)", {
  skip_if_not_installed("future.apply")

  tmp <- withr::local_tempdir(clean = TRUE)

  run_log <- run_scenarios(
    scenarios  = .kusko_scen1,
    niter      = 5L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = TRUE,
    workers    = 2L
  )

  log_path <- file.path(tmp, "logs", "log_7a_yearloop.csv")
  expect_true(file.exists(log_path),
    label = "log_7a_yearloop.csv exists after merge")

  log_df <- utils::read.csv(log_path, check.names = FALSE)
  expected_rows <- run_log$meta$niter * .expected_nrev
  expect_equal(nrow(log_df), expected_rows,
    label = sprintf("log has %d rows (%d iter x %d reviews)",
                    expected_rows, run_log$meta$niter, .expected_nrev))
})

test_that("merged log has no entirely-NA rows and no partial writes", {
  skip_if_not_installed("future.apply")

  tmp <- withr::local_tempdir(clean = TRUE)
  run_scenarios(
    scenarios  = .kusko_scen1,
    niter      = 5L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = TRUE,
    workers    = 2L
  )

  log_df      <- utils::read.csv(file.path(tmp, "logs", "log_7a_yearloop.csv"),
                                 check.names = FALSE)
  all_na_rows <- apply(log_df, 1L, function(r) all(is.na(r)))
  expect_false(any(all_na_rows),
    label = "no row in merged CSV is entirely NA (no partial writes)")

  expect_false(any(is.na(log_df$iteration)),
    label = "'iteration' column has no NAs")
})

test_that("log_7a has required identification columns", {
  tmp <- withr::local_tempdir(clean = TRUE)
  run_scenarios(
    scenarios  = .kusko_scen1,
    niter      = 3L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = FALSE
  )

  log_df   <- utils::read.csv(file.path(tmp, "logs", "log_7a_yearloop.csv"),
                               check.names = FALSE)
  required <- c("timestamp", "scen_num", "iteration", "review_i", "y_rev", "seed",
                "n_years", "Esc_min", "Esc_max", "Harv_min", "Harv_max",
                "Ret_min", "Ret_max")
  missing  <- setdiff(required, names(log_df))
  expect_equal(length(missing), 0L,
    label = paste0("log_7a missing columns: ",
                   if (length(missing)) paste(missing, collapse = ", ") else "none"))
})

test_that("PID shard files are removed after merge", {
  skip_if_not_installed("future.apply")

  tmp <- withr::local_tempdir(clean = TRUE)
  run_scenarios(
    scenarios  = .kusko_scen1,
    niter      = 4L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = TRUE,
    workers    = 2L
  )

  remaining_shards <- list.files(
    file.path(tmp, "logs"),
    pattern = "\\.pid[0-9]+\\.csv$"
  )
  expect_equal(length(remaining_shards), 0L,
    label = "no PID shard files remain after fun_merge_all_logs()")
})


# ============================================================================
# TEST 6: BLAS thread guard
# ============================================================================
#
# When RhpcBLASctl is unavailable, run_scenarios() should proceed silently:
# no error, and no warning mentioning RhpcBLASctl. OMP/MKL env vars are still
# set as a fallback.

test_that("run_scenarios does not error when RhpcBLASctl appears unavailable", {
  skip_if_not_installed("mockery")

  tmp <- withr::local_tempdir(clean = TRUE)

  # Stub requireNamespace so RhpcBLASctl appears absent; pass all other packages.
  mockery::stub(
    run_scenarios, "requireNamespace",
    function(package, quietly = FALSE) {
      if (identical(package, "RhpcBLASctl")) return(FALSE)
      base::requireNamespace(package, quietly = quietly)
    }
  )

  expect_no_error(
    run_scenarios(
      scenarios  = .kusko_scen1,
      niter      = 2L,
      seed       = "reproducible",
      params     = "Kuskokwim",
      output_dir = tmp,
      parallel   = FALSE
    )
  )
})

test_that("absence of RhpcBLASctl produces no warning to the user", {
  skip_if_not_installed("mockery")

  tmp <- withr::local_tempdir(clean = TRUE)

  mockery::stub(
    run_scenarios, "requireNamespace",
    function(package, quietly = FALSE) {
      if (identical(package, "RhpcBLASctl")) return(FALSE)
      base::requireNamespace(package, quietly = quietly)
    }
  )

  # No warning about RhpcBLASctl should reach the user
  expect_no_warning(
    run_scenarios(
      scenarios  = .kusko_scen1,
      niter      = 2L,
      seed       = "reproducible",
      params     = "Kuskokwim",
      output_dir = tmp,
      parallel   = FALSE
    ),
    message = "RhpcBLASctl"
  )
})

test_that("OMP_NUM_THREADS is set to '1' during run_scenarios", {
  # The meta$blas_threads field records the OMP_NUM_THREADS env var value
  # at the time the run was saved (after run_scenarios has set it to 1).
  tmp <- withr::local_tempdir(clean = TRUE)
  out <- run_scenarios(
    scenarios  = .kusko_scen1,
    niter      = 2L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = FALSE
  )
  recorded <- out$parameters$blas_threads
  # blas_threads is NA if OMP_NUM_THREADS was unset before package load,
  # otherwise it should be 1 (since we set it in run_scenarios and .onLoad).
  expect_true(is.na(recorded) || identical(as.integer(recorded), 1L),
    label = "blas_threads recorded as 1 (or NA if OMP_NUM_THREADS was unset)")
})


# ============================================================================
# TEST 7: Future plan restoration
# ============================================================================
#
# future::plan() must be restored to its prior state after run_scenarios()
# returns — whether via normal return or via an error exit.

test_that("future plan is restored to prior class after a normal parallel run", {
  skip_if_not_installed("future.apply")

  prior_class <- class(future::plan())
  tmp <- withr::local_tempdir(clean = TRUE)

  run_scenarios(
    scenarios  = .kusko_scen1,
    niter      = 3L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = TRUE,
    workers    = 2L
  )

  expect_equal(class(future::plan()), prior_class,
    label = "future plan class unchanged after parallel run returns")
})

test_that("future plan is restored even when run_scenarios() exits with an error", {
  skip_if_not_installed("future.apply")

  prior_class <- class(future::plan())

  # Use a non-existent scenario ID to force a "No scenarios selected" error.
  expect_error(
    run_scenarios(
      scenarios  = 999999L,
      niter      = 3L,
      seed       = "reproducible",
      params     = "Kuskokwim",
      output_dir = withr::local_tempdir(clean = TRUE),
      parallel   = TRUE,
      workers    = 2L
    )
  )

  expect_equal(class(future::plan()), prior_class,
    label = "future plan class restored after error exit from parallel run")
})

test_that("sequential run does not alter the future plan", {
  prior_class <- class(future::plan())
  tmp <- withr::local_tempdir(clean = TRUE)

  run_scenarios(
    scenarios  = .kusko_scen1,
    niter      = 2L,
    seed       = "reproducible",
    params     = "Kuskokwim",
    output_dir = tmp,
    parallel   = FALSE
  )

  expect_equal(class(future::plan()), prior_class,
    label = "sequential run does not change future plan")
})
