# End-to-end smoke test: run_scenarios() -> ssi_run structure + biological sanity.
#
# Uses 3 iterations of a single scenario to keep runtime short (~10–20 s).
# Exercises the full pipeline: scenario selection -> config -> model -> output.

# Shared fixture: run once, reuse across all tests in this file.
# Each test session gets its own timestamped subfolder under outputs/test_runs/
# (gitignored), so every run is preserved and can be loaded for inspection.
local({
  .test_out_dir <<- file.path(
    here::here(), "outputs", "test_runs",
    format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  )
  dir.create(.test_out_dir, recursive = TRUE, showWarnings = FALSE)

  .run <<- run_scenarios(
    scenarios  = 13L,
    niter      = 3L,
    seed       = "reproducible",
    params     = "Ohlberger",
    output_dir = .test_out_dir,
    parallel   = FALSE
  )
})


# ── ssi_run object structure ────────────────────────────────────────────────

test_that("run_scenarios() returns an ssi_run object", {
  expect_s3_class(.run, "ssi_run")
})

test_that("ssi_run has required top-level fields", {
  expect_named(.run, c("meta", "parameters", "scenarios", "results"),
               ignore.order = TRUE)
})

test_that("meta records correct niter and nscen", {
  expect_equal(.run$meta$niter, 3L)
  expect_equal(.run$meta$nscen, 1L)
})

test_that("parameters records timeline fields", {
  p <- .run$parameters
  expect_true(all(c("nyi", "nyh", "ny", "goalfreq", "review_years") %in% names(p)))
  expect_equal(p$nyi, 10L)
  expect_equal(p$nyh, 50L)
  expect_equal(p$ny,  110L)
})

test_that("scenarios tibble has one row matching scen_num 13", {
  expect_equal(nrow(.run$scenarios), 1L)
  expect_equal(.run$scenarios$scen_num, 13)
})


# ── results completeness ────────────────────────────────────────────────────

test_that("results list has all expected sublists", {
  expected <- c("para", "sr_sim", "fec", "egg", "S_msy", "data",
                "obs", "ret_by_age", "meanSaA", "propfemale",
                "selectivity", "MSY_Goals", "impl_errors")
  expect_true(all(expected %in% names(.run$results)))
})

test_that("each iteration produced a non-NULL data result", {
  iters <- .run$results$data[[1]]
  for (k in seq_len(.run$meta$niter)) {
    expect_false(is.null(iters[[k]]),
                 label = sprintf("data[[1]][[%d]] should not be NULL", k))
  }
})

test_that("no iteration result is a try-error", {
  for (field in c("data", "obs", "para")) {
    for (k in seq_len(.run$meta$niter)) {
      expect_false(inherits(.run$results[[field]][[1]][[k]], "try-error"),
                   label = sprintf("%s iter %d should not be try-error", field, k))
    }
  }
})


# ── biological sanity ────────────────────────────────────────────────────────

test_that("annual data has expected population columns", {
  d <- .run$results$data[[1]][[1]]
  expect_true(all(c("Ret", "Rec", "Esc", "Harv") %in% names(d)))
})

test_that("escapement and harvest are non-negative in all iterations", {
  for (k in seq_len(.run$meta$niter)) {
    d <- .run$results$data[[1]][[k]]
    expect_true(all(d$Esc  >= 0, na.rm = TRUE),
                label = sprintf("iter %d: Esc must be non-negative", k))
    expect_true(all(d$Harv >= 0, na.rm = TRUE),
                label = sprintf("iter %d: Harv must be non-negative", k))
  }
})

test_that("Ret = Esc + Harv (mass balance) in all iterations", {
  for (k in seq_len(.run$meta$niter)) {
    d <- .run$results$data[[1]][[k]]
    residual <- abs(d$Ret - (d$Esc + d$Harv))
    expect_true(all(residual <= 1, na.rm = TRUE),
                label = sprintf("iter %d: Ret should equal Esc + Harv (±1)", k))
  }
})

test_that("SR parameters are positive and finite in all iterations", {
  for (k in seq_len(.run$meta$niter)) {
    p <- .run$results$para[[1]][[k]]
    expect_true(all(is.finite(p) & p > 0),
                label = sprintf("iter %d: all SR params should be finite and positive", k))
  }
})

test_that("meanSaA has correct column dimension (nage = 9) in all iterations", {
  for (k in seq_len(.run$meta$niter)) {
    m <- .run$results$meanSaA[[1]][[k]]
    expect_equal(ncol(m), 9L,
                 label = sprintf("iter %d: meanSaA should have 9 age columns", k))
  }
})

test_that("meanSaA columns are named age1..age9 (dimnames preserved after + operator)", {
  m <- .run$results$meanSaA[[1]][[1]]
  expect_equal(colnames(m), paste0("age", 1:9))
})

test_that("long-run mean escapement is positive across iterations", {
  mean_esc <- sapply(seq_len(.run$meta$niter), function(k) {
    mean(.run$results$data[[1]][[k]]$Esc, na.rm = TRUE)
  })
  expect_true(all(mean_esc > 0))
})

test_that("results are reproducible: same seed gives identical recruits", {
  # Re-run with the same seed and compare iter 1 to the fixture
  run2 <- run_scenarios(
    scenarios  = 13L,
    niter      = 1L,
    seed       = "reproducible",
    params     = "Ohlberger",
    output_dir = .test_out_dir,
    parallel   = FALSE
  )
  d1 <- .run$results$data[[1]][[1]]$Rec
  d2 <- run2$results$data[[1]][[1]]$Rec
  expect_equal(d1, d2)
})
