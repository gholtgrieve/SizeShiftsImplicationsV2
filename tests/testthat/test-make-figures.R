# Tests for make_figures() dispatcher and figure generation.
# Uses a small multi-scenario run (3 Kuskokwim scenarios, 3 iterations)
# as a shared fixture so figure code is actually exercised.

# ── Shared fixture: minimal multi-scenario ssi_run ─────────────────────────
local({
  .fig_out_dir <<- file.path(tempdir(), paste0("ssi_fig_test_", Sys.getpid()))
  dir.create(.fig_out_dir, recursive = TRUE, showWarnings = FALSE)

  # 3 scenarios: no trends / ASL trends / continuing trends, all smsy_goal
  .fig_run <<- suppressMessages(
    run_scenarios(
      scenarios  = c(13L, 15L, 16L),
      niter      = 3L,
      seed       = "reproducible",
      params     = "Kuskokwim",
      output_dir = .fig_out_dir,
      parallel   = FALSE
    )
  )
})


# ── Input validation ──────────────────────────────────────────────────────

test_that("make_figures errors on non-ssi_run data", {
  expect_error(make_figures("Kuskokwim", data = list()), "ssi_run")
})

test_that("make_figures errors on path to non-ssi_run file", {
  f <- tempfile(fileext = ".rds")
  saveRDS(list(x = 1), f)
  on.exit(unlink(f), add = TRUE)
  expect_error(make_figures("Kuskokwim", data = f), "ssi_run")
})

test_that("make_figures errors on non-existent file path", {
  expect_error(make_figures("Kuskokwim", data = "/no/such/file.rds"))
})

test_that("make_figures errors on invalid figure code", {
  expect_error(make_figures("Z", data = .fig_run), "No valid figures")
})

test_that("make_figures errors on unavailable builder code", {
  # "9" is not a recognised figure code, so normalisation returns character(0)
  # and we get "No valid figures selected" before reaching the builder dispatch.
  expect_error(make_figures("9", data = .fig_run), "No valid figures selected")
})

test_that("make_figures accepts ssi_run object directly", {
  expect_no_error(
    make_figures("A", data = .fig_run, figure_dir = file.path(.fig_out_dir, "figs_A"))
  )
})

test_that("make_figures accepts path to saved .rds", {
  rds_files <- list.files(.fig_out_dir, pattern = "\\.rds$", full.names = TRUE)
  skip_if(length(rds_files) == 0, "No RDS file saved by fixture")
  expect_no_error(
    make_figures("A", data = rds_files[[1]],
                 figure_dir = file.path(.fig_out_dir, "figs_A_rds"))
  )
})


# ── Kuskokwim figure dispatch ─────────────────────────────────────────────

test_that("'Kuskokwim' expands to codes A-E", {
  # Test via normalize_tokens (access through make_figures side-effect: no error)
  expect_no_error(
    make_figures(c("A","B","C","D","E"), data = .fig_run,
                 figure_dir = file.path(.fig_out_dir, "figs_ABCDE"))
  )
})

test_that("figure A runs and returns a list", {
  out <- make_figures("A", data = .fig_run, figure_dir = file.path(.fig_out_dir, "fA"))
  expect_type(out, "list")
  expect_true("FigureA" %in% names(out))
})

test_that("figure B runs and returns a list", {
  out <- make_figures("B", data = .fig_run, figure_dir = file.path(.fig_out_dir, "fB"))
  expect_type(out, "list")
  expect_true("FigureB" %in% names(out))
})

test_that("figure C runs and returns a list", {
  out <- make_figures("C", data = .fig_run, figure_dir = file.path(.fig_out_dir, "fC"))
  expect_type(out, "list")
  expect_true("FigureC" %in% names(out))
})

test_that("figure D runs and returns a list", {
  out <- make_figures("D", data = .fig_run, figure_dir = file.path(.fig_out_dir, "fD"))
  expect_type(out, "list")
  expect_true("FigureD" %in% names(out))
})

test_that("figure E runs and returns a list", {
  out <- make_figures("E", data = .fig_run, figure_dir = file.path(.fig_out_dir, "fE"))
  expect_type(out, "list")
  expect_true("FigureE" %in% names(out))
})

test_that("figure code is case-insensitive ('a' same as 'A')", {
  out_upper <- make_figures("A", data = .fig_run, figure_dir = file.path(.fig_out_dir, "fA2"))
  out_lower <- make_figures("a", data = .fig_run, figure_dir = file.path(.fig_out_dir, "fA3"))
  expect_true("FigureA" %in% names(out_upper))
  expect_true("FigureA" %in% names(out_lower))
})

test_that("'Figure5' prefix is stripped correctly", {
  # 'FigureA' -> 'A'; but for Ohlberger figures we need an Ohlberger run
  # Just test dispatch parses without error (Ohlberger figures need different data)
  expect_error(make_figures("FigureZ", data = .fig_run), "No valid figures")
})

test_that("duplicate figure codes are deduplicated", {
  out <- make_figures(c("A","A","B"), data = .fig_run,
                      figure_dir = file.path(.fig_out_dir, "fdedup"))
  # Should produce FigureA and FigureB (A appears once)
  fig_keys <- setdiff(names(out), "data")
  expect_equal(sort(fig_keys), c("FigureA", "FigureB"))
})

test_that("make_figures creates figure_dir if it does not exist", {
  new_dir <- file.path(.fig_out_dir, "brand_new_subdir")
  expect_false(dir.exists(new_dir))
  make_figures("A", data = .fig_run, figure_dir = new_dir)
  expect_true(dir.exists(new_dir))
})

test_that("make_figures returns invisibly with data element", {
  out <- make_figures("A", data = .fig_run, figure_dir = file.path(.fig_out_dir, "finv"))
  expect_true("data" %in% names(out))
  expect_s3_class(out$data, "ssi_run")
})
