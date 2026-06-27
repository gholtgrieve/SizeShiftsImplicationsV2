# Tests for log helper utilities in log_helpers.R

pkg <- "SizeShiftsImplicationsV2"


# ── fun_reason ────────────────────────────────────────────────────────────

test_that("fun_reason returns NA when all inputs are NULL", {
  expect_true(is.na(SizeShiftsImplicationsV2:::fun_reason(NULL, NULL)))
})

test_that("fun_reason returns NA when all inputs are empty strings", {
  expect_true(is.na(SizeShiftsImplicationsV2:::fun_reason("", "")))
})

test_that("fun_reason concatenates non-empty strings with semicolons", {
  out <- SizeShiftsImplicationsV2:::fun_reason("a", NULL, "b")
  expect_equal(out, "a;b")
})

test_that("fun_reason deduplicates repeated strings", {
  out <- SizeShiftsImplicationsV2:::fun_reason("foo", "foo", "bar")
  expect_equal(out, "foo;bar")
})

test_that("fun_reason filters empty strings but keeps non-empty", {
  out <- SizeShiftsImplicationsV2:::fun_reason("", "x", NULL, "y")
  expect_equal(out, "x;y")
})

test_that("fun_reason returns a single string with no semicolon for one input", {
  out <- SizeShiftsImplicationsV2:::fun_reason("only_one")
  expect_equal(out, "only_one")
  expect_false(grepl(";", out))
})


# ── fun_rng_range ────────────────────────────────────────────────────────

test_that("fun_rng_range returns min and max for numeric vector", {
  out <- SizeShiftsImplicationsV2:::fun_rng_range(c(3, 1, 4, 1, 5, 9))
  expect_equal(out$min, 1)
  expect_equal(out$max, 9)
})

test_that("fun_rng_range returns NA for all-NA input", {
  out <- SizeShiftsImplicationsV2:::fun_rng_range(c(NA_real_, NA_real_))
  expect_true(is.na(out$min))
  expect_true(is.na(out$max))
})

test_that("fun_rng_range ignores non-finite values", {
  out <- SizeShiftsImplicationsV2:::fun_rng_range(c(1, Inf, 5, -Inf, NA))
  expect_equal(out$min, 1)
  expect_equal(out$max, 5)
})


# ── fun_ypr_converged ────────────────────────────────────────────────────

test_that("fun_ypr_converged returns TRUE for valid interior solution", {
  ok <- SizeShiftsImplicationsV2:::fun_ypr_converged(
    alpha_rep_out = 2.0, beta_rep_out = 1e-5,
    zPR0 = 1e6, R0_val = 1e5,
    F_eq = 0.0, fit_value = -50, H_eq = 5000, S_eq = 10000, U_eq = 0.33,
    fit_convergence = 0L, lower = -10, upper = 2
  )
  expect_true(ok)
})

test_that("fun_ypr_converged returns FALSE for non-finite alpha", {
  ok <- SizeShiftsImplicationsV2:::fun_ypr_converged(
    alpha_rep_out = NA, beta_rep_out = 1e-5,
    zPR0 = 1e6, R0_val = 1e5,
    F_eq = 0.0, fit_value = -50, H_eq = 5000, S_eq = 10000, U_eq = 0.33,
    fit_convergence = 0L
  )
  expect_false(ok)
})

test_that("fun_ypr_converged returns FALSE when F_eq is at lower boundary", {
  ok <- SizeShiftsImplicationsV2:::fun_ypr_converged(
    alpha_rep_out = 2.0, beta_rep_out = 1e-5,
    zPR0 = 1e6, R0_val = 1e5,
    F_eq = -10.0, fit_value = -50, H_eq = 5000, S_eq = 10000, U_eq = 0.33,
    fit_convergence = 0L, lower = -10, upper = 2
  )
  expect_false(ok)
})

test_that("fun_ypr_converged returns FALSE when F_eq is at upper boundary", {
  ok <- SizeShiftsImplicationsV2:::fun_ypr_converged(
    alpha_rep_out = 2.0, beta_rep_out = 1e-5,
    zPR0 = 1e6, R0_val = 1e5,
    F_eq = 2.0, fit_value = -50, H_eq = 5000, S_eq = 10000, U_eq = 0.33,
    fit_convergence = 0L, lower = -10, upper = 2
  )
  expect_false(ok)
})

test_that("fun_ypr_converged returns FALSE when fit_convergence != 0", {
  ok <- SizeShiftsImplicationsV2:::fun_ypr_converged(
    alpha_rep_out = 2.0, beta_rep_out = 1e-5,
    zPR0 = 1e6, R0_val = 1e5,
    F_eq = 0.0, fit_value = -50, H_eq = 5000, S_eq = 10000, U_eq = 0.33,
    fit_convergence = 1L
  )
  expect_false(ok)
})

test_that("fun_ypr_converged returns FALSE when H_eq + S_eq <= 0", {
  ok <- SizeShiftsImplicationsV2:::fun_ypr_converged(
    alpha_rep_out = 2.0, beta_rep_out = 1e-5,
    zPR0 = 1e6, R0_val = 1e5,
    F_eq = 0.0, fit_value = -50, H_eq = 0, S_eq = 0, U_eq = 0,
    fit_convergence = 0L
  )
  expect_false(ok)
})

test_that("fun_ypr_converged allows boundary when allow_boundary=TRUE", {
  ok <- SizeShiftsImplicationsV2:::fun_ypr_converged(
    alpha_rep_out = 2.0, beta_rep_out = 1e-5,
    zPR0 = 1e6, R0_val = 1e5,
    F_eq = -10.0, fit_value = -50, H_eq = 5000, S_eq = 10000, U_eq = 0.33,
    fit_convergence = 0L, lower = -10, upper = 2,
    allow_boundary = TRUE
  )
  expect_true(ok)
})


# ── fun_log_dir ──────────────────────────────────────────────────────────

test_that("fun_log_dir returns the config log_dir when set", {
  cfg <- list(log_dir = tempdir())
  out <- SizeShiftsImplicationsV2:::fun_log_dir(cfg)
  expect_equal(out, tempdir())
})

test_that("fun_log_dir falls back to option when config is NULL", {
  old <- getOption("ssi.log_dir")
  on.exit(options(ssi.log_dir = old), add = TRUE)
  options(ssi.log_dir = tempdir())
  out <- SizeShiftsImplicationsV2:::fun_log_dir(NULL)
  expect_equal(out, tempdir())
})

test_that("fun_log_dir creates the directory if it does not exist", {
  new_dir <- file.path(tempdir(), paste0("ssi_test_", Sys.getpid()))
  on.exit(unlink(new_dir, recursive = TRUE), add = TRUE)
  SizeShiftsImplicationsV2:::fun_log_dir(list(log_dir = new_dir))
  expect_true(dir.exists(new_dir))
})
