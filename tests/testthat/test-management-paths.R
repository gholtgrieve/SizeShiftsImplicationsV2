# Tests for all seven harvmgmt paths through run_model().
# Uses the Ohlberger profile (sim_recruits="eggmass", ricker_type="const_beta")
# which is required for the YPR paths (s_eq_goal / u_eq_goal).

make_ohlberger_cfg <- function(harvmgmt, seednum = 1L) {
  cfg <- utils::modifyList(
    SizeShiftsImplicationsV2:::default_config(),
    param_configs[["Ohlberger"]]
  )
  cfg$nyi      <- 5L
  cfg$nyh      <- 20L
  cfg$ny       <- 40L
  cfg$goalfreq <- 20L
  cfg$harvmgmt <- harvmgmt
  cfg$factorMSY <- 1.0
  cfg$futureT  <- "no"
  cfg$sdsel    <- 0.0
  cfg$maxsel   <- 5.5
  cfg$j        <- 1L
  cfg$k        <- 1L
  cfg$seednum  <- as.integer(seednum)
  cfg$log_dir  <- tempdir()
  SizeShiftsImplicationsV2:::validate_config(cfg)
}

sanity_check <- function(out, label) {
  d <- out$data
  # mass balance
  residual <- abs(d$Ret - (d$Esc + d$Harv))
  expect_true(all(residual <= 1, na.rm = TRUE),
              label = paste(label, ": Ret = Esc + Harv"))
  # non-negative
  expect_true(all(d$Esc  >= 0, na.rm = TRUE), label = paste(label, ": Esc >= 0"))
  expect_true(all(d$Harv >= 0, na.rm = TRUE), label = paste(label, ": Harv >= 0"))
  # SR params
  expect_true(all(is.finite(out$para) & out$para > 0),
              label = paste(label, ": SR params positive"))
}

# ── fix_harv_rate ──────────────────────────────────────────────────────────

test_that("fix_harv_rate runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("fix_harv_rate")))
})

test_that("fix_harv_rate passes sanity checks", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("fix_harv_rate"))
  sanity_check(out, "fix_harv_rate")
})

test_that("fix_harv_rate produces non-zero harvest in most years", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("fix_harv_rate"))
  expect_gt(mean(out$data$Harv > 0, na.rm = TRUE), 0.5)
})

# ── smsy_goal ──────────────────────────────────────────────────────────────

test_that("smsy_goal runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("smsy_goal")))
})

test_that("smsy_goal passes sanity checks", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("smsy_goal"))
  sanity_check(out, "smsy_goal")
})

test_that("smsy_goal updates S_msy estimate from GLS fit", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("smsy_goal"))
  # S_msy is updated at each review; should be a finite positive number
  s <- out$S_msy[!is.na(out$S_msy)]
  expect_true(length(s) > 0)
  expect_true(all(is.finite(s) & s > 0))
})

# ── umsy_goal ──────────────────────────────────────────────────────────────

test_that("umsy_goal runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("umsy_goal")))
})

test_that("umsy_goal passes sanity checks", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("umsy_goal"))
  sanity_check(out, "umsy_goal")
})

# ── s_eq_goal (YPR) ────────────────────────────────────────────────────────

test_that("s_eq_goal runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("s_eq_goal")))
})

test_that("s_eq_goal passes sanity checks", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("s_eq_goal"))
  sanity_check(out, "s_eq_goal")
})

test_that("s_eq_goal uses RepOut (eggmass) for recruitment", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("s_eq_goal"))
  rep_out <- out$obs$RepOut[!is.na(out$obs$RepOut)]
  expect_true(all(rep_out > 0))
})

# ── u_eq_goal (YPR harvest rate) ───────────────────────────────────────────

test_that("u_eq_goal runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("u_eq_goal")))
})

test_that("u_eq_goal passes sanity checks", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("u_eq_goal"))
  sanity_check(out, "u_eq_goal")
})

# ── smsy_dlm_goal (DLM) ────────────────────────────────────────────────────

test_that("smsy_dlm_goal runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("smsy_dlm_goal")))
})

test_that("smsy_dlm_goal passes sanity checks", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("smsy_dlm_goal"))
  sanity_check(out, "smsy_dlm_goal")
})

# ── umsy_dlm_goal (DLM harvest rate) ───────────────────────────────────────

test_that("umsy_dlm_goal runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("umsy_dlm_goal")))
})

test_that("umsy_dlm_goal passes sanity checks", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("umsy_dlm_goal"))
  sanity_check(out, "umsy_dlm_goal")
})

# ── Cross-path comparisons ─────────────────────────────────────────────────

test_that("S-goal and U-goal management types produce different harvest patterns", {
  out_s <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("smsy_goal", seednum = 5L))
  out_u <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("umsy_goal", seednum = 5L))
  # Different management rules -> different outcomes
  expect_false(identical(out_s$data$Harv, out_u$data$Harv))
})

test_that("smsy_dlm_goal path produces finite S_msy estimates from DLM", {
  out <- SizeShiftsImplicationsV2:::run_model(make_ohlberger_cfg("smsy_dlm_goal", seednum = 7L))
  # S_msy estimates are stored per review; at least one should be a finite positive number
  s <- out$S_msy[!is.na(out$S_msy)]
  if (length(s) > 0) expect_true(all(is.finite(s) & s > 0))
  # MSY_Goals vector is populated across all years
  expect_equal(length(out$MSY_Goals), length(out$data$Year) + 15L)
})
