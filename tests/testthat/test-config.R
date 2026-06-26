# Tests for config building, validation, and scenario override pipeline
# (config_helpers.R)

pkg <- "SizeShiftsImplicationsV2"


# ============================================================
# resolve_profile
# ============================================================

test_that("resolve_profile accepts 'Ohlberger' by name", {
  p <- SizeShiftsImplicationsV2:::resolve_profile("Ohlberger")
  expect_type(p, "list")
  expect_true("alpha_mean" %in% names(p))
})

test_that("resolve_profile accepts 'Kuskokwim' by name", {
  p <- SizeShiftsImplicationsV2:::resolve_profile("Kuskokwim")
  expect_type(p, "list")
  expect_true("alpha_mean" %in% names(p))
})

test_that("resolve_profile accepts a list directly", {
  custom <- list(alpha_mean = 3, beta_mean = 1e-5)
  p <- SizeShiftsImplicationsV2:::resolve_profile(custom)
  expect_identical(p, custom)
})

test_that("resolve_profile errors on unknown name", {
  expect_error(SizeShiftsImplicationsV2:::resolve_profile("Unknown"), "must be one of")
})

test_that("resolve_profile errors on non-string, non-list input", {
  expect_error(SizeShiftsImplicationsV2:::resolve_profile(42), "must be a profile name")
})


# ============================================================
# param_configs (exported dataset)
# ============================================================

test_that("param_configs is a named list with Ohlberger and Kuskokwim", {
  expect_type(param_configs, "list")
  expect_true("Ohlberger" %in% names(param_configs))
  expect_true("Kuskokwim" %in% names(param_configs))
})

test_that("Ohlberger profile has required SR parameters", {
  p <- param_configs[["Ohlberger"]]
  expect_true(all(c("alpha_mean", "beta_mean", "procerr", "rho") %in% names(p)))
})

test_that("Kuskokwim profile has higher rho than Ohlberger", {
  expect_gt(param_configs[["Kuskokwim"]]$rho, param_configs[["Ohlberger"]]$rho)
})

test_that("Kuskokwim profile has higher alpha_mean than Ohlberger", {
  expect_gt(param_configs[["Kuskokwim"]]$alpha_mean, param_configs[["Ohlberger"]]$alpha_mean)
})

test_that("profiles have matching allometry and alt_sr_param", {
  # Both profiles use the same allometry and alt SR params (Staton et al. 2021)
  expect_equal(param_configs[["Ohlberger"]]$allometry,
               param_configs[["Kuskokwim"]]$allometry)
  expect_equal(param_configs[["Ohlberger"]]$alt_sr_param,
               param_configs[["Kuskokwim"]]$alt_sr_param)
})


# ============================================================
# apply_scenario_overrides
# ============================================================

# Helper to make a minimal scenario row
make_scen_row <- function(...) {
  as.data.frame(list(...), stringsAsFactors = FALSE)
}

test_that("apply_scenario_overrides sets harvmgmt from 'mgmt' column", {
  cfg <- SizeShiftsImplicationsV2:::resolve_profile("Ohlberger")
  scen <- make_scen_row(mgmt = "smsy_goal", factorMSY = 1.0,
                        futureT = "no", sdsel = 0, maxsel = 5.5,
                        ageT = 1, sizeT = 1, sexT = 1)
  out <- SizeShiftsImplicationsV2:::apply_scenario_overrides(cfg, scen)
  expect_equal(out$harvmgmt, "smsy_goal")
})

test_that("apply_scenario_overrides multiplies agetrend by ageT flag", {
  cfg <- SizeShiftsImplicationsV2:::resolve_profile("Ohlberger")
  base_agetrend <- cfg$agetrend

  scen0 <- make_scen_row(mgmt = "fix_harv_rate", factorMSY = 1.0,
                         futureT = "no", sdsel = 0, maxsel = 5.5,
                         ageT = 0, sizeT = 1, sexT = 1)
  out0 <- SizeShiftsImplicationsV2:::apply_scenario_overrides(cfg, scen0)
  expect_equal(out0$agetrend, 0)  # ageT=0 zeroes the trend

  scen1 <- make_scen_row(mgmt = "fix_harv_rate", factorMSY = 1.0,
                         futureT = "no", sdsel = 0, maxsel = 5.5,
                         ageT = 1, sizeT = 1, sexT = 1)
  out1 <- SizeShiftsImplicationsV2:::apply_scenario_overrides(cfg, scen1)
  expect_equal(out1$agetrend, base_agetrend)  # ageT=1 preserves the trend
})

test_that("apply_scenario_overrides multiplies sizetrends by sizeT flag", {
  cfg <- SizeShiftsImplicationsV2:::resolve_profile("Ohlberger")
  scen <- make_scen_row(mgmt = "fix_harv_rate", factorMSY = 1.0,
                        futureT = "no", sdsel = 0, maxsel = 5.5,
                        ageT = 1, sizeT = 0, sexT = 1)
  out <- SizeShiftsImplicationsV2:::apply_scenario_overrides(cfg, scen)
  expect_true(all(out$sizetrends == 0))
})

test_that("apply_scenario_overrides errors if scen_row is not a single row", {
  cfg <- SizeShiftsImplicationsV2:::resolve_profile("Ohlberger")
  bad <- data.frame(mgmt = c("fix_harv_rate", "smsy_goal"))
  expect_error(SizeShiftsImplicationsV2:::apply_scenario_overrides(cfg, bad),
               "one-row")
})


# ============================================================
# validate_config
# ============================================================

test_that("validate_config passes for a complete config", {
  scen <- make_scen_row(mgmt = "fix_harv_rate", factorMSY = 1.0,
                        futureT = "no", sdsel = 0, maxsel = 5.5,
                        ageT = 0, sizeT = 0, sexT = 0)
  cfg <- SizeShiftsImplicationsV2:::build_config(
    params   = "Ohlberger",
    scen_row = scen,
    j = 1L, k = 1L, seednum = 1L,
    log_dir  = tempdir()
  )
  expect_type(cfg, "list")
})

test_that("validate_config errors when required fields are missing", {
  bad <- list(alpha_mean = 5)  # far too few fields
  expect_error(SizeShiftsImplicationsV2:::validate_config(bad),
               "missing required field")
})


# ============================================================
# build_config
# ============================================================

test_that("build_config returns a list with j, k, seednum set", {
  scen <- make_scen_row(mgmt = "fix_harv_rate", factorMSY = 1.0,
                        futureT = "no", sdsel = 0, maxsel = 5.5,
                        ageT = 0, sizeT = 0, sexT = 0)
  cfg <- SizeShiftsImplicationsV2:::build_config(
    params = "Ohlberger", scen_row = scen,
    j = 3L, k = 7L, seednum = 42L, log_dir = tempdir()
  )
  expect_equal(cfg$j, 3L)
  expect_equal(cfg$k, 7L)
  expect_equal(cfg$seednum, 42L)
})

test_that("build_config applies scenario overrides before returning", {
  scen <- make_scen_row(mgmt = "smsy_goal", factorMSY = 0.75,
                        futureT = "yes", sdsel = 0, maxsel = 5.5,
                        ageT = 1, sizeT = 1, sexT = 1)
  cfg <- SizeShiftsImplicationsV2:::build_config(
    params = "Ohlberger", scen_row = scen,
    j = 1L, k = 1L, seednum = 1L, log_dir = tempdir()
  )
  expect_equal(cfg$harvmgmt, "smsy_goal")
  expect_equal(cfg$factorMSY, 0.75)
  expect_equal(cfg$futureT, "yes")
})
