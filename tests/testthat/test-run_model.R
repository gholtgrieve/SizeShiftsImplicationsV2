# Tests for run_model() output structure, reproducibility, and known-bug regressions.
# These tests use a minimal config (nyi=5, nyh=20, ny=40) to keep runtime short.

pkg <- "SizeShiftsImplicationsV2"

# ---------------------------------------------------------------------------
# Shared minimal config factory
# ---------------------------------------------------------------------------
make_minimal_config <- function(harvmgmt = "fix_harv_rate",
                                sim_recruits = "spawners",
                                reglength = 0,
                                regstr = 1,
                                seednum = 1L,
                                futureT = "no") {
  base <- utils::modifyList(
    SizeShiftsImplicationsV2:::default_config(),
    param_configs[["Ohlberger"]]
  )
  base$nyi          <- 5L
  base$nyh          <- 20L
  base$ny           <- 40L
  base$goalfreq     <- 20L
  base$harvmgmt     <- harvmgmt
  base$factorMSY    <- 1.0
  base$futureT      <- futureT
  base$sdsel        <- 0.0
  base$maxsel       <- 5.5
  base$sim_recruits <- sim_recruits
  base$reglength    <- reglength
  base$regstr       <- regstr
  base$j            <- 1L
  base$k            <- 1L
  base$seednum      <- as.integer(seednum)
  base$log_dir      <- tempdir()
  SizeShiftsImplicationsV2:::validate_config(base)
}


# ============================================================
# Output structure
# ============================================================

test_that("run_model returns a named list with expected top-level fields", {
  cfg <- make_minimal_config()
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  expected <- c("para", "sr_sim", "fec", "egg", "S_msy", "data", "obs",
                "ret_by_age", "meanSaA", "propfemale",
                "selectivities_by_age", "MSY_Goals", "impl_errors")
  expect_named(out, expected, ignore.order = TRUE)
})

test_that("run_model $data contains population dynamics columns", {
  cfg <- make_minimal_config()
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  expect_true(all(c("Ret", "Rec", "Esc", "Harv") %in% names(out$data)))
})

test_that("run_model $data has no negative Esc or Harv", {
  cfg <- make_minimal_config()
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  expect_true(all(out$data$Esc  >= 0, na.rm = TRUE))
  expect_true(all(out$data$Harv >= 0, na.rm = TRUE))
})

test_that("run_model Ret = Esc + Harv (mass balance) for each year", {
  cfg <- make_minimal_config()
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  d <- out$data
  residual <- abs(d$Ret - (d$Esc + d$Harv))
  # Allow small integer rounding discrepancy (<=1 fish per year)
  expect_true(all(residual <= 1, na.rm = TRUE))
})

test_that("run_model $obs has obsEsc, obsHarv, obsRet columns", {
  cfg <- make_minimal_config()
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  expect_true(all(c("obsEsc", "obsHarv", "obsRet") %in% names(out$obs)))
})

test_that("run_model $para returns named SR parameter vector", {
  cfg <- make_minimal_config()
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  expect_named(out$para, c("alpha.low", "alpha.high", "beta", "maxr"))
  expect_true(all(out$para > 0))
})

test_that("run_model $meanSaA has correct dimensions (nyr x nage)", {
  cfg <- make_minimal_config()
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  # meanSaA should have nage = length(ages) = 9 columns
  expect_equal(ncol(out$meanSaA), length(cfg$ages))
})


# ============================================================
# Reproducibility (CRN)
# ============================================================

test_that("run_model is reproducible with the same seed", {
  cfg <- make_minimal_config(seednum = 42L)
  out1 <- SizeShiftsImplicationsV2:::run_model(cfg)
  out2 <- SizeShiftsImplicationsV2:::run_model(cfg)
  expect_equal(out1$data$Rec, out2$data$Rec)
  expect_equal(out1$data$Esc, out2$data$Esc)
})

test_that("run_model gives different results with different seeds", {
  out1 <- SizeShiftsImplicationsV2:::run_model(make_minimal_config(seednum = 1L))
  out2 <- SizeShiftsImplicationsV2:::run_model(make_minimal_config(seednum = 2L))
  expect_false(identical(out1$data$Rec, out2$data$Rec))
})


# ============================================================
# Management modes
# ============================================================

test_that("fix_harv_rate management produces non-zero harvest in most years", {
  cfg <- make_minimal_config(harvmgmt = "fix_harv_rate")
  out <- SizeShiftsImplicationsV2:::run_model(cfg)
  frac_harvested <- mean(out$data$Harv > 0, na.rm = TRUE)
  expect_gt(frac_harvested, 0.5)
})

test_that("smsy_goal management runs without error", {
  cfg <- make_minimal_config(harvmgmt = "smsy_goal")
  expect_no_error(SizeShiftsImplicationsV2:::run_model(cfg))
})

test_that("umsy_goal management runs without error", {
  cfg <- make_minimal_config(harvmgmt = "umsy_goal")
  expect_no_error(SizeShiftsImplicationsV2:::run_model(cfg))
})


# ============================================================
# sim_recruits modes
# ============================================================

test_that("sim_recruits = 'eggmass' runs without error", {
  cfg <- make_minimal_config(sim_recruits = "eggmass")
  expect_no_error(SizeShiftsImplicationsV2:::run_model(cfg))
})

test_that("sim_recruits = 'fecundity' runs without error", {
  cfg <- make_minimal_config(sim_recruits = "fecundity")
  expect_no_error(SizeShiftsImplicationsV2:::run_model(cfg))
})

test_that("sim_recruits = 'eggmass' and 'spawners' give different recruits", {
  out_s <- SizeShiftsImplicationsV2:::run_model(make_minimal_config(sim_recruits = "spawners",  seednum = 10L))
  out_e <- SizeShiftsImplicationsV2:::run_model(make_minimal_config(sim_recruits = "eggmass",   seednum = 10L))
  expect_false(identical(out_s$data$Rec, out_e$data$Rec))
})


# ============================================================
# KNOWN BUG REGRESSION TESTS
#
# These tests are marked skip() to document the desired (correct)
# behavior without blocking the rest of the test suite.
# To activate a test after fixing the bug:
#   1. Remove the skip() line.
#   2. Run the test and confirm it passes.
# ============================================================

test_that("BUG1 (regime shift): high-productivity regime produces more recruitment", {
  # Bug location: run_model.R:347-349
  # When a regime shift occurs, alpha is set to the OLD regime before aset
  # flips. This means alpha never changes in the shift year -- the productivity
  # change is always delayed by one shift cycle.
  #
  # Correct behavior: when a shift is triggered, alpha should immediately
  # reflect the NEW regime (alpha.high if switching to high, alpha.low if low).
  #
  # Test: with reglength=1 (shift guaranteed every year), regstr=10 makes
  # alpha.high = 10 * alpha.low. Starting in the LOW regime, the FIRST shift
  # should put the model in the HIGH regime immediately. Over many seeds, the
  # HIGH-regime years should produce noticeably more recruitment than a
  # no-shift baseline (reglength=0).

  make_cfg <- function(reglength, regstr, seed) {
    make_minimal_config(reglength = reglength, regstr = regstr, seednum = seed)
  }

  n <- 30L
  seeds <- seq_len(n) * 7L

  mean_rec_noshift <- mean(vapply(seeds, function(s) {
    mean(SizeShiftsImplicationsV2:::run_model(make_cfg(0, 1, s))$data$Rec, na.rm = TRUE)
  }, numeric(1)))

  mean_rec_shift <- mean(vapply(seeds, function(s) {
    mean(SizeShiftsImplicationsV2:::run_model(make_cfg(1, 10, s))$data$Rec, na.rm = TRUE)
  }, numeric(1)))

  # With very strong regimes alternating every year, mean productivity should
  # be substantially higher than the no-shift baseline.
  expect_gt(mean_rec_shift, mean_rec_noshift * 1.5)
})


test_that("BUG2 (init fecundity): RepOut is consistent between init and main sim years", {
  # Bug location: run_model.R:235-239
  # In the initialization loop (y in 1:nyi), sizes_y is built from MALE age
  # composition and proportions, then scaled to TOTAL escapement. This is used
  # to compute fecundity/eggmass for recruits. The main simulation loop
  # (section 7, lines 460-468) correctly uses female-only escapement.
  #
  # Correct behavior: sizes_y in init years should reflect female fish only,
  # matching the approach in section 7.
  #
  # This test checks that RepOut per-spawner (a proxy for fecundity consistency)
  # does not show a discontinuity between the last init year and first main year.
  # A large jump would indicate the two sections are computing RepOut differently.
  skip("BUG2 not yet fixed -- run_model.R:235-239: init-year sizes_y uses male age composition scaled to total escapement. See code_review.txt for details.")

  cfg <- make_minimal_config(sim_recruits = "eggmass")
  out <- SizeShiftsImplicationsV2:::run_model(cfg)

  obs <- out$obs
  obs <- obs[!is.na(obs$RepOut) & obs$RepOut > 0, ]

  nyi <- cfg$nyi
  # RepOut per escapement as a size proxy
  rep_per_esc <- obs$RepOut / obs$obsEsc

  # Last init year vs first main year: ratio should be within 2x
  last_init  <- rep_per_esc[nyi]
  first_main <- rep_per_esc[nyi + 1L]
  ratio <- first_main / last_init
  expect_lt(ratio, 2)
  expect_gt(ratio, 0.5)
})
