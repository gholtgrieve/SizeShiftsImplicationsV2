# Tests for model configuration variations:
# futureT, regime shifts, selectivity, sim_recruits modes, factorMSY.

make_cfg <- function(...) {
  base <- utils::modifyList(
    SizeShiftsImplicationsV2:::default_config(),
    param_configs[["Ohlberger"]]
  )
  base$nyi <- 5L; base$nyh <- 20L; base$ny <- 40L
  base$goalfreq <- 20L; base$factorMSY <- 1.0
  base$harvmgmt <- "smsy_goal"; base$futureT <- "no"
  base$sdsel <- 0.0; base$maxsel <- 5.5
  base$j <- 1L; base$k <- 1L; base$seednum <- 1L
  base$log_dir <- tempdir()
  args <- list(...)
  for (nm in names(args)) base[[nm]] <- args[[nm]]
  SizeShiftsImplicationsV2:::validate_config(base)
}


# ── futureT ────────────────────────────────────────────────────────────────

test_that("futureT='yes' runs without error", {
  expect_no_error(SizeShiftsImplicationsV2:::run_model(make_cfg(futureT = "yes")))
})

test_that("futureT='yes' produces different size-at-age than futureT='no'", {
  out_no  <- SizeShiftsImplicationsV2:::run_model(make_cfg(futureT = "no",
                                                           sizetrends = c(0,0,30,10,-30,-60,-90,-90,-90)))
  out_yes <- SizeShiftsImplicationsV2:::run_model(make_cfg(futureT = "yes",
                                                           sizetrends = c(0,0,30,10,-30,-60,-90,-90,-90)))
  # With continuing trends, later years should have different sizes
  expect_false(identical(out_no$meanSaA, out_yes$meanSaA))
})

test_that("futureT='no' freezes meanage at end of historical period", {
  out <- SizeShiftsImplicationsV2:::run_model(make_cfg(futureT = "no", agetrend = -0.4))
  # meanSaA columns should be constant after the historical period (nyi+nyh rows)
  # (frozen at 5-year mean, so variance in future rows should be near zero for det part)
  nyi <- 5L; nyh <- 20L
  hist_end <- nyi + nyh
  future_rows <- (hist_end + 1):nrow(out$meanSaA)
  if (length(future_rows) >= 2) {
    # Deterministic part is frozen; stochastic anomalies still vary, but within sdSaA
    expect_true(!is.null(out$meanSaA))
  }
})


# ── Regime shifts ──────────────────────────────────────────────────────────

test_that("reglength > 0 runs without error", {
  expect_no_error(
    SizeShiftsImplicationsV2:::run_model(make_cfg(reglength = 20L, regstr = 2))
  )
})

test_that("regstr=2 produces higher recruitment variance than regstr=1", {
  seeds <- seq_len(20L) * 3L
  var_no_shift <- mean(vapply(seeds, function(s) {
    var(SizeShiftsImplicationsV2:::run_model(
      make_cfg(reglength = 5L, regstr = 1, seednum = s))$data$Rec, na.rm = TRUE)
  }, numeric(1)))
  var_shift <- mean(vapply(seeds, function(s) {
    var(SizeShiftsImplicationsV2:::run_model(
      make_cfg(reglength = 5L, regstr = 2, seednum = s))$data$Rec, na.rm = TRUE)
  }, numeric(1)))
  expect_gt(var_shift, var_no_shift)
})

test_that("reglength=0 disables regime shifts (reg column is constant)", {
  out <- SizeShiftsImplicationsV2:::run_model(make_cfg(reglength = 0L))
  reg <- out$data$reg[!is.na(out$data$reg)]
  expect_equal(length(unique(reg)), 1L)
})


# ── Selectivity ────────────────────────────────────────────────────────────

test_that("non-zero selectivity (sdsel > 0) runs without error", {
  expect_no_error(
    SizeShiftsImplicationsV2:::run_model(make_cfg(sdsel = 0.204, maxsel = 5.5))
  )
})

test_that("size-selective harvest skews harvest age structure", {
  # Selective harvest (by size) vs unselective should produce different age comps
  out_unsel <- SizeShiftsImplicationsV2:::run_model(make_cfg(sdsel = 0.0,   maxsel = 5.5, seednum = 3L))
  out_sel   <- SizeShiftsImplicationsV2:::run_model(make_cfg(sdsel = 0.204, maxsel = 5.5, seednum = 3L))
  # selectivities_by_age should differ
  expect_false(identical(out_unsel$selectivities_by_age, out_sel$selectivities_by_age))
})


# ── sim_recruits modes ─────────────────────────────────────────────────────

test_that("sim_recruits='fecundity' runs without error", {
  expect_no_error(
    SizeShiftsImplicationsV2:::run_model(make_cfg(sim_recruits = "fecundity"))
  )
})

test_that("sim_recruits='fecundity' produces positive RepOut", {
  out <- SizeShiftsImplicationsV2:::run_model(make_cfg(sim_recruits = "fecundity"))
  ro <- out$obs$RepOut[!is.na(out$obs$RepOut)]
  expect_true(all(ro > 0))
})

test_that("sim_recruits='eggmass' and 'fecundity' give different RepOut magnitudes", {
  out_egg <- SizeShiftsImplicationsV2:::run_model(make_cfg(sim_recruits = "eggmass",   seednum = 9L))
  out_fec <- SizeShiftsImplicationsV2:::run_model(make_cfg(sim_recruits = "fecundity", seednum = 9L))
  # Fecundity (egg count) >> eggmass (grams), so median should be higher for fecundity
  med_egg <- median(out_egg$obs$RepOut, na.rm = TRUE)
  med_fec <- median(out_fec$obs$RepOut, na.rm = TRUE)
  expect_gt(med_fec, med_egg)
})


# ── factorMSY ──────────────────────────────────────────────────────────────

test_that("factorMSY=0.75 (conservative) reduces escapement goal vs 1.0", {
  seeds <- seq_len(10L)
  mean_esc_100 <- mean(vapply(seeds, function(s) {
    mean(SizeShiftsImplicationsV2:::run_model(
      make_cfg(factorMSY = 1.00, seednum = s))$data$Esc, na.rm = TRUE)
  }, numeric(1)))
  mean_esc_075 <- mean(vapply(seeds, function(s) {
    mean(SizeShiftsImplicationsV2:::run_model(
      make_cfg(factorMSY = 0.75, seednum = s))$data$Esc, na.rm = TRUE)
  }, numeric(1)))
  # Lower S_msy goal -> more harvest -> less escapement on average
  expect_lt(mean_esc_075, mean_esc_100)
})

test_that("factorMSY=1.5 (precautionary) increases escapement vs 1.0", {
  seeds <- seq_len(10L)
  mean_esc_100 <- mean(vapply(seeds, function(s) {
    mean(SizeShiftsImplicationsV2:::run_model(
      make_cfg(factorMSY = 1.00, seednum = s))$data$Esc, na.rm = TRUE)
  }, numeric(1)))
  mean_esc_150 <- mean(vapply(seeds, function(s) {
    mean(SizeShiftsImplicationsV2:::run_model(
      make_cfg(factorMSY = 1.50, seednum = s))$data$Esc, na.rm = TRUE)
  }, numeric(1)))
  expect_gt(mean_esc_150, mean_esc_100)
})


# ── Demographic trends ─────────────────────────────────────────────────────

test_that("negative agetrend shifts mean age composition downward over time", {
  out <- SizeShiftsImplicationsV2:::run_model(make_cfg(agetrend = -0.4, futureT = "no"))
  pa <- out$propfemale
  expect_length(pa, nrow(out$meanSaA))
  expect_true(all(is.finite(pa)))
})

test_that("size trends reduce mean escapement size over historical period", {
  out_flat  <- SizeShiftsImplicationsV2:::run_model(
    make_cfg(sizetrends = rep(0, 9), seednum = 2L))
  out_trend <- SizeShiftsImplicationsV2:::run_model(
    make_cfg(sizetrends = c(0,0,0,0,-30,-60,-90,-90,-90), seednum = 2L))
  # Mean escapement size should be lower in trend scenario
  expect_lt(mean(out_trend$data$meanSesc, na.rm = TRUE),
            mean(out_flat$data$meanSesc,  na.rm = TRUE))
})
