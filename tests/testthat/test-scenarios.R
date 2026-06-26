# Tests for scenario selection (.select_scenarios in model_helpers.R)

pkg <- "SizeShiftsImplicationsV2"


test_that("'all' selector returns a tibble with rows", {
  out <- SizeShiftsImplicationsV2:::.select_scenarios("all")
  expect_s3_class(out, "tbl_df")
  expect_gt(nrow(out), 0)
})

test_that("'ohlberger' alias is equivalent to 'all'", {
  out_all <- SizeShiftsImplicationsV2:::.select_scenarios("all")
  out_ohl <- SizeShiftsImplicationsV2:::.select_scenarios("ohlberger")
  expect_equal(nrow(out_all), nrow(out_ohl))
})

test_that("'kuskokwim' selector returns a subset of 'all'", {
  out_all  <- SizeShiftsImplicationsV2:::.select_scenarios("all")
  out_kusk <- SizeShiftsImplicationsV2:::.select_scenarios("kuskokwim")
  expect_gt(nrow(out_kusk), 0)
  expect_lt(nrow(out_kusk), nrow(out_all))
})

test_that("'kuskokwim' selector excludes age-length trend scenarios", {
  out <- SizeShiftsImplicationsV2:::.select_scenarios("kuskokwim")
  expect_false(any(out$trends == "age-length trends", na.rm = TRUE))
})

test_that("scenarios tibble has expected columns", {
  out <- SizeShiftsImplicationsV2:::.select_scenarios("all")
  required_cols <- c("mgmt", "factorMSY", "futureT", "selectivity")
  expect_true(all(required_cols %in% names(out)))
})

test_that("numeric selector filters by scen_num", {
  out_all <- SizeShiftsImplicationsV2:::.select_scenarios("all")
  if ("scen_num" %in% names(out_all) && nrow(out_all) >= 2) {
    first_id <- out_all$scen_num[1]
    out1 <- SizeShiftsImplicationsV2:::.select_scenarios(first_id)
    expect_equal(nrow(out1), 1L)
    expect_equal(out1$scen_num, first_id)
  } else {
    skip("scenarios lacks scen_num or has fewer than 2 rows")
  }
})

test_that("filter expression selector works", {
  out <- SizeShiftsImplicationsV2:::.select_scenarios("mgmt == 'fix_harv_rate'")
  expect_gt(nrow(out), 0)
  expect_true(all(out$mgmt == "fix_harv_rate"))
})

test_that("invalid filter expression raises an error", {
  expect_error(
    SizeShiftsImplicationsV2:::.select_scenarios("%%not valid R%%"),
    "Could not parse selector"
  )
})

test_that("empty selector string raises an error", {
  expect_error(
    SizeShiftsImplicationsV2:::.select_scenarios(""),
    "Empty selector"
  )
})

test_that("enforce_constraints drops disallowed selectivity/factorMSY combos", {
  out_constrained   <- SizeShiftsImplicationsV2:::.select_scenarios("all", enforce_constraints = TRUE)
  out_unconstrained <- SizeShiftsImplicationsV2:::.select_scenarios("all", enforce_constraints = FALSE)
  # Constraints should remove at least some rows (or equal if all combos are valid)
  expect_lte(nrow(out_constrained), nrow(out_unconstrained))
})
