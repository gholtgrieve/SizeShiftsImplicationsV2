# Tests for internal model helper functions in model_helpers.R
# Access via ::: because these are not exported.

pkg <- "SizeShiftsImplicationsV2"


# ============================================================
# .calc_ricker
# ============================================================

test_that(".calc_ricker returns named list with rec and eps", {
  out <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0.4, last.eps = 0, seed = 1L
  )
  expect_type(out, "list")
  expect_named(out, c("rec", "eps"))
  expect_true(is.numeric(out$rec))
  expect_true(is.numeric(out$eps))
})

test_that(".calc_ricker produces positive recruits for positive spawners", {
  out <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0.4, last.eps = 0, seed = 42L
  )
  expect_gt(out$rec, 0)
})

test_that(".calc_ricker is reproducible with a seed", {
  args <- list(spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
               rho = 0.4, last.eps = 0, seed = 99L)
  out1 <- do.call(SizeShiftsImplicationsV2:::.calc_ricker, args)
  out2 <- do.call(SizeShiftsImplicationsV2:::.calc_ricker, args)
  expect_equal(out1$rec, out2$rec)
  expect_equal(out1$eps, out2$eps)
})

test_that(".calc_ricker gives different results with different seeds", {
  out1 <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0.4, last.eps = 0, seed = 1L
  )
  out2 <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0.4, last.eps = 0, seed = 2L
  )
  expect_false(out1$rec == out2$rec)
})

test_that(".calc_ricker AR(1): eps autocorrelates with last.eps", {
  # With rho=1 and no noise variance, eps should track last.eps exactly.
  # Here we just verify that a non-zero last.eps shifts eps relative to zero.
  out0 <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0.4, last.eps = 0, seed = 7L
  )
  out1 <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0.4, last.eps = 1.0, seed = 7L
  )
  # eps should differ because last.eps differs
  expect_false(out0$eps == out1$eps)
})

test_that(".calc_ricker with rho=0 ignores last.eps", {
  out_a <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0, last.eps = 0, seed = 5L
  )
  out_b <- SizeShiftsImplicationsV2:::.calc_ricker(
    spawn = 5000, sigma = 0.3, alpha = 5, beta = 5e-5,
    rho = 0, last.eps = 999, seed = 5L
  )
  expect_equal(out_a$rec, out_b$rec)
})


# ============================================================
# .calc_agecomp
# ============================================================

test_that(".calc_agecomp returns a vector with length equal to ages", {
  ages <- 1:9
  out <- SizeShiftsImplicationsV2:::.calc_agecomp(
    ages = ages, recruits = 1000, meanage = 5, sdage = 0.6, seed = 1L
  )
  expect_length(out, length(ages))
})

test_that(".calc_agecomp counts sum to approximately recruits", {
  out <- SizeShiftsImplicationsV2:::.calc_agecomp(
    ages = 1:9, recruits = 1000, meanage = 5, sdage = 0.6, seed = 1L
  )
  # Rounding means exact equality isn't guaranteed, but should be very close
  expect_lte(abs(sum(out) - 1000), 9)  # at most one fish off per age class
})

test_that(".calc_agecomp returns non-negative integers", {
  out <- SizeShiftsImplicationsV2:::.calc_agecomp(
    ages = 1:9, recruits = 500, meanage = 5, sdage = 0.6, seed = 1L
  )
  expect_true(all(out >= 0))
  expect_true(all(out == floor(out)))
})

test_that(".calc_agecomp is reproducible with a seed", {
  args <- list(ages = 1:9, recruits = 800, meanage = 5, sdage = 0.6, seed = 77L)
  out1 <- do.call(SizeShiftsImplicationsV2:::.calc_agecomp, args)
  out2 <- do.call(SizeShiftsImplicationsV2:::.calc_agecomp, args)
  expect_equal(out1, out2)
})

test_that(".calc_agecomp concentrates fish near meanage", {
  # With tight sd, most fish should be at ages near meanage = 5
  out <- SizeShiftsImplicationsV2:::.calc_agecomp(
    ages = 1:9, recruits = 10000, meanage = 5, sdage = 0.5, seed = 1L
  )
  # Ages 4, 5, 6 should hold the majority of fish
  expect_gt(sum(out[4:6]), sum(out) * 0.7)
})

test_that(".calc_agecomp returns zeros for recruits = 0", {
  out <- SizeShiftsImplicationsV2:::.calc_agecomp(
    ages = 1:9, recruits = 0, meanage = 5, sdage = 0.6, seed = 1L
  )
  expect_true(all(out == 0))
})


# ============================================================
# .calc_selectivity
# ============================================================

test_that(".calc_selectivity returns a numeric vector", {
  sizes <- c(400, 500, 600, 700, 800)
  out <- SizeShiftsImplicationsV2:::.calc_selectivity(
    size = sizes, meshsize = 5.5, s = 0.204
  )
  expect_type(out, "double")
  expect_length(out, length(sizes))
})

test_that(".calc_selectivity returns non-negative values", {
  sizes <- seq(300, 900, by = 100)
  out <- SizeShiftsImplicationsV2:::.calc_selectivity(
    size = sizes, meshsize = 5.5, s = 0.204
  )
  expect_true(all(out >= 0, na.rm = TRUE))
})

test_that(".calc_selectivity is unimodal (peak near mesh-matched size)", {
  sizes <- seq(300, 1000, by = 10)
  out <- SizeShiftsImplicationsV2:::.calc_selectivity(
    size = sizes, meshsize = 5.5, s = 0.204
  )
  peak_idx <- which.max(out)
  # Peak should be somewhere in the middle, not at the extremes
  expect_gt(peak_idx, 5)
  expect_lt(peak_idx, length(sizes) - 5)
})


# ============================================================
# .calc_reprod_output
# ============================================================

test_that(".calc_reprod_output returns list with fecundity and eggmass", {
  allom <- c(size_ref = 800, fec_ref = 6600, b_fec = 2.4,
             egg_ref = 916, b_eggs = 4.8)
  out <- SizeShiftsImplicationsV2:::.calc_reprod_output(size = 700, allometry = allom)
  expect_type(out, "list")
  expect_named(out, c("fecundity", "eggmass"))
})

test_that(".calc_reprod_output returns log-scale values (fecundity > 0 after exp)", {
  allom <- c(size_ref = 800, fec_ref = 6600, b_fec = 2.4,
             egg_ref = 916, b_eggs = 4.8)
  out <- SizeShiftsImplicationsV2:::.calc_reprod_output(size = 700, allometry = allom)
  expect_gt(exp(out$fecundity), 0)
  expect_gt(exp(out$eggmass), 0)
})

test_that(".calc_reprod_output is monotonically increasing in size", {
  allom <- c(size_ref = 800, fec_ref = 6600, b_fec = 2.4,
             egg_ref = 916, b_eggs = 4.8)
  sizes <- c(500, 600, 700, 800, 900)
  fecs <- vapply(sizes, function(s) {
    SizeShiftsImplicationsV2:::.calc_reprod_output(s, allom)$fecundity
  }, numeric(1))
  expect_true(all(diff(fecs) > 0))
})

test_that(".calc_reprod_output at reference size matches reference values", {
  allom <- c(size_ref = 800, fec_ref = 6600, b_fec = 2.4,
             egg_ref = 916, b_eggs = 4.8)
  out <- SizeShiftsImplicationsV2:::.calc_reprod_output(size = 800, allometry = allom)
  # unname() strips names inherited from the allometry vector before comparing
  expect_equal(unname(exp(out$fecundity)), 6600, tolerance = 1e-6)
  expect_equal(unname(exp(out$eggmass)),   916,  tolerance = 1e-6)
})
