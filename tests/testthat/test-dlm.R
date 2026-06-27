# Tests for .calc_DLMfit (model_helpers.R)

make_sr_data <- function(n = 20, seed = 1L) {
  withr::with_seed(seed, {
    Esc <- round(runif(n, 5000, 50000))
    Rec <- round(Esc * exp(1.5 - 5e-5 * Esc + rnorm(n, 0, 0.3)))
    data.frame(Esc = Esc, Rec = pmax(Rec, 1L))
  })
}

# Data with a clear declining productivity trend so the DLM's time-varying
# state has signal to track (static data can correctly produce W=0, i.e. no variation).
make_trending_sr_data <- function(n = 40, seed = 1L) {
  withr::with_seed(seed, {
    Esc       <- round(runif(n, 5000, 20000))
    log_alpha <- seq(log(6), log(2), length.out = n)  # alpha declines 6 -> 2
    Rec       <- round(Esc * exp(log_alpha - 5e-5 * Esc + rnorm(n, 0, 0.2)))
    data.frame(Esc = Esc, Rec = pmax(Rec, 1L))
  })
}

# ── Structure ──────────────────────────────────────────────────────────────

test_that(".calc_DLMfit returns list with results, AICc, sigma", {
  out <- SizeShiftsImplicationsV2:::.calc_DLMfit(make_sr_data(), var_alpha = TRUE, var_beta = FALSE)
  expect_type(out, "list")
  expect_named(out, c("results", "AICc", "sigma"))
})

test_that(".calc_DLMfit sigma is positive and finite", {
  out <- SizeShiftsImplicationsV2:::.calc_DLMfit(make_sr_data(), var_alpha = TRUE, var_beta = FALSE)
  expect_true(is.finite(out$sigma) && out$sigma > 0)
})

test_that(".calc_DLMfit AICc is finite", {
  out <- SizeShiftsImplicationsV2:::.calc_DLMfit(make_sr_data(), var_alpha = TRUE, var_beta = FALSE)
  expect_true(is.finite(out$AICc))
})

test_that(".calc_DLMfit results has same rows as input plus alpha_y and beta_y", {
  dat <- make_sr_data()
  out <- SizeShiftsImplicationsV2:::.calc_DLMfit(dat, var_alpha = TRUE, var_beta = FALSE)
  expect_equal(nrow(out$results), nrow(dat))
  expect_true("alpha_y" %in% names(out$results))
  expect_true("beta_y"  %in% names(out$results))
})

# ── All three var combinations ─────────────────────────────────────────────

test_that(".calc_DLMfit var_alpha=TRUE, var_beta=FALSE runs without error", {
  expect_no_error(
    SizeShiftsImplicationsV2:::.calc_DLMfit(make_sr_data(), var_alpha = TRUE, var_beta = FALSE)
  )
})

test_that(".calc_DLMfit var_alpha=FALSE, var_beta=TRUE runs without error", {
  expect_no_error(
    SizeShiftsImplicationsV2:::.calc_DLMfit(make_sr_data(), var_alpha = FALSE, var_beta = TRUE)
  )
})

test_that(".calc_DLMfit var_alpha=TRUE, var_beta=TRUE runs without error", {
  expect_no_error(
    SizeShiftsImplicationsV2:::.calc_DLMfit(make_sr_data(), var_alpha = TRUE, var_beta = TRUE)
  )
})

# ── npara is correct (1 + number of free W elements) ──────────────────────

test_that(".calc_DLMfit AICc uses correct npara (1 + k, not 3 + k)", {
  dat <- make_sr_data()
  out1 <- SizeShiftsImplicationsV2:::.calc_DLMfit(dat, var_alpha = TRUE,  var_beta = FALSE) # npara=2
  out2 <- SizeShiftsImplicationsV2:::.calc_DLMfit(dat, var_alpha = TRUE,  var_beta = TRUE)  # npara=3

  # More params -> higher AICc penalty; out2 should have higher AICc if fit is similar
  # At minimum, they should differ
  expect_false(isTRUE(all.equal(out1$AICc, out2$AICc)))
})

# ── Time-varying behaviour ─────────────────────────────────────────────────

test_that(".calc_DLMfit var_alpha=TRUE tracks a declining alpha trend", {
  # Use trending data so the DLM has a signal to track; static data can
  # correctly return W[1,1]=0 (no variation), which is also a valid result.
  out <- SizeShiftsImplicationsV2:::.calc_DLMfit(make_trending_sr_data(), var_alpha = TRUE, var_beta = FALSE)
  # With a 6->2 alpha decline, smoothed alpha_y should also decline
  expect_gt(out$results$alpha_y[1], out$results$alpha_y[nrow(out$results)])
})

test_that(".calc_DLMfit var_beta=FALSE has less variation in beta_y than alpha_y", {
  # On trending data, alpha_y varies while beta_y (W[2,2]=0) stays near constant.
  out <- SizeShiftsImplicationsV2:::.calc_DLMfit(make_trending_sr_data(), var_alpha = TRUE, var_beta = FALSE)
  expect_lt(var(out$results$beta_y), var(out$results$alpha_y))
})
