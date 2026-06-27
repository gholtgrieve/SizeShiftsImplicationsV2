# Tests for convert_legacy_to_ssi()

make_legacy_rdata <- function(path, nscen = 2L, niter = 3L) {
  # Build minimal synthetic legacy objects matching the expected structure
  make_df <- function() data.frame(
    Year = seq_len(50L),
    Esc  = round(runif(50, 5000, 20000)),
    Rec  = round(runif(50, 10000, 40000)),
    Harv = round(runif(50, 2000, 10000))
  )

  data.list  <- lapply(seq_len(nscen), function(j) lapply(seq_len(niter), function(k) make_df()))
  parms.list <- list(ny = 110L, nyh = 50L, nyi = 10L,
                     goalfreq = 10L, firstrev = 20L, seednum = 1L)
  scen       <- data.frame(scen_num = seq_len(nscen), mgmt = rep("smsy_goal", nscen))

  # save() uses bare variable names; backtick names work
  e <- new.env(parent = emptyenv())
  e[["data.list"]]  <- data.list
  e[["parms.list"]] <- parms.list
  e[["scen"]]       <- scen
  save(list = ls(e), envir = e, file = path)
  invisible(path)
}


# ── Error handling ────────────────────────────────────────────────────────

test_that("convert_legacy_to_ssi errors on missing file", {
  expect_error(
    convert_legacy_to_ssi("/nonexistent/path/file.RData"),
    "File not found"
  )
})

test_that("convert_legacy_to_ssi errors when data.list is missing", {
  bad_file <- tempfile(fileext = ".RData")
  on.exit(unlink(bad_file), add = TRUE)
  e <- new.env(parent = emptyenv())
  e[["parms.list"]] <- list(ny = 110L, nyh = 50L)
  save(list = "parms.list", envir = e, file = bad_file)
  expect_error(convert_legacy_to_ssi(bad_file), "data.list")
})

test_that("convert_legacy_to_ssi errors when parms.list is missing", {
  bad_file <- tempfile(fileext = ".RData")
  on.exit(unlink(bad_file), add = TRUE)
  e <- new.env(parent = emptyenv())
  e[["data.list"]] <- list(list(data.frame(Year = 1:5)))
  save(list = "data.list", envir = e, file = bad_file)
  expect_error(convert_legacy_to_ssi(bad_file), "parms.list")
})


# ── Output structure ──────────────────────────────────────────────────────

local({
  .legacy_file <<- tempfile(fileext = ".RData")
  make_legacy_rdata(.legacy_file, nscen = 2L, niter = 3L)
  .legacy_run  <<- convert_legacy_to_ssi(.legacy_file)
})

test_that("convert_legacy_to_ssi returns an ssi_run object", {
  expect_s3_class(.legacy_run, "ssi_run")
})

test_that("convert_legacy_to_ssi has required top-level fields", {
  expect_true(all(c("meta", "scenarios", "results", "parameters") %in% names(.legacy_run)))
})

test_that("convert_legacy_to_ssi meta has correct nscen and niter", {
  expect_equal(.legacy_run$meta$nscen, 2L)
  expect_equal(.legacy_run$meta$niter, 3L)
})

test_that("convert_legacy_to_ssi scenarios has scen_num column", {
  expect_true("scen_num" %in% names(.legacy_run$scenarios))
})

test_that("convert_legacy_to_ssi scenarios are named scen_1, scen_2, ...", {
  expect_equal(names(.legacy_run$results$data.list), c("scen_1", "scen_2"))
})

test_that("convert_legacy_to_ssi parameters has ny, nyh, nyi, goalfreq", {
  p <- .legacy_run$parameters
  expect_true(all(c("ny", "nyh", "nyi", "goalfreq") %in% names(p)))
  expect_equal(p$ny,  110L)
  expect_equal(p$nyh,  50L)
  expect_equal(p$nyi,  10L)
})

test_that("convert_legacy_to_ssi data.list has correct nscen x niter structure", {
  dl <- .legacy_run$results$data.list
  expect_length(dl, 2L)          # nscen = 2
  expect_length(dl[[1]], 3L)     # niter = 3
  expect_s3_class(dl[[1]][[1]], "data.frame")
})

test_that("convert_legacy_to_ssi drops No column if present", {
  file2 <- tempfile(fileext = ".RData")
  on.exit(unlink(file2), add = TRUE)
  # add a No column to the scenario table
  e <- new.env(parent = emptyenv())
  e[["data.list"]]  <- list(list(data.frame(Year = 1:5)))
  e[["parms.list"]] <- list(ny = 110L, nyh = 50L, nyi = 10L,
                             goalfreq = 10L, firstrev = 20L, seednum = 1L)
  e[["scen"]]       <- data.frame(No = 1L, mgmt = "smsy_goal")
  save(list = ls(e), envir = e, file = file2)
  run2 <- convert_legacy_to_ssi(file2)
  expect_false("No" %in% names(run2$scenarios))
})
