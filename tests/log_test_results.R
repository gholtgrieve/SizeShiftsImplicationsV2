# Run this after devtools::test() to append the current results to the log.
# Usage (from project root):
#   source("tests/log_test_results.R")
#   log_test_results()
#
# Or call with explicit counts if logging a known-result run:
#   log_test_results(pass=90, fail=0, skip=1, warn=0, label="Bug 1 fix")

log_test_results <- function(pass = NULL, fail = NULL, skip = NULL, warn = NULL,
                             label = NULL,
                             log_path = "tests/test_history.csv") {

  if (is.null(pass) || is.null(fail) || is.null(skip)) {
    # Run the tests live and capture results
    res <- testthat::test_local(reporter = testthat::SilentReporter$new())
    s   <- as.data.frame(res)
    pass <- sum(s$passed,  na.rm = TRUE)
    fail <- sum(s$failed,  na.rm = TRUE)
    skip <- sum(s$skipped, na.rm = TRUE)
    warn <- sum(s$warning, na.rm = TRUE)
  }
  warn <- if (is.null(warn)) 0L else as.integer(warn)

  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    commit    = tryCatch(
      trimws(system("git rev-parse --short HEAD", intern = TRUE)),
      error = function(e) NA_character_
    ),
    label     = if (is.null(label)) NA_character_ else as.character(label),
    pass      = as.integer(pass),
    fail      = as.integer(fail),
    skip      = as.integer(skip),
    warn      = as.integer(warn),
    stringsAsFactors = FALSE
  )

  has_file <- file.exists(log_path)
  write.table(row, file = log_path, sep = ",",
              row.names = FALSE, col.names = !has_file,
              append = has_file, qmethod = "double")

  cat(sprintf("Logged: PASS=%d  FAIL=%d  SKIP=%d  WARN=%d  [%s]\n",
              pass, fail, skip, warn, row$commit))
  invisible(row)
}
