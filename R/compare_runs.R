#' Compare Two SSI Run Outputs
#'
#' Compares all components of two `ssi_run` objects, including nested results per scenario and iteration.
#' Assumes same number of scenarios and iterations. Reports mismatches and saves them to a CSV.
#'
#' @param run1 First `ssi_run` object.
#' @param run2 Second `ssi_run` object.
#' @param tolerance Numeric tolerance for numeric comparisons.
#' @param verbose Logical. Print detailed mismatch info.
#' @param csv_path File path to save mismatches as CSV (default: "ssi_run_mismatches.csv").
#' @keywords internal
#'
.compare_ssi_runs <- function(run1, run2, tolerance = 1e-8, verbose = TRUE, csv_path = "ssi_run_mismatches.csv") {
  stopifnot(inherits(run1, "ssi_run"), inherits(run2, "ssi_run"))
  
  # Compare metadata
  if (!identical(run1$meta$niter, run2$meta$niter) ||
      !identical(run1$meta$nscen, run2$meta$nscen)) {
    stop("Mismatch in number of iterations or scenarios.")
  }
  
  nscen <- run1$meta$nscen
  niter <- run1$meta$niter
  all_passed <- TRUE
  mismatch_log <- data.frame(
    scenario = integer(),
    iteration = integer(),
    field = character(),
    value1 = character(),
    value2 = character(),
    stringsAsFactors = FALSE
  )
  
  # Helper to compare two objects with tolerance
  compare_obj <- function(a, b, label, j, k, field) {
    tryCatch({
      if (is.numeric(a) && is.numeric(b)) {
        if (!isTRUE(all.equal(a, b, tolerance = tolerance))) {
          if (verbose) cat(sprintf("Mismatch in %s\n", label))
          mismatch_log <<- rbind(mismatch_log, data.frame(
            scenario = j, iteration = k, field = field,
            value1 = toString(round(a, 5)[1:5]),
            value2 = toString(round(b, 5)[1:5]),
            stringsAsFactors = FALSE
          ))
          return(FALSE)
        }
      } else if (!identical(a, b)) {
        if (verbose) cat(sprintf("Mismatch in %s\n", label))
        mismatch_log <<- rbind(mismatch_log, data.frame(
          scenario = j, iteration = k, field = field,
          value1 = toString(a)[1],
          value2 = toString(b)[1],
          stringsAsFactors = FALSE
        ))
        return(FALSE)
      }
      TRUE
    }, error = function(e) {
      if (verbose) cat(sprintf("Error comparing %s: %s\n", label, e$message))
      mismatch_log <<- rbind(mismatch_log, data.frame(
        scenario = j, iteration = k, field = field,
        value1 = "<error>",
        value2 = "<error>",
        stringsAsFactors = FALSE
      ))
      FALSE
    })
  }
  
  for (j in seq_len(nscen)) {
    for (k in seq_len(niter)) {
      for (field in names(run1$results)) {
        x1 <- run1$results[[field]][[j]][[k]]
        x2 <- run2$results[[field]][[j]][[k]]
        label <- sprintf("Scenario %d - Iteration %d - %s", j, k, field)
        if (!compare_obj(x1, x2, label, j, k, field)) all_passed <- FALSE
      }
    }
  }
  
  if (nrow(mismatch_log) > 0) {
    utils::write.csv(mismatch_log, file = csv_path, row.names = FALSE)
    if (verbose) cat(sprintf("\nSaved mismatch log to: %s\n", csv_path))
  }
  
  if (all_passed && verbose) cat("\nAll components matched successfully.\n")
  invisible(all_passed)
}
