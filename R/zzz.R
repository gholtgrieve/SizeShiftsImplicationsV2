#' @noRd
#' Package startup configuration
#'
#' This function is automatically called by R when the package is loaded.
#' It configures numerical libraries to use a single thread to ensure
#' reproducibility and to avoid conflicts when running parallel simulations.

.onLoad <- function(libname, pkgname) {
  # ---------------------------------------------------------------------------
  # Limit numerical libraries to a single thread
  #
  # Many R packages (or underlying C/Fortran libraries like OpenBLAS or MKL)
  # use multi-threading by default. This can lead to:
  #   - non-reproducible floating-point behavior,
  #   - over-subscription of CPU cores during parallel execution,
  #   - inconsistent runtime performance.
  #
  # To avoid these issues, we explicitly set OMP and MKL thread limits
  # and—if available—also limit the number of BLAS threads via RhpcBLASctl.
  # ---------------------------------------------------------------------------

  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    # Sets the number of BLAS threads to 1 (if using OpenBLAS, Accelerate, etc.)
    RhpcBLASctl::blas_set_num_threads(1)
  }
}
