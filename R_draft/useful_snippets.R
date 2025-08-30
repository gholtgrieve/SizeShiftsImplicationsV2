


# Find a replace text string in package
files <- c(
  list.files("R", "\\.R$", full.names = TRUE),
  list.files("tests", "\\.R$", full.names = TRUE),
  list.files("vignettes", "\\.(Rmd|qmd|Rnw)$", full.names = TRUE)
)
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  txt <- gsub("\\bresolve_profile\\b", ".resolve_profile", txt, perl = TRUE)
  writeLines(txt, f)
}




#Find text string in package and return file name and line.
find_in_pkg <- function(pattern, dirs = c("R","tests","vignettes","inst","src"),
                        fixed = TRUE) {
  files <- unlist(lapply(dirs, function(d) if (dir.exists(d))
    list.files(d, pattern = "\\.(R|Rmd|qmd|Rnw|Rd|c|cpp|h|hpp)$",
               recursive = TRUE, full.names = TRUE) else character()))
  for (f in files) {
    x <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL)
    if (is.null(x)) next
    hits <- grep(pattern, x, fixed = fixed)
    if (length(hits)) {
      for (i in hits) cat(sprintf("%s:%d: %s\n", f, i, x[i]))
    }
  }
}

find_in_pkg("resolve_profile")  # change pattern as needed

#small test run
out <- run_scenarios("Ohlberger", params = "Ohlberger",
                     niter = 10, seed = "reproducible",
                     parallel = FALSE)
