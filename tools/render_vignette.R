args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript tools/render_vignette.R <package-root> [library]")
}

package_root <- normalizePath(args[[1L]], mustWork = TRUE)
if (length(args) >= 2L) {
  .libPaths(c(normalizePath(args[[2L]], mustWork = TRUE), .libPaths()))
}

source_file <- file.path(package_root, "vignettes", "Clonotrace.Rmd")
output_file <- file.path(package_root, "vignettes", "Clonotrace.html")

rmarkdown::render(
  input = source_file,
  output_file = output_file,
  clean = TRUE,
  envir = new.env(parent = globalenv())
)
