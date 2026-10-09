args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript tools/build_pkgdown.R <package-root> [library]")
}

package_root <- normalizePath(args[[1L]], mustWork = TRUE)
system_libraries <- c(
  "/usr/local/lib/R/site-library",
  "/usr/lib/R/site-library",
  "/usr/lib/R/library"
)
system_libraries <- system_libraries[dir.exists(system_libraries)]
.libPaths(c(system_libraries, .libPaths()))
if (length(args) >= 2L) {
  .libPaths(c(normalizePath(args[[2L]], mustWork = TRUE), .libPaths()))
}

old_working_directory <- setwd(package_root)
on.exit(setwd(old_working_directory), add = TRUE)

pkgdown::build_article(
  name = "Clonotrace",
  pkg = ".",
  data = list(source = NULL),
  lazy = FALSE
)
