if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(Clonotrace)
  test_check("Clonotrace")
} else {
  message("Skipping testthat suite because the suggested package `testthat` is unavailable.")
}
