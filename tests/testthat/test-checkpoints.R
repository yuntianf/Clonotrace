test_that("checkpoint metadata rejects stale inputs", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  metadata <- list(input = "first")
  Clonotrace:::.clonotrace_atomic_save_rds(list(metadata = metadata, value = 42), path)

  expect_identical(Clonotrace:::.clonotrace_read_cache(path, metadata, "test"), 42)
  expect_error(
    Clonotrace:::.clonotrace_read_cache(path, list(input = "second"), "test"),
    "does not match"
  )
})

test_that("legacy checkpoints remain readable with a warning", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(data.frame(group1 = 1, group2 = 2, dis = 3), path)
  expect_warning(
    value <- Clonotrace:::.clonotrace_read_cache(path, list(), "test"),
    "legacy"
  )
  expect_s3_class(value, "data.frame")
})

test_that("clone distance validates paths and mode-specific arguments", {
  embedding <- matrix(rnorm(20), nrow = 10, dimnames = list(paste0("cell", 1:10), NULL))
  membership <- Matrix::sparseMatrix(
    i = 1:10, j = rep(1:2, each = 5), x = 1,
    dimnames = list(rownames(embedding), c("clone1", "clone2"))
  )
  expect_error(
    clone_disance(embedding, membership, tempfile(), graph_k = 2, exact = FALSE, cores = 2),
    "inapplicable"
  )
  expect_error(
    clone_disance(embedding, membership, tempfile(), graph_k = 2, exact = TRUE,
                   exact_cell_limit = 5),
    "safety limits"
  )
})
