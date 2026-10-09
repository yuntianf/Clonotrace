test_that("label_spreading_bootstrap uses the blocked backend by default", {
  adjacency <- Matrix::bandSparse(
    12, 12, k = c(-1, 0, 1),
    diagonals = list(rep(1, 11), rep(1, 12), rep(1, 11))
  )
  labels <- c("a", rep(NA_character_, 4), "b", rep(NA_character_, 5), "a")
  h5_dir <- tempfile("clonotrace-blocked-")
  dir.create(h5_dir)
  on.exit(unlink(h5_dir, recursive = TRUE), add = TRUE)

  set.seed(11)
  default <- label_spreading_bootstrap(
    adjacency, labels, sample_n = 3, block_size = 1,
    tmpdir = h5_dir, verbose = FALSE
  )
  set.seed(11)
  explicit <- label_spreading_bootstrap_blocked(
    adjacency, labels, sample_n = 3, block_size = 1,
    tmpdir = h5_dir, verbose = FALSE
  )

  expect_equal(default$prob, explicit$prob, tolerance = 1e-12)
  expect_equal(default$deviance, explicit$deviance, tolerance = 1e-12)
  expect_equal(dim(default$prob), c(12L, 2L))
  expect_true(all(is.finite(default$prob)))
  expect_true(all(is.finite(default$deviance)))
})

test_that("the legacy dense bootstrap remains an explicit option", {
  adjacency <- Matrix::bandSparse(
    6, 6, k = c(-1, 0, 1),
    diagonals = list(rep(1, 5), rep(1, 6), rep(1, 5))
  )
  labels <- c("a", NA, "b", NA, "a", NA)

  set.seed(7)
  result <- label_spreading_bootstrap(
    adjacency, labels, sample_n = 2, backend = "dense", verbose = FALSE
  )

  expect_equal(dim(result$prob), c(6L, 2L))
  expect_true(all(is.finite(result$prob)))
  expect_true(all(is.finite(result$deviance)))
})

test_that("a caller-owned blocked reference cache is preserved", {
  adjacency <- Matrix::bandSparse(
    8, 8, k = c(-1, 0, 1),
    diagonals = list(rep(1, 7), rep(1, 8), rep(1, 7))
  )
  labels <- c("a", NA, NA, "b", NA, NA, "a", NA)
  preindexed <- as.integer(as.factor(labels))
  h5_dir <- tempfile("clonotrace-reference-")
  dir.create(h5_dir)
  on.exit(unlink(h5_dir, recursive = TRUE), add = TRUE)
  reference_path <- file.path(h5_dir, "reference.h5")

  label_spreading_blocked(
    adjacency, preindexed, label_n = 2L, block_size = 1,
    outfile = reference_path, verbose = FALSE
  )
  normalize_hdf5_rows(reference_path, "prob", block_cols = 1)
  set.seed(12)
  label_spreading_bootstrap(
    adjacency, labels, sample_n = 2, block_size = 1,
    refer_h5 = reference_path, tmpdir = h5_dir, verbose = FALSE
  )

  expect_true(file.exists(reference_path))
})
