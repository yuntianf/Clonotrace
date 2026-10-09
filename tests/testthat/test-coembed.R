test_that("cell_clone_coembed does not rely on a global membership object", {
  set.seed(2)
  cells <- paste0("cell", 1:12)
  clones <- paste0("clone", 1:4)
  cell_embedding <- matrix(rnorm(24), nrow = 12, dimnames = list(cells, NULL))
  clone_embedding <- matrix(rnorm(8), nrow = 4, dimnames = list(clones, NULL))
  membership <- Matrix::sparseMatrix(
    i = 1:12, j = rep(1:4, each = 3), x = 1,
    dimnames = list(cells, clones)
  )

  expect_error(
    cell_clone_coembed(cell_embedding, clone_embedding),
    "cell_clone_prob"
  )
  expect_error(
    cell_clone_coembed(cell_embedding, clone_embedding, membership[, -1], cell_k = 2, clone_k = 2),
    "does not contain"
  )
})
