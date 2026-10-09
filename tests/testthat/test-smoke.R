test_that("small clone-distance workflow produces a square clustering input", {
  skip_on_cran()
  set.seed(3)
  cells <- paste0("cell", 1:18)
  clones <- paste0("clone", 1:3)
  embedding <- rbind(
    matrix(rnorm(12, -2), ncol = 2),
    matrix(rnorm(12, 0), ncol = 2),
    matrix(rnorm(12, 2), ncol = 2)
  )
  rownames(embedding) <- cells
  membership <- Matrix::sparseMatrix(
    i = 1:18, j = rep(1:3, each = 6), x = 1,
    dimnames = list(cells, clones)
  )
  output <- tempfile("clonotrace-smoke-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)

  long <- clone_disance(
    embedding, membership, output,
    graph_k = 10, exact = FALSE, k = 1, revised = FALSE
  )
  square <- long2square(long, "group1", "group2", "dis", symmetric = TRUE,
                        nodes = as.character(seq_along(clones)))
  diag(square) <- 0

  expect_true(file.exists(file.path(output, "cell_graph.rds")))
  expect_true(file.exists(file.path(output, "clone_graph_dis.rds")))
  expect_true(all(long$group1 %in% clones), all(long$group2 %in% clones))
  expect_true(isSymmetric(square))
  expect_identical(dim(square), c(3L, 3L))
})
