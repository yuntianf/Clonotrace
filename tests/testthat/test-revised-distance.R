test_that("revised distance accepts named and positional long inputs", {
  set.seed(21)
  n_cells <- 240L
  clone_names <- paste0("clone", seq_len(8L))
  cells <- paste0("cell", seq_len(n_cells))
  clone <- rep(clone_names, each = n_cells / length(clone_names))
  embedding <- cbind(
    PC1 = rep(seq(-2, 2, length.out = 8L), each = 30L) + rnorm(n_cells, sd = 0.5),
    PC2 = rnorm(n_cells),
    PC3 = rnorm(n_cells)
  )
  rownames(embedding) <- cells
  clone_labels <- data.frame(cell = cells, clone = clone)

  raw <- as.matrix(stats::dist(matrix(rnorm(8L * 4L), nrow = 8L)))
  dimnames(raw) <- list(clone_names, clone_names)
  upper <- which(upper.tri(raw), arr.ind = TRUE)
  named_long <- data.frame(
    group1 = clone_names[upper[, 1L]],
    group2 = clone_names[upper[, 2L]],
    dis = raw[upper]
  )
  positional_long <- transform(
    named_long,
    group1 = match(group1, clone_names),
    group2 = match(group2, clone_names)
  )

  named <- revise_clone_distance(
    named_long, embedding, clone_labels = clone_labels,
    clone_names = clone_names, local_k = 2L,
    local_rescale = TRUE,
    mean_weight = 0.25, verbose = FALSE
  )
  positional <- revise_clone_distance(
    positional_long, embedding, clone_labels = clone_labels,
    clone_names = clone_names, local_k = 2L,
    local_rescale = TRUE,
    mean_weight = 0.25, verbose = FALSE
  )

  expect_equal(unname(named), unname(positional), tolerance = 1e-12)
  expect_true(isSymmetric(unname(named)))
  expect_equal(unname(diag(named)), rep(0, 8L))
  expect_equal(attr(named, "mean_weight"), 0.25)

  local_scale <- vapply(seq_len(8L), function(i) sort(raw[i, -i])[2L], numeric(1))
  local <- raw / sqrt(outer(local_scale, local_scale))
  centroids <- rowsum(embedding, clone) / as.numeric(table(clone)[clone_names])
  centroid <- as.matrix(stats::dist(centroids[clone_names, , drop = FALSE]))
  expected <- 0.75 * local / stats::median(local[upper.tri(local)]) +
    0.25 * centroid / stats::median(centroid[upper.tri(centroid)])
  diag(expected) <- 0
  expect_equal(as.vector(named), as.vector(expected), tolerance = 1e-10)
})

test_that("local clone-distance rescaling is optional and off by default", {
  set.seed(2101)
  clone_names <- paste0("clone", seq_len(8L))
  cells <- paste0("cell", seq_len(160L))
  clone <- rep(clone_names, each = 20L)
  embedding <- matrix(
    rnorm(160L * 3L), ncol = 3L,
    dimnames = list(cells, NULL)
  )
  raw <- as.matrix(stats::dist(matrix(rnorm(8L * 3L), nrow = 8L)))
  dimnames(raw) <- list(clone_names, clone_names)

  unscaled <- revise_clone_distance(
    raw, embedding,
    clone_labels = stats::setNames(clone, cells),
    mean_weight = 0, verbose = FALSE
  )
  scaled <- revise_clone_distance(
    raw, embedding,
    clone_labels = stats::setNames(clone, cells),
    local_k = 2L, local_rescale = TRUE,
    mean_weight = 0, verbose = FALSE
  )
  expected <- raw / stats::median(raw[upper.tri(raw)])

  expect_equal(as.vector(unscaled), as.vector(expected), tolerance = 1e-12)
  expect_false(attr(unscaled, "local_rescale"))
  expect_true(attr(scaled, "local_rescale"))
  expect_false(isTRUE(all.equal(as.vector(unscaled), as.vector(scaled))))
})

test_that("automatic revised-distance selection is fate-blind and auditable", {
  set.seed(22)
  cells <- paste0("cell", seq_len(160L))
  clone_names <- paste0("clone", seq_len(8L))
  clone <- rep(clone_names, each = 20L)
  embedding <- matrix(rnorm(160L * 3L), ncol = 3L, dimnames = list(cells, NULL))
  raw <- as.matrix(stats::dist(matrix(rnorm(8L * 3L), nrow = 8L)))
  dimnames(raw) <- list(clone_names, clone_names)
  clusters <- data.frame(
    resolution_a = rep(1:4, length.out = 160L),
    resolution_b = rep(1:5, each = 32L),
    row.names = cells
  )
  revised <- revise_clone_distance(
    raw, embedding,
    clone_labels = stats::setNames(clone, cells),
    clusters = clusters, local_k = 2L,
    s_grid = c(0, 0.5, 1), verbose = FALSE
  )

  expect_true(attr(revised, "mean_weight") %in% c(0, 0.25, 0.5, 0.75, 1))
  expect_length(attr(revised, "picks"), 2L)
  expect_equal(
    dim(attr(revised, "selection_correlations")),
    c(3L, 2L)
  )
  expect_equal(attr(revised, "clonotrace"), raw)
})

test_that("automatic revised-distance selection builds an internal Leiden ladder", {
  set.seed(2201)
  cells <- paste0("cell", seq_len(160L))
  clone_names <- paste0("clone", seq_len(8L))
  clone <- rep(clone_names, each = 20L)
  embedding <- matrix(
    rnorm(160L * 3L), ncol = 3L,
    dimnames = list(cells, NULL)
  )
  raw <- as.matrix(stats::dist(matrix(rnorm(8L * 3L), nrow = 8L)))
  dimnames(raw) <- list(clone_names, clone_names)

  revised <- revise_clone_distance(
    raw, embedding,
    clone_labels = stats::setNames(clone, cells),
    local_k = 2L,
    resolutions = c(0.3, 0.5),
    cluster_k = 10L,
    s_grid = c(0, 0.5, 1),
    verbose = FALSE
  )

  expect_true(isSymmetric(revised))
  expect_length(attr(revised, "picks"), 2L)
  expect_equal(
    dim(attr(revised, "selection_correlations")),
    c(3L, 2L)
  )
})

test_that("clone_disance returns revised and raw distances together", {
  set.seed(23)
  cells <- paste0("cell", seq_len(240L))
  clone_names <- paste0("clone", seq_len(8L))
  clone <- sample(rep(clone_names, each = 30L))
  embedding <- matrix(rnorm(240L * 4L), ncol = 4L, dimnames = list(cells, NULL))
  membership <- Matrix::sparseMatrix(
    i = seq_along(cells), j = match(clone, clone_names), x = 1,
    dims = c(length(cells), length(clone_names)),
    dimnames = list(cells, clone_names)
  )
  output <- tempfile("clonotrace-revised-distance-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)

  revised <- clone_disance(
    embedding, membership, output,
    graph_k = 10L, exact = FALSE, revised = TRUE,
    clone_labels = stats::setNames(clone, cells),
    local_k = 2L, local_rescale = TRUE, mean_weight = 0.25,
    k = 3L, clone_block_size = 4L, verbose = FALSE
  )
  raw <- clone_disance(
    embedding, membership, output,
    graph_k = 10L, exact = FALSE, revised = FALSE,
    k = 3L, clone_block_size = 4L, verbose = FALSE
  )

  expect_true(all(c("group1", "group2", "dis", "dis_clonotrace") %in% names(revised)))
  expect_equal(revised$dis_clonotrace, raw$dis)
  expect_equal(attr(revised, "mean_weight"), 0.25)
  expect_true(all(is.finite(revised$dis)))
})
