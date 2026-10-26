# Independent reference: the original statistic evaluated one clone pair at
# a time, without the new support indexing or source-distance reuse.
nn_reference <- function(graph, probability, threshold, k) {
  membership <- probability
  membership[membership < threshold] <- 0
  membership[membership >= threshold] <- 1
  rows <- lapply(seq_len(ncol(membership) - 1L), function(i) {
    do.call(rbind, lapply(seq.int(i + 1L, ncol(membership)), function(j) {
      x <- igraph::distances(
        graph, v = which(membership[, i] > 0),
        to = which(membership[, j] > 0)
      )
      # Preserve legacy singleton and k-larger-than-support behavior as well.
      x <- x[seq_len(nrow(x)), seq_len(ncol(x))]
      value <- if (is.null(nrow(x))) mean(x) else {
        mean(c(
          mean(apply(x, 1, function(z) mean(sort(z)[seq_len(k)]))),
          mean(apply(x, 2, function(z) mean(sort(z)[seq_len(k)])))
        ))
      }
      c(i, j, value)
    }))
  })
  out <- do.call(rbind, rows)
  colnames(out) <- c("group1", "group2", "dis")
  out
}

test_that("cached and fallback NN preserve the original distance statistic", {
  graph <- igraph::make_ring(16)
  graph <- igraph::set_edge_attr(graph, "weight", value = seq(0.2, 1, length.out = 16))
  probability <- matrix(0, 16, 5)
  probability[c(1, 4, 7), 1] <- c(0.1, 0.3, 1)
  probability[c(2, 4, 9, 10), 2] <- c(0.09, 0.2, 0.8, 1)
  probability[12, 3] <- 1
  probability[c(3, 5, 11, 16), 4] <- 0.5
  probability[c(6, 8, 14), 5] <- 0.7
  for (input in list(probability, Matrix::Matrix(probability, sparse = TRUE))) {
    for (k in c(1, 2, 5)) {
      reference <- nn_reference(graph, input, 0.1, k)
      for (budget in c(0, 0.000001, 64)) {
        for (block in c(1L, 3L)) {
          observed <- graph_clone_nn(
            graph, input, k = k, clone_block_size = block,
            distance_cache_mb = budget
          )
          expect_equal(observed, reference, tolerance = 1e-12)
        }
      }
    }
  }
})

test_that("disconnected graph and two-clone cases retain their results", {
  graph <- igraph::make_empty_graph(6, directed = FALSE)
  graph <- igraph::add_edges(graph, c(1, 2, 2, 3, 4, 5, 5, 6))
  probability <- matrix(c(1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0), 6, 2)
  for (budget in c(0, 64)) {
    expect_equal(
      graph_clone_nn(graph, probability, distance_cache_mb = budget),
      nn_reference(graph, probability, 0.1, 2)
    )
  }
})

test_that("invalid cache budgets and empty clone supports fail clearly", {
  graph <- igraph::make_ring(4)
  probability <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), 4, 2)
  for (budget in list(-1, Inf, NA_real_, numeric(), c(1, 2), "64")) {
    expect_error(
      graph_clone_nn(graph, probability, distance_cache_mb = budget),
      "finite non-negative"
    )
  }
  probability[, 2] <- 0
  expect_error(graph_clone_nn(graph, probability), "no cells")
})

test_that("clone_disance forwards the cache budget without changing outputs", {
  skip_on_cran()
  set.seed(26)
  embedding <- matrix(rnorm(60), 20, 3,
                      dimnames = list(paste0("cell", 1:20), NULL))
  probability <- Matrix::sparseMatrix(
    i = 1:20, j = rep(1:4, each = 5), x = 1,
    dimnames = list(rownames(embedding), paste0("clone", 1:4))
  )
  output <- tempfile("nn-cache-test-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  cached <- clone_disance(
    embedding, probability, output, graph_k = 5, revised = FALSE,
    distance_cache_mb = 64
  )
  fallback <- clone_disance(
    embedding, probability, output, graph_k = 5, revised = FALSE,
    distance_cache_mb = 0, overwrite = TRUE
  )
  expect_equal(cached, fallback)
})

test_that("NN distance reuse works in installed-package future workers", {
  skip_on_cran()
  skip_if_not_installed("future")
  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  graph <- igraph::make_ring(12)
  probability <- matrix(0, 12, 4)
  probability[cbind(1:12, rep(1:4, each = 3))] <- 1
  expected <- nn_reference(graph, probability, 0.1, 2)
  future::plan(future::multisession, workers = 2)
  for (budget in c(0, 64)) {
    expect_equal(
      graph_clone_nn(graph, probability, distance_cache_mb = budget),
      expected
    )
  }
})
