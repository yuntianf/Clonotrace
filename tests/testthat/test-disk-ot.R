test_that("disk-backed exact OT equals dense shortest-path OT", {
  graph <- igraph::make_graph(c(
    1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 1, 2, 5
  ), directed = FALSE)
  igraph::E(graph)$weight <- c(0.7, 1.1, 0.4, 0.9, 1.3, 0.8, 0.6)
  cell_ids <- paste0("cell", seq_len(igraph::vcount(graph)))
  membership <- matrix(
    c(
      0.8, 0.0, 0.2,
      0.2, 0.4, 0.0,
      0.0, 0.6, 0.1,
      0.0, 0.0, 0.7,
      0.3, 0.2, 0.0,
      0.0, 0.1, 0.4
    ), nrow = 6, byrow = TRUE,
    dimnames = list(cell_ids, c("cloneA", "cloneB", "cloneC"))
  )
  path <- tempfile(fileext = ".bin")
  checkpoint_dir <- tempfile("disk-ot-parts-")
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)

  store <- build_graph_distance_store(
    graph, path, cell_ids = cell_ids, block_size = 2L,
    verbose = FALSE
  )
  result <- graph_clone_ot_disk(
    store, membership, target_cell_limit = 3L,
    checkpoint_dir = checkpoint_dir, verbose = FALSE
  )

  full_distance <- igraph::distances(graph, weights = igraph::E(graph)$weight)
  pairs <- utils::combn(seq_len(ncol(membership)), 2L)
  expected <- data.frame(
    group1 = colnames(membership)[pairs[1L, ]],
    group2 = colnames(membership)[pairs[2L, ]],
    dis = apply(pairs, 2L, function(pair) {
      source <- which(membership[, pair[[1L]]] >= 0.05)
      target <- which(membership[, pair[[2L]]] >= 0.05)
      clone_2_ot(
        full_distance[source, target, drop = FALSE],
        membership[source, pair[[1L]]], membership[target, pair[[2L]]]
      )
    }),
    stringsAsFactors = FALSE
  )
  expect_equal(result, expected, tolerance = 1e-12)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2L)
  parallel_checkpoint_dir <- tempfile("disk-ot-parallel-parts-")
  on.exit(unlink(parallel_checkpoint_dir, recursive = TRUE), add = TRUE)
  parallel_result <- graph_clone_ot_disk(
    store, membership, target_cell_limit = 3L, cores = 2L,
    checkpoint_dir = parallel_checkpoint_dir, verbose = FALSE
  )
  expect_equal(parallel_result, expected, tolerance = 1e-12)
  expect_identical(attr(parallel_result, "diagnostics"), NULL)
  future::plan(future::sequential)

  store_time <- file.info(path)$mtime
  checkpoint_files <- list.files(checkpoint_dir, full.names = TRUE)
  checkpoint_times <- file.info(checkpoint_files)$mtime
  resumed_store <- build_graph_distance_store(
    graph, path, cell_ids = cell_ids, block_size = 2L,
    verbose = FALSE
  )
  resumed_result <- graph_clone_ot_disk(
    resumed_store, membership, target_cell_limit = 3L,
    checkpoint_dir = checkpoint_dir, verbose = FALSE
  )
  expect_equal(result, resumed_result)
  expect_identical(file.info(path)$mtime, store_time)
  expect_identical(file.info(checkpoint_files)$mtime, checkpoint_times)
})

test_that("distance stores reject stale and disconnected inputs", {
  graph <- igraph::make_ring(6)
  igraph::E(graph)$weight <- rep(1, igraph::ecount(graph))
  path <- tempfile(fileext = ".bin")
  on.exit(unlink(c(path, paste0(path, ".metadata.rds"))), add = TRUE)
  build_graph_distance_store(graph, path, block_size = 2L, verbose = FALSE)
  changed <- graph
  igraph::E(changed)$weight[1] <- 2
  expect_error(
    build_graph_distance_store(changed, path, block_size = 2L, verbose = FALSE),
    "does not match"
  )

  disconnected <- igraph::make_graph(c(1, 2, 3, 4), directed = FALSE)
  expect_error(
    build_graph_distance_store(disconnected, tempfile(fileext = ".bin"), verbose = FALSE),
    "connected"
  )
})
