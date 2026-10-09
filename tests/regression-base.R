# Base-R regression checks for development environments where the suggested
# testthat package is unavailable. The installed package is loaded by
# R CMD check before this file is sourced.

library(Clonotrace)

graph <- igraph::make_ring(4)
graph <- igraph::set_edge_attr(graph, "weight", value = c(1, 2, 4, 8))
corrected_graph <- Clonotrace:::.clonotrace_snn_similarity_to_cost(graph)
stopifnot(
  isTRUE(all.equal(
    as.numeric(igraph::edge_attr(corrected_graph, "weight")),
    exp(-c(1, 2, 4, 8))
  )),
  isTRUE(all.equal(
    as.numeric(igraph::edge_attr(graph, "weight")),
    c(1, 2, 4, 8)
  ))
)

long <- data.frame(
  group1 = c("clone1", "clone1", "clone2"),
  group2 = c("clone2", "clone3", "clone3"),
  dis = c(1, 2, 3)
)
square <- long2square(long, "group1", "group2", "dis", symmetric = TRUE)
stopifnot(
  is.double(square),
  identical(rownames(square), colnames(square)),
  square["clone1", "clone2"] == 1,
  square["clone2", "clone1"] == 1
)

duplicate <- data.frame(from = c("a", "a"), to = c("b", "b"), value = c(1, 1))
duplicate_error <- try(long2square(duplicate, "from", "to", "value"), silent = TRUE)
stopifnot(inherits(duplicate_error, "try-error"))

cache_path <- tempfile(fileext = ".rds")
on.exit(unlink(cache_path), add = TRUE)
metadata <- list(input = "first")
Clonotrace:::.clonotrace_atomic_save_rds(list(metadata = metadata, value = 42), cache_path)
stopifnot(identical(Clonotrace:::.clonotrace_read_cache(cache_path, metadata, "test"), 42))
stale_error <- try(
  Clonotrace:::.clonotrace_read_cache(cache_path, list(input = "second"), "test"),
  silent = TRUE
)
stopifnot(inherits(stale_error, "try-error"))

set.seed(1)
coordinates <- matrix(rnorm(40), nrow = 10)
distances <- as.matrix(stats::dist(coordinates))
dimnames(distances) <- list(paste0("clone", 1:10), paste0("clone", 1:10))
membership <- leiden_dis(distances, k = 3, if_umap = FALSE)
stopifnot(is.factor(membership), identical(names(membership), rownames(distances)))
set.seed(2)
layout <- leiden_dis(distances, k = 3, if_umap = TRUE)
stopifnot(
  is.data.frame(layout),
  all(c("umap1", "umap2", "cluster", "clone") %in% colnames(layout)),
  identical(layout$clone, rownames(distances))
)

set.seed(3)
cells <- paste0("cell", 1:18)
clones <- paste0("clone", 1:3)
embedding <- rbind(
  matrix(rnorm(12, -2), ncol = 2),
  matrix(rnorm(12, 0), ncol = 2),
  matrix(rnorm(12, 2), ncol = 2)
)
rownames(embedding) <- cells
clone_prob <- Matrix::sparseMatrix(
  i = 1:18, j = rep(1:3, each = 6), x = 1,
  dimnames = list(cells, clones)
)
output_dir <- tempfile("clonotrace-smoke-")
dir.create(output_dir)
on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
long_distance <- clone_disance(
  embedding, clone_prob, output_dir,
  graph_k = 10, exact = FALSE, k = 1, clone_block_size = 2,
  revised = FALSE
)
reloaded_distance <- clone_disance(
  embedding, clone_prob, output_dir,
  graph_k = 10, exact = FALSE, k = 1, clone_block_size = 2,
  revised = FALSE
)
stopifnot(
  identical(long_distance, reloaded_distance),
  all(long_distance$group1 %in% clones),
  all(long_distance$group2 %in% clones),
  file.exists(file.path(output_dir, "cell_graph.rds")),
  file.exists(file.path(output_dir, "clone_graph_dis.rds")),
  !file.exists(file.path(output_dir, "groups_graph_dis.rds"))
)
changed_parameter_error <- try(
  clone_disance(
    embedding, clone_prob, output_dir,
    graph_k = 10, exact = FALSE, k = 2, clone_block_size = 2,
    revised = FALSE
  ),
  silent = TRUE
)
stopifnot(inherits(changed_parameter_error, "try-error"))

exact_output_dir <- tempfile("clonotrace-exact-smoke-")
dir.create(exact_output_dir)
on.exit(unlink(exact_output_dir, recursive = TRUE), add = TRUE)
exact_distance <- clone_disance(
  embedding[1:12, , drop = FALSE],
  clone_prob[1:12, 1:2, drop = FALSE],
  exact_output_dir,
  graph_k = 8,
  exact = TRUE,
  cores = 1,
  distance_block_size = 4,
  target_cell_limit = 8,
  verbose = FALSE,
  revised = FALSE
)
stopifnot(
  nrow(exact_distance) == 1L,
  all(exact_distance$group1 %in% clones[1:2]),
  all(exact_distance$group2 %in% clones[1:2]),
  is.finite(exact_distance$dis),
  file.exists(file.path(exact_output_dir, "clone_graph_dis_exact.rds")),
  file.exists(file.path(exact_output_dir, "cell_graph_distance_exact.bin")),
  file.exists(file.path(exact_output_dir, "cell_graph_distance_exact.bin.metadata.rds")),
  dir.exists(file.path(exact_output_dir, "exact_ot_checkpoints"))
)

legacy_exact_output_dir <- tempfile("clonotrace-memory-exact-smoke-")
dir.create(legacy_exact_output_dir)
on.exit(unlink(legacy_exact_output_dir, recursive = TRUE), add = TRUE)
legacy_exact_distance <- clone_disance(
  embedding[1:12, , drop = FALSE],
  clone_prob[1:12, 1:2, drop = FALSE],
  legacy_exact_output_dir,
  graph_k = 8,
  exact = TRUE,
  exact_backend = "memory",
  cores = 1,
  cache = 100,
  verbose = FALSE,
  revised = FALSE
)
stopifnot(isTRUE(all.equal(exact_distance$dis, legacy_exact_distance$dis)))

# Revised-distance regression: preserve the raw distance, accept named observed
# clone identities, and keep the former local-rescaling behavior available as
# an explicit opt-in.
set.seed(4)
revision_cells <- paste0("revision_cell", seq_len(160L))
revision_clones <- paste0("revision_clone", seq_len(8L))
revision_labels <- rep(revision_clones, each = 20L)
revision_embedding <- cbind(
  PC1 = rep(seq(-2, 2, length.out = 8L), each = 20L) + rnorm(160L, sd = 0.4),
  PC2 = rnorm(160L),
  PC3 = rnorm(160L)
)
rownames(revision_embedding) <- revision_cells
raw_square <- as.matrix(stats::dist(matrix(rnorm(8L * 3L), nrow = 8L)))
dimnames(raw_square) <- list(revision_clones, revision_clones)
revised_square <- revise_clone_distance(
  raw_square,
  revision_embedding,
  clone_labels = stats::setNames(revision_labels, revision_cells),
  local_k = 2L,
  local_rescale = TRUE,
  mean_weight = 0.25,
  verbose = FALSE
)
stopifnot(
  isTRUE(isSymmetric(unname(revised_square))),
  all(diag(revised_square) == 0),
  identical(dim(revised_square), c(8L, 8L)),
  identical(attr(revised_square, "mean_weight"), 0.25),
  isTRUE(all.equal(attr(revised_square, "clonotrace"), raw_square))
)

revision_clusters <- data.frame(
  resolution_a = rep(1:4, length.out = length(revision_cells)),
  resolution_b = rep(1:5, each = length(revision_cells) / 5L),
  row.names = revision_cells
)
automatic_square <- revise_clone_distance(
  raw_square,
  revision_embedding,
  clone_labels = stats::setNames(revision_labels, revision_cells),
  clusters = revision_clusters,
  local_k = 2L,
  s_grid = c(0, 0.5, 1),
  verbose = FALSE
)
stopifnot(
  is.finite(attr(automatic_square, "mean_weight")),
  length(attr(automatic_square, "picks")) == 2L,
  identical(
    dim(attr(automatic_square, "selection_correlations")),
    c(3L, 2L)
  )
)

internal_ladder_square <- revise_clone_distance(
  raw_square,
  revision_embedding,
  clone_labels = stats::setNames(revision_labels, revision_cells),
  local_k = 2L,
  resolutions = c(0.3, 0.5),
  cluster_k = 10L,
  s_grid = c(0, 0.5, 1),
  verbose = FALSE
)
stopifnot(
  isTRUE(isSymmetric(internal_ladder_square)),
  length(attr(internal_ladder_square, "picks")) == 2L,
  identical(
    dim(attr(internal_ladder_square, "selection_correlations")),
    c(3L, 2L)
  )
)

# The blocked/HDF5 implementation is now the public default, while the dense
# backend remains available explicitly for legacy reproduction.
bootstrap_adj <- Matrix::bandSparse(
  8, 8, k = c(-1, 0, 1),
  diagonals = list(rep(1, 7), rep(1, 8), rep(1, 7))
)
bootstrap_labels <- c("a", NA, NA, "b", NA, NA, "a", NA)
bootstrap_tmp <- tempfile("clonotrace-blocked-regression-")
dir.create(bootstrap_tmp)
on.exit(unlink(bootstrap_tmp, recursive = TRUE), add = TRUE)
set.seed(5)
blocked_default <- label_spreading_bootstrap(
  bootstrap_adj, bootstrap_labels,
  sample_n = 2L, block_size = 1L,
  tmpdir = bootstrap_tmp, verbose = FALSE
)
stopifnot(
  identical(dim(blocked_default$prob), c(8L, 2L)),
  all(is.finite(blocked_default$prob)),
  all(is.finite(blocked_default$deviance))
)

message("Base Clonotrace regression checks passed.")
