# Reproducible clone-distance resource benchmark.
#
# Run from an installed Clonotrace source tree with:
# Rscript inst/benchmarks/clone_distance_benchmark.R output.csv

output_path <- if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[[1L]]
} else {
  "clone-distance-benchmark.csv"
}

library(Clonotrace)
library(Matrix)

make_fixture <- function(n_cells, n_clones, dimensions = 10L, seed = 1L) {
  set.seed(seed)
  embedding <- matrix(rnorm(n_cells * dimensions), nrow = n_cells)
  rownames(embedding) <- paste0("cell", seq_len(n_cells))
  clone_id <- rep(seq_len(n_clones), length.out = n_cells)
  membership <- sparseMatrix(
    i = seq_len(n_cells), j = clone_id, x = 1,
    dims = c(n_cells, n_clones),
    dimnames = list(rownames(embedding), paste0("clone", seq_len(n_clones)))
  )
  list(embedding = embedding, membership = membership)
}

run_case <- function(n_cells, n_clones, exact, block_size = 25L, seed = 1L) {
  fixture <- make_fixture(n_cells, n_clones, seed = seed)
  output_dir <- tempfile("clonotrace-benchmark-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  gc(reset = TRUE)
  mode_arguments <- if (exact) {
    list(force_exact = TRUE, cores = 1L, verbose = FALSE)
  } else {
    list(clone_block_size = block_size, verbose = FALSE)
  }
  elapsed <- system.time({
    result <- do.call(
      clone_disance,
      c(
        list(
          embedding = fixture$embedding,
          cell_clone_prob = fixture$membership,
          outpath = output_dir,
          graph_k = min(10L, n_cells - 1L),
          exact = exact
        ),
        mode_arguments
      )
    )
  })
  gc_stats <- gc()
  data.frame(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = as.character(packageVersion("Clonotrace")),
    r_version = R.version.string,
    platform = R.version$platform,
    n_cells = n_cells,
    n_clones = n_clones,
    exact = exact,
    clone_block_size = if (exact) NA_integer_ else block_size,
    elapsed_seconds = unname(elapsed[["elapsed"]]),
    max_vcells_mb = gc_stats["Vcells", "max used"] * 8 / 1024^2,
    result_rows = nrow(result)
  )
}

cases <- list(
  list(n_cells = 200L, n_clones = 10L, exact = FALSE, block_size = 5L),
  list(n_cells = 1000L, n_clones = 50L, exact = FALSE, block_size = 10L),
  list(n_cells = 200L, n_clones = 10L, exact = TRUE, block_size = NA_integer_)
)

results <- do.call(rbind, lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  do.call(run_case, c(case, seed = i))
}))
write.csv(results, output_path, row.names = FALSE)
message("Wrote benchmark results to ", normalizePath(output_path, mustWork = FALSE))
