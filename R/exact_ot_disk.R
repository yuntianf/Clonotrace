.clonotrace_distance_store_metadata_path <- function(path) {
  paste0(path, ".metadata.rds")
}

.clonotrace_read_distance_store <- function(store) {
  path <- if (inherits(store, "clonotrace_distance_store")) store$path else store
  if (!is.character(path) || length(path) != 1L || is.na(path) || path == "") {
    stop("`store` must be a distance-store object or one non-empty binary-store path.")
  }
  path <- normalizePath(path, mustWork = FALSE)
  metadata_path <- .clonotrace_distance_store_metadata_path(path)
  if (!file.exists(path) || !file.exists(metadata_path)) {
    stop("The binary distance store and its metadata sidecar must both exist.")
  }
  metadata <- tryCatch(readRDS(metadata_path), error = function(error) {
    stop("Unable to read distance-store metadata: ", conditionMessage(error))
  })
  required <- c("format_version", "specification", "completed_blocks", "complete")
  if (!is.list(metadata) || !all(required %in% names(metadata)) ||
      !identical(metadata$format_version, 2L)) {
    stop("The distance-store metadata is missing or unsupported.")
  }
  if (!isTRUE(metadata$complete)) {
    stop("The distance store is incomplete; resume `build_graph_distance_store()` first.")
  }
  expected_n <- length(metadata$specification$cell_ids)
  expected_bytes <- as.double(expected_n) * as.double(expected_n) * 8
  observed_bytes <- file.info(path)$size
  if (!is.finite(observed_bytes) || observed_bytes != expected_bytes) {
    stop("The binary distance-store size does not match the metadata dimensions.")
  }
  structure(
    list(path = path, metadata_path = metadata_path, metadata = metadata),
    class = "clonotrace_distance_store"
  )
}

.clonotrace_read_binary_rows <- function(connection, rows, column_count) {
  rows <- as.integer(rows)
  if (!length(rows) || anyNA(rows) || is.unsorted(rows, strictly = TRUE)) {
    stop("Internal error: distance-store row indices must be non-empty, unique, and sorted.")
  }
  result <- matrix(NA_real_, nrow = length(rows), ncol = column_count)
  run_id <- cumsum(c(TRUE, diff(rows) != 1L))
  runs <- split(rows, run_id)
  for (i in seq_along(runs)) {
    run <- runs[[i]]
    offset <- as.double(run[[1L]] - 1L) * as.double(column_count) * 8
    seek(connection, where = offset, origin = "start")
    values <- readBin(
      connection, what = double(), n = length(run) * column_count,
      size = 8L, endian = "little"
    )
    if (length(values) != length(run) * column_count) {
      stop("Unexpected end of the binary distance store.")
    }
    destination <- which(run_id == i)
    result[destination, ] <- matrix(
      values, nrow = length(run), ncol = column_count, byrow = TRUE
    )
  }
  result
}

#' @title Build a Disk-Backed Exact Cell-Graph Distance Store
#' @description Calculate cell-graph shortest-path distances once and write them
#'   to a resumable, row-major binary matrix for exact clone-pair OT calculations.
#' @details
#' Distances are calculated only between the selected endpoint cells, but paths
#' may traverse every vertex in `graph`. Values are stored as double precision
#' in a full symmetric matrix. Each source-cell row is contiguous on disk so
#' clone supports can be fetched with direct, bounded reads. The row-block
#' computation is restartable through an atomic metadata sidecar.
#'
#' @param graph An undirected `igraph` cell graph with finite, non-negative
#'   `weight` edge costs. An unweighted graph uses unit costs.
#' @param path Output path ending in `.bin`.
#' @param cells Optional endpoint cells to store, supplied as graph indices or
#'   identifiers from `cell_ids`. Default: every graph vertex.
#' @param cell_ids Optional unique identifiers in graph vertex order. If omitted,
#'   graph vertex names are used; unnamed graphs use character vertex indices.
#' @param block_size Number of source rows calculated and written per resumable
#'   block (default: 256).
#' @param overwrite Logical. Replace an existing store and metadata sidecar.
#' @param verbose Logical. Report completed row blocks.
#'
#' @return A `clonotrace_distance_store` object containing the paths and
#'   validated metadata.
#' @export
build_graph_distance_store <- function(graph, path, cells = NULL,
                                       cell_ids = NULL, block_size = 256L,
                                       overwrite = FALSE, verbose = TRUE) {
  if (!inherits(graph, "igraph")) stop("`graph` must be an igraph object.")
  if (igraph::is_directed(graph)) {
    stop("Disk-backed exact graph distances currently require an undirected graph.")
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) || path == "") {
    stop("`path` must be one non-empty binary-store path.")
  }
  if (!grepl("[.]bin$", path, ignore.case = TRUE)) {
    stop("`path` must end in `.bin`.")
  }
  validate_positive_integer <- function(value, label) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        value < 1L || value != as.integer(value)) {
      stop("`", label, "` must be one positive integer.")
    }
    as.integer(value)
  }
  block_size <- validate_positive_integer(block_size, "block_size")
  overwrite <- isTRUE(overwrite)
  verbose <- isTRUE(verbose)

  graph_n <- igraph::vcount(graph)
  if (is.null(cell_ids)) cell_ids <- igraph::vertex_attr(graph, "name")
  if (is.null(cell_ids)) cell_ids <- as.character(seq_len(graph_n))
  cell_ids <- as.character(cell_ids)
  if (length(cell_ids) != graph_n || anyNA(cell_ids) || any(cell_ids == "") ||
      anyDuplicated(cell_ids)) {
    stop("`cell_ids` must provide one unique, non-missing identifier per graph vertex.")
  }
  if (is.null(cells)) {
    cell_index <- seq_len(graph_n)
  } else if (is.character(cells)) {
    if (anyNA(cells) || anyDuplicated(cells) || !all(cells %in% cell_ids)) {
      stop("Character `cells` must be unique identifiers present in `cell_ids`.")
    }
    cell_index <- match(cells, cell_ids)
  } else {
    cell_index <- as.integer(cells)
    if (anyNA(cell_index) || anyDuplicated(cell_index) ||
        any(cell_index < 1L | cell_index > graph_n)) {
      stop("Numeric `cells` must be unique graph vertex indices.")
    }
  }
  endpoint_ids <- cell_ids[cell_index]
  endpoint_n <- length(cell_index)
  if (endpoint_n < 2L) stop("At least two endpoint cells are required.")

  weights <- igraph::edge_attr(graph, "weight")
  if (is.null(weights)) weights <- rep(1, igraph::ecount(graph))
  weights <- as.numeric(weights)
  if (length(weights) != igraph::ecount(graph) || any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("Graph edge costs must be finite and non-negative.")
  }
  components <- igraph::components(graph)$membership[cell_index]
  if (length(unique(components)) != 1L) {
    stop("All stored endpoint cells must lie in one connected graph component.")
  }

  path <- normalizePath(path, mustWork = FALSE)
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create the distance-store directory: ", parent)
  }
  metadata_path <- .clonotrace_distance_store_metadata_path(path)
  specification <- list(
    graph = .clonotrace_graph_signature(graph, weights),
    graph_vertex_count = graph_n,
    cell_index = as.integer(cell_index),
    cell_ids = endpoint_ids,
    block_size = block_size,
    storage_backend = "row-major-binary",
    storage_mode = "double",
    byte_order = "little",
    value_bytes = 8L
  )
  block_starts <- seq.int(1L, endpoint_n, by = block_size)

  if (overwrite) {
    if (file.exists(path) && !file.remove(path)) stop("Unable to replace existing store: ", path)
    if (file.exists(metadata_path) && !file.remove(metadata_path)) {
      stop("Unable to replace existing store metadata: ", metadata_path)
    }
  }
  if (file.exists(path) || file.exists(metadata_path)) {
    if (!file.exists(path) || !file.exists(metadata_path)) {
      stop("Existing distance-store files are incomplete; use `overwrite = TRUE`.")
    }
    metadata <- tryCatch(readRDS(metadata_path), error = function(error) {
      stop("Unable to read existing distance-store metadata: ", conditionMessage(error))
    })
    if (!is.list(metadata) || !identical(metadata$specification, specification)) {
      stop("The existing distance store does not match the current graph, cells, or parameters.")
    }
    if (isTRUE(metadata$complete)) return(.clonotrace_read_distance_store(path))
  } else {
    connection <- file(path, open = "wb")
    close(connection)
    metadata <- list(
      format_version = 2L,
      specification = specification,
      completed_blocks = rep(FALSE, length(block_starts)),
      complete = FALSE,
      timing = list(
        shortest_path_seconds = 0,
        disk_write_seconds = 0,
        completed_blocks = 0L
      )
    )
    .clonotrace_atomic_save_rds(metadata, metadata_path)
  }

  for (block_id in seq_along(block_starts)) {
    if (isTRUE(metadata$completed_blocks[[block_id]])) next
    first <- block_starts[[block_id]]
    rows <- seq.int(first, min(endpoint_n, first + block_size - 1L))
    if (verbose) message("Calculating distance-store block ", block_id, "/", length(block_starts))
    shortest_start <- proc.time()[["elapsed"]]
    values <- igraph::distances(
      graph, v = cell_index[rows], to = cell_index,
      weights = weights, algorithm = "automatic"
    )
    shortest_elapsed <- proc.time()[["elapsed"]] - shortest_start
    if (any(!is.finite(values))) {
      stop("A shortest-path block contains non-finite distances.")
    }
    write_start <- proc.time()[["elapsed"]]
    connection <- file(path, open = "r+b")
    write_result <- tryCatch({
      offset <- as.double(rows[[1L]] - 1L) * as.double(endpoint_n) * 8
      seek(connection, where = offset, origin = "start", rw = "write")
      writeBin(
        as.double(t(values)), connection, size = 8L,
        endian = "little", useBytes = TRUE
      )
      flush(connection)
      TRUE
    }, error = function(error) error)
    close(connection)
    if (inherits(write_result, "error")) {
      stop("Unable to write the binary distance store: ", conditionMessage(write_result))
    }
    write_elapsed <- proc.time()[["elapsed"]] - write_start
    metadata$completed_blocks[[block_id]] <- TRUE
    metadata$timing$shortest_path_seconds <-
      metadata$timing$shortest_path_seconds + unname(shortest_elapsed)
    metadata$timing$disk_write_seconds <-
      metadata$timing$disk_write_seconds + unname(write_elapsed)
    metadata$timing$completed_blocks <- metadata$timing$completed_blocks + 1L
    .clonotrace_atomic_save_rds(metadata, metadata_path)
    rm(values)
  }
  expected_bytes <- as.double(endpoint_n) * as.double(endpoint_n) * 8
  observed_bytes <- file.info(path)$size
  if (observed_bytes != expected_bytes) {
    stop(
      "The completed binary distance store has an unexpected file size (expected ",
      expected_bytes, " bytes; observed ", observed_bytes, ")."
    )
  }
  metadata$complete <- TRUE
  .clonotrace_atomic_save_rds(metadata, metadata_path)
  .clonotrace_read_distance_store(path)
}

#' @title Exact Clone OT from a Disk-Backed Graph-Distance Store
#' @description Read bounded cell-distance blocks from a binary store and apply
#'   the existing exact `clone_2_ot()` calculation to every clone pair.
#' @details
#' For each source clone, all distance-store columns for its positive-mass cells
#' are read once. This row-oriented access uses direct file seeks and bounded
#' buffers. Target clones are then grouped so the union
#' of target cells in one in-memory cost block does not exceed
#' `target_cell_limit` where possible. A single clone larger than the limit is
#' processed alone. This changes storage and data access only: the full
#' source-by-target cost matrix and original probability masses are supplied to
#' `clone_2_ot()` for each pair.
#'
#' @param store A completed `clonotrace_distance_store` object or its binary path.
#' @param cell_clone_prob Numeric or sparse cell-by-clone membership matrix with
#'   unique cell row names and clone column names.
#' @param prob_thresh Values below this threshold are set to zero, matching
#'   `graph_clone_ot()` (default: 0.05).
#' @param target_cell_limit Approximate maximum number of unique target cells in
#'   one in-memory cost block (default: 5000). Each disk source-row read contains
#'   all stored cells and uses approximately `8 * source_support * stored_cells`
#'   bytes.
#' @param cores Positive integer number of source-clone partitions to process.
#'   With `cores > 1`, configure a parallel \pkg{future} plan before calling;
#'   each worker opens its own read-only distance-store connection. Per-source
#'   checkpoints remain compatible between serial and parallel runs.
#' @param checkpoint_dir Optional directory for atomic, validated
#'   per-source-clone checkpoints.
#' @param overwrite_checkpoint Logical. Recompute existing checkpoints.
#' @param monitor Logical. Attach accumulated disk-read and OT-solver elapsed
#'   timings and maximum loaded cost-block size to the returned data frame.
#' @param verbose Logical. Report source-clone and target-block progress.
#'
#' @return A data frame with `group1`, `group2`, and exact OT `dis` columns.
#' @export
graph_clone_ot_disk <- function(store, cell_clone_prob, prob_thresh = 0.05,
                                target_cell_limit = 5000L,
                                cores = 1L,
                                checkpoint_dir = NULL,
                                overwrite_checkpoint = FALSE,
                                monitor = FALSE,
                                verbose = TRUE) {
  store <- .clonotrace_read_distance_store(store)
  cell_ids <- store$metadata$specification$cell_ids
  if (is.null(rownames(cell_clone_prob)) || is.null(colnames(cell_clone_prob))) {
    stop("`cell_clone_prob` must have cell row names and clone column names.")
  }
  if (anyNA(rownames(cell_clone_prob)) || anyDuplicated(rownames(cell_clone_prob)) ||
      anyNA(colnames(cell_clone_prob)) || anyDuplicated(colnames(cell_clone_prob)) ||
      any(rownames(cell_clone_prob) == "") || any(colnames(cell_clone_prob) == "")) {
    stop("Cell and clone identifiers must be unique, non-missing, and non-empty.")
  }
  values <- if (inherits(cell_clone_prob, "Matrix")) cell_clone_prob@x else as.numeric(cell_clone_prob)
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0)) {
    stop("`cell_clone_prob` must contain finite, non-negative numeric values.")
  }
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L ||
      is.na(prob_thresh) || prob_thresh < 0) {
    stop("`prob_thresh` must be one non-negative numeric value.")
  }
  if (!is.numeric(target_cell_limit) || length(target_cell_limit) != 1L ||
      is.na(target_cell_limit) || target_cell_limit < 1L ||
      target_cell_limit != as.integer(target_cell_limit)) {
    stop("`target_cell_limit` must be one positive integer.")
  }
  target_cell_limit <- as.integer(target_cell_limit)
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores) ||
      cores < 1L || cores != as.integer(cores)) {
    stop("`cores` must be one positive integer.")
  }
  cores <- as.integer(cores)
  overwrite_checkpoint <- isTRUE(overwrite_checkpoint)
  monitor <- isTRUE(monitor)
  verbose <- isTRUE(verbose)

  if (inherits(cell_clone_prob, "Matrix")) {
    cell_clone_prob@x[cell_clone_prob@x < prob_thresh] <- 0
    cell_clone_prob <- Matrix::drop0(cell_clone_prob)
  } else {
    cell_clone_prob[cell_clone_prob < prob_thresh] <- 0
  }
  clone_sums <- Matrix::colSums(cell_clone_prob)
  if (ncol(cell_clone_prob) < 2L || any(!is.finite(clone_sums)) || any(clone_sums <= 0)) {
    stop("At least two clones must retain positive total mass after thresholding.")
  }

  supports <- lapply(seq_len(ncol(cell_clone_prob)), function(i) {
    column <- cell_clone_prob[, i]
    membership_index <- which(column > 0)
    identifiers <- rownames(cell_clone_prob)[membership_index]
    store_index <- match(identifiers, cell_ids)
    if (anyNA(store_index)) {
      missing <- identifiers[is.na(store_index)]
      stop(
        "The distance store omits positive-mass cell(s), including: ",
        paste(utils::head(missing, 3L), collapse = ", "), "."
      )
    }
    mass <- as.numeric(column[membership_index])
    order_index <- order(store_index)
    list(
      store_index = as.integer(store_index[order_index]),
      mass = mass[order_index] / sum(mass)
    )
  })
  clone_names <- colnames(cell_clone_prob)

  if (!is.null(checkpoint_dir)) {
    if (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L ||
        is.na(checkpoint_dir) || checkpoint_dir == "") {
      stop("`checkpoint_dir` must be NULL or one non-empty directory path.")
    }
    if (!dir.exists(checkpoint_dir) &&
        !dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Unable to create checkpoint directory: ", checkpoint_dir)
    }
  }
  base_metadata <- list(
    format_version = 1L,
    method = "exact-disk-backed-ot",
    store_specification = store$metadata$specification,
    membership = .clonotrace_input_signature(cell_clone_prob),
    prob_thresh = prob_thresh,
    target_cell_limit = target_cell_limit
  )

  make_target_batches <- function(target_ids) {
    batches <- list()
    current <- integer()
    current_cells <- integer()
    for (target in target_ids) {
      candidate_cells <- union(current_cells, supports[[target]]$store_index)
      if (length(current) > 0L && length(candidate_cells) > target_cell_limit) {
        batches[[length(batches) + 1L]] <- current
        current <- target
        current_cells <- supports[[target]]$store_index
      } else {
        current <- c(current, target)
        current_cells <- candidate_cells
      }
    }
    if (length(current) > 0L) batches[[length(batches) + 1L]] <- current
    batches
  }

  # Bind package-internal helpers into the closure submitted to future workers.
  # This is necessary when the package is loaded from an installed namespace:
  # automatic globals discovery does not reliably resolve dotted internals
  # referenced only inside a nested worker function.
  read_binary_rows <- .clonotrace_read_binary_rows
  read_validated_cache <- .clonotrace_read_cache
  atomic_save_rds <- .clonotrace_atomic_save_rds
  solve_clone_ot <- clone_2_ot

  source_ids <- seq_len(ncol(cell_clone_prob) - 1L)
  cores <- min(cores, length(source_ids))
  support_sizes <- vapply(supports, function(x) length(x$store_index), integer(1))
  suffix_support <- rev(cumsum(rev(as.double(support_sizes))))
  estimated_work <- as.double(support_sizes[source_ids]) * suffix_support[source_ids + 1L]
  partition_id <- integer(length(source_ids))
  partition_load <- numeric(cores)
  for (i in order(estimated_work, decreasing = TRUE)) {
    worker <- which.min(partition_load)
    partition_id[[i]] <- worker
    partition_load[[worker]] <- partition_load[[worker]] + estimated_work[[i]]
  }
  source_partitions <- lapply(seq_len(cores), function(worker) {
    sort(source_ids[partition_id == worker])
  })

  solve_partition <- function(partition) {
    connection <- file(store$path, open = "rb")
    on.exit(close(connection), add = TRUE)
    local_diagnostics <- list(
      disk_read_seconds = 0,
      ot_solve_seconds = 0,
      disk_reads = 0L,
      max_cost_block_bytes = 0
    )
    values <- lapply(partition, function(i) {
      metadata <- c(base_metadata, list(source_clone = clone_names[[i]]))
      checkpoint_path <- if (is.null(checkpoint_dir)) NULL else {
        file.path(checkpoint_dir, sprintf("source-clone-%06d.rds", i))
      }
      if (!is.null(checkpoint_path) && file.exists(checkpoint_path) &&
          !overwrite_checkpoint) {
        if (verbose) message("Reusing disk-backed OT checkpoint ", i, "/", length(source_ids))
        return(read_validated_cache(checkpoint_path, metadata, "disk-backed exact OT"))
      }
      if (verbose) {
        message("Solving disk-backed exact OT for source clone ", i, "/",
                length(source_ids), " on PID ", Sys.getpid())
      }
      target_ids <- seq.int(i + 1L, ncol(cell_clone_prob))
      batches <- make_target_batches(target_ids)
      if (monitor) {
        read_start <- proc.time()[["elapsed"]]
        source_rows <- read_binary_rows(
          connection, supports[[i]]$store_index, length(cell_ids)
        )
        read_elapsed <- proc.time()[["elapsed"]] - read_start
        local_diagnostics$disk_read_seconds <<-
          local_diagnostics$disk_read_seconds + unname(read_elapsed)
        local_diagnostics$disk_reads <<- local_diagnostics$disk_reads + 1L
      } else {
        source_rows <- read_binary_rows(
          connection, supports[[i]]$store_index, length(cell_ids)
        )
      }
      source_rows <- matrix(
        as.numeric(source_rows), nrow = length(supports[[i]]$store_index),
        ncol = length(cell_ids)
      )
      if (monitor) {
        local_diagnostics$max_cost_block_bytes <<- max(
          local_diagnostics$max_cost_block_bytes,
          as.numeric(utils::object.size(source_rows))
        )
      }
      batch_results <- lapply(seq_along(batches), function(batch_id) {
        targets <- batches[[batch_id]]
        target_union <- sort(unique(unlist(lapply(
          targets, function(j) supports[[j]]$store_index
        ))))
        if (verbose && length(batches) > 1L) {
          message("  Reading target block ", batch_id, "/", length(batches),
                  " (", length(target_union), " unique cells)")
        }
        cost_block <- source_rows[, target_union, drop = FALSE]
        if (monitor) {
          local_diagnostics$max_cost_block_bytes <<- max(
            local_diagnostics$max_cost_block_bytes,
            as.numeric(utils::object.size(cost_block))
          )
        }
        distances <- vapply(targets, function(j) {
          columns <- match(supports[[j]]$store_index, target_union)
          if (monitor) {
            solve_start <- proc.time()[["elapsed"]]
            value <- solve_clone_ot(
              cost_block[, columns, drop = FALSE],
              supports[[i]]$mass, supports[[j]]$mass
            )
            solve_elapsed <- proc.time()[["elapsed"]] - solve_start
            local_diagnostics$ot_solve_seconds <<-
              local_diagnostics$ot_solve_seconds + unname(solve_elapsed)
            value
          } else {
            solve_clone_ot(
              cost_block[, columns, drop = FALSE],
              supports[[i]]$mass, supports[[j]]$mass
            )
          }
        }, numeric(1))
        data.frame(
          group1 = rep(clone_names[[i]], length(targets)),
          group2 = clone_names[targets],
          dis = distances,
          stringsAsFactors = FALSE
        )
      })
      value <- do.call(rbind, batch_results)
      rownames(value) <- NULL
      if (!is.null(checkpoint_path)) {
        atomic_save_rds(list(metadata = metadata, value = value), checkpoint_path)
      }
      value
    })
    list(values = values, diagnostics = local_diagnostics)
  }

  if (verbose && cores > 1L) {
    message("Processing ", length(source_ids), " source clones in ", cores,
            " balanced future partitions")
  }
  partition_results <- if (cores > 1L) {
    future.apply::future_lapply(
      source_partitions, solve_partition,
      future.seed = TRUE, future.scheduling = 1
    )
  } else {
    lapply(source_partitions, solve_partition)
  }
  source_results <- unlist(lapply(partition_results, `[[`, "values"), recursive = FALSE)
  result <- do.call(rbind, source_results)
  rownames(result) <- NULL
  result <- .clonotrace_normalize_distance_result(result, clone_names)
  result <- result[order(match(result$group1, clone_names),
                         match(result$group2, clone_names)), , drop = FALSE]
  rownames(result) <- NULL
  if (monitor) {
    worker_diagnostics <- lapply(partition_results, `[[`, "diagnostics")
    attr(result, "diagnostics") <- list(
      disk_read_seconds = sum(vapply(
        worker_diagnostics, `[[`, numeric(1), "disk_read_seconds"
      )),
      ot_solve_seconds = sum(vapply(
        worker_diagnostics, `[[`, numeric(1), "ot_solve_seconds"
      )),
      disk_reads = sum(vapply(
        worker_diagnostics, `[[`, integer(1), "disk_reads"
      )),
      max_cost_block_bytes = max(vapply(
        worker_diagnostics, `[[`, numeric(1), "max_cost_block_bytes"
      )),
      source_partitions = cores
    )
  }
  result
}

.clonotrace_graph_signature <- function(graph, weights) {
  edges <- igraph::as_edgelist(graph, names = FALSE)
  list(
    vertices = igraph::vcount(graph),
    edges = igraph::ecount(graph),
    directed = igraph::is_directed(graph),
    edge_from = .clonotrace_numeric_signature(edges[, 1L]),
    edge_to = .clonotrace_numeric_signature(edges[, 2L]),
    weights = .clonotrace_numeric_signature(weights),
    vertex_names = igraph::vertex_attr(graph, "name")
  )
}
