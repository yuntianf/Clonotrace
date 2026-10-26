.clonotrace_numeric_signature <- function(x) {
  values <- if (inherits(x, "Matrix")) x@x else as.numeric(x)
  finite <- values[is.finite(values)]
  c(
    length = length(values),
    finite = length(finite),
    sum = if (length(finite)) sum(finite) else 0,
    abs_sum = if (length(finite)) sum(abs(finite)) else 0,
    sum_squares = if (length(finite)) sum(finite^2) else 0,
    min = if (length(finite)) min(finite) else NA_real_,
    max = if (length(finite)) max(finite) else NA_real_
  )
}

.clonotrace_input_signature <- function(x) {
  list(
    class = class(x),
    dim = dim(x),
    rownames = rownames(x),
    colnames = colnames(x),
    numeric = .clonotrace_numeric_signature(x)
  )
}

.clonotrace_atomic_save_rds <- function(object, path) {
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (file.exists(path) && !file.remove(path)) {
    stop("Unable to replace existing checkpoint: ", path)
  }
  if (!file.rename(temporary, path)) {
    stop("Unable to atomically write checkpoint: ", path)
  }
  invisible(path)
}

.clonotrace_read_cache <- function(path, expected, label) {
  cached <- tryCatch(readRDS(path), error = function(error) {
    stop("Unable to read ", label, " checkpoint `", path, "`: ", conditionMessage(error))
  })
  if (!is.list(cached) || !all(c("metadata", "value") %in% names(cached))) {
    warning(
      "Loading legacy ", label, " checkpoint without validation metadata. ",
      "Use `overwrite = TRUE` once to replace it with a validated checkpoint.",
      call. = FALSE
    )
    return(cached)
  }
  if (!identical(cached$metadata, expected)) {
    stop(
      "The existing ", label, " checkpoint does not match the current inputs or parameters. ",
      "Use a different `outpath` or set `overwrite = TRUE` to recompute it."
    )
  }
  cached$value
}

.clonotrace_snn_similarity_to_cost <- function(graph) {
  similarity <- as.numeric(igraph::edge_attr(graph, "weight"))
  if (length(similarity) != igraph::ecount(graph) ||
      any(!is.finite(similarity)) || any(similarity < 0)) {
    stop("SNN edge weights must be finite, non-negative similarities.")
  }

  # `set_edge_attr()` returns the modified graph; it does not mutate `graph`.
  graph <- igraph::set_edge_attr(
    graph,
    "weight",
    value = exp(-similarity)
  )
  graph
}

.clonotrace_normalize_distance_result <- function(result, clone_names) {
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  required <- c("group1", "group2", "dis")
  if (!all(required %in% colnames(result))) {
    stop("Clone-distance output must contain columns: group1, group2, and dis.")
  }
  result <- result[, required, drop = FALSE]
  normalize_ids <- function(ids) {
    ids_character <- as.character(ids)
    ids_integer <- suppressWarnings(as.integer(ids_character))
    positional <- !anyNA(ids_integer) && all(ids_integer >= 1L & ids_integer <= length(clone_names))
    if (positional) ids_character <- clone_names[ids_integer]
    if (!all(ids_character %in% clone_names)) {
      stop("Clone-distance output contains identifiers absent from `cell_clone_prob`.")
    }
    ids_character
  }
  result$group1 <- normalize_ids(result$group1)
  result$group2 <- normalize_ids(result$group2)
  result$dis <- as.numeric(result$dis)
  if (any(is.infinite(result$dis))) {
    stop(
      "Clone-distance output contains infinite distances because the cell graph has disconnected clone-bearing components. ",
      "Increase `graph_k`, inspect/filter disconnected components, or provide an embedding that yields a connected graph."
    )
  }
  if (anyNA(result$dis) || any(result$dis < 0)) {
    stop("Clone-distance output contains missing or negative distance values.")
  }
  result
}

#' @title clone_disance
#' @description Compute Clone-to-Clone Distances from Cell Embedding and Clone Assignments
#' @details
#'
#' Constructs a cell-cell shared nearest neighbor (SNN) graph from an embedding and uses it to estimate
#' pairwise distances between clones. Users can choose between exact OT-based distances or fast
#' approximate k-nearest neighbor methods.
#'
#' @param embedding A numeric matrix of cell embeddings (cells × dimensions).
#' @param cell_clone_prob A numeric matrix of clone membership probabilities (cells × clones).
#' @param outpath File path to save and load intermediate graph and distance results.
#' @param graph_k Integer. Number of neighbors to use in the SNN graph (default: 10).
#' @param overwrite Logical. Whether to overwrite previously saved graphs and results (default: FALSE).
#' @param exact Logical. If `TRUE`, computes exact OT-based clone distances with
#'   the disk-backed backend by default; otherwise, uses approximate NN
#'   (default: FALSE).
#' @param exact_backend Exact-OT implementation. `"disk"` (default) writes a
#'   resumable cell-distance store and reads bounded cost blocks. `"memory"`
#'   retains the legacy dense in-memory implementation for small reproductions.
#' @param force_exact Logical. Allow exact mode to proceed when the input exceeds
#'   the conservative cell or clone safety limits (default: `FALSE`).
#' @param exact_cell_limit,exact_clone_limit Conservative safety limits used to
#'   warn before an exact calculation that is likely to require substantial
#'   runtime or memory.
#' @param revised Logical. If `TRUE` (the default), optionally locally rescales
#'   the raw Clonotrace distance and blends it with clone-centroid distance
#'   using a fate-blind automatically selected weight; see
#'   `revise_clone_distance()`. Local rescaling is off by default and can be
#'   enabled with `local_rescale = TRUE`.
#'   Set `FALSE` to return the original exact-OT or approximate distance.
#' @param clone_labels Optional observed clone membership used for clone
#'   centroids, as a data frame with columns `cell` and `clone`, or a vector
#'   named by cell. If `NULL`, `cell_clone_prob` supplies probability weights.
#' @param sample Optional clone-level grouping, such as patient, passed to
#'   `revise_clone_distance()` so normalization and weight selection use
#'   within-group clone pairs.
#' @param ... Named arguments passed to `graph_clone_ot()` or
#'   `graph_clone_nn()`, and—when `revised = TRUE`—to
#'   `revise_clone_distance()` (`local_k`, `local_rescale`, `resolutions`, `clusters`,
#'   `cluster_k`, `mean_weight`, `s_grid`, `seed`, and `verbose`).
#'
#' @return A data frame with columns `group1`, `group2`, and `dis`. With
#'   `revised = TRUE`, `dis` is the revised distance, `dis_clonotrace` retains
#'   the raw distance, and attributes `mean_weight` and `picks` record the
#'   automatic blend. With `revised = FALSE`, `dis` is the raw distance.
#' @importFrom  bluster makeSNNGraph
#' @importFrom igraph edge_attr
#' @importFrom igraph set_edge_attr
#' @importFrom utils packageVersion
#' @export
clone_disance = function(embedding,cell_clone_prob,outpath,graph_k = 10,
                         overwrite = FALSE, exact = FALSE, force_exact = FALSE,
                         exact_cell_limit = 20000L, exact_clone_limit = 500L,
                         exact_backend = c("disk", "memory"),
                         revised = TRUE, clone_labels = NULL, sample = NULL, ...){
  if (!is.character(outpath) || length(outpath) != 1L || is.na(outpath) || outpath == "") {
    stop("`outpath` must be one non-empty directory path.")
  }
  if (!dir.exists(outpath) && !dir.create(outpath, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create `outpath`: ", outpath)
  }
  write_probe <- tempfile(pattern = ".clonotrace-write-test-", tmpdir = outpath)
  if (!file.create(write_probe)) stop("`outpath` is not writable: ", outpath)
  unlink(write_probe)
  if (!is.numeric(graph_k) || length(graph_k) != 1L || is.na(graph_k) || graph_k < 1L) {
    stop("`graph_k` must be one positive integer.")
  }
  graph_k <- as.integer(graph_k)
  exact <- isTRUE(exact)
  exact_backend <- match.arg(exact_backend)
  overwrite <- isTRUE(overwrite)
  revised <- isTRUE(revised)

  args <- list(...)
  distance_args_allowed <- if (exact && identical(exact_backend, "disk")) {
    c(
      "prob_thresh", "cores", "distance_block_size", "target_cell_limit",
      "checkpoint_dir", "monitor", "verbose"
    )
  } else if (exact) {
    c("prob_thresh", "cores", "cache", "verbose")
  } else {
    c("prob_thresh", "k", "clone_block_size", "verbose", "distance_cache_mb")
  }
  revision_args_allowed <- c(
    "local_k", "local_rescale", "resolutions", "clusters", "cluster_k", "mean_weight",
    "s_grid", "seed", "verbose"
  )
  allowed_args <- union(
    distance_args_allowed,
    if (revised) revision_args_allowed else character()
  )
  unknown_args <- setdiff(names(args), allowed_args)
  if (length(args) > 0L && (is.null(names(args)) || any(names(args) == ""))) {
    stop("All additional arguments in `...` must be named.")
  }
  if (length(unknown_args) > 0L) {
    stop(
      "Unknown or inapplicable argument(s) for ", if (exact) "exact" else "approximate",
      if (revised) " revised" else " unrevised", " mode: ",
      paste(unknown_args, collapse = ", "), "."
    )
  }
  distance_args <- args[names(args) %in% distance_args_allowed]
  revision_args <- args[names(args) %in% revision_args_allowed]

  if (is.null(rownames(embedding)) || is.null(rownames(cell_clone_prob)) ||
      is.null(colnames(cell_clone_prob))) {
    stop("`embedding` and `cell_clone_prob` must use cell row names; clone columns must also be named.")
  }
  if (anyDuplicated(rownames(embedding)) || anyDuplicated(rownames(cell_clone_prob)) ||
      anyDuplicated(colnames(cell_clone_prob))) {
    stop("Cell and clone identifiers must be unique.")
  }
  cells = intersect(rownames(embedding),rownames(cell_clone_prob))
  missing_membership <- setdiff(rownames(embedding), rownames(cell_clone_prob))
  if (length(missing_membership) > 0L) {
    stop(
      "`cell_clone_prob` is missing ", length(missing_membership),
      " cell(s) present in `embedding`; supply aligned inputs."
    )
  }
  if(length(setdiff(rownames(cell_clone_prob),cells)) > 0){
    warning("The cell_clone_prob has cells which are not in the cell embedding, those cells would be removed!")
  }

  embedding = as.matrix(embedding)
  if (!is.numeric(embedding) || any(!is.finite(embedding))) {
    stop("`embedding` must be a finite numeric matrix.")
  }
  cell_clone_prob = sync_sparse_rows(cell_clone_prob,rownames(embedding))
  membership_values <- if (inherits(cell_clone_prob, "Matrix")) cell_clone_prob@x else as.numeric(cell_clone_prob)
  if (any(!is.finite(membership_values)) || any(membership_values < 0)) {
    stop("`cell_clone_prob` must contain finite, non-negative values.")
  }
  if (any(Matrix::colSums(cell_clone_prob) <= 0)) {
    stop("Every clone column in `cell_clone_prob` must have positive total mass.")
  }
  if (nrow(embedding) <= graph_k) {
    stop("`graph_k` must be smaller than the number of cells.")
  }
  if (exact && !isTRUE(force_exact) &&
      (nrow(embedding) > exact_cell_limit || ncol(cell_clone_prob) > exact_clone_limit)) {
    stop(
      "Exact clone distance was requested for ", nrow(embedding), " cells and ",
      ncol(cell_clone_prob), " clones, above the conservative safety limits (",
      exact_cell_limit, " cells / ", exact_clone_limit, " clones). ",
      "Use approximate mode (`exact = FALSE`) or set `force_exact = TRUE` after confirming resources."
    )
  }

  graph_metadata <- list(
    format_version = 2L,
    package_version = tryCatch(as.character(utils::packageVersion("Clonotrace")), error = function(error) NA_character_),
    embedding = .clonotrace_input_signature(embedding),
    graph_k = graph_k,
    edge_cost_transform = "exp(-shared_neighbor_count)"
  )
  graph_path <- file.path(outpath, "cell_graph.rds")

  if(file.exists(graph_path) && !overwrite){
    message("Reusing cell-graph checkpoint: ", graph_path)
    cell_graph = .clonotrace_read_cache(graph_path, graph_metadata, "cell graph")
  }
  else{
    message("Computing cell graph (k = ", graph_k, ").")
    cell_graph = bluster::makeSNNGraph(embedding,k = graph_k,type="number")
    cell_graph <- .clonotrace_snn_similarity_to_cost(cell_graph)
    .clonotrace_atomic_save_rds(list(metadata = graph_metadata, value = cell_graph), graph_path)
  }

  result_metadata <- list(
    format_version = 1L,
    package_version = graph_metadata$package_version,
    graph = graph_metadata,
    cell_clone_prob = .clonotrace_input_signature(cell_clone_prob),
    method = if (exact) paste0("exact_", exact_backend) else "approximate",
    arguments = distance_args
  )
  result_filename <- if (exact) "clone_graph_dis_exact.rds" else "clone_graph_dis.rds"
  result_path <- file.path(outpath, result_filename)

  if(file.exists(result_path) && !overwrite){
    message("Reusing clone-distance checkpoint: ", result_path)
    dis_result = .clonotrace_read_cache(result_path, result_metadata, "clone distance")
    dis_result = .clonotrace_normalize_distance_result(dis_result, colnames(cell_clone_prob))
  }
  else{
    message("Computing ", if (exact) "exact" else "approximate", " clone distances.")
    if(exact == TRUE && identical(exact_backend, "disk")){
      prob_thresh <- if (is.null(distance_args$prob_thresh)) 0.05 else distance_args$prob_thresh
      if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L ||
          is.na(prob_thresh) || prob_thresh < 0) {
        stop("`prob_thresh` must be one non-negative numeric value.")
      }
      positive <- if (prob_thresh == 0) cell_clone_prob > 0 else cell_clone_prob >= prob_thresh
      clone_cells <- rownames(cell_clone_prob)[Matrix::rowSums(positive) > 0]
      if (length(clone_cells) < 2L) {
        stop("At least two cells must retain clone mass after thresholding.")
      }

      distance_block_size <- if (is.null(distance_args$distance_block_size)) {
        256L
      } else {
        distance_args$distance_block_size
      }
      target_cell_limit <- if (is.null(distance_args$target_cell_limit)) {
        5000L
      } else {
        distance_args$target_cell_limit
      }
      checkpoint_dir <- if (is.null(distance_args$checkpoint_dir)) {
        file.path(outpath, "exact_ot_checkpoints")
      } else {
        distance_args$checkpoint_dir
      }
      cores <- if (is.null(distance_args$cores)) 1L else distance_args$cores
      monitor <- isTRUE(distance_args$monitor)
      verbose <- if (is.null(distance_args$verbose)) TRUE else isTRUE(distance_args$verbose)

      distance_store <- build_graph_distance_store(
        graph = cell_graph,
        path = file.path(outpath, "cell_graph_distance_exact.bin"),
        cells = clone_cells,
        cell_ids = rownames(embedding),
        block_size = distance_block_size,
        overwrite = overwrite,
        verbose = verbose
      )
      dis_result <- graph_clone_ot_disk(
        store = distance_store,
        cell_clone_prob = cell_clone_prob[clone_cells, , drop = FALSE],
        prob_thresh = prob_thresh,
        target_cell_limit = target_cell_limit,
        cores = cores,
        checkpoint_dir = checkpoint_dir,
        overwrite_checkpoint = overwrite,
        monitor = monitor,
        verbose = verbose
      )
    }
    else if(exact == TRUE){
      ot_args <- c("prob_thresh","cores", "cache","verbose")
      call_args <- distance_args[names(distance_args) %in% ot_args]
      dis_result <- do.call(graph_clone_ot, c(list(graph = cell_graph, cell_clone_prob = cell_clone_prob), call_args))
    }
    else{
      nn_args <- c("prob_thresh", "k", "clone_block_size", "verbose", "distance_cache_mb")
      call_args <- distance_args[names(distance_args) %in% nn_args]
      dis_result <- do.call(graph_clone_nn, c(list(graph = cell_graph, cell_clone_prob = cell_clone_prob), call_args))
    }

    dis_result <- .clonotrace_normalize_distance_result(dis_result, colnames(cell_clone_prob))
    .clonotrace_atomic_save_rds(list(metadata = result_metadata, value = dis_result), result_path)
  }

  if (revised) {
    revised_distance <- do.call(
      revise_clone_distance,
      c(
        list(
          clone_dis = dis_result,
          embedding = embedding,
          clone_labels = clone_labels,
          cell_clone_prob = cell_clone_prob,
          clone_names = colnames(cell_clone_prob),
          sample = sample
        ),
        revision_args
      )
    )
    result <- as.data.frame(dis_result, stringsAsFactors = FALSE)
    first <- match(as.character(result$group1), colnames(cell_clone_prob))
    second <- match(as.character(result$group2), colnames(cell_clone_prob))
    if (anyNA(first) || anyNA(second)) {
      stop("Raw clone-distance identifiers could not be matched during revision.")
    }
    result$dis_clonotrace <- result$dis
    result$dis <- revised_distance[cbind(first, second)]
    attr(result, "mean_weight") <- attr(revised_distance, "mean_weight")
    attr(result, "picks") <- attr(revised_distance, "picks")
    attr(result, "selection_correlations") <-
      attr(revised_distance, "selection_correlations")
    dis_result <- result
  }

  return(dis_result)
}

#' @title clone_partition
#' @description Partition Clones Based on Containment Similarity With Centers From Farthest Point Sampling
#' @details
#'
#' Groups clones into `k` roughly balanced partitions based on containment similarity.
#' Useful for parallelizing clone-level computations (e.g., OT) by minimizing inter-group redundancy.
#'
#' @param clone_matrix A cells × clones matrix, where non-zero entries indicate clone membership.
#' @param k Integer. Number of clone partitions (default: 10).
#' @param similarity_threshold Numeric threshold for assigning a clone to a group based on similarity (default: 0).
#'
#' @return A named list of clone ID vectors, one per partition.
#'
#' @importFrom Matrix colSums t
clone_partition <- function(clone_matrix, k = 10, similarity_threshold = 0) {
  bin <- clone_matrix > 0
  n_clones <- ncol(bin)
  sizes <- Matrix::colSums(bin)
  clone_names <- colnames(bin)

  # Compute containment similarity: sim(i, j) = |i ∩ j| / |j|
  intersection <- as.matrix(Matrix::t(bin) %*% bin)
  denom <- matrix(rep(sizes, each = n_clones), nrow = n_clones)
  sim <- intersection / denom
  diag(sim) <- 1

  # Step 1: Farthest Point Sampling (FPS) for seed selection
  seeds <- c(which.max(sizes))  # start with largest clone
  while (length(seeds) < k) {
    remaining <- setdiff(1:n_clones, seeds)
    min_sim_to_seeds <- apply(sim[remaining, seeds, drop = FALSE], 1, max)
    next_seed <- remaining[which.min(min_sim_to_seeds)]
    seeds <- c(seeds, next_seed)
  }

  # Step 2: Greedy assignment with balancing
  group_ids <- rep(NA_integer_, n_clones)
  group_sizes <- rep(0, k)
  target_size <- ceiling(n_clones / k)

  # Order clones by max similarity to any seed
  clone_sim_to_seeds <- sim[seeds, , drop = FALSE]
  assign_order <- order(apply(clone_sim_to_seeds, 2, max), decreasing = TRUE)

  for (i in assign_order) {
    clone_sims <- sim[seeds, i]
    clone_sims[is.na(clone_sims)] <- 0
    ranked <- order(clone_sims, decreasing = TRUE)

    for (g in ranked) {
      if (clone_sims[g] >= similarity_threshold && group_sizes[g] < target_size) {
        group_ids[i] <- g
        group_sizes[g] <- group_sizes[g] + 1
        break
      }
    }

    # Fallback if no high-similar group is under capacity
    if (is.na(group_ids[i])) {
      g <- which.min(group_sizes)
      group_ids[i] <- g
      group_sizes[g] <- group_sizes[g] + 1
    }
  }

  names(group_ids) <- clone_names
  return(split(names(group_ids), group_ids))
}

#' @title graph_clone_ot_sub
#' @description Compute Pairwise OT Distances Between Clones Using Graph-Based Subset Strategy
#' @details
#'
#' Computes clone-to-clone optimal transport distances by expanding a local graph-based cell neighborhood
#' using a subset pool strategy to reduce computational burden.
#'
#' @param graph An `igraph` object representing cell-cell distances.
#' @param cell_clone_prob A matrix of clone membership probabilities (cells × clones).
#' @param target_clone Optional. Vector of clone indices to compute pairwise OT from (default: `NULL` = all).
#' @param cache Integer. Maximum number of pooled cells retained in memory during iteration (default: 5000).
#' @param verbose Bool flag to indicate if print the intermediate log for each pair of clones.
#'
#' @return A matrix with columns: `group1`, `group2`, `dis` for each clone pair and their OT distance.
#'
#' @importFrom igraph distances
#'
#' @export
graph_clone_ot_sub = function(graph,cell_clone_prob,target_clone = NULL,cache = 5000,verbose = TRUE){
  if(is.null(target_clone)){
    target_clone = 1:ncol(cell_clone_prob)
  }

  target_clone_ident = cell_clone_prob[,target_clone,drop = FALSE] > 0
  flag = rep(0,ncol(target_clone_ident))

  pool = c()
  cell_dis = c()

  full_ot = c()

  target_id = which((colSums(target_clone_ident) == max(colSums(target_clone_ident))))[1]
  while(sum(flag) < length(target_clone)){
    flag[target_id] = 1
    global_id = target_clone[target_id]
    if(verbose){
      cat("Start to process the clone ",global_id," on PID ", Sys.getpid(), ".\n")
    }

    if(global_id == ncol(cell_clone_prob)){
      next
    }

    cell_id = which(cell_clone_prob[,global_id] > 0)

    if(length(cell_id) > cache){
      message = paste(c("The size of the clone is too large(",length(cell_id),"),
                        please consider to remove this clone or increase the cache to be over the size of the clone"),sep = "")
      stop(message)
    }

    append_cells = setdiff(cell_id,pool)
    append_dis = distances(graph,v = append_cells)

    if(length(append_cells)+length(pool) > cache){
      remove_n = length(append_cells)+length(pool)-cache

      cell_freq = rowSums(target_clone_ident[setdiff(pool,cell_id),!flag,drop = FALSE])

      remove_id = which(rank(cell_freq) <= remove_n)
      pool = pool[-remove_id]
      cell_dis = cell_dis[-remove_id,]
    }

    pool = c(pool,append_cells)
    cell_dis = rbind(cell_dis,append_dis)

    clone_mass = cell_clone_prob[cell_id,global_id]
    row_id = match(cell_id,pool)

    clone_ot = lapply((global_id+1):ncol(cell_clone_prob),function(i){
      clone2_cells = cell_clone_prob[,i] > 0
      sub_dis = cell_dis[row_id,clone2_cells]
      ot = clone_2_ot(sub_dis,clone_mass,cell_clone_prob[clone2_cells,i])
      return(c(global_id,i,ot))
    })
    clone_ot = do.call(rbind,clone_ot)
    full_ot = rbind(full_ot,clone_ot)

    pool_ident = matrix(rep(0,nrow(cell_clone_prob)),nrow = 1)
    pool_ident[,pool] = 1
    pool_subset_dis = colSums(target_clone_ident[,!flag,drop = FALSE]) - pool_ident %*% target_clone_ident[,!flag,drop = FALSE]
    pool_subset_dis = as.numeric(pool_subset_dis)

    target_id = which(!flag)[which.min(pool_subset_dis)]
  }

  return(full_ot)
}

#' @title graph_clone_ot
#' @description Parallel Clone-to-Clone OT Computation on Cell Graph
#' @details
#'
#' Computes clone-to-clone optimal transport (OT) distances in parallel by first partitioning clones
#' and applying `graph_clone_ot_sub()` on each group. Supports memory-efficient computation with subset caching.
#'
#' @param graph An `igraph` object representing the cell-cell distance graph.
#' @param cell_clone_prob A matrix of clone membership probabilities (cells × clones).
#' @param prob_thresh Float. The minimum weight of a cell in a clone to be preserved, otherwise it would be sparsified.
#' @param cache Integer. Maximum number of pooled cells retained in memory during each subgraph operation (default: 5000).
#' @param cores Integer. Number of parallel threads to use (default: 1).
#' @param verbose Bool flag to indicate if print the intermediate log for each pair of clones.
#'
#' @return A data frame with columns `group1`, `group2`, and `dis` representing OT distances between clone pairs.
#'
#' @importFrom future.apply future_lapply
#'
#' @export
graph_clone_ot = function(graph,cell_clone_prob,prob_thresh = 0.05,cache = 5000,cores = 1,verbose = TRUE){
  cell_clone_prob[cell_clone_prob < prob_thresh] = 0

  partition = clone_partition(cell_clone_prob,k = cores)
  partition = lapply(partition,function(x) return(match(x,colnames(cell_clone_prob))))
  if(verbose){
    cat("There are ",length(partition)," clone partitions\n")
  }

  dis = future_lapply(1:length(partition),function(i){
    cat("Running clone partition", i, "on PID", Sys.getpid(), "\n")
    x = partition[[i]]
    if(verbose){
      cat("This partition includes clone",x,"\n")
    }
    result = graph_clone_ot_sub(graph,cell_clone_prob,x,cache,verbose = verbose)
  },future.seed=TRUE)

  dis = do.call(rbind,dis)
  colnames(dis) = c("group1","group2","dis")
  return(dis)
}

#' @title clone_2_ot
#' @description Compute Optimal Transport Distance Between Two Clones
#' @details
#'
#' Calculates the optimal transport (OT) distance between two sets of cells (e.g., clones) based on
#' a cost matrix and their mass distributions. Uses the `transport` package for solving the OT plan.
#'
#' @param distance A numeric cost matrix of pairwise distances (rows = group 1 cells, cols = group 2 cells).
#' @param group1_mass A numeric vector of mass/probabilities for group 1 (length must equal number of rows in `distance`).
#' @param group2_mass A numeric vector of mass/probabilities for group 2 (length must equal number of columns in `distance`).
#'
#' @return A single numeric value representing the OT distance between the two groups.
#'
#' @importFrom transport transport
#'
#' @export
clone_2_ot = function(distance,group1_mass,group2_mass){
  if(nrow(distance) != length(group1_mass) | ncol(distance) != length(group2_mass)){
    stop("The size of mass doesn't match with the dimension of the distance!")
  }

  group1_mass = group1_mass/sum(group1_mass)
  group2_mass = group2_mass/sum(group2_mass)

  plan = transport::transport(group1_mass, group2_mass, costm = distance)
  # plan$cost = diag(distance[plan$from,plan$to])
  plan$cost <- distance[cbind(plan$from, plan$to)]

  return(sum(plan$cost*plan$mass))
}

#' @title graph_clone_nn
#' @description Compute Clone-to-Clone Nearest Neighbor Distances on a Cell Graph
#' @details
#'
#' Calculates clone-level nearest neighbor distances based on a shared cell-cell graph. Uses binarized
#' clone membership to find inter-clone distances, and applies a greedy k-nearest neighbor approximation.
#'
#' @param graph An `igraph` object representing the cell-cell distance graph.
#' @param cell_clone_prob A numeric matrix of clone membership probabilities (cells × clones).
#' @param prob_thresh Threshold for binarizing clone assignment probabilities (default: 0.1).
#' @param k Integer. Number of nearest neighbors to select for each clone comparison (default: 2).
#' @param clone_block_size Integer. Number of target clones whose cells are
#'   included in each fallback graph-distance block when the source-clone
#'   distance cache would exceed `distance_cache_mb`. Smaller values reduce
#'   peak memory at the cost of additional graph-distance calls (default: 25).
#' @param verbose Logical. Whether to print progress during processing (default: FALSE).
#' @param distance_cache_mb Non-negative finite number. Maximum size in MiB
#'   of the numeric source-clone distance cache per worker (default: 64).
#'   Distances to all remaining target cells are computed once and reused
#'   when they fit; otherwise target-clone blocks are used. Set to 0 to
#'   disable this cache. This is not a total RAM limit: the graph, temporary
#'   matrices, fallback blocks, and other workers require additional memory.
#'
#' @return A data frame with columns `group1`, `group2`, and `dis` indicating the NN distance between clones.
#'
#' @importFrom future.apply future_lapply
#' @importFrom igraph distances
#'
#' @export
graph_clone_nn = function(graph, cell_clone_prob, prob_thresh = 0.1, k = 2,
                          clone_block_size = 25L, verbose = FALSE,
                          distance_cache_mb = 64){
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1L) {
    stop("`k` must be one positive integer.")
  }
  if (!is.numeric(clone_block_size) || length(clone_block_size) != 1L ||
      is.na(clone_block_size) || clone_block_size < 1L) {
    stop("`clone_block_size` must be one positive integer.")
  }
  if (!is.numeric(distance_cache_mb) || length(distance_cache_mb) != 1L ||
      !is.finite(distance_cache_mb) || distance_cache_mb < 0) {
    stop("`distance_cache_mb` must be one finite non-negative number.")
  }
  k <- as.integer(k)
  clone_block_size <- as.integer(clone_block_size)
  cell_group_mat = cell_clone_prob
  cell_group_mat[cell_group_mat < prob_thresh] = 0
  cell_group_mat[cell_group_mat >= prob_thresh] = 1

  group_n = ncol(cell_group_mat)
  retained_sizes <- Matrix::colSums(cell_group_mat > 0)
  if (any(retained_sizes == 0L)) {
    stop(
      "The following clones have no cells after applying `prob_thresh`: ",
      paste(colnames(cell_group_mat)[retained_sizes == 0L], collapse = ", "), "."
    )
  }

  if (group_n < 2L) stop("At least two clones are required to calculate clone distances.")
  # Resolve supports once rather than scanning the cell-by-clone matrix for
  # every clone pair. Keep only these index vectors in the worker closure.
  supports <- lapply(seq_len(group_n), function(i) which(cell_group_mat[, i] > 0))
  rm(cell_group_mat)
  cache_bytes <- distance_cache_mb * 1024^2
  pair_distance <- group_2_min
  dis = future_lapply(1:(group_n-1),function(i){
    if(verbose){
      cat("Starting to proceed clone ",i," on PID ", Sys.getpid(),"\n")
    }

    from <- supports[[i]]
    target_clones <- (i + 1L):group_n
    all_targets <- sort(unique(unlist(supports[target_clones], use.names = FALSE)))
    estimated_bytes <- 8 * as.double(length(from)) * length(all_targets)
    if (cache_bytes > 0 && estimated_bytes <= cache_bytes) {
      if (verbose) {
        message("Reusing source-clone distances across ",
                length(target_clones), " target clones")
      }
      source_distances <- igraph::distances(graph, v = from, to = all_targets)
      dis_out <- vapply(target_clones, function(j) {
        x <- source_distances[, match(supports[[j]], all_targets), drop = FALSE]
        pair_distance(x, seq_len(nrow(x)), seq_len(ncol(x)), k = k)
      }, numeric(1))
      return(cbind(i, target_clones, dis_out))
    }
    target_blocks <- split(target_clones, ceiling(seq_along(target_clones) / clone_block_size))
    block_output <- lapply(target_blocks, function(block) {
      to <- sort(unique(unlist(supports[block], use.names = FALSE)))
      block_distances <- igraph::distances(graph, v = from, to = to)
      dis_out <- vapply(block, function(j) {
        x <- block_distances[, match(supports[[j]], to), drop = FALSE]
        pair_distance(x, seq_len(nrow(x)), seq_len(ncol(x)), k = k)
      }, numeric(1))
      cbind(i, block, dis_out)
    })
    do.call(rbind, block_output)
  },future.seed=TRUE)
  dis = do.call(rbind,dis)
  colnames(dis) = c("group1","group2","dis")

  return(dis)
}

#' @title group_2_min
#' @description Compute Average Minimal Pairwise Distance Between Two Groups
#' @details
#'
#' Estimates a symmetric inter-group distance based on averaging the top-k smallest pairwise distances between two groups.
#' This is useful for comparing cell groups or clone sets using a local-neighborhood approximation of intergroup proximity.
#' If the submatrix between `group1` and `group2` has only one value, that value is returned.
#' Otherwise, the function returns the mean of the average top-k distances across both rows and columns.
#'
#' @param distance A numeric matrix of pairwise distances (e.g., from `dist()` or `igraph::distances()`).
#' @param group1 A vector of row indices corresponding to group 1.
#' @param group2 A vector of column indices corresponding to group 2.
#' @param k Integer. Number of nearest neighbors to average over (default: 3).
#'
#' @return A single numeric value representing the symmetric average distance between the two groups.

group_2_min = function(distance,group1,group2, k = 3){
  sub_dis = distance[group1,group2]
  if(is.null(nrow(sub_dis))){
    out = mean(sub_dis)
  }
  else{
    out = mean(c(mean(apply(sub_dis,1,function(x) mean(top_k(x,k)))),
                 mean(apply(sub_dis,2,function(x) mean(top_k(x,k))))))
  }
  return(out)
}
