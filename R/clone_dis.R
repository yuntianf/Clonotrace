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
#' @param exact Logical. If `TRUE`, computes exact OT-based clone distances; otherwise, uses approximate NN (default: FALSE).
#' @param ... Additional arguments passed to `graph_clone_ot()` or `graph_clone_nn()`.
#'
#' @return A data frame with columns `group1`, `group2`, and `dis`, representing pairwise clone distances.
#' @importFrom  bluster makeSNNGraph
#' @importFrom igraph edge_attr
#' @export
clone_disance = function(embedding,cell_clone_prob,outpath,graph_k = 10,
                         overwrite = FALSE,exact = FALSE,...){
  cells = intersect(rownames(embedding),rownames(cell_clone_prob))
  if(length(setdiff(rownames(cell_clone_prob),cells)) > 0){
    warning("The cell_clone_prob has cells which are not in the cell embedding, those cells would be removed!")
  }

  embedding = as.matrix(embedding)
  cell_clone_prob = sync_sparse_rows(cell_clone_prob,rownames(embedding))

  if(file.exists(file.path(outpath,"cell_graph.rds")) & !overwrite){
    cat("The cell_graph is already calculated and overwrite is not allowed, will use the exsting one!\n")
    cell_graph = readRDS(file.path(outpath,"cell_graph.rds"))
  }
  else{
    cell_graph = makeSNNGraph(embedding,k = graph_k,type="number")
    # cell_graph = delete.edges(cell_graph, E(cell_graph)[edge_attr(cell_graph)$weight < weight])
    edge_attr(cell_graph)$weight = exp(-edge_attr(cell_graph)$weight)
    saveRDS(cell_graph,file.path(outpath,"cell_graph.rds"))
  }

  if(file.exists(file.path(outpath,"clone_graph_dis.rds")) & !overwrite){
    cat("The groups_graph is already calculated and overwrite is not allowed, will use the exsting one!\n")
    dis_result = readRDS(file.path(outpath,"groups_graph_dis.rds"))
  }
  else{
    args <- list(...)
    if(exact == TRUE){
      ot_args <- c("prob_thresh","cores", "cache","verbose")
      call_args <- args[names(args) %in% ot_args]
      dis_result <- do.call(graph_clone_ot, c(list(graph = cell_graph, cell_clone_prob = cell_clone_prob), call_args))
    }
    else{
      nn_args <- c("prob_thresh","k","verbose")
      call_args <- args[names(args) %in% nn_args]
      dis_result <- do.call(graph_clone_nn, c(list(graph = cell_graph, cell_clone_prob = cell_clone_prob), call_args))
    }

    saveRDS(dis_result,file.path(outpath,"clone_graph_dis.rds"))
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
#' @examples
#' mat <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#' colnames(mat) <- paste0("Clone", 1:10)
#' partitions <- clone_partition(mat, k = 3)
#' str(partitions)
clone_partition <- function(clone_matrix, k = 10, similarity_threshold = 0) {
  library(Matrix)

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
#' @param cache Integer. Maximum number of pooled cells retained in memory during each subgraph operation (default: 5000).
#' @param cores Integer. Number of parallel threads to use (default: 1).
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

  dis = future_lapply(partition,function(x){
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
#' @param nn_k Integer. Number of nearest neighbors to select for each clone comparison (default: 2).
#' @param verbose Logical. Whether to print progress during processing (default: FALSE).
#'
#' @return A data frame with columns `group1`, `group2`, and `dis` indicating the NN distance between clones.
#'
#' @importFrom future.apply future_lapply
#' @importFrom igraph distances
#'
#' @export
graph_clone_nn = function(graph,cell_clone_prob,prob_thresh = 0.1,k = 2,verbose = FALSE){
  cell_group_mat = cell_clone_prob
  cell_group_mat[cell_group_mat < prob_thresh] = 0
  cell_group_mat[cell_group_mat >= prob_thresh] = 1

  group_n = ncol(cell_group_mat)

  dis = future_lapply(1:(group_n-1),function(i){
    if(verbose){
      cat("Starting to proceed clone ",i," on PID ", Sys.getpid(),"\n")
    }

    from = which(cell_group_mat[,i] > 0)
    mass_1 = cell_group_mat[from,i]
    # cat("The size of mass 1:",length(mass_1),"\n")
    to = which(rowSums(cell_group_mat[,(i+1):group_n,drop = FALSE]) > 0)

    dis_i = distances(graph, v = from, to = to)

    dis_out = sapply((i+1):group_n,function(j){
      id = which(cell_group_mat[,j] > 0)

      x = dis_i[,to %in% id,drop = FALSE]
      # cat("The size of the distance matrix:",dim(x),"\n")

      out = group_2_min(x,group1 = 1:nrow(x),group2 = 1:ncol(x), k = k)

      return(out)
    })

    out = cbind(i,(i+1):group_n,dis_out)
    return(out)
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
#' @examples
#' dist_mat <- matrix(runif(100), nrow = 10)
#' d <- group_2_min(dist_mat, group1 = 1:3, group2 = 4:6, k = 2)
#' print(d)
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
