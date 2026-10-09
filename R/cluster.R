#' @title leiden_embedding
#' @description Run Leiden Clustering from Embedding via Shared Nearest Neighbors
#' @details
#'
#' Constructs a shared nearest neighbor (SNN) graph from an embedding matrix (e.g., PCA, UMAP),
#' optionally prunes weak connections, and performs Louvain or Leiden clustering.
#'
#' @param data A numeric matrix or data frame where rows are observations and columns are features.
#' @param k Integer. Number of nearest neighbors for SNN construction (default: 30).
#' @param prune.snn Numeric. Threshold below which SNN edges are removed (default: 0).
#' @param weight Character. Column to use as edge weight (`"jaccard"` or `"dis"`) (default: `"jaccard"`).
#' @param resolution Numeric resolution parameter for clustering (default: 1).
#'
#' @return A factor vector of cluster memberships for each row in `data`.
#'
#' @importFrom dbscan sNN
#' @importFrom igraph graph_from_data_frame cluster_louvain
#' @importFrom dplyr filter
#' @examples
#' mat <- matrix(rnorm(500), nrow = 100)
#' clusters <- leiden_embedding(mat, k = 20)
#' table(clusters)
#'
#' @export
leiden_embedding = function(data,k = 30,prune.snn = 0,weight = "jaccard",resolution = 1){
  snn = dbscan::sNN(data,k = k)
  snn$jaccard = snn$shared/(2*k-snn$shared)

  snn_edge = as.data.frame(cbind(rep(1:nrow(data),each = k),c(t(snn$id)),
                                     c(t(snn$dis)),c(t(snn$jaccard))))
  colnames(snn_edge) = c("start","end","dis","jaccard")
  snn_edge = dplyr::filter(snn_edge, jaccard > prune.snn)

  colnames(snn_edge)[which(colnames(snn_edge) == weight)] = "weight"
  snn_graph = igraph::graph_from_data_frame(snn_edge,directed = FALSE,vertices = 1:nrow(data))

  # cluster = cluster_leiden(dis_snn_graph,resolution = resolution)
  cluster = igraph::cluster_louvain(snn_graph,resolution = resolution)
  # cluster = leiden::leiden(snn_graph,resolution_parameter = resolution)

  return(as.factor(cluster$membership))
  # return(as.factor(cluster))
}

#' @title leiden_embedding_fast
#' @description Fast Leiden Clustering Using Sparse Shared Nearest Neighbors
#' @details
#'
#' Efficient implementation of SNN graph construction and Leiden clustering.
#' Uses matrix algebra for fast computation of shared neighbors and Jaccard similarity.
#'
#' @param data A numeric matrix of embeddings (e.g., PCA, UMAP).
#' @param k Integer. Number of nearest neighbors (default: 30).
#' @param prune.snn Numeric threshold to prune weak SNN edges (default: 0).
#' @param weight Character. Weighting scheme (currently unused; default: `"jaccard"`).
#' @param resolution Resolution parameter for Leiden clustering (default: 1).
#'
#' @return A factor vector of cluster labels.
#'
#' @importFrom RANN nn2
#' @importFrom Matrix sparseMatrix
#' @importFrom igraph graph_from_adjacency_matrix cluster_leiden
leiden_embedding_fast <- function(data, k = 30, prune.snn = 0, weight = "jaccard", resolution = 1) {
  # 1. Compute k-nearest neighbors (Faster than dbscan::sNN)
  knn <- RANN::nn2(data, k = k + 1)  # +1 to exclude self

  # 2. Compute shared neighbors and Jaccard weights (sparse matrix)
  n_cells <- nrow(data)
  idx <- knn$nn.idx[, -1]  # Remove self (1st column)

  # Sparse adjacency matrix (binary)
  adj <- Matrix::sparseMatrix(
    i = rep(1:n_cells, each = k),
    j = as.vector(idx),
    dims = c(n_cells, n_cells)
  )

  # Shared neighbors: A %*% t(A) [faster than dbscan::sNN]
  shared <- adj %*% Matrix::t(adj)

  # Jaccard weights: shared / (2k - shared)
  jaccard <- shared / (2*k - shared)

  # 3. Prune edges below threshold
  jaccard[jaccard <= prune.snn] <- 0
  diag(jaccard) <- 0
  jaccard <- Matrix::drop0(jaccard)
  jaccard <- methods::as(jaccard, "dgCMatrix")

  # 4. Convert to igraph (weighted)
  snn_graph <- igraph::graph_from_adjacency_matrix(
    jaccard,
    mode = "undirected",
    weighted = TRUE
  )

  # 5. Run Leiden (faster than Louvain)
  cluster <- igraph::cluster_leiden(
    snn_graph,
    resolution_parameter = resolution,
    objective_function = "modularity"
  )

  return(as.factor(cluster$membership))
}

#' @title leiden_dis
#' @description Leiden Clustering from Precomputed Distance Matrix
#' @details
#'
#' Constructs an SNN graph from a distance matrix, applies Leiden clustering,
#' and optionally returns a UMAP layout annotated with clusters.
#'
#' @param dismat A symmetric distance matrix (e.g., from `dist()`).
#' @param k Integer. Number of neighbors for SNN construction (default: 10).
#' @param prune.snn Numeric threshold to prune SNN edges (default: 0).
#' @param weight Character. Column to use as edge weight (`"jaccard"` or `"dis"`) (default: `"jaccard"`).
#' @param resolution Numeric resolution parameter for Leiden clustering (default: 1).
#' @param if_umap Logical. If `TRUE`, returns UMAP coordinates with cluster annotations (default: TRUE).
#'
#' @return If `if_umap = TRUE`, returns a data frame with UMAP coordinates and cluster labels.
#' Otherwise, returns a factor vector of cluster memberships.
#'
#' @importFrom igraph graph_from_data_frame cluster_leiden
#' @importFrom uwot umap
#' @importFrom dplyr filter
#' @importFrom stats as.dist
#' @examples
#' mat <- matrix(rnorm(500), nrow = 100)
#' dmat <- dist(mat)
#' clusters <- leiden_dis(as.matrix(dmat), k = 15, if_umap = FALSE)
#'
#' @export
leiden_dis = function(dismat,k = 10,prune.snn = 0,weight = "jaccard",resolution = 1,if_umap = TRUE){
    node_names <- if (inherits(dismat, "dist")) attr(dismat, "Labels") else rownames(dismat)
    dismat <- as.matrix(dismat)
    storage.mode(dismat) <- "double"
    if (nrow(dismat) != ncol(dismat) || nrow(dismat) < 2L) {
      stop("`dismat` must be a square distance matrix with at least two rows.")
    }
    if (any(!is.finite(dismat)) || any(dismat < 0)) {
      stop("`dismat` must contain finite, non-negative distances.")
    }
    if (!isTRUE(all.equal(dismat, t(dismat), tolerance = sqrt(.Machine$double.eps), check.attributes = FALSE))) {
      stop("`dismat` must be symmetric.")
    }
    diag(dismat) <- 0
    if (is.null(node_names)) node_names <- as.character(seq_len(nrow(dismat)))
    if (length(node_names) != nrow(dismat) || anyNA(node_names) || anyDuplicated(node_names)) {
      stop("Distance labels must be unique and match the matrix dimensions.")
    }
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1L || k >= nrow(dismat)) {
      stop("`k` must be a positive integer smaller than the number of rows in `dismat`.")
    }
    k <- as.integer(k)
    weight <- match.arg(weight, c("jaccard", "dis"))
    if (!is.numeric(prune.snn) || length(prune.snn) != 1L || is.na(prune.snn)) {
      stop("`prune.snn` must be one numeric value.")
    }
    if(if_umap){
      umap_config <- umap::umap.defaults
      umap_config$n_neighbors <- k
      group_umap = umap::umap(dismat, input = "dist", config = umap_config)
      group_umap = as.data.frame(group_umap$layout)
      colnames(group_umap) = c("umap1","umap2")
      rownames(group_umap) <- node_names
    }

    neighbor_id <- t(vapply(seq_len(nrow(dismat)), function(i) {
      candidates <- setdiff(seq_len(nrow(dismat)), i)
      candidates[order(dismat[i, candidates], na.last = NA)[seq_len(k)]]
    }, integer(k)))
    neighbor_dis <- matrix(
      dismat[cbind(rep(seq_len(nrow(dismat)), each = k), as.vector(t(neighbor_id)))],
      nrow = nrow(dismat), byrow = TRUE
    )
    shared <- matrix(0L, nrow = nrow(dismat), ncol = k)
    for (i in seq_len(nrow(dismat))) {
      for (j in seq_len(k)) {
        shared[i, j] <- length(intersect(neighbor_id[i, ], neighbor_id[neighbor_id[i, j], ]))
      }
    }
    jaccard <- shared / (2 * k - shared)

    dis_snn_edge = as.data.frame(cbind(rep(seq_len(nrow(dismat)), each = k),
                           as.vector(t(neighbor_id)), as.vector(t(neighbor_dis)), as.vector(t(jaccard))))
    colnames(dis_snn_edge) = c("start","end","dis","jaccard")
    dis_snn_edge = dplyr::filter(dis_snn_edge, jaccard > prune.snn)

    colnames(dis_snn_edge)[which(colnames(dis_snn_edge) == weight)] = "weight"
    vertices <- data.frame(name = seq_len(nrow(dismat)), label = node_names)
    dis_snn_graph = igraph::graph_from_data_frame(dis_snn_edge,directed = FALSE,vertices = vertices)

    cluster = igraph::cluster_leiden(dis_snn_graph,resolution = resolution,objective_function = "modularity")
    # cluster = cluster_louvain(dis_snn_graph,resolution = resolution)
    # cluster = leiden(dis_snn_graph,resolution_parameter = resolution)

    if(if_umap){
      group_umap$cluster = as.factor(cluster$membership)
      group_umap$clone = node_names
      # group_umap$cluster = as.factor(cluster)
      return(group_umap)
    }
    membership <- as.factor(cluster$membership)
    names(membership) <- node_names
    return(membership)
    # return(as.factor(cluster))
}
