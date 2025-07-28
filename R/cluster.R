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
  snn = sNN(data,k = k)
  snn$jaccard = snn$shared/(2*k-snn$shared)

  snn_edge = as.data.frame(cbind(rep(1:nrow(data),each = k),c(t(snn$id)),
                                     c(t(snn$dis)),c(t(snn$jaccard))))
  colnames(snn_edge) = c("start","end","dis","jaccard")
  snn_edge = snn_edge %>% filter(jaccard > prune.snn)

  colnames(snn_edge)[which(colnames(snn_edge) == weight)] = "weight"
  snn_graph = graph_from_data_frame(snn_edge,directed = FALSE,vertices = 1:nrow(data))

  # cluster = cluster_leiden(dis_snn_graph,resolution = resolution)
  cluster = igraph::cluster_louvain(dis_snn_graph,resolution = resolution)
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
#' @examples
#' mat <- matrix(rnorm(500), nrow = 100)
#' clusters <- leiden_embedding_fast(mat, k = 20)
#' table(clusters)
#'
#' @export
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
#' @importFrom dbscan sNN
#' @importFrom igraph graph_from_data_frame cluster_leiden
#' @importFrom uwot umap
#' @importFrom dplyr filter
#' @examples
#' mat <- matrix(rnorm(500), nrow = 100)
#' dmat <- dist(mat)
#' clusters <- leiden_dis(as.matrix(dmat), k = 15, if_umap = FALSE)
#'
#' @export
leiden_dis = function(dismat,k = 10,prune.snn = 0,weight = "jaccard",resolution = 1,if_umap = TRUE){
    if(if_umap){
      group_umap = umap::umap(dismat,input="dist")
      group_umap = as.data.frame(group_umap$layout)
      colnames(group_umap) = c("umap1","umap2")
    }

    dis_snn = sNN(as.dist(dismat),k = k)
    dis_snn$jaccard = dis_snn$shared/(2*k-dis_snn$shared)

    dis_snn_edge = as.data.frame(cbind(rep(1:nrow(dismat),each = k),c(t(dis_snn$id)),
                           c(t(dis_snn$dis)),c(t(dis_snn$jaccard))))
    colnames(dis_snn_edge) = c("start","end","dis","jaccard")
    dis_snn_edge = dis_snn_edge %>% filter(jaccard > prune.snn)

    colnames(dis_snn_edge)[which(colnames(dis_snn_edge) == weight)] = "weight"
    dis_snn_graph = graph_from_data_frame(dis_snn_edge,directed = FALSE,vertices = 1:nrow(dismat))

    cluster = cluster_leiden(dis_snn_graph,resolution = resolution,objective_function = "modularity")
    # cluster = cluster_louvain(dis_snn_graph,resolution = resolution)
    # cluster = leiden(dis_snn_graph,resolution_parameter = resolution)

    if(if_umap){
      group_umap$cluster = as.factor(cluster$membership)
      # group_umap$cluster = as.factor(cluster)
      return(group_umap)
    }
    return(as.factor(cluster$membership))
    # return(as.factor(cluster))
}
