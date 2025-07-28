#' @title cell_clone_coembed
#'
#' @description Construct Cell-Cell Distance Matrix via Clone-Aware Coembedding
#' @details
#'
#' Integrates cell and clone embeddings by using clone-level proximity to reweight cell-cell distances.
#' Generates a sparse distance matrix suitable for downstream embedding or trajectory inference.
#'
#' @param cell_embedding A numeric matrix of cell embeddings (rows = cells).
#' @param clone_embedding A numeric matrix of clone embeddings (rows = clones).
#' @param cell_k Integer. Number of neighbors to compute for cell-level kNN (default: 30).
#' @param clone_k Integer. Number of neighbors to compute for clone-level kNN (default: 15).
#'
#' @return A symmetric sparse distance matrix (class `dgCMatrix`) encoding clone-informed cell-cell distances.
#'
#' @importFrom dplyr filter mutate group_by slice_min
#' @importFrom Matrix drop0
#'
#' @export
cell_clone_coembed = function(cell_embedding,clone_embedding,cell_k = 30,clone_k = 15){
  cell_knn_flat = knn_flat(cell_embedding,k = cell_k,if_dedup = TRUE,symmetric = FALSE,if_self = TRUE)
  clone_knn = embedding2knn(clone_embedding,clone_k,if_self = TRUE)
  clone_knn@x = rep(1,length(clone_knn@x))

  cell_knn_flat$weight = cell_knn_matrix_multiplication_parallel(cell_knn_flat,cell_clone_prob,clone_knn,chunk_size = 5000)

  cell_knn_mat = long_symmetry(cell_knn_flat,row_names_from = "node1",col_names_from = "node2")

  distance = cell_knn_mat %>% filter(weight > 0.1) %>% mutate(dis = dist/weight)
  distance = distance %>%  group_by(node1) %>%  slice_min(dis,n = 20)

  distance = drop0(long2sparse(distance,
                               row_names_from = "node1",col_names_from = "node2",
                               values_from = "dis",
                               unique_rows = 1:nrow(cell_embedding),
                               unique_cols = 1:nrow(cell_embedding),symmetric = TRUE))
  rownames(distance) = colnames(distance) = rownames(cell_embedding)
  diag(distance) = 0
  distance = drop0(distance)

  return(distance)
}

#' @title cell_knn_matrix_mutiplication
#'
#' @description Matrix Multiplication for Cell-Feature Weighted Edges
#' @details
#'
#' Computes a score for each kNN edge using the dot product between a cell’s features and a transformed clone matrix.
#' Used to integrate clone-to-clone and cell-to-clone relationships into edge-level scores.
#'
#' @param knn A matrix or data frame with two columns: `i` (source cell) and `j` (neighbor cell).
#' @param cell_feature_mat A numeric matrix of cell-by-feature values.
#' @param feature_feature_mat A matrix representing relationships between features (e.g., clone similarity).
#'
#' @return A numeric vector of scores, one per kNN edge.
cell_knn_matrix_mutiplication = function(knn,cell_feature_mat,feature_feature_mat){
  i = knn[,1]
  j = knn[,2]

  intermediate <- feature_feature_mat %*% t(cell_feature_mat)

  cell_subset <- cell_feature_mat[i, , drop = FALSE]
  intermediate_subset <- intermediate[, j, drop = FALSE]

  result_vals <- rowSums(cell_subset * t(intermediate_subset))
  return(result_vals)
}

#' @title cell_knn_matrix_multiplication_parallel
#'
#' @description Parallelized Matrix Multiplication for kNN Edge Weighting
#' @details
#'
#' Computes edge-level scores via parallelized matrix multiplication between cell and feature matrices,
#' allowing scalable construction of clone-informed kNN graphs.
#'
#' @param knn A matrix or data frame of kNN edges (columns: `i`, `j`).
#' @param cell_feature_mat A numeric matrix of cell-by-feature weights.
#' @param feature_feature_mat A sparse feature similarity matrix (e.g., clone-by-clone).
#' @param chunk_size Integer. Number of edges to process per parallel chunk (default: 5000).
#'
#' @return A numeric vector of computed scores, one per kNN edge.
#'
#' @importFrom future.apply future_lapply
#' @importFrom Matrix rowSums
cell_knn_matrix_multiplication_parallel <- function(knn, cell_feature_mat, feature_feature_mat, chunk_size = 5000) {
  # Convert to optimal formats
  feature_feature_mat <- as(feature_feature_mat, "CsparseMatrix")
  t_cell_feature <- t(cell_feature_mat)

  # Split indices into chunks
  n_edges <- nrow(knn)
  chunks <- split(1:n_edges, ceiling(seq_along(1:n_edges)/chunk_size))

  # Parallel process chunks
  result <- future_lapply(chunks, function(chunk_idx) {
    i <- knn[chunk_idx, 1]
    j <- knn[chunk_idx, 2]

    intermediate_chunk <- feature_feature_mat %*% t_cell_feature[, j, drop = FALSE]
    cell_subset <- cell_feature_mat[i, , drop = FALSE]

    Matrix::rowSums(cell_subset * Matrix::t(intermediate_chunk))
  }, future.seed = TRUE,future.packages = c("Matrix"))

  return(unlist(result))
}
