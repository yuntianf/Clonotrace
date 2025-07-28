#' @title acct
#'
#' @description Compute Accumulated Commute Time Matrix from Transition Matrix
#' @details
#'
#' Computes the accumulated commute time (ACT) matrix from a normalized transition matrix.
#' This is used as a basis for diffusion-based distances and pseudotime inference.
#'
#' @param T_mat A square transition matrix (symmetric, row-normalized).
#'
#' @return A matrix of accumulated commute times between all pairs of nodes.
#'
#' @importFrom Matrix Diagonal tcrossprod solve
#' @importFrom RSpectra eigs
#' @examples
#' mat <- matrix(runif(100), nrow = 10)
#' mat <- mat / rowSums(mat)
#' act <- acct(mat)
#' print(dim(act))
#'
acct = function(T_mat) {
  # Step 1: Compute the dominant eigenvalue and eigenvector using RSpectra
  ev = RSpectra::eigs(T_mat, k = 1)
  phi0 = as.numeric(ev$vectors)

  # Step 2: Compute T_bar using sparse matrix operations
  # T_bar = T_mat - phi0 %*% t(phi0) can be computed more efficiently
  I = Matrix::Diagonal(nrow(T_mat))  # Identity matrix as a sparse matrix
  # P = matrix(as.numeric(phi0 %*% t(phi0)),nrow = nrow(I))
  A = I + Matrix::tcrossprod(phi0, phi0) - T_mat

  # Use iterative solver for solving (I - T_bar)^{-1} * something
  # This avoids the dense inversion of the matrix
  M = Matrix::solve(A,I,method = "CG") - I  # Solve sparse system and subtract identity matrix

  return(M)
}

#' @title DPT_T
#'
#' @description Compute Diffusion Pseudotime from Accumulated Commute Time Matrix
#' @details
#'
#' Computes diffusion pseudotime (DPT) distances from a specified root node using the accumulated commute time matrix.
#'
#' @param T_mat A row-normalized transition matrix (symmetric or stochastic).
#' @param start Integer. Index of the root cell for pseudotime computation.
#'
#' @return A numeric vector of diffusion pseudotime values for all nodes.
#'
#' @examples
#' T <- matrix(runif(100), 10, 10)
#' T <- T / rowSums(T)
#' pt <- DPT_T(T, start = 1)
#'
DPT_T = function(T_mat,start){
  M = acct(T_mat)

  dpt = sqrt(rowSums((t(t(M)-M[start,]))^2))
  return(dpt)
}

#' @title dpt
#'
#' @description Compute Diffusion Pseudotime from Transition Matrix (Eigen Decomposition)
#' @details
#'
#' Calculates pseudotime by spectral decomposition of a diffusion transition matrix,
#' using a root cell as a reference point in the reduced diffusion space.
#'
#' @param T_mat A square transition matrix (e.g., from kNN graph).
#' @param root Integer. Index of the root cell.
#' @param k Integer. Number of eigenvectors to use (default: 30).
#'
#' @return A numeric vector of pseudotime values normalized to [0, 1].
#'
#' @importFrom RSpectra eigs
#' @examples
#' T <- matrix(runif(100), 10, 10)
#' T <- T / rowSums(T)
#' pt <- dpt(T, root = 1)
#'
dpt = function(T_mat,root,k = 30){
  T_eigen = RSpectra::eigs(T_mat,k = k)

  eigen_value = as.numeric(T_eigen$values)
  eigen_vector = matrix(as.numeric(T_eigen$vectors),nrow(T_mat),k)

  flag = eigen_value < 0.9999

  eigen_value = eigen_value[flag]
  eigen_vector = eigen_vector[,flag]

  eigen_value = eigen_value/(1-eigen_value)
  eigen_vector = t(t(eigen_vector)*eigen_value)

  pseudotime = sqrt(rowSums((t(t(eigen_vector) - eigen_vector[root,]))^2))
  pseudotime = as.numeric(pseudotime)
  pseudotime = pseudotime/max(pseudotime)

  # return(list(pseudotime,eigen_vector))
  return(pseudotime)
}


#' @title embedding2dpt
#'
#' @description Compute Diffusion Pseudotime from Embedding
#' @details
#'
#' Builds a kNN graph from an embedding matrix and computes diffusion pseudotime starting from a root cell.
#'
#' @param embedding A numeric matrix of embeddings (e.g., UMAP or PCA).
#' @param nn_k Integer. Number of neighbors in the kNN graph.
#' @param root Integer. Index of the root cell.
#' @param dpt_k Integer. Number of eigenvectors used in DPT (default: 30).
#'
#' @return A numeric vector of pseudotime values for each row in `embedding`.
#'
#' @examples
#' mat <- matrix(rnorm(200), ncol = 5)
#' pt <- embedding2dpt(mat, nn_k = 10, root = 1)
#'
#' @export
embedding2dpt = function(embedding,nn_k,root,dpt_k = 30){
  knn = embedding2knn(embedding,nn_k)
  T_mat = knn/rowSums(knn)

  pseudotime = dpt(T_mat = T_mat,root = root,k = dpt_k)
  return(pseudotime)
}

#' @title clone_root
#'
#' @description Identify Root Clone from Cluster Enrichment
#' @details
#'
#' Finds the most enriched clone in a given starting cluster, to use as a pseudotime root.
#'
#' @param clones A character vector of clone identifiers.
#' @param cell_meta A data frame with clone and cluster information per cell.
#' @param clone_col Name of the column containing clone labels.
#' @param cluster_col Name of the column containing cluster labels.
#' @param start_cluster Name or value of the cluster to define the root state.
#'
#' @return The name of the most enriched clone in the start cluster.
#'
#' @importFrom dplyr filter_at group_by_at summarise arrange
#' @examples
#' root_clone <- clone_root(clones = c("A", "B", "C"),
#'                          cell_meta = cell_metadata,
#'                          clone_col = "clone_id",
#'                          cluster_col = "cluster",
#'                          start_cluster = "Naive")
#'
clone_root = function(clones,cell_meta,clone_col,cluster_col,start_cluster){
  root = cell_meta %>% filter_at(clone_col,~. %in% clones) %>%
    group_by_at(clone_col) %>% summarise(ratio = sum(!!sym(cluster_col) == start_cluster)/n(),count = n()) %>%
    arrange(-ratio,-count)

  return(unlist(root[,clone_col])[1])
}

#' @title clone_root
#'
#' @description Compute Clone-Level Diffusion Pseudotime
#' @details
#'
#' Estimates pseudotime on the clone embedding by building a kNN graph and applying diffusion pseudotime from a root clone.
#'
#' @param clone_embedding A numeric matrix where rows are clones and columns are features.
#' @param cell_meta A data frame of metadata with clone and cluster annotations.
#' @param clone_col Name of the column containing clone identifiers.
#' @param cluster_col Name of the column containing cluster labels.
#' @param start_cluster The cluster label used to select the root clone.
#' @param k Integer. Number of neighbors for kNN graph construction (default: 10).
#' @param dpt_k Integer. Number of eigenvectors for DPT (default: 30).
#'
#' @return A numeric vector of clone-level pseudotime values.
#'
#' @examples
#' clone_pt <- clone_dpt(clone_embedding, cell_meta, "clone", "cluster", "Naive")
#' hist(clone_pt)
#'
#' @export
clone_dpt = function(clone_embedding,cell_meta,
                     clone_col,cluster_col,start_cluster,
                     k = 10,dpt_k = 30){
  clone_knn = embedding2knn(clone_embedding,k)

  T_mat = clone_knn/rowSums(clone_knn)

  clones = rownames(clone_embedding)
  root = clone_root(clones,cell_meta,clone_col,cluster_col,start_cluster)

  pseudotime = dpt(T_mat = T_mat,root = which(clones == root),k = dpt_k)

  return(pseudotime)
}


