#' @title label_spreading
#' @description Label Propagation via Iterative Graph-Based Spreading
#' @details
#'
#' Performs semi-supervised label propagation over a graph using a sparse or dense adjacency matrix.
#' The algorithm propagates known labels across the graph structure, allowing soft label assignment
#' for unlabeled nodes.
#'
#' @param adj A square adjacency matrix (preferably sparse) representing the graph.
#' @param labels An integer vector of length equal to the number of nodes. Use `NA` for unlabeled entries.
#' @param label_n Optional. The number of label classes. If `NULL`, inferred as `max(labels, na.rm = TRUE)`.
#' @param alpha Float in (0, 1). The propagation coefficient controlling the balance between prior and propagated labels (default: 0.9).
#' @param max_iter Maximum number of iterations for propagation (default: 100).
#' @param tol Convergence threshold (default: 1e-3).
#' @param epsilon Small prior assigned to unlabeled entries (default: 0).
#' @param verbose Logical. Whether to print progress and convergence status (default: TRUE).
#'
#' @return A matrix of size N x C where each row contains soft label probabilities for a node across `C` classes.
#'
#' @importFrom Matrix rowSums Diagonal
#' @examples
#' \dontrun{
#' set.seed(1)
#' adj <- Matrix::rsparsematrix(100, 100, density = 0.05)
#' labels <- rep(NA, 100)
#' labels[1:10] <- sample(1:3, 10, replace = TRUE)
#' prob_matrix <- label_spreading(adj, labels)
#' }
label_spreading <- function(adj, labels, label_n = NULL, alpha = 0.9, max_iter = 100, tol = 1e-3, epsilon = 0, verbose = TRUE) {
  # Initial setup
  labels = as.numeric(as.factor(labels))
  N <- length(labels)
  C <- if (is.null(label_n)) max(labels, na.rm = TRUE) else label_n

  # Precompute label weights (log2(count)/count)
  label_counts <- tabulate(na.omit(labels), nbins = C)
  label_weights <- ifelse(label_counts > 0, log2(label_counts)/label_counts, 0)

  # Initialize Y matrix more efficiently
  Y <- matrix(epsilon/C, nrow = N, ncol = C)
  labeled_idx <- which(!is.na(labels))
  Y[labeled_idx, ] <- 0
  Y[cbind(labeled_idx, labels[labeled_idx])] <- label_weights[labels[labeled_idx]]

  # Row-normalize adjacency matrix (sparse-safe)
  degrees <- Matrix::rowSums(adj)
  Dinv <- Matrix::Diagonal(x = 1 / pmax(degrees, 1e-6))
  P <- Dinv %*% adj

  # Iterative propagation with optimized matrix operations
  F_current <- Y
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    F_new <- alpha * (P %*% F_current) + (1 - alpha) * Y

    # Faster convergence check using L1 norm
    diff <- max(abs(F_new - F_current))
    if (verbose && (iter %% 10 == 0 || iter == 1 || diff < tol)) {
      message(sprintf("Iter %3d: loss = %.4e", iter, diff))
    }

    if (diff < tol) {
      converged <- TRUE
      break
    }

    F_current <- F_new
  }

  if (verbose) {
    if (converged) {
      message(sprintf("Converged in %d iterations", iter))
    } else {
      message(sprintf("Reached maximum iterations (%d)", max_iter))
    }
  }

  return(F_current)
}


#' @title label_spreading_bootstrap
#' @description Bootstrap-Based Stability Estimation for Label Propagation
#' @details
#'
#' Repeatedly applies `label_spreading()` on subsampled label sets to assess the stability or
#' uncertainty of label propagation results. Returns a node-level deviance score indicating variability across bootstrap runs.
#'
#' @inheritParams label_spreading
#' @param refer Optional. A reference soft label matrix (N x C). If not provided, computed from full label set.
#' @param sample_rate Fraction of labeled nodes used in each bootstrap sample (default: 0.8).
#' @param sample_n Number of bootstrap replicates (default: 50).
#' @param ... Additional arguments passed to `label_spreading()`.
#'
#' @return A numeric vector of length N giving the deviance between each node's label probabilities and the reference across bootstraps.
#'
#' @importFrom future.apply future_lapply
#' @importFrom Matrix rowSums
#' @examples
#' \dontrun{
#' adj <- Matrix::rsparsematrix(100, 100, density = 0.05)
#' labels <- rep(NA, 100)
#' labels[1:10] <- sample(1:3, 10, replace = TRUE)
#' deviance_scores <- label_spreading_bootstrap(adj, labels)
#' }
#' @export
label_spreading_bootstrap <- function(adj, labels,refer = NULL, alpha = 0.8, sample_rate = 0.8, sample_n = 50, ...) {
  # estimate with all labels
  labels = as.numeric(as.factor(labels))
  if(is.null(refer)){
    refer <- label_spreading(
      adj = adj,
      labels = labels,
      alpha = alpha,
      verbose = FALSE,
      ...
    )
  }

  refer = as.matrix(refer / rowSums(refer))

  # Initial setup
  labeled_idx <- which(!is.na(labels))
  label_sample_count <- round(length(labeled_idx) * sample_rate)
  label_n <- max(labels, na.rm = TRUE)

  full_label_flag = labels
  full_label_flag[is.na(full_label_flag)] = 0
  # Main bootstrap loop
  result_list <- future_lapply(1:sample_n, function(i) {
    # Create subsampled labels
    subsampled_labels <- rep(NA, length(labels))
    sampled_idx <- sample(labeled_idx, size = label_sample_count, replace = FALSE)
    subsampled_labels[sampled_idx] <- labels[sampled_idx]

    sample_label_flag = subsampled_labels
    sample_label_flag[is.na(sample_label_flag)] = 0

    # Perform label spreading
    prob_matrix <- label_spreading(
      adj = adj,
      labels = subsampled_labels,
      label_n = label_n,
      alpha = alpha,
      epsilon = 0,
      verbose = FALSE,
      ...
    )

    # Normalize probabilities
    prob_matrix = prob_matrix + 1e-12
    prob_matrix = prob_matrix / rowSums(prob_matrix)
    prob_matrix = as.matrix(prob_matrix)

    L1 = rowSums(abs(prob_matrix-refer))
    flag = !xor(full_label_flag,sample_label_flag)

    return(cbind(L1,flag))
  }, future.seed = TRUE)


  deviance <- do.call(cbind,result_list)

  flag = deviance[,seq(2,ncol(deviance),2)]
  deviance = deviance[,seq(1,ncol(deviance),2)]
  # deviance = rowMeans(deviance)

  deviance_norm = rowSums(deviance*flag)/rowSums(flag)
  return(list(prob = refer,deviance = deviance_norm))
}
