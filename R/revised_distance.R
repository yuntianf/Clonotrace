#' @title revise_clone_distance
#' @description Combine the Clonotrace Clone Distance With the Clone-Mean Distance, Weighted Without Labels
#' @details
#'
#' The revised Clonotrace distance, which `clone_disance()` returns by default. It is built from the
#' graph-based clone distance in three steps:
#'
#' 1. **Local rescaling.** Each clone distance is divided by `sqrt(sigma_i * sigma_j)`, where `sigma_i`
#'    is the distance from clone `i` to its `local_k`-th nearest clone. Clones in dense and sparse
#'    regions of clone space then sit on a common scale.
#' 2. **Clone-mean distance.** The Euclidean distance between the clones' mean positions in the cell
#'    embedding. It is not rescaled.
#' 3. **Blend.** Both distances are divided by their median and combined as
#'    `(1 - s) * clonotrace + s * mean`. The weight `s` is chosen without cell-type labels. For each
#'    Leiden resolution in `resolutions`, the cells are clustered, and `s` is set to the grid value
#'    whose blended distance has the highest Spearman correlation with the clones' cluster-composition
#'    difference (total variation distance). The final `s` is the median across resolutions.
#'
#' In the benchmark (14 CD8 T-cell cohorts), the revised distance matched a weight fitted on
#' cell-type labels (mean Spearman 0.766 against 0.763) and beat the unrevised distance (0.725).
#' Its cost is small next to the clone distance itself: one Leiden ladder over the cells and some
#' arithmetic over the clone pairs.
#'
#' @param clone_dis Clone distances, either as a square matrix (clones x clones) or as the long output
#'   of `graph_clone_nn()`/`graph_clone_ot()` (columns `group1`, `group2`, `dis`, with integer indices
#'   into `clone_names`).
#' @param embedding A numeric matrix of cell embeddings (cells x dimensions) with cell names as row names;
#'   normally the same embedding passed to `clone_disance()`.
#' @param clone_labels Observed clone membership, used for the clone means and the cluster
#'   compositions. Either a data frame with columns `cell` and `clone` (the `cell_clone` table built
#'   before label spreading) or a named vector (names are cells). Cells of other clones and `NA` are
#'   ignored. If `NULL`, `cell_clone_prob` is used instead, weighting cells by their clone probability.
#' @param cell_clone_prob A cells x clones probability matrix. Only used when `clone_labels` is `NULL`.
#' @param clone_names Clone names, in the order of `clone_dis`. Taken from the matrix dimnames or from
#'   `colnames(cell_clone_prob)` when missing.
#' @param sample Optional clone-level grouping, such as the patient or sample of each clone (a vector
#'   named by clone, or in the order of `clone_names`). When given, medians and the weight `s` are
#'   computed on pairs of clones from the same group only. The benchmark used the patient.
#' @param local_k Integer. Neighbour rank for the local scale (default: 7).
#' @param resolutions Numeric vector of Leiden resolutions for choosing `s` (default: 0.3, 0.5, 0.8, 1.2).
#' @param clusters Optional precomputed cell clusters, a data frame with one column per resolution and
#'   cell names as row names. Replaces the Leiden step.
#' @param cluster_k Integer. Number of neighbours for the cell graph used by Leiden (default: 15).
#' @param mean_weight Optional fixed weight `s` in `[0, 1]`. Skips the selection step.
#' @param s_grid The weights tried (default: 0 to 1 by 0.05).
#' @param seed Integer seed for the Leiden step (default: 0).
#' @param verbose Logical. Whether to print the chosen weight (default: TRUE).
#'
#' @return A symmetric clones x clones matrix of revised distances. Attributes: `mean_weight` (the
#'   weight `s`), `picks` (the `s` chosen at each resolution) and `clonotrace` (the input distance
#'   as a square matrix).
#'
#' @importFrom stats median cor
#' @importFrom Matrix sparseMatrix
#' @importFrom uwot similarity_graph
#' @importFrom igraph graph_from_adjacency_matrix cluster_leiden
#' @export
revise_clone_distance = function(clone_dis, embedding, clone_labels = NULL, cell_clone_prob = NULL,
                                 clone_names = NULL, sample = NULL, local_k = 7,
                                 resolutions = c(0.3, 0.5, 0.8, 1.2), clusters = NULL,
                                 cluster_k = 15, mean_weight = NULL,
                                 s_grid = seq(0, 1, by = 0.05), seed = 0, verbose = TRUE){
  is_long = all(c("group1", "group2", "dis") %in% colnames(clone_dis))
  if(is.null(clone_names)){
    clone_names = if(!is_long) colnames(clone_dis) else colnames(cell_clone_prob)
  }
  if(is.null(clone_names)){
    stop("clone_names is required when neither clone_dis nor cell_clone_prob carries clone names.")
  }
  D = clone_dis_square(clone_dis, clone_names, is_long)
  n = nrow(D)
  if(n <= local_k){
    stop("There are ", n, " clones; the local scale needs more than local_k = ", local_k, ".")
  }
  embedding = as.matrix(embedding)
  if(is.null(rownames(embedding))){
    stop("The embedding needs cell names as row names.")
  }

  # cell -> clone index, or a cells x clones weight matrix when only probabilities are given
  W = clone_cell_weights(embedding, clone_labels, cell_clone_prob, clone_names)
  mass = Matrix::colSums(W)
  if(any(mass == 0)){
    stop("No cells found for clone(s): ", paste(clone_names[mass == 0][seq_len(min(5, sum(mass == 0)))], collapse = ", "))
  }
  M = as.matrix(Matrix::t(W) %*% embedding) / mass
  Dm = sqrt(pmax(outer(rowSums(M^2), rowSums(M^2), "+") - 2 * tcrossprod(M), 0))

  # pairs used for the medians and the weight: all pairs, or within-sample pairs
  pairs = which(upper.tri(D), arr.ind = TRUE)
  if(!is.null(sample)){
    if(!is.null(names(sample))){
      sample = sample[clone_names]
    }
    if(length(sample) != n || anyNA(sample)){
      stop("sample must give a group for every clone.")
    }
    pairs = pairs[sample[pairs[, 1]] == sample[pairs[, 2]], , drop = FALSE]
  }
  ii = pairs[, 1]
  jj = pairs[, 2]

  sigma = pmax(apply(D, 1, function(x) sort(x, partial = local_k + 1)[local_k + 1]), 1e-12)
  D_loc = D / sqrt(outer(sigma, sigma))
  zc = D_loc / stats::median(D_loc[pairs])
  zm = Dm / stats::median(Dm[pairs])

  picks = NULL
  if(is.null(mean_weight)){
    if(is.null(clusters)){
      clusters = leiden_ladder(embedding, resolutions, k = cluster_k, seed = seed)
    }
    zc_p = zc[pairs]
    zm_p = zm[pairs]
    picks = sapply(colnames(clusters), function(r){
      lab = clusters[rownames(embedding), r]
      ok = !is.na(lab)
      H = Matrix::sparseMatrix(i = which(ok), j = as.integer(factor(lab[ok])), x = 1,
                               dims = c(nrow(embedding), length(unique(lab[ok]))))
      C = as.matrix(Matrix::t(W) %*% H)
      keep = rowSums(C) > 0
      m = keep[ii] & keep[jj]
      ry = rank(composition_tv(C, ii[m], jj[m], integer = is.null(W@x) || all(W@x == 1)))
      curve = sapply(s_grid, function(s) stats::cor(rank((1 - s) * zc_p[m] + s * zm_p[m]), ry))
      s_grid[which.max(curve)]
    })
    mean_weight = stats::median(picks)
  }
  if(verbose){
    cat("Revised clone distance: weight on the clone mean s =", round(mean_weight, 3),
        if(!is.null(picks)) paste0("(per resolution: ", paste(picks, collapse = ", "), ")") else "", "\n")
  }

  out = (1 - mean_weight) * zc + mean_weight * zm
  diag(out) = 0
  dimnames(out) = list(clone_names, clone_names)
  attr(out, "mean_weight") = mean_weight
  attr(out, "picks") = picks
  attr(out, "clonotrace") = D
  return(out)
}

#' @title leiden_ladder
#' @description Leiden Clusters of Cells at Several Resolutions
#' @details
#'
#' Builds a UMAP-style fuzzy k-nearest-neighbour graph over the cell embedding (as scanpy's
#' `pp.neighbors` does) and runs modularity Leiden on it at each resolution.
#'
#' @param embedding A numeric matrix of cell embeddings (cells x dimensions).
#' @param resolutions Numeric vector of Leiden resolutions.
#' @param k Integer. Number of neighbours, counting the cell itself (default: 15).
#' @param seed Integer random seed (default: 0).
#'
#' @return A data frame of cluster labels, one column per resolution, cells as row names.
#'
#' @importFrom uwot similarity_graph
#' @importFrom igraph graph_from_adjacency_matrix cluster_leiden
leiden_ladder = function(embedding, resolutions = c(0.3, 0.5, 0.8, 1.2), k = 15, seed = 0){
  set.seed(seed)
  A = uwot::similarity_graph(as.matrix(embedding), n_neighbors = k)
  g = igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  out = sapply(resolutions, function(r){
    set.seed(seed)
    igraph::cluster_leiden(g, objective_function = "modularity", resolution = r,
                           n_iterations = -1)$membership
  })
  out = as.data.frame(matrix(out, nrow = nrow(embedding)))
  colnames(out) = paste0("leiden_r", resolutions)
  rownames(out) = rownames(embedding)
  return(out)
}

# Square, symmetric clone distance matrix from either input form.
clone_dis_square = function(clone_dis, clone_names, is_long){
  n = length(clone_names)
  if(!is_long){
    D = as.matrix(clone_dis)
    if(nrow(D) != n || ncol(D) != n){
      stop("clone_dis is ", nrow(D), " x ", ncol(D), " but there are ", n, " clone names.")
    }
    if(!is.null(rownames(D))){
      D = D[clone_names, clone_names]
    }
  }
  else{
    cd = as.data.frame(clone_dis)
    D = matrix(0, n, n)
    gi = as.integer(cd$group1)
    gj = as.integer(cd$group2)
    if(max(gi, gj) > n){
      stop("clone_dis indexes clone ", max(gi, gj), " but there are ", n, " clone names.")
    }
    D[cbind(gi, gj)] = as.numeric(cd$dis)
    D[cbind(gj, gi)] = as.numeric(cd$dis)
  }
  storage.mode(D) = "double"
  D = (D + t(D)) / 2
  diag(D) = 0
  if(anyNA(D)){
    stop("The clone distance has missing pairs.")
  }
  dimnames(D) = list(clone_names, clone_names)
  return(D)
}

# cells x clones weights: 0/1 observed membership, or clone probabilities.
clone_cell_weights = function(embedding, clone_labels, cell_clone_prob, clone_names){
  cells = rownames(embedding)
  if(!is.null(clone_labels)){
    if(is.data.frame(clone_labels)){
      clone_labels = stats::setNames(as.character(clone_labels$clone), clone_labels$cell)
    }
    if(is.null(names(clone_labels))){
      stop("clone_labels must be named by cell, or a data frame with columns cell and clone.")
    }
    lab = as.character(clone_labels[cells])
    j = match(lab, clone_names)
    ok = !is.na(j)
    return(Matrix::sparseMatrix(i = which(ok), j = j[ok], x = 1,
                                dims = c(length(cells), length(clone_names))))
  }
  if(is.null(cell_clone_prob)){
    stop("Either clone_labels or cell_clone_prob is needed for the clone means.")
  }
  P = methods::as(methods::as(Matrix::Matrix(cell_clone_prob, sparse = TRUE), "generalMatrix"), "CsparseMatrix")
  colnames(P) = colnames(cell_clone_prob)
  if(is.null(colnames(P))){
    colnames(P) = clone_names
  }
  P = sync_sparse_rows(P, cells)
  return(P[, clone_names, drop = FALSE])
}

# Total variation distance between the rows of a composition matrix, for the pairs (ii, jj).
# With integer counts it is computed as sum|Ka*Nb - Kb*Na| / (2*Na*Nb), so pairs sharing no cluster
# come out as exactly 1 and tie, as they should.
composition_tv = function(C, ii, jj, integer = TRUE, chunk = 200000){
  if(integer){
    C = round(C)
  }
  N = rowSums(C)
  out = numeric(length(ii))
  for(s in seq(1, length(ii), by = chunk)){
    e = min(s + chunk - 1, length(ii))
    a = ii[s:e]
    b = jj[s:e]
    out[s:e] = rowSums(abs(C[a, , drop = FALSE] * N[b] - C[b, , drop = FALSE] * N[a])) / (2 * N[a] * N[b])
  }
  if(!integer){
    out = round(out, 12)
  }
  return(out)
}
