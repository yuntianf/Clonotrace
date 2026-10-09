#' Revise a Clonotrace clone distance
#'
#' Locally rescales a Clonotrace clone-distance matrix and optionally blends it
#' with Euclidean distance between clone centroids in the supplied cell
#' embedding. The centroid weight can be fixed by the user or selected without
#' fate/cell-type labels by agreement with multiresolution RNA-cluster
#' composition differences.
#'
#' @details
#' The revision has three steps:
#' \enumerate{
#'   \item Optional local rescaling: when \code{local_rescale = TRUE}, use
#'     \eqn{D_{ij}/\sqrt{\sigma_i\sigma_j}}, where \eqn{\sigma_i} is the
#'     distance to clone \code{i}'s \code{local_k}-th nearest clone. This is
#'     disabled by default because it changes clone-neighbour rankings rather
#'     than merely changing distance units.
#'   \item Clone-centroid distance: Euclidean distance between observed clone
#'     means in the cell embedding. If \code{clone_labels} is unavailable,
#'     \code{cell_clone_prob} supplies probability weights.
#'   \item Median normalization and blending:
#'     \eqn{(1-s)D_{local} + sD_{centroid}}. When \code{mean_weight} is
#'     \code{NULL}, \code{s} is selected at each cell-clustering resolution by
#'     maximizing Spearman agreement with clone cluster-composition total
#'     variation, and the median selected value is used.
#' }
#'
#' The automatic procedure is fate-blind but RNA-cluster-informed. For method
#' comparisons, retain the unrevised distance as an ablation and avoid treating
#' the selected weight as independent of transcriptomic cluster structure.
#'
#' @param clone_dis Clone distances as either a named square matrix or a long
#'   table with columns \code{group1}, \code{group2}, and \code{dis}. Long-form
#'   clone identifiers may be names or one-based positions in
#'   \code{clone_names}.
#' @param embedding Finite numeric cells-by-dimensions embedding with unique
#'   cell names as row names.
#' @param clone_labels Observed clone membership as a data frame with columns
#'   \code{cell} and \code{clone}, or as a vector named by cell. This is
#'   preferred for centroid calculation when exact observed membership is
#'   available.
#' @param cell_clone_prob Optional non-negative cells-by-clones membership or
#'   probability matrix. Used only when \code{clone_labels} is \code{NULL}.
#' @param clone_names Clone identifiers in distance-matrix order. Inferred from
#'   named square distances or \code{cell_clone_prob} when possible.
#' @param sample Optional clone-level grouping, such as patient. Supply a vector
#'   named by clone or in \code{clone_names} order. If given, normalization and
#'   automatic weight selection use only within-group clone pairs.
#' @param local_k Positive integer nearest-neighbor rank used when
#'   \code{local_rescale = TRUE}.
#' @param local_rescale Logical. Apply clone-specific local distance rescaling.
#'   Defaults to \code{FALSE}; the input Clonotrace geometry is otherwise
#'   retained up to a single global median normalization.
#' @param resolutions Positive Leiden resolutions used to choose the automatic
#'   centroid weight.
#' @param clusters Optional precomputed cell-cluster data frame or matrix with
#'   cell row names and one column per resolution. This replaces internal
#'   Leiden clustering.
#' @param cluster_k Positive integer number of neighbors in the cell graph used
#'   by the internal Leiden ladder.
#' @param mean_weight Optional fixed centroid weight in \code{[0, 1]}. If
#'   supplied, automatic selection is skipped.
#' @param s_grid Candidate centroid weights in \code{[0, 1]}.
#' @param seed Integer random seed for the Leiden ladder.
#' @param verbose Logical; report the selected centroid weight.
#'
#' @return A symmetric named clone-distance matrix. Attributes include
#'   \code{mean_weight}, per-resolution \code{picks}, the original
#'   \code{clonotrace} distance, and clone-specific \code{local_scale}.
#'
#' @export
revise_clone_distance <- function(
    clone_dis, embedding, clone_labels = NULL, cell_clone_prob = NULL,
    clone_names = NULL, sample = NULL, local_k = 7L,
    local_rescale = FALSE,
    resolutions = c(0.3, 0.5, 0.8, 1.2), clusters = NULL,
    cluster_k = 15L, mean_weight = NULL,
    s_grid = seq(0, 1, by = 0.05), seed = 0L, verbose = TRUE
) {
  is_long <- !is.null(colnames(clone_dis)) &&
    all(c("group1", "group2", "dis") %in% colnames(clone_dis))
  if (is.null(clone_names)) {
    clone_names <- if (!is_long) colnames(clone_dis) else colnames(cell_clone_prob)
  }
  if (is.null(clone_names)) {
    stop(
      "`clone_names` is required when neither `clone_dis` nor ",
      "`cell_clone_prob` carries clone names."
    )
  }
  clone_names <- as.character(clone_names)
  if (anyNA(clone_names) || any(clone_names == "") || anyDuplicated(clone_names)) {
    stop("`clone_names` must contain unique, non-missing, non-empty identifiers.")
  }

  distance <- .clonotrace_clone_distance_square(clone_dis, clone_names, is_long)
  n_clones <- nrow(distance)
  if (!is.logical(local_rescale) || length(local_rescale) != 1L || is.na(local_rescale)) {
    stop("`local_rescale` must be one non-missing logical value.")
  }
  if (!is.numeric(local_k) || length(local_k) != 1L || is.na(local_k) || local_k < 1L) {
    stop("`local_k` must be one positive integer.")
  }
  local_k <- as.integer(local_k)
  if (isTRUE(local_rescale) && n_clones <= local_k) {
    stop(
      "There are ", n_clones, " clones; local rescaling requires more than ",
      "`local_k = ", local_k, "` clones."
    )
  }

  embedding <- as.matrix(embedding)
  if (!is.numeric(embedding) || any(!is.finite(embedding))) {
    stop("`embedding` must be a finite numeric matrix.")
  }
  if (is.null(rownames(embedding)) || anyDuplicated(rownames(embedding))) {
    stop("`embedding` must have unique cell names as row names.")
  }

  weights <- .clonotrace_clone_cell_weights(
    embedding, clone_labels, cell_clone_prob, clone_names
  )
  mass <- Matrix::colSums(weights)
  if (any(!is.finite(mass)) || any(mass <= 0)) {
    bad <- clone_names[!is.finite(mass) | mass <= 0]
    stop(
      "No positive cell mass was found for clone(s): ",
      paste(utils::head(bad, 5L), collapse = ", "), "."
    )
  }
  clone_centroids <- as.matrix(Matrix::t(weights) %*% embedding) / mass
  centroid_distance <- sqrt(pmax(
    outer(rowSums(clone_centroids^2), rowSums(clone_centroids^2), "+") -
      2 * tcrossprod(clone_centroids),
    0
  ))

  pairs <- which(upper.tri(distance), arr.ind = TRUE)
  if (!is.null(sample)) {
    if (!is.null(names(sample))) sample <- sample[clone_names]
    if (length(sample) != n_clones || anyNA(sample)) {
      stop("`sample` must give one non-missing group for every clone.")
    }
    pairs <- pairs[sample[pairs[, 1L]] == sample[pairs[, 2L]], , drop = FALSE]
  }
  if (nrow(pairs) == 0L) {
    stop("No clone pairs are available for distance normalization.")
  }
  pair_i <- pairs[, 1L]
  pair_j <- pairs[, 2L]

  if (isTRUE(local_rescale)) {
    local_scale <- vapply(seq_len(n_clones), function(i) {
      values <- distance[i, -i]
      sort(values, partial = local_k)[local_k]
    }, numeric(1))
    local_scale <- pmax(local_scale, .Machine$double.eps)
    local_distance <- distance / sqrt(outer(local_scale, local_scale))
  } else {
    local_scale <- rep(1, n_clones)
    local_distance <- distance
  }

  local_median <- stats::median(local_distance[pairs])
  centroid_median <- stats::median(centroid_distance[pairs])
  if (!is.finite(local_median) || local_median <= 0 ||
      !is.finite(centroid_median) || centroid_median <= 0) {
    stop("Distance normalization medians must be finite and positive.")
  }
  local_normalized <- local_distance / local_median
  centroid_normalized <- centroid_distance / centroid_median

  if (!is.null(mean_weight)) {
    if (!is.numeric(mean_weight) || length(mean_weight) != 1L ||
        is.na(mean_weight) || mean_weight < 0 || mean_weight > 1) {
      stop("`mean_weight` must be one number in [0, 1].")
    }
    picks <- NULL
    selection_correlations <- NULL
  } else {
    if (!is.numeric(s_grid) || !length(s_grid) || anyNA(s_grid) ||
        any(s_grid < 0 | s_grid > 1)) {
      stop("`s_grid` must contain finite candidate weights in [0, 1].")
    }
    s_grid <- sort(unique(as.numeric(s_grid)))
    if (is.null(clusters)) {
      clusters <- .clonotrace_leiden_ladder(
        embedding, resolutions = resolutions, k = cluster_k, seed = seed
      )
    } else {
      clusters <- as.data.frame(clusters)
      if (is.null(rownames(clusters)) || anyDuplicated(rownames(clusters))) {
        stop("`clusters` must have unique cell names as row names.")
      }
      missing_cells <- setdiff(rownames(embedding), rownames(clusters))
      if (length(missing_cells)) {
        stop("`clusters` is missing cells present in `embedding`.")
      }
      clusters <- clusters[rownames(embedding), , drop = FALSE]
    }
    if (!ncol(clusters)) stop("`clusters` must contain at least one resolution.")
    if (is.null(colnames(clusters))) {
      colnames(clusters) <- paste0("resolution_", seq_len(ncol(clusters)))
    }

    local_pairs <- local_normalized[pairs]
    centroid_pairs <- centroid_normalized[pairs]
    selection_correlations <- matrix(
      NA_real_, nrow = length(s_grid), ncol = ncol(clusters),
      dimnames = list(as.character(s_grid), colnames(clusters))
    )
    picks <- vapply(seq_len(ncol(clusters)), function(column) {
      labels <- clusters[[column]]
      observed <- !is.na(labels)
      if (length(unique(labels[observed])) < 2L) {
        warning(
          "Cell clustering `", colnames(clusters)[column],
          "` has fewer than two observed groups; selecting weight 0 for it.",
          call. = FALSE
        )
        return(0)
      }
      cluster_index <- as.integer(factor(labels[observed]))
      indicator <- Matrix::sparseMatrix(
        i = which(observed), j = cluster_index, x = 1,
        dims = c(nrow(embedding), max(cluster_index))
      )
      composition <- as.matrix(Matrix::t(weights) %*% indicator)
      keep <- rowSums(composition) > 0
      eligible <- keep[pair_i] & keep[pair_j]
      if (sum(eligible) < 3L) {
        warning(
          "Too few eligible clone pairs at `", colnames(clusters)[column],
          "`; selecting weight 0 for it.", call. = FALSE
        )
        return(0)
      }
      integer_weights <- !length(weights@x) ||
        all(abs(weights@x - round(weights@x)) < 1e-12)
      composition_difference <- .clonotrace_composition_tv(
        composition, pair_i[eligible], pair_j[eligible],
        integer = integer_weights
      )
      target_rank <- rank(composition_difference, ties.method = "average")
      curve <- vapply(s_grid, function(weight) {
        blended <- (1 - weight) * local_pairs[eligible] +
          weight * centroid_pairs[eligible]
        suppressWarnings(stats::cor(
          rank(blended, ties.method = "average"), target_rank,
          method = "pearson", use = "complete.obs"
        ))
      }, numeric(1))
      selection_correlations[, column] <<- curve
      finite <- which(is.finite(curve))
      if (!length(finite)) return(0)
      s_grid[finite[which.max(curve[finite])]]
    }, numeric(1))
    names(picks) <- colnames(clusters)
    mean_weight <- stats::median(picks)
  }

  if (verbose) {
    detail <- if (is.null(picks)) "" else paste0(
      " (per resolution: ", paste(round(picks, 3), collapse = ", "), ")"
    )
    message(
      "Revised clone distance: local rescaling = ", isTRUE(local_rescale),
      "; centroid weight = ", round(mean_weight, 3), detail
    )
  }

  revised <- (1 - mean_weight) * local_normalized +
    mean_weight * centroid_normalized
  revised <- (revised + t(revised)) / 2
  diag(revised) <- 0
  dimnames(revised) <- list(clone_names, clone_names)
  attr(revised, "mean_weight") <- mean_weight
  attr(revised, "picks") <- picks
  attr(revised, "selection_correlations") <- selection_correlations
  attr(revised, "clonotrace") <- distance
  attr(revised, "local_scale") <- stats::setNames(local_scale, clone_names)
  attr(revised, "local_rescale") <- isTRUE(local_rescale)
  revised
}


.clonotrace_clone_distance_square <- function(clone_dis, clone_names, is_long) {
  n_clones <- length(clone_names)
  if (!is_long) {
    distance <- as.matrix(clone_dis)
    if (!is.numeric(distance) || nrow(distance) != n_clones ||
        ncol(distance) != n_clones) {
      stop("`clone_dis` must be a numeric square matrix matching `clone_names`.")
    }
    if (!is.null(rownames(distance)) || !is.null(colnames(distance))) {
      if (is.null(rownames(distance)) || is.null(colnames(distance)) ||
          anyDuplicated(rownames(distance)) || anyDuplicated(colnames(distance)) ||
          !all(clone_names %in% rownames(distance)) ||
          !all(clone_names %in% colnames(distance))) {
        stop("Named square `clone_dis` must contain every clone name exactly once.")
      }
      distance <- distance[clone_names, clone_names, drop = FALSE]
    } else {
      dimnames(distance) <- list(clone_names, clone_names)
    }
    if (!isTRUE(all.equal(distance, t(distance), tolerance = 1e-8, check.attributes = FALSE))) {
      stop("Square `clone_dis` must be symmetric.")
    }
  } else {
    long <- as.data.frame(clone_dis, stringsAsFactors = FALSE)
    normalize_ids <- function(ids) {
      values <- as.character(ids)
      matched <- match(values, clone_names)
      missing <- is.na(matched)
      if (any(missing)) {
        positions <- suppressWarnings(as.integer(values[missing]))
        valid <- !is.na(positions) & positions >= 1L & positions <= n_clones
        matched[which(missing)[valid]] <- positions[valid]
      }
      if (anyNA(matched)) {
        stop("Long-form `clone_dis` contains clone identifiers absent from `clone_names`.")
      }
      matched
    }
    first <- normalize_ids(long$group1)
    second <- normalize_ids(long$group2)
    values <- suppressWarnings(as.numeric(as.character(long$dis)))
    if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
      stop("Long-form clone distances must be finite, non-negative numbers.")
    }
    non_diagonal <- first != second
    first <- first[non_diagonal]
    second <- second[non_diagonal]
    values <- values[non_diagonal]
    low <- pmin(first, second)
    high <- pmax(first, second)
    key <- paste(low, high, sep = "\r")
    if (anyDuplicated(key)) {
      duplicates <- unique(key[duplicated(key) | duplicated(key, fromLast = TRUE)])
      conflicting <- vapply(duplicates, function(current) {
        length(unique(values[key == current])) > 1L
      }, logical(1))
      if (any(conflicting)) stop("Conflicting distances were supplied for a clone pair.")
      keep <- !duplicated(key)
      low <- low[keep]
      high <- high[keep]
      values <- values[keep]
    }
    distance <- matrix(
      NA_real_, nrow = n_clones, ncol = n_clones,
      dimnames = list(clone_names, clone_names)
    )
    distance[cbind(low, high)] <- values
    distance[cbind(high, low)] <- values
    diag(distance) <- 0
    if (anyNA(distance)) {
      stop("Long-form `clone_dis` is missing one or more clone pairs.")
    }
  }
  storage.mode(distance) <- "double"
  if (anyNA(distance) || any(!is.finite(distance)) || any(distance < 0)) {
    stop("Clone distances must be finite and non-negative.")
  }
  diag(distance) <- 0
  dimnames(distance) <- list(clone_names, clone_names)
  distance
}


.clonotrace_clone_cell_weights <- function(
    embedding, clone_labels, cell_clone_prob, clone_names
) {
  cells <- rownames(embedding)
  if (!is.null(clone_labels)) {
    if (is.data.frame(clone_labels)) {
      if (!all(c("cell", "clone") %in% colnames(clone_labels))) {
        stop("Data-frame `clone_labels` must contain `cell` and `clone` columns.")
      }
      if (anyDuplicated(clone_labels$cell)) {
        stop("Data-frame `clone_labels` must contain at most one row per cell.")
      }
      clone_labels <- stats::setNames(
        as.character(clone_labels$clone), as.character(clone_labels$cell)
      )
    }
    if (is.null(names(clone_labels)) || anyDuplicated(names(clone_labels))) {
      stop("Vector `clone_labels` must have unique cell names.")
    }
    labels <- as.character(clone_labels[cells])
    clone_index <- match(labels, clone_names)
    observed <- !is.na(clone_index)
    return(Matrix::sparseMatrix(
      i = which(observed), j = clone_index[observed], x = 1,
      dims = c(length(cells), length(clone_names)),
      dimnames = list(cells, clone_names)
    ))
  }
  if (is.null(cell_clone_prob)) {
    stop("Supply either `clone_labels` or `cell_clone_prob`.")
  }
  if (is.null(rownames(cell_clone_prob)) || is.null(colnames(cell_clone_prob)) ||
      anyDuplicated(rownames(cell_clone_prob)) ||
      anyDuplicated(colnames(cell_clone_prob))) {
    stop("`cell_clone_prob` must have unique cell and clone names.")
  }
  missing_cells <- setdiff(cells, rownames(cell_clone_prob))
  missing_clones <- setdiff(clone_names, colnames(cell_clone_prob))
  if (length(missing_cells) || length(missing_clones)) {
    stop("`cell_clone_prob` is missing required cells or clone columns.")
  }
  weights <- methods::as(
    methods::as(Matrix::Matrix(
      cell_clone_prob[cells, clone_names, drop = FALSE], sparse = TRUE
    ), "generalMatrix"),
    "CsparseMatrix"
  )
  if (any(!is.finite(weights@x)) || any(weights@x < 0)) {
    stop("`cell_clone_prob` must contain finite, non-negative values.")
  }
  weights
}


.clonotrace_composition_tv <- function(
    composition, first, second, integer = TRUE, chunk = 200000L
) {
  if (integer) composition <- round(composition)
  mass <- rowSums(composition)
  result <- numeric(length(first))
  for (start in seq.int(1L, length(first), by = chunk)) {
    end <- min(start + chunk - 1L, length(first))
    a <- first[start:end]
    b <- second[start:end]
    result[start:end] <- rowSums(abs(
      composition[a, , drop = FALSE] * mass[b] -
        composition[b, , drop = FALSE] * mass[a]
    )) / (2 * mass[a] * mass[b])
  }
  if (!integer) result <- round(result, 12)
  result
}


.clonotrace_leiden_ladder <- function(
    embedding, resolutions = c(0.3, 0.5, 0.8, 1.2), k = 15L, seed = 0L
) {
  if (!is.numeric(resolutions) || !length(resolutions) || anyNA(resolutions) ||
      any(resolutions <= 0)) {
    stop("`resolutions` must contain positive finite values.")
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2L ||
      k >= nrow(embedding)) {
    stop("`cluster_k` must be at least 2 and smaller than the number of cells.")
  }
  set.seed(seed)
  adjacency <- uwot::similarity_graph(as.matrix(embedding), n_neighbors = as.integer(k))
  graph <- igraph::graph_from_adjacency_matrix(
    adjacency, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  membership <- vapply(resolutions, function(resolution) {
    set.seed(seed)
    as.integer(igraph::cluster_leiden(
      graph, objective_function = "modularity",
      resolution = resolution, n_iterations = -1L
    )$membership)
  }, integer(nrow(embedding)))
  membership <- as.data.frame(membership)
  colnames(membership) <- paste0("leiden_r", resolutions)
  rownames(membership) <- rownames(embedding)
  membership
}
