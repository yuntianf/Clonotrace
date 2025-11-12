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

#' @title label_spreading_blocked
#' @description Memory-efficient label spreading (blocked, HDF5 optional)
#' Runs label spreading on a sparse graph in **class blocks** to control memory,
#' optionally streaming class probabilities to an on-disk HDF5 dataset.
#'
#' @details
#' - Normalizes the adjacency to a row-stochastic matrix.
#' - Splits the label space into \code{block_size}-sized groups to limit RAM use.
#' - Optionally writes each probability block directly into an HDF5 dataset.
#'
#'
#' @param adj \code{dgCMatrix} sparse adjacency; will be row-normalized
#'   in-place to act as a transition matrix.
#' @param labels Integer or numeric vector of length \eqn{N} with class labels in
#'   \code{1..C} and \code{NA} for unlabeled cells. Do not refactor.
#' @param label_n Optional total number of classes; inferred if \code{NULL}.
#' @param alpha Damping parameter in label propagation update.
#' @param max_iter Maximum iterations per block.
#' @param tol Convergence tolerance for updates.
#' @param block_size Number of classes per block to process at once.
#' @param outfile Optional HDF5 file path; if \code{NULL}, returns dense matrix.
#' @param verbose Logical; show progress messages.
#' @param epsilon Small uniform prior mass per class to avoid zeros.
#'
#' @return If \code{outfile = NULL}, a dense \eqn{N \times C} numeric matrix of
#'   class probabilities. Otherwise, invisibly returns the HDF5 file path with
#'   dataset \code{"prob"} written.
#'
#'
#' @seealso \code{\link{create_hdf5_matrix}}, \code{\link{normalize_hdf5_rows}},
#'   \code{\link{label_spreading_bootstrap_blocked}}
#'
#' @importFrom Matrix rowSums Diagonal
#' @importFrom rhdf5 h5createFile h5createDataset h5write
#' @examples
#' \dontrun{
#' library(Matrix)
#' set.seed(1)
#' N <- 1000; C <- 5
#' A <- rsparsematrix(N, N, 0.001); A <- abs(A) + Diagonal(N)
#' labs <- sample(c(1:C, NA), N, TRUE, prob = c(rep(0.05, C), 0.75))
#' P <- label_spreading_blocked(A, labs, label_n = C, block_size = 2, max_iter = 50)
#' }
#' @export
label_spreading_blocked <- function(
    adj, labels, label_n = NULL, alpha = 0.9,
    max_iter = 100, tol = 1e-3, block_size = 128,
    outfile = NULL, verbose = TRUE, epsilon = 0
) {
  # labels are assumed integer in 1..C with NAs (DO NOT re-factor)
  stopifnot(is.integer(labels) || is.numeric(labels))
  labs <- as.integer(labels)
  C <- if (is.null(label_n)) max(labels, na.rm = TRUE) else label_n
  stopifnot(C >= 1L)
  N <- length(labs)

  # weights log2(count)/count computed over present classes in 1..C
  cnt <- tabulate(na.omit(labs), nbins = C)
  wts <- numeric(C); nz <- which(cnt > 0L)
  if (length(nz)) wts[nz] <- log2(cnt[nz]) / cnt[nz]

  labeled_idx <- which(!is.na(labs))

  # Row-normalize adj in-place so it acts as P
  adj <- as(adj, "dgCMatrix")
  deg <- Matrix::rowSums(adj); deg[deg < 1e-12] <- 1e-12
  adj@x <- adj@x / deg[adj@i + 1L]

  # Prepare sink
  if (is.null(outfile)) {
    Prob <- matrix(NA_real_, nrow = N, ncol = C)
  } else {
    if (file.exists(outfile)) file.remove(outfile)
    h5createFile(outfile)
    h5createDataset(outfile, "prob",
                    dims = c(N, C),
                    chunk = c(min(N, 4096L), min(C, block_size)),
                    storage.mode = "double", level = 7)
  }

  cls_seq <- split(seq_len(C), ceiling(seq_len(C) / block_size))

  for (b in seq_along(cls_seq)) {
    cols <- cls_seq[[b]]
    B <- length(cols)

    # Build Y block EXACTLY like your dense Y semantics
    Yb <- matrix(epsilon / C, nrow = N, ncol = B)
    if (length(labeled_idx)) {
      Yb[labeled_idx, ] <- 0
      rows_in_block <- which(!is.na(labs) & labs %in% cols)
      if (length(rows_in_block)) {
        col_pos <- match(labs[rows_in_block], cols)
        Yb[cbind(rows_in_block, col_pos)] <- wts[labs[rows_in_block]]
      }
    }

    F <- matrix(0, nrow = N, ncol = B)
    conv <- FALSE
    for (it in seq_len(max_iter)) {
      F_new <- alpha * (adj %*% F) + (1 - alpha) * Yb
      F_new <- as.matrix(F_new)
      diff <- max(abs(F_new - F))
      if (verbose && (it %% 10L == 0L || it == 1L || diff < tol)) {
        message(sprintf("Block %d/%d, iter %d: diff=%.3e", b, length(cls_seq), it, diff))
      }
      F <- F_new
      if (diff < tol) { conv <- TRUE; break }
    }
    if (verbose && !conv) message(sprintf("Block %d reached max_iter", b))

    if (is.null(outfile)) {
      if (B == 1L) Prob[, cols] <- as.numeric(F) else Prob[, cols] <- F
    } else {
      h5write(F, outfile, "prob", index = list(seq_len(N), cols))
    }
    rm(F, Yb); gc()
  }

  if (is.null(outfile)) return(Prob) else return(invisible(outfile))
}

#' @title create_hdf5_matrix
#' @description Create a chunked HDF5 matrix for streaming label-spreading output
#'
#' Preallocates a compressed HDF5 dataset (default name \code{"prob"})
#' for storing class-probability blocks streamed from
#' \code{\link{label_spreading_blocked}}.
#'
#' @details Removes any existing file before creation. Uses gzip compression
#' (\code{level = 7}) and double precision storage.
#'
#'
#' @param h5file Output HDF5 file path.
#' @param dataset Dataset name to create.
#' @param nrow,ncol Dimensions of the dataset.
#' @param chunk_rows,chunk_cols Chunk sizes for row/column blocks.
#'
#' @return Invisibly returns \code{TRUE} after successful creation.
#'
#'
#' @importFrom rhdf5 h5createFile h5createDataset
#' @examples
#' \dontrun{
#' create_hdf5_matrix("prob_out.h5", "prob", nrow = 5000, ncol = 256)
#' }
create_hdf5_matrix <- function(h5file, dataset = "prob", nrow, ncol,
                               chunk_rows = 4096, chunk_cols = 128) {
  if (file.exists(h5file)) file.remove(h5file)
  h5createFile(h5file)
  h5createDataset(h5file, dataset,
                  dims = c(nrow, ncol),
                  chunk = c(min(nrow, chunk_rows), min(ncol, chunk_cols)),
                  storage.mode = "double", level = 7)
  invisible(TRUE)
}

#' @title normalize_hdf5_rows
#' @description Row-normalize an HDF5 probability matrix in column blocks
#' Iteratively reads an \eqn{N \times C} HDF5 dataset in column blocks,
#' adds a small epsilon if requested, computes row sums, and writes back
#' normalized values so each row sums to 1.
#'
#' @param h5file Path to HDF5 file.
#' @param dataset Dataset name (default \code{"prob"}).
#' @param block_cols Number of columns per I/O block.
#' @param add_eps Logical; if \code{TRUE}, adds \code{1e-12} before summing.
#'
#' @return Numeric vector of length \eqn{N} containing the pre-normalization row sums.
#'
#' @importFrom rhdf5 h5ls h5read h5write
#' @examples
#' \dontrun{
#' normalize_hdf5_rows("prob_out.h5", "prob", block_cols = 128)
#' }
normalize_hdf5_rows <- function(h5file, dataset = "prob", block_cols = 256, add_eps = TRUE) {
  info <- h5ls(h5file); ds <- subset(info, name == dataset)
  dims <- as.integer(strsplit(ds$dim, " x ")[[1]])
  N <- dims[1]; C <- dims[2]
  rowsum <- numeric(N)
  col_blocks <- split(seq_len(C), ceiling(seq_len(C)/block_cols))
  for (cols in col_blocks) {
    X <- h5read(h5file, dataset, index = list(seq_len(N), cols))
    if (add_eps) X <- X + 1e-12  # mirror your original tweak
    rowsum <- rowsum + rowSums(X)
    h5write(X, h5file, dataset, index = list(seq_len(N), cols))
    rm(X); gc()
  }
  rowsum[rowsum == 0] <- 1
  for (cols in col_blocks) {
    X <- h5read(h5file, dataset, index = list(seq_len(N), cols))
    X <- X / rowsum
    h5write(X, h5file, dataset, index = list(seq_len(N), cols))
    rm(X); gc()
  }
  invisible(rowsum)
}

#' @title label_spreading_bootstrap_blocked
#' @description Bootstrap label spreading with HDF5 reference and per-cell deviance
#'
#' Performs multiple label-spreading runs on random subsets of labeled cells,
#' compares each to a full reference solution, and returns the reference
#' probabilities plus mean L1 deviance per cell.
#'
#' @details
#' - Computes an exact reference with all labeled cells, normalizes it, and stores it in HDF5.
#' - Each bootstrap uses a random subset of labeled cells, writes its results to a temporary HDF5,
#'   and accumulates the L1 deviation vs the reference in streaming mode.
#' - Parallelizable via \pkg{future.apply}: set a plan before calling.
#'
#' @param adj Sparse adjacency (\code{dgCMatrix}); will be row-normalized.
#' @param labels Vector of original labels (length \eqn{N}); mapped once to \code{1..C}.
#' @param alpha Damping parameter for label spreading.
#' @param sample_rate Fraction of labeled cells retained in each bootstrap.
#' @param sample_n Number of bootstrap replicates.
#' @param block_size Classes per block when streaming to HDF5.
#' @param tol,max_iter Convergence parameters.
#' @param epsilon Uniform prior mass added to avoid zeros.
#' @param refer_h5 Optional existing reference HDF5 file; if \code{NULL}, computed internally.
#' @param tmpdir Directory for temporary HDF5 files.
#' @param verbose Logical; print progress.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{prob}}{Dense reference \eqn{N \times C} probability matrix.}
#'   \item{\code{deviance}}{Per-cell mean L1 deviation vs reference.}
#' }
#'
#'
#' @seealso \code{\link{label_spreading_blocked}}, \code{\link{normalize_hdf5_rows}}
#'
#' @importFrom Matrix rowSums Diagonal
#' @importFrom rhdf5 h5read h5write h5createFile h5createDataset
#' @importFrom future.apply future_lapply
#'
#' @examples
#' \dontrun{
#' library(Matrix); library(future.apply)
#' plan(multisession, workers = 4)
#' N <- 1000; C <- 4
#' A <- rsparsematrix(N, N, 0.001); A <- abs(A) + Diagonal(N)
#' labs <- sample(c(1:C, NA), N, TRUE, prob = c(rep(0.05, C), 0.8))
#' res <- label_spreading_bootstrap_blocked(A, labs, sample_rate = 0.8, sample_n = 10)
#' str(res)
#' }
#' @export
label_spreading_bootstrap_blocked <- function(
    adj, labels,                    # labels are preindexed 1..C or will be preindexed below
    alpha = 0.8,
    sample_rate = 0.8,
    sample_n = 50,
    block_size = 128,
    tol = 1e-3,
    max_iter = 100,
    epsilon = 0,
    refer_h5 = NULL,
    tmpdir = tempdir(),
    verbose = TRUE
) {
  # Preindex ONCE from the full labels; keep this mapping fixed
  labs_full <- as.integer(as.factor(labels))
  N <- length(labs_full)
  C <- max(labs_full, na.rm = TRUE)

  # Build or reuse reference (normalized, aligned to 1..C)
  if (is.null(refer_h5)) {
    refer_h5 <- file.path(tmpdir, sprintf("refer_%s.h5", paste0(sample(letters, 6), collapse = "")))
    if (verbose) message("Computing reference (exact, preindexed), then global normalize…")
    label_spreading_blocked(
      adj, labs_full, label_n = C, alpha = alpha,
      max_iter = max_iter, tol = tol,
      block_size = block_size, outfile = refer_h5,
      verbose = FALSE, epsilon = epsilon
    )
    normalize_hdf5_rows(refer_h5, "prob", block_cols = block_size, add_eps = TRUE)
  } else if (verbose) {
    message(sprintf("Using provided reference: %s", refer_h5))
  }

  labeled_idx_all <- which(!is.na(labs_full))
  full_flag <- labs_full;  full_flag[is.na(full_flag)] <- 0L

  worker_fun <- function(i) {
    subs <- rep(NA_integer_, N)
    if (length(labeled_idx_all)) {
      take <- max(1L, round(length(labeled_idx_all) * sample_rate))
      idx  <- sample(labeled_idx_all, size = take, replace = FALSE)
      subs[idx] <- labs_full[idx]     # IMPORTANT: keep original indices (no refactor)
    }

    samp_flag <- subs; samp_flag[is.na(samp_flag)] <- 0L
    same_flag <- as.integer(!xor(full_flag, samp_flag))

    tmp_h5 <- file.path(tmpdir, sprintf("boot_%04d_%s.h5", i, paste0(sample(letters, 4), collapse = "")))
    create_hdf5_matrix(tmp_h5, "prob", nrow = N, ncol = C, chunk_rows = 4096, chunk_cols = block_size)

    # Run exact solver with preindexed labels (no refactor)
    label_spreading_blocked(
      adj, subs, label_n = C, alpha = alpha,
      max_iter = max_iter, tol = tol,
      block_size = block_size, outfile = tmp_h5,
      verbose = FALSE, epsilon = epsilon
    )

    # Normalize with +1e-12 like your original
    normalize_hdf5_rows(tmp_h5, "prob", block_cols = block_size, add_eps = TRUE)

    # Stream L1 deviance vs reference
    dev <- numeric(N)
    cls_seq <- split(seq_len(C), ceiling(seq_len(C)/block_size))
    for (cols in cls_seq) {
      Fb <- rhdf5::h5read(tmp_h5, "prob", index = list(seq_len(N), cols))
      Rb <- rhdf5::h5read(refer_h5, "prob", index = list(seq_len(N), cols))
      dev <- dev + rowSums(abs(Fb - Rb))
      rm(Fb, Rb); gc()
    }

    try(unlink(tmp_h5), silent = TRUE)
    dev <- dev * same_flag
    list(dev = dev, flag = same_flag)
  }

  # Run in parallel (set future::plan() before calling)
  res_list <- future.apply::future_lapply(seq_len(sample_n), worker_fun, future.seed = TRUE)

  dev_sum  <- Reduce(`+`, lapply(res_list, `[[`, "dev"))
  flag_sum <- Reduce(`+`, lapply(res_list, `[[`, "flag"))
  deviance_norm <- dev_sum / pmax(1, flag_sum)

  prob = h5read(refer_h5,"prob")

  try(unlink(refer_h5), silent = TRUE)

  list(prob = prob, deviance = deviance_norm)
}
