# cell_knn_matrix_multiplication_parallel

Parallelized Matrix Multiplication for kNN Edge Weighting

## Usage

``` r
cell_knn_matrix_multiplication_parallel(
  knn,
  cell_feature_mat,
  feature_feature_mat,
  chunk_size = 5000
)
```

## Arguments

- knn:

  A matrix or data frame of kNN edges (columns: \`i\`, \`j\`).

- cell_feature_mat:

  A numeric matrix of cell-by-feature weights.

- feature_feature_mat:

  A sparse feature similarity matrix (e.g., clone-by-clone).

- chunk_size:

  Integer. Number of edges to process per parallel chunk (default:
  5000).

## Value

A numeric vector of computed scores, one per kNN edge.

## Details

Computes edge-level scores via parallelized matrix multiplication
between cell and feature matrices, allowing scalable construction of
clone-informed kNN graphs.
