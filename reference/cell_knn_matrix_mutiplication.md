# cell_knn_matrix_mutiplication

Matrix Multiplication for Cell-Feature Weighted Edges

## Usage

``` r
cell_knn_matrix_mutiplication(knn, cell_feature_mat, feature_feature_mat)
```

## Arguments

- knn:

  A matrix or data frame with two columns: \`i\` (source cell) and \`j\`
  (neighbor cell).

- cell_feature_mat:

  A numeric matrix of cell-by-feature values.

- feature_feature_mat:

  A matrix representing relationships between features (e.g., clone
  similarity).

## Value

A numeric vector of scores, one per kNN edge.

## Details

Computes a score for each kNN edge using the dot product between a
cell’s features and a transformed clone matrix. Used to integrate
clone-to-clone and cell-to-clone relationships into edge-level scores.
