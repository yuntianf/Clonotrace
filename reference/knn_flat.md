# knn_flat

Compute Flattened k-Nearest Neighbor Edges

## Usage

``` r
knn_flat(
  x,
  k,
  input = "matrix",
  symmetric = FALSE,
  if_dedup = FALSE,
  if_self = FALSE
)
```

## Arguments

- x:

  A numeric matrix (rows = points) or a distance matrix.

- k:

  Integer. Number of neighbors to compute per point.

- input:

  Character. Either \`"matrix"\` or \`"dist"\` (default: \`"matrix"\`).

- symmetric:

  Logical. If \`TRUE\`, ensures bidirectional edges (default:
  \`FALSE\`).

- if_dedup:

  Logical. If \`TRUE\`, removes redundant (i,j)/(j,i) pairs (default:
  \`FALSE\`).

- if_self:

  Logical. If \`TRUE\`, includes self-edges (i,i) with distance 0
  (default: \`FALSE\`).

## Value

A data frame with columns \`node1\`, \`node2\`, and \`dist\`.

## Details

Computes a flattened list of k-nearest neighbor edges for each point in
a matrix or distance matrix. Returns a long-format data frame of
neighbor pairs and their distances, with optional symmetry, self-edges,
or deduplication.

## Examples

``` r
mat <- matrix(rnorm(100), ncol = 2)
edges <- knn_flat(mat, k = 5, symmetric = TRUE)
head(edges)
#>   node1 node2      dist
#> 1     1     9 0.7116306
#> 2     1    11 1.0831579
#> 3     1    20 0.6729854
#> 4     1    23 1.1161600
#> 5     1    26 0.2518130
#> 6     1    28 0.3707401
```
