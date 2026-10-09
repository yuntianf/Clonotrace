# mnn_dist

Identify Mutual Nearest Neighbors (MNN) from Distance Matrix

## Usage

``` r
mnn_dist(dis, k)
```

## Arguments

- dis:

  A symmetric distance matrix (as \`dist\` or matrix) representing
  pairwise distances.

- k:

  Integer. Number of nearest neighbors to consider.

## Value

A data frame with columns \`i\` and \`j\`, representing mutual nearest
neighbor pairs with \`i \< j\`.

## Details

Computes mutual nearest neighbors between items based on a symmetric
distance matrix. A pair (i, j) is considered mutual if i is among the k
nearest neighbors of j, and j is among the k nearest neighbors of i.

## Examples

``` r
dist_mat <- dist(matrix(rnorm(100), ncol = 5))
mnn_pairs <- mnn_dist(dist_mat, k = 5)
head(mnn_pairs)
#>   i  j
#> 1 1  3
#> 2 2 20
#> 3 3 11
#> 4 4 13
#> 5 5  9
#> 6 6  9
```
