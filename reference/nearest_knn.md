# nearest_knn

Extract Top k Nearest Neighbor Pairs by Distance

## Usage

``` r
nearest_knn(dis, k, top = 3)
```

## Arguments

- dis:

  A symmetric distance matrix (as \`dist\` or matrix).

- k:

  Integer. Number of nearest neighbors to compute per item.

- top:

  Integer. Number of top-ranked unique nearest neighbor pairs to return
  (default: 3).

## Value

A data frame with columns \`i\`, \`j\`, and \`dis\`, listing the top
nearest neighbor pairs and their distances.

## Details

Returns the top closest k-nearest neighbor pairs based on pairwise
distances, excluding duplicates.

## Examples

``` r
dist_mat <- dist(matrix(rnorm(100), ncol = 5))
top_pairs <- nearest_knn(dist_mat, k = 5, top = 10)
head(top_pairs)
#>     i  j       dis
#> 1   4 17 0.7213296
#> 3  12 20 0.9302488
#> 5  10 15 1.0678299
#> 7   4 13 1.1722862
#> 9   2 12 1.2086909
#> 11 14 16 1.2303938
```
