# knn_between_groups

Compute k-Nearest Neighbors Between Groups

## Usage

``` r
knn_between_groups(distance, k)
```

## Arguments

- distance:

  A numeric symmetric distance matrix.

- k:

  Integer. Number of nearest neighbors to find for each group.

## Value

A data frame with columns:

- group1:

  Row index of the source group

- group2:

  Index of the neighbor group

- dist:

  Distance between group1 and group2

## Details

For each row in a distance matrix, returns the indices and distances of
the k nearest neighbors (smallest distances). Returns a long-format data
frame of all (i, j) pairs.

## Examples

``` r
dist_mat <- as.matrix(dist(matrix(rnorm(100), nrow = 10)))
knn_df <- knn_between_groups(dist_mat, k = 3)
head(knn_df)
#>   group1 group2     dist
#> 1      1      1 0.000000
#> 2      1      2 2.995119
#> 3      1      6 3.253930
#> 4      2      2 0.000000
#> 5      2      9 2.410759
#> 6      2      4 2.463784
```
