# find_mutual_nn

Identify Mutual Nearest Neighbors Between Groups

## Usage

``` r
find_mutual_nn(distance, k, dis_thresh = NULL)
```

## Arguments

- distance:

  A numeric symmetric distance matrix.

- k:

  Integer. Number of neighbors to consider per group.

- dis_thresh:

  Optional numeric threshold. If set, only mutual neighbors with
  distance \< threshold are returned.

## Value

A data frame of mutual nearest neighbor pairs with columns:

- group1:

  Index of first group

- group2:

  Index of second group

- dist:

  Distance between groups

## Details

Finds mutually shared k-nearest neighbors from both row and column
perspectives of a distance matrix. Optionally filters by a distance
threshold.

## Examples

``` r
dist_mat <- as.matrix(dist(matrix(rnorm(100), nrow = 10)))
mnn_df <- find_mutual_nn(dist_mat, k = 3, dis_thresh = 0.5)
head(mnn_df)
#>   group1 group2 dist
#> 1      1      1    0
#> 2      2      2    0
#> 3      3      3    0
#> 4      4      4    0
#> 5      5      5    0
#> 6      6      6    0
```
