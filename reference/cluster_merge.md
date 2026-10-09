# cluster_merge

Merge Elements from a List According to Cluster Assignments

## Usage

``` r
cluster_merge(input_list, cluster)
```

## Arguments

- input_list:

  A list of vectors or objects to be grouped.

- cluster:

  A vector of cluster labels of the same length as \`input_list\`.

## Value

A list where each element contains the merged values from the same
cluster.

## Details

Given a list and a vector of cluster assignments, merges list elements
that belong to the same cluster into grouped lists.

## Examples

``` r
lst <- list(1:2, 3:4, 5, 6)
clust <- c(1, 1, 2, 2)
merged <- cluster_merge(lst, clust)
print(merged)
#> [[1]]
#> [1] 1 2 3 4
#> 
#> [[2]]
#> [1] 5 6
#> 
```
