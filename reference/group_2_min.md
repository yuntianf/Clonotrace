# group_2_min

Compute Average Minimal Pairwise Distance Between Two Groups

## Usage

``` r
group_2_min(distance, group1, group2, k = 3)
```

## Arguments

- distance:

  A numeric matrix of pairwise distances (e.g., from \`dist()\` or
  \`igraph::distances()\`).

- group1:

  A vector of row indices corresponding to group 1.

- group2:

  A vector of column indices corresponding to group 2.

- k:

  Integer. Number of nearest neighbors to average over (default: 3).

## Value

A single numeric value representing the symmetric average distance
between the two groups.

## Details

Estimates a symmetric inter-group distance based on averaging the top-k
smallest pairwise distances between two groups. This is useful for
comparing cell groups or clone sets using a local-neighborhood
approximation of intergroup proximity. If the submatrix between
\`group1\` and \`group2\` has only one value, that value is returned.
Otherwise, the function returns the mean of the average top-k distances
across both rows and columns.
