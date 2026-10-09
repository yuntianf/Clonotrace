# leiden_embedding

Run Leiden Clustering from Embedding via Shared Nearest Neighbors

## Usage

``` r
leiden_embedding(
  data,
  k = 30,
  prune.snn = 0,
  weight = "jaccard",
  resolution = 1
)
```

## Arguments

- data:

  A numeric matrix or data frame where rows are observations and columns
  are features.

- k:

  Integer. Number of nearest neighbors for SNN construction (default:
  30).

- prune.snn:

  Numeric. Threshold below which SNN edges are removed (default: 0).

- weight:

  Character. Column to use as edge weight (\`"jaccard"\` or \`"dis"\`)
  (default: \`"jaccard"\`).

- resolution:

  Numeric resolution parameter for clustering (default: 1).

## Value

A factor vector of cluster memberships for each row in \`data\`.

## Details

Constructs a shared nearest neighbor (SNN) graph from an embedding
matrix (e.g., PCA, UMAP), optionally prunes weak connections, and
performs Louvain or Leiden clustering.

## Examples

``` r
mat <- matrix(rnorm(500), nrow = 100)
clusters <- leiden_embedding(mat, k = 20)
table(clusters)
#> clusters
#>  1  2  3  4 
#> 31 23 21 25 
```
