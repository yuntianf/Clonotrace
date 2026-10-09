# leiden_dis

Leiden Clustering from Precomputed Distance Matrix

## Usage

``` r
leiden_dis(
  dismat,
  k = 10,
  prune.snn = 0,
  weight = "jaccard",
  resolution = 1,
  if_umap = TRUE
)
```

## Arguments

- dismat:

  A symmetric distance matrix (e.g., from \`dist()\`).

- k:

  Integer. Number of neighbors for SNN construction (default: 10).

- prune.snn:

  Numeric threshold to prune SNN edges (default: 0).

- weight:

  Character. Column to use as edge weight (\`"jaccard"\` or \`"dis"\`)
  (default: \`"jaccard"\`).

- resolution:

  Numeric resolution parameter for Leiden clustering (default: 1).

- if_umap:

  Logical. If \`TRUE\`, returns UMAP coordinates with cluster
  annotations (default: TRUE).

## Value

If \`if_umap = TRUE\`, returns a data frame with UMAP coordinates and
cluster labels. Otherwise, returns a factor vector of cluster
memberships.

## Details

Constructs an SNN graph from a distance matrix, applies Leiden
clustering, and optionally returns a UMAP layout annotated with
clusters.

## Examples

``` r
mat <- matrix(rnorm(500), nrow = 100)
dmat <- dist(mat)
clusters <- leiden_dis(as.matrix(dmat), k = 15, if_umap = FALSE)
```
