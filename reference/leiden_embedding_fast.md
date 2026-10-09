# leiden_embedding_fast

Fast Leiden Clustering Using Sparse Shared Nearest Neighbors

## Usage

``` r
leiden_embedding_fast(
  data,
  k = 30,
  prune.snn = 0,
  weight = "jaccard",
  resolution = 1
)
```

## Arguments

- data:

  A numeric matrix of embeddings (e.g., PCA, UMAP).

- k:

  Integer. Number of nearest neighbors (default: 30).

- prune.snn:

  Numeric threshold to prune weak SNN edges (default: 0).

- weight:

  Character. Weighting scheme (currently unused; default:
  \`"jaccard"\`).

- resolution:

  Resolution parameter for Leiden clustering (default: 1).

## Value

A factor vector of cluster labels.

## Details

Efficient implementation of SNN graph construction and Leiden
clustering. Uses matrix algebra for fast computation of shared neighbors
and Jaccard similarity.
