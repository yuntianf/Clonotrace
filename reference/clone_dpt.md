# clone_dpt

Compute Clone-Level Diffusion Pseudotime

## Usage

``` r
clone_dpt(
  clone_embedding,
  cell_meta,
  clone_col,
  cluster_col,
  start_cluster,
  k = 10,
  dpt_k = 30
)
```

## Arguments

- clone_embedding:

  A numeric matrix where rows are clones and columns are features.

- cell_meta:

  A data frame of metadata with clone and cluster annotations.

- clone_col:

  Name of the column containing clone identifiers.

- cluster_col:

  Name of the column containing cluster labels.

- start_cluster:

  The cluster label used to select the root clone.

- k:

  Integer. Number of neighbors for kNN graph construction (default: 10).

- dpt_k:

  Integer. Number of eigenvectors for DPT (default: 30).

## Value

A numeric vector of clone-level pseudotime values.

## Details

Estimates pseudotime on the clone embedding by building a kNN graph and
applying diffusion pseudotime from a root clone.

## Examples

``` r
if (FALSE) { # \dontrun{
clone_pt <- clone_dpt(clone_embedding, cell_meta, "clone", "cluster", "Naive")
hist(clone_pt)
} # }
```
