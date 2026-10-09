# graph_clone_ot_sub

Compute Pairwise OT Distances Between Clones Using Graph-Based Subset
Strategy

## Usage

``` r
graph_clone_ot_sub(
  graph,
  cell_clone_prob,
  target_clone = NULL,
  cache = 5000,
  verbose = TRUE
)
```

## Arguments

- graph:

  An \`igraph\` object representing cell-cell distances.

- cell_clone_prob:

  A matrix of clone membership probabilities (cells × clones).

- target_clone:

  Optional. Vector of clone indices to compute pairwise OT from
  (default: \`NULL\` = all).

- cache:

  Integer. Maximum number of pooled cells retained in memory during
  iteration (default: 5000).

- verbose:

  Bool flag to indicate if print the intermediate log for each pair of
  clones.

## Value

A matrix with columns: \`group1\`, \`group2\`, \`dis\` for each clone
pair and their OT distance.

## Details

Computes clone-to-clone optimal transport distances by expanding a local
graph-based cell neighborhood using a subset pool strategy to reduce
computational burden.
