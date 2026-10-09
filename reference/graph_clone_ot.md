# graph_clone_ot

Parallel Clone-to-Clone OT Computation on Cell Graph

## Usage

``` r
graph_clone_ot(
  graph,
  cell_clone_prob,
  prob_thresh = 0.05,
  cache = 5000,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- graph:

  An \`igraph\` object representing the cell-cell distance graph.

- cell_clone_prob:

  A matrix of clone membership probabilities (cells × clones).

- prob_thresh:

  Float. The minimum weight of a cell in a clone to be preserved,
  otherwise it would be sparsified.

- cache:

  Integer. Maximum number of pooled cells retained in memory during each
  subgraph operation (default: 5000).

- cores:

  Integer. Number of parallel threads to use (default: 1).

- verbose:

  Bool flag to indicate if print the intermediate log for each pair of
  clones.

## Value

A data frame with columns \`group1\`, \`group2\`, and \`dis\`
representing OT distances between clone pairs.

## Details

Computes clone-to-clone optimal transport (OT) distances in parallel by
first partitioning clones and applying \`graph_clone_ot_sub()\` on each
group. Supports memory-efficient computation with subset caching.
