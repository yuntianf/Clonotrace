# graph_clone_nn

Compute Clone-to-Clone Nearest Neighbor Distances on a Cell Graph

## Usage

``` r
graph_clone_nn(
  graph,
  cell_clone_prob,
  prob_thresh = 0.1,
  k = 2,
  clone_block_size = 25L,
  verbose = FALSE,
  distance_cache_mb = 64
)
```

## Arguments

- graph:

  An \`igraph\` object representing the cell-cell distance graph.

- cell_clone_prob:

  A numeric matrix of clone membership probabilities (cells × clones).

- prob_thresh:

  Threshold for binarizing clone assignment probabilities (default:
  0.1).

- k:

  Integer. Number of nearest neighbors to select for each clone
  comparison (default: 2).

- clone_block_size:

  Integer. Number of target clones whose cells are included in each
  fallback graph-distance block when the source-clone distance cache
  would exceed \`distance_cache_mb\`. Smaller values reduce peak memory
  at the cost of additional graph-distance calls (default: 25).

- verbose:

  Logical. Whether to print progress during processing (default: FALSE).

- distance_cache_mb:

  Non-negative finite number. Maximum size in MiB of the numeric
  source-clone distance cache per worker (default: 64). Distances to all
  remaining target cells are computed once and reused when they fit;
  otherwise target-clone blocks are used. Set to 0 to disable this
  cache. This is not a total RAM limit: the graph, temporary matrices,
  fallback blocks, and other workers require additional memory.

## Value

A data frame with columns \`group1\`, \`group2\`, and \`dis\` indicating
the NN distance between clones.

## Details

Calculates clone-level nearest neighbor distances based on a shared
cell-cell graph. Uses binarized clone membership to find inter-clone
distances, and applies a greedy k-nearest neighbor approximation.
