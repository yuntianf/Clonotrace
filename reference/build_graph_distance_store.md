# Build a Disk-Backed Exact Cell-Graph Distance Store

Computes exact cell-graph shortest-path distances in resumable row
blocks and writes a row-major, double-precision binary matrix for
disk-backed clone OT.

## Usage

``` r
build_graph_distance_store(
  graph,
  path,
  cells = NULL,
  cell_ids = NULL,
  block_size = 256L,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- graph:

  An undirected `igraph` cell graph with finite, non-negative edge costs
  in its `weight` attribute.

- path:

  Output path ending in `.bin`.

- cells:

  Optional endpoint-cell identifiers or graph indices. Paths may still
  traverse every graph vertex.

- cell_ids:

  Optional unique identifiers in graph-vertex order.

- block_size:

  Number of source rows written per resumable block.

- overwrite:

  Whether to replace an existing store and metadata sidecar.

- verbose:

  Whether to report completed row blocks.

## Value

A validated `clonotrace_distance_store` object.

## Details

The store requires `8 * n_endpoint^2` bytes. It lowers peak memory and
permits reuse, but does not change or approximate the distance
objective.
