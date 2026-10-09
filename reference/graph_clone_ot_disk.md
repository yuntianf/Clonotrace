# Exact Clone OT from a Disk-Backed Graph-Distance Store

Reads bounded cell-distance blocks from a binary store and applies the
same unregularized
[`clone_2_ot()`](https://yuntianf.github.io/Clonotrace/reference/clone_2_ot.md)
calculation to every clone pair.

## Usage

``` r
graph_clone_ot_disk(
  store,
  cell_clone_prob,
  prob_thresh = 0.05,
  target_cell_limit = 5000L,
  cores = 1L,
  checkpoint_dir = NULL,
  overwrite_checkpoint = FALSE,
  monitor = FALSE,
  verbose = TRUE
)
```

## Arguments

- store:

  A completed `clonotrace_distance_store` object or its binary path.

- cell_clone_prob:

  A numeric or sparse cell-by-clone membership matrix with named rows
  and columns.

- prob_thresh:

  Membership values below this threshold are set to zero.

- target_cell_limit:

  Approximate maximum number of unique target cells loaded in one cost
  block.

- cores:

  Number of balanced source-clone partitions. Configure a `future` plan
  before using more than one.

- checkpoint_dir:

  Optional directory for validated per-source checkpoints.

- overwrite_checkpoint:

  Whether to recompute existing checkpoints.

- monitor:

  Whether to attach timing and memory diagnostics.

- verbose:

  Whether to report progress.

## Value

A data frame with `group1`, `group2`, and exact OT `dis`.

## Details

Parallel workers open independent read-only connections. Storage and
data access differ from
[`graph_clone_ot()`](https://yuntianf.github.io/Clonotrace/reference/graph_clone_ot.md),
but the graph-Wasserstein objective and clone probability masses are
unchanged.
