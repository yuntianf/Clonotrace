# clone_disance

Compute Clone-to-Clone Distances from Cell Embedding and Clone
Assignments

## Usage

``` r
clone_disance(
  embedding,
  cell_clone_prob,
  outpath,
  graph_k = 10,
  overwrite = FALSE,
  exact = FALSE,
  force_exact = FALSE,
  exact_cell_limit = 20000L,
  exact_clone_limit = 500L,
  exact_backend = c("disk", "memory"),
  revised = TRUE,
  clone_labels = NULL,
  sample = NULL,
  ...
)
```

## Arguments

- embedding:

  A numeric matrix of cell embeddings (cells × dimensions).

- cell_clone_prob:

  A numeric matrix of clone membership probabilities (cells × clones).

- outpath:

  File path to save and load intermediate graph and distance results.

- graph_k:

  Integer. Number of neighbors to use in the SNN graph (default: 10).

- overwrite:

  Logical. Whether to overwrite previously saved graphs and results
  (default: FALSE).

- exact:

  Logical. If \`TRUE\`, computes exact OT-based clone distances with the
  disk-backed implementation by default; otherwise, uses approximate NN
  (default: FALSE).

- exact_backend:

  Exact-OT implementation. `"disk"` (default) writes a resumable
  cell-distance store and reads bounded cost blocks. `"memory"` retains
  the legacy dense in-memory implementation for small reproductions.

- force_exact:

  Logical. Allow exact mode to proceed when the input exceeds the
  conservative cell or clone safety limits (default: \`FALSE\`).

- exact_cell_limit, exact_clone_limit:

  Conservative safety limits used to warn before an exact calculation
  that is likely to require substantial runtime or memory.

- revised:

  Logical. If \`TRUE\` (the default), optionally locally rescales the
  raw Clonotrace distance and blends it with clone-centroid distance
  using a fate-blind automatically selected weight; see
  \`revise_clone_distance()\`. Local rescaling is off by default and can
  be enabled with \`local_rescale = TRUE\`. Set \`FALSE\` to return the
  original exact-OT or approximate distance.

- clone_labels:

  Optional observed clone membership used for clone centroids, as a data
  frame with columns \`cell\` and \`clone\`, or a vector named by cell.
  If \`NULL\`, \`cell_clone_prob\` supplies probability weights.

- sample:

  Optional clone-level grouping, such as patient, passed to
  \`revise_clone_distance()\` so normalization and weight selection use
  within-group clone pairs.

- ...:

  Named arguments passed to
  [`graph_clone_ot_disk()`](https://yuntianf.github.io/Clonotrace/reference/graph_clone_ot_disk.md)
  (disk exact),
  [`graph_clone_ot()`](https://yuntianf.github.io/Clonotrace/reference/graph_clone_ot.md)
  (legacy memory exact), or
  [`graph_clone_nn()`](https://yuntianf.github.io/Clonotrace/reference/graph_clone_nn.md),
  and—when `revised = TRUE`—to \`revise_clone_distance()\` (\`local_k\`,
  \`local_rescale\`, \`resolutions\`, \`clusters\`, \`cluster_k\`,
  \`mean_weight\`, \`s_grid\`, \`seed\`, and \`verbose\`).

## Value

A data frame with columns \`group1\`, \`group2\`, and \`dis\`. With
\`revised = TRUE\`, \`dis\` is the revised distance, \`dis_clonotrace\`
retains the raw distance, and attributes \`mean_weight\` and \`picks\`
record the automatic blend. With \`revised = FALSE\`, \`dis\` is the raw
distance.

## Details

Constructs a cell-cell shared nearest neighbor (SNN) graph from an
embedding and uses it to estimate pairwise distances between clones.
Exact mode uses a resumable disk-backed shortest-path store by default;
the legacy dense in-memory backend remains available explicitly. Users
can alternatively select the faster approximate k-nearest-neighbor
distance.
