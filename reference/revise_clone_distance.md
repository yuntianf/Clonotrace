# Revise a Clonotrace clone distance

Locally rescales a Clonotrace clone-distance matrix and optionally
blends it with Euclidean distance between clone centroids in the
supplied cell embedding. The centroid weight can be fixed by the user or
selected without fate/cell-type labels by agreement with multiresolution
RNA-cluster composition differences.

## Usage

``` r
revise_clone_distance(
  clone_dis,
  embedding,
  clone_labels = NULL,
  cell_clone_prob = NULL,
  clone_names = NULL,
  sample = NULL,
  local_k = 7L,
  local_rescale = FALSE,
  resolutions = c(0.3, 0.5, 0.8, 1.2),
  clusters = NULL,
  cluster_k = 15L,
  mean_weight = NULL,
  s_grid = seq(0, 1, by = 0.05),
  seed = 0L,
  verbose = TRUE
)
```

## Arguments

- clone_dis:

  Clone distances as either a named square matrix or a long table with
  columns `group1`, `group2`, and `dis`. Long-form clone identifiers may
  be names or one-based positions in `clone_names`.

- embedding:

  Finite numeric cells-by-dimensions embedding with unique cell names as
  row names.

- clone_labels:

  Observed clone membership as a data frame with columns `cell` and
  `clone`, or as a vector named by cell. This is preferred for centroid
  calculation when exact observed membership is available.

- cell_clone_prob:

  Optional non-negative cells-by-clones membership or probability
  matrix. Used only when `clone_labels` is `NULL`.

- clone_names:

  Clone identifiers in distance-matrix order. Inferred from named square
  distances or `cell_clone_prob` when possible.

- sample:

  Optional clone-level grouping, such as patient. Supply a vector named
  by clone or in `clone_names` order. If given, normalization and
  automatic weight selection use only within-group clone pairs.

- local_k:

  Positive integer nearest-neighbor rank used when
  `local_rescale = TRUE`.

- local_rescale:

  Logical. Apply clone-specific local distance rescaling. Defaults to
  `FALSE`; the input Clonotrace geometry is otherwise retained up to a
  single global median normalization.

- resolutions:

  Positive Leiden resolutions used to choose the automatic centroid
  weight.

- clusters:

  Optional precomputed cell-cluster data frame or matrix with cell row
  names and one column per resolution. This replaces internal Leiden
  clustering.

- cluster_k:

  Positive integer number of neighbors in the cell graph used by the
  internal Leiden ladder.

- mean_weight:

  Optional fixed centroid weight in `[0, 1]`. If supplied, automatic
  selection is skipped.

- s_grid:

  Candidate centroid weights in `[0, 1]`.

- seed:

  Integer random seed for the Leiden ladder.

- verbose:

  Logical; report the selected centroid weight.

## Value

A symmetric named clone-distance matrix. Attributes include
`mean_weight`, per-resolution `picks`, the original `clonotrace`
distance, and clone-specific `local_scale`.

## Details

The revision has three steps:

1.  Optional local rescaling: when `local_rescale = TRUE`, use
    \\D\_{ij}/\sqrt{\sigma_i\sigma_j}\\, where \\\sigma_i\\ is the
    distance to clone `i`'s `local_k`-th nearest clone. This is disabled
    by default because it changes clone-neighbour rankings rather than
    merely changing distance units.

2.  Clone-centroid distance: Euclidean distance between observed clone
    means in the cell embedding. If `clone_labels` is unavailable,
    `cell_clone_prob` supplies probability weights.

3.  Median normalization and blending: \\(1-s)D\_{local} +
    sD\_{centroid}\\. When `mean_weight` is `NULL`, `s` is selected at
    each cell-clustering resolution by maximizing Spearman agreement
    with clone cluster-composition total variation, and the median
    selected value is used.

The automatic procedure is fate-blind but RNA-cluster-informed. For
method comparisons, retain the unrevised distance as an ablation and
avoid treating the selected weight as independent of transcriptomic
cluster structure.
