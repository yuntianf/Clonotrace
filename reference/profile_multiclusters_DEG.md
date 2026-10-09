# profile_multiclusters_DEG

Profile-Specific Differential Expression Across Multiple Clusters

## Usage

``` r
profile_multiclusters_DEG(
  profile,
  exprs,
  cell_meta,
  cell_profile_prob,
  clusters = NULL,
  cluster_col = "cluster",
  pseudotime_col = "cell_t",
  mass_thresh = 100
)
```

## Arguments

- profile:

  Character. Profile name or index to test.

- exprs:

  A gene expression matrix (genes × cells).

- cell_meta:

  A data frame of cell-level metadata (includes pseudotime and cluster
  labels).

- cell_profile_prob:

  Matrix of soft profile probabilities (cells × profiles).

- clusters:

  Optional. A vector of cluster names to test. If \`NULL\`, all clusters
  are used.

- cluster_col:

  Name of the cluster label column in \`cell_meta\` (default:
  \`"cluster"\`).

- pseudotime_col:

  Name of the pseudotime column in \`cell_meta\` (default:
  \`"cell_t"\`).

- mass_thresh:

  Minimum profile mass in a cluster required to perform DE test
  (default: 100).

## Value

A named list of DE test results, one for each qualifying cluster.

## Details

Applies \`profile_cluster_DEG()\` to a set of clusters to identify genes
associated with a specific profile across the transcriptional landscape.
Skips clusters where the profile mass is below a threshold.
