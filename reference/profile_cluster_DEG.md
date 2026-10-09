# profile_cluster_DEG

Differential Expression Between Profile and Background in a Single
Cluster

## Usage

``` r
profile_cluster_DEG(
  profile,
  cluster,
  exprs,
  cell_meta,
  cell_profile_prob,
  cluster_col = "cluster",
  pseudotime_col = "cell_t",
  permute_n = 50
)
```

## Arguments

- profile:

  Character or numeric vector indicating which profiles to test.

- cluster:

  Character or factor. Name or ID of the cluster to analyze.

- exprs:

  A gene expression matrix (genes × cells).

- cell_meta:

  A data frame of cell-level metadata (must include pseudotime and
  cluster labels).

- cell_profile_prob:

  A numeric matrix of profile probabilities (cells × profiles).

- cluster_col:

  Name of the column in \`cell_meta\` for cluster labels (default:
  \`"cluster"\`).

- pseudotime_col:

  Name of the column in \`cell_meta\` representing pseudotime (default:
  \`"cell_t"\`).

- permute_n:

  Integer. Number of permutations for null model generation (default:
  50).

## Value

A list containing:

- stat:

  Data frame with statistics, empirical p-values, and FDR-adjusted
  p-values

- cell:

  Character vector of cell IDs used in the test

- design_null, design_full, coef, df:

  Model components from GAM fitting

## Details

Performs pseudotime-aware differential expression testing between a
target profile and control background within a specific transcriptional
cluster. Combines soft cluster-weighted regression with
permutation-derived null distributions to assess gene-level
significance.
