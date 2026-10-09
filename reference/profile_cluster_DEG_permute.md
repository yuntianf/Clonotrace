# profile_cluster_DEG_permute

Permutation Null Distribution for Profile-Based Differential Expression

## Usage

``` r
profile_cluster_DEG_permute(P, G, dpt, n = 50)
```

## Arguments

- P:

  A matrix of soft profile assignments for each cell (cells × 2,
  typically: target vs. control).

- G:

  A gene expression matrix (genes × cells).

- dpt:

  A numeric pseudotime vector (length = number of columns in \`G\`).

- n:

  Integer. Number of permutations (default: 50).

## Value

A list of p-value vectors (one per permutation) representing the null
distribution.

## Details

Computes a null distribution of p-values for profile-based differential
gene expression by repeatedly permuting the soft profile assignment
matrix. This supports empirical calibration of DE gene significance in
pseudotime analysis.
