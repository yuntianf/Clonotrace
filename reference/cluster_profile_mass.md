# cluster_profile_mass

Aggregate Profile Probabilities by Cluster

## Usage

``` r
cluster_profile_mass(cell_profile_prob, cluster_label)
```

## Arguments

- cell_profile_prob:

  A numeric matrix of profile probabilities (cells × profiles).

- cluster_label:

  A vector of cluster identities (length must match number of rows in
  \`cell_profile_prob\`).

## Value

A numeric matrix where rows represent clusters and columns represent
profile sums.

## Details

Computes the total probability mass for each profile across a set of
clusters. Typically used to summarize soft profile assignments (e.g.,
from Clonotrace) across transcriptional clusters.
