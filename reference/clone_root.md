# clone_root

Identify Root Clone from Cluster Enrichment

## Usage

``` r
clone_root(clones, cell_meta, clone_col, cluster_col, start_cluster)
```

## Arguments

- clones:

  A character vector of clone identifiers.

- cell_meta:

  A data frame with clone and cluster information per cell.

- clone_col:

  Name of the column containing clone labels.

- cluster_col:

  Name of the column containing cluster labels.

- start_cluster:

  Name or value of the cluster to define the root state.

## Value

The name of the most enriched clone in the start cluster.

## Details

Finds the most enriched clone in a given starting cluster, to use as a
pseudotime root.
