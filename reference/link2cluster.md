# link2cluster

Cluster Nodes Based on Link Connectivity

## Usage

``` r
link2cluster(link, nodes)
```

## Arguments

- link:

  A data frame with two columns (i, j) indicating node pairs with links.

- nodes:

  A character or factor vector of all node names to be clustered.

## Value

A numeric vector of cluster labels for each node, as assigned by DBSCAN.

## Details

Given a set of links between node pairs and a list of all nodes, this
function constructs a symmetric adjacency matrix, applies diffusion to
expand connectivity, and clusters the result using DBSCAN.
