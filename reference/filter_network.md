# filter_network

Filter Sparse Graph to Retain Nodes with Minimum Degree

## Usage

``` r
filter_network(adj, n_neighbors = 5)
```

## Arguments

- adj:

  A symmetric adjacency matrix (numeric or sparse).

- n_neighbors:

  Integer. The minimum number of neighbors each node must have to remain
  in the network (default: 5).

## Value

A filtered adjacency matrix containing only nodes with at least
\`n_neighbors\` connections.

## Details

Iteratively removes nodes from a symmetric adjacency matrix until all
remaining nodes have at least \`n_neighbors\` edges.
