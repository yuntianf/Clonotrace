# umap_from_knn

Compute UMAP Embedding from a kNN Graph

## Usage

``` r
umap_from_knn(adj, n_neighbors = 5, seed = 1024)
```

## Arguments

- adj:

  A sparse adjacency matrix (class \`dgCMatrix\`) representing a kNN
  graph.

- n_neighbors:

  Integer. Minimum number of neighbors required per node (default: 5).

- seed:

  Integer. Random seed for UMAP initialization (default: 1024).

## Value

A data frame with columns \`umap_1\` and \`umap_2\` and rownames
matching input adjacency matrix.

## Details

Runs UMAP on a filtered kNN adjacency matrix, ensuring that all nodes
have at least \`n_neighbors\`. Assumes the embedding is computed via a
Python-based backend (e.g., \`umap_from_knn_py()\`).

This function filters the graph to remove low-degree nodes before
computing UMAP.
