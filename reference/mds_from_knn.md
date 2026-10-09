# mds_from_knn

Compute MDS Embedding from a kNN Graph

## Usage

``` r
mds_from_knn(adj, n_components = 15)
```

## Arguments

- adj:

  A sparse adjacency matrix (class \`dgCMatrix\`) representing a kNN
  graph.

- n_components:

  Integer. Number of MDS dimensions to return (default: 15).

## Value

A data frame of size (nodes × \`n_components\`) with column names
\`mds_1\`, \`mds_2\`, ..., and rownames from \`adj\`.

## Details

Computes a multi-dimensional scaling (MDS) embedding from a filtered
adjacency matrix representing a kNN graph. Requires a Python-based
backend function (e.g., \`mds_from_knn_py()\`) to perform the embedding.
