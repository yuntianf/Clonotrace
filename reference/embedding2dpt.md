# embedding2dpt

Compute Diffusion Pseudotime from Embedding

## Usage

``` r
embedding2dpt(embedding, nn_k, root, dpt_k = 30)
```

## Arguments

- embedding:

  A numeric matrix of embeddings (e.g., UMAP or PCA).

- nn_k:

  Integer. Number of neighbors in the kNN graph.

- root:

  Integer. Index of the root cell.

- dpt_k:

  Integer. Number of eigenvectors used in DPT (default: 30).

## Value

A numeric vector of pseudotime values for each row in \`embedding\`.

## Details

Builds a kNN graph from an embedding matrix and computes diffusion
pseudotime starting from a root cell.

## Examples

``` r
mat <- matrix(rnorm(200), ncol = 5)
pt <- embedding2dpt(mat, nn_k = 10, root = 1)
```
