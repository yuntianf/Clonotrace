# dpt

Compute Diffusion Pseudotime from Transition Matrix (Eigen
Decomposition)

## Usage

``` r
dpt(T_mat, root, k = 30)
```

## Arguments

- T_mat:

  A square transition matrix (e.g., from kNN graph).

- root:

  Integer. Index of the root cell.

- k:

  Integer. Number of eigenvectors to use (default: 30).

## Value

A numeric vector of pseudotime values normalized to \[0, 1\].

## Details

Calculates pseudotime by spectral decomposition of a diffusion
transition matrix, using a root cell as a reference point in the reduced
diffusion space.
