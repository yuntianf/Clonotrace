# DPT_T

Compute Diffusion Pseudotime from Accumulated Commute Time Matrix

## Usage

``` r
DPT_T(T_mat, start)
```

## Arguments

- T_mat:

  A row-normalized transition matrix (symmetric or stochastic).

- start:

  Integer. Index of the root cell for pseudotime computation.

## Value

A numeric vector of diffusion pseudotime values for all nodes.

## Details

Computes diffusion pseudotime (DPT) distances from a specified root node
using the accumulated commute time matrix.
