# acct

Compute Accumulated Commute Time Matrix from Transition Matrix

## Usage

``` r
acct(T_mat)
```

## Arguments

- T_mat:

  A square transition matrix (symmetric, row-normalized).

## Value

A matrix of accumulated commute times between all pairs of nodes.

## Details

Computes the accumulated commute time (ACT) matrix from a normalized
transition matrix. This is used as a basis for diffusion-based distances
and pseudotime inference.
