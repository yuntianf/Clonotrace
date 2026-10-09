# sparse_manupulation

Apply a Transformation Function to a Sparse Matrix

## Usage

``` r
sparse_manupulation(mat, func)
```

## Arguments

- mat:

  A sparse matrix (class \`dgCMatrix\`).

- func:

  A function to apply to the values (e.g., \`log\`, \`exp\`, or custom).

## Value

A sparse matrix with the transformed values.

## Details

Applies a user-defined function to the non-zero entries of a sparse
matrix and returns a new sparse matrix.

## Examples

``` r
library(Matrix)
mat <- rsparsematrix(5, 5, density = 0.5)
transformed <- sparse_manupulation(mat, log1p)
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `dis = func(x)`.
#> Caused by warning in `func()`:
#> ! NaNs produced
print(transformed)
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>             1         2         3         4           5
#> 1         NaN .         .         0.1133287  0.03149867
#> 2  .          .         0.8754687 .          .         
#> 3 -0.02839947 .         .         1.0986123 -0.17435339
#> 4  0.45107562 .         0.5709795 0.4317824  .         
#> 5  .          0.6205765 .         0.1906204  .         
```
