# is_symmetric

Check if a Matrix is Symmetric

## Usage

``` r
is_symmetric(matrix)
```

## Arguments

- matrix:

  A square numeric matrix.

## Value

Logical value: \`TRUE\` if the matrix is symmetric, \`FALSE\` otherwise.

## Details

Determines whether a matrix is equal to its transpose within numerical
tolerance.

## Examples

``` r
mat <- matrix(c(1, 2, 2, 1), nrow = 2)
is_symmetric(mat)  # TRUE
#> [1] TRUE
```
