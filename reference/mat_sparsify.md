# mat_sparsify

Sparsify a Dense Matrix by Cumulative Mass Filtering

## Usage

``` r
mat_sparsify(mat, row_mass = 0.9, col_mass = 0.9)
```

## Arguments

- mat:

  A numeric matrix.

- row_mass:

  Proportion of row-wise cumulative mass to retain (default: 0.9).

- col_mass:

  Proportion of column-wise cumulative mass to retain (default: 0.9).

## Value

A sparsified matrix of the same shape with many zero entries.

## Details

Applies mass-based thresholding to rows and columns to keep only the
most informative values. For each row/column, entries are retained until
the cumulative sum reaches a specified mass (e.g., 90

## Examples

``` r
mat <- matrix(runif(100), nrow = 10)
mat_sparse <- mat_sparsify(mat, row_mass = 0.8, col_mass = 0.8)
print(mat_sparse)
#>            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,] 0.0000000 0.6772402 0.0000000 0.0000000 0.0000000 0.0000000 0.9751419
#>  [2,] 0.0000000 0.0000000 0.9216622 0.0000000 0.0000000 0.0000000 0.9349492
#>  [3,] 0.0000000 0.8943790 0.0000000 0.0000000 0.9522314 0.8714917 0.0000000
#>  [4,] 0.0000000 0.0000000 0.0000000 0.0000000 0.8590911 0.0000000 0.0000000
#>  [5,] 0.7442099 0.0000000 0.0000000 0.8301855 0.9042145 0.0000000 0.0000000
#>  [6,] 0.0000000 0.6855531 0.8367631 0.0000000 0.6255544 0.0000000 0.0000000
#>  [7,] 0.0000000 0.9749299 0.0000000 0.9771000 0.0000000 0.0000000 0.9187957
#>  [8,] 0.8580543 0.0000000 0.0000000 0.0000000 0.0000000 0.9072732 0.0000000
#>  [9,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.6479407 0.9948942
#> [10,] 0.0000000 0.0000000 0.8261784 0.0000000 0.0000000 0.0000000 0.0000000
#>            [,8]      [,9]     [,10]
#>  [1,] 0.0000000 0.0000000 0.8473921
#>  [2,] 0.0000000 0.7456405 0.0000000
#>  [3,] 0.0000000 0.0000000 0.0000000
#>  [4,] 0.7183564 0.0000000 0.9100310
#>  [5,] 0.0000000 0.0000000 0.0000000
#>  [6,] 0.7244172 0.0000000 0.0000000
#>  [7,] 0.7781644 0.7985098 0.9165301
#>  [8,] 0.0000000 0.8074137 0.0000000
#>  [9,] 0.0000000 0.0000000 0.0000000
#> [10,] 0.9295000 0.0000000 0.0000000
```
