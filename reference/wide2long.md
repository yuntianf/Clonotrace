# wide2long

Convert a Matrix to Long Format

## Usage

``` r
wide2long(mat)
```

## Arguments

- mat:

  A numeric matrix or data frame.

## Value

A data frame in long format with columns:

- i:

  Row index

- j:

  Column index (as numeric)

- value:

  Matrix entry at (i, j)

## Details

Converts a numeric matrix or data frame into a long-format data frame
with integer-based \`i\`, \`j\`, and \`value\` columns.

## Examples

``` r
mat <- matrix(1:9, nrow = 3)
long_df <- wide2long(mat)
head(long_df)
#>   i j value
#> 1 1 1     1
#> 2 1 2     4
#> 3 1 3     7
#> 4 2 1     2
#> 5 2 2     5
#> 6 2 3     8
```
