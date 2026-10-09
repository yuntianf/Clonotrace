# mat_split

Split a Matrix into Row or Column Chunks

## Usage

``` r
mat_split(mat, size, index = 0)
```

## Arguments

- mat:

  A matrix to be split.

- size:

  An integer vector specifying the size of each chunk (rows or columns).

- index:

  Integer. If \`0\` (default), splits by rows; if \`1\`, splits by
  columns.

## Value

A list of matrices representing the chunks.

## Details

Splits a matrix into smaller matrices based on user-defined row or
column sizes.

## Examples

``` r
mat <- matrix(1:20, nrow = 5)
chunks <- mat_split(mat, size = c(2, 3), index = 0)
str(chunks)
#> List of 2
#>  $ : int [1:2, 1:4] 1 2 6 7 11 12 16 17
#>  $ : int [1:3, 1:4] 3 4 5 8 9 10 13 14 15 18 ...
```
