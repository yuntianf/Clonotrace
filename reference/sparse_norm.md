# sparse_norm

Normalize a Sparse Matrix by Rows or Columns

## Usage

``` r
sparse_norm(mat, dim = 1)
```

## Arguments

- mat:

  A sparse matrix (class \`dgCMatrix\`).

- dim:

  Integer. Dimension to normalize over: \`1\` = rows (default), \`2\` =
  columns.

## Value

A sparse matrix with normalized values.

## Details

Applies row-wise or column-wise L1 normalization to a sparse matrix.
Returns a new sparse matrix with normalized entries, preserving
sparsity.

## Examples

``` r
library(Matrix)
m <- rsparsematrix(5, 5, density = 0.4)
norm_mat <- sparse_norm(m, dim = 1)
print(norm_mat)
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                                                       
#> [1,] .         0.3953488 0.6046512 .         .        
#> [2,] 0.2794118 .         .         0.7205882 .        
#> [3,] 0.3596491 .         .         .         0.6403509
#> [4,] .         .         .         1.0000000 .        
#> [5,] 0.3575130 .         0.5181347 .         0.1243523
```
