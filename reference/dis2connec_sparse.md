# dis2connec_sparse

Convert Sparse Distance Matrix to Sparse Connectivity Matrix

## Usage

``` r
dis2connec_sparse(D)
```

## Arguments

- D:

  A symmetric sparse distance matrix (class \`dgCMatrix\`).

## Value

A sparse matrix of connectivities (class \`dgCMatrix\`) where weights
represent local affinities.

## Details

Transforms a sparse distance matrix into a Gaussian kernel-based
connectivity matrix using a locally adaptive bandwidth. The resulting
matrix can be interpreted as a similarity graph.

## Examples

``` r
library(Matrix)
D <- rsparsematrix(100, 100, density = 0.05)
D <- 0.5 * (D + t(D))  # ensure symmetry
C <- dis2connec_sparse(D)
print(C[1:5, 1:5])
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>   1 2 3 4 5
#> 1 . . . . .
#> 2 . . . . .
#> 3 . . . . .
#> 4 . . . . .
#> 5 . . . . .
```
