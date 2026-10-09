# compute_transition

Compute Normalized Transition Matrix from Connectivity Matrix

## Usage

``` r
compute_transition(connectivity)
```

## Arguments

- connectivity:

  A symmetric numeric or sparse matrix representing pairwise affinities
  (e.g., k-NN graph).

## Value

A symmetric matrix representing the normalized transition matrix.

## Details

Converts a symmetric connectivity matrix into a double-normalized
transition matrix suitable for diffusion-based analyses. The output
matrix is scaled such that it preserves local affinities and is
symmetric with unit row and column scaling.

The transition matrix is computed as: \$\$T = Z D^{-1} A D^{-1} Z\$\$
where \\A\\ is the input matrix, \\D\\ is the row sum diagonal matrix,
and \\Z\\ rescales the result for symmetry.

## Examples

``` r
library(Matrix)
mat <- Matrix::rsparsematrix(5, 5, density = 0.6)
mat <- 0.5 * (mat + t(mat))  # ensure symmetry
trans <- compute_transition(mat)
#> Warning: NaNs produced
print(trans)
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                              
#> [1,]   . .        NaN NaN NaN
#> [2,]   . 0.710773 NaN NaN NaN
#> [3,] NaN      NaN NaN NaN NaN
#> [4,] NaN      NaN NaN NaN   .
#> [5,] NaN      NaN NaN   . NaN
```
