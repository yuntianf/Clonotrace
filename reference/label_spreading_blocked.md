# label_spreading_blocked

Memory-efficient label spreading (blocked, HDF5 optional) Runs label
spreading on a sparse graph in \*\*class blocks\*\* to control memory,
optionally streaming class probabilities to an on-disk HDF5 dataset.

## Usage

``` r
label_spreading_blocked(
  adj,
  labels,
  label_n = NULL,
  alpha = 0.9,
  max_iter = 100,
  tol = 0.001,
  block_size = 128,
  outfile = NULL,
  verbose = TRUE,
  epsilon = 0
)
```

## Arguments

- adj:

  `dgCMatrix` sparse adjacency; will be row-normalized in-place to act
  as a transition matrix.

- labels:

  Integer or numeric vector of length \\N\\ with class labels in `1..C`
  and `NA` for unlabeled cells. Do not refactor.

- label_n:

  Optional total number of classes; inferred if `NULL`.

- alpha:

  Damping parameter in label propagation update.

- max_iter:

  Maximum iterations per block.

- tol:

  Convergence tolerance for updates.

- block_size:

  Number of classes per block to process at once.

- outfile:

  Optional HDF5 file path; if `NULL`, returns dense matrix.

- verbose:

  Logical; show progress messages.

- epsilon:

  Small uniform prior mass per class to avoid zeros.

## Value

If `outfile = NULL`, a dense \\N \times C\\ numeric matrix of class
probabilities. Otherwise, invisibly returns the HDF5 file path with
dataset `"prob"` written.

## Details

\- Normalizes the adjacency to a row-stochastic matrix. - Splits the
label space into `block_size`-sized groups to limit RAM use. -
Optionally writes each probability block directly into an HDF5 dataset.

## See also

[`create_hdf5_matrix`](https://yuntianf.github.io/Clonotrace/reference/create_hdf5_matrix.md),
[`normalize_hdf5_rows`](https://yuntianf.github.io/Clonotrace/reference/normalize_hdf5_rows.md),
[`label_spreading_bootstrap_blocked`](https://yuntianf.github.io/Clonotrace/reference/label_spreading_bootstrap_blocked.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Matrix)
set.seed(1)
N <- 1000; C <- 5
A <- rsparsematrix(N, N, 0.001); A <- abs(A) + Diagonal(N)
labs <- sample(c(1:C, NA), N, TRUE, prob = c(rep(0.05, C), 0.75))
P <- label_spreading_blocked(A, labs, label_n = C, block_size = 2, max_iter = 50)
} # }
```
