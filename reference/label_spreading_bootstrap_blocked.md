# label_spreading_bootstrap_blocked

Bootstrap label spreading with HDF5 reference and per-cell deviance

Performs multiple label-spreading runs on random subsets of labeled
cells, compares each to a full reference solution, and returns the
reference probabilities plus mean L1 deviance per cell.

## Usage

``` r
label_spreading_bootstrap_blocked(
  adj,
  labels,
  alpha = 0.8,
  sample_rate = 0.8,
  sample_n = 50,
  block_size = 128,
  tol = 0.001,
  max_iter = 100,
  epsilon = 0,
  refer_h5 = NULL,
  tmpdir = tempdir(),
  verbose = TRUE
)
```

## Arguments

- adj:

  Sparse adjacency (`dgCMatrix`); will be row-normalized.

- labels:

  Vector of original labels (length \\N\\); mapped once to `1..C`.

- alpha:

  Damping parameter for label spreading.

- sample_rate:

  Fraction of labeled cells retained in each bootstrap.

- sample_n:

  Number of bootstrap replicates.

- block_size:

  Classes per block when streaming to HDF5.

- tol, max_iter:

  Convergence parameters.

- epsilon:

  Uniform prior mass added to avoid zeros.

- refer_h5:

  Optional existing reference HDF5 file; if `NULL`, computed internally.

- tmpdir:

  Directory for temporary HDF5 files.

- verbose:

  Logical; print progress.

## Value

A list with:

- `prob`:

  Dense reference \\N \times C\\ probability matrix.

- `deviance`:

  Per-cell mean L1 deviation vs reference.

## Details

\- Computes an exact reference with all labeled cells, normalizes it,
and stores it in HDF5. - Each bootstrap uses a random subset of labeled
cells, writes its results to a temporary HDF5, and accumulates the L1
deviation vs the reference in streaming mode. - Parallelizable via
future.apply: set a plan before calling. - A reference created
internally is removed on exit. A caller-supplied `refer_h5` is treated
as a cache and is preserved.

Only intermediate class blocks are disk-backed. The returned `prob`
matrix is dense for compatibility, so callers should still budget
roughly \\8NC\\ bytes for one final double-precision cells-by-label
matrix.

## See also

[`label_spreading_blocked`](https://yuntianf.github.io/Clonotrace/reference/label_spreading_blocked.md),
[`normalize_hdf5_rows`](https://yuntianf.github.io/Clonotrace/reference/normalize_hdf5_rows.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Matrix); library(future.apply)
plan(multisession, workers = 4)
N <- 1000; C <- 4
A <- rsparsematrix(N, N, 0.001); A <- abs(A) + Diagonal(N)
labs <- sample(c(1:C, NA), N, TRUE, prob = c(rep(0.05, C), 0.8))
res <- label_spreading_bootstrap_blocked(A, labs, sample_rate = 0.8, sample_n = 10)
str(res)
} # }
```
