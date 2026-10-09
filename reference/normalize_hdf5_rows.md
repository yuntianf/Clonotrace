# normalize_hdf5_rows

Row-normalize an HDF5 probability matrix in column blocks Iteratively
reads an \\N \times C\\ HDF5 dataset in column blocks, adds a small
epsilon if requested, computes row sums, and writes back normalized
values so each row sums to 1.

## Usage

``` r
normalize_hdf5_rows(h5file, dataset = "prob", block_cols = 256, add_eps = TRUE)
```

## Arguments

- h5file:

  Path to HDF5 file.

- dataset:

  Dataset name (default `"prob"`).

- block_cols:

  Number of columns per I/O block.

- add_eps:

  Logical; if `TRUE`, adds `1e-12` before summing.

## Value

Numeric vector of length \\N\\ containing the pre-normalization row
sums.

## Examples

``` r
if (FALSE) { # \dontrun{
normalize_hdf5_rows("prob_out.h5", "prob", block_cols = 128)
} # }
```
