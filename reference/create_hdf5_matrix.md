# create_hdf5_matrix

Create a chunked HDF5 matrix for streaming label-spreading output

Preallocates a compressed HDF5 dataset (default name `"prob"`) for
storing class-probability blocks streamed from
[`label_spreading_blocked`](https://yuntianf.github.io/Clonotrace/reference/label_spreading_blocked.md).

## Usage

``` r
create_hdf5_matrix(
  h5file,
  dataset = "prob",
  nrow,
  ncol,
  chunk_rows = 4096,
  chunk_cols = 128
)
```

## Arguments

- h5file:

  Output HDF5 file path.

- dataset:

  Dataset name to create.

- nrow, ncol:

  Dimensions of the dataset.

- chunk_rows, chunk_cols:

  Chunk sizes for row/column blocks.

## Value

Invisibly returns `TRUE` after successful creation.

## Details

Removes any existing file before creation. Uses gzip compression
(`level = 7`) and double precision storage.

## Examples

``` r
if (FALSE) { # \dontrun{
create_hdf5_matrix("prob_out.h5", "prob", nrow = 5000, ncol = 256)
} # }
```
