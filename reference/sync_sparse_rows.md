# sync_sparse_rows

Synchronize Rows of a Sparse Matrix to Match Target Row Names

## Usage

``` r
sync_sparse_rows(sparse_mat, target_rownames)
```

## Arguments

- sparse_mat:

  A sparse matrix of class \`dgCMatrix\`, with rownames defined.

- target_rownames:

  A character vector of desired rownames in the final matrix.

## Value

A sparse matrix with the same number of columns and rownames exactly
matching \`target_rownames\`.

## Details

Ensures that a sparse matrix of class \`dgCMatrix\` has rows that
exactly match a target set of row names. Adds missing rows filled with
zeros if necessary, and reorders existing rows to match the specified
order.

This is useful when aligning matrices from different sources (e.g.,
profiles, probabilities, annotations) in cases where missing rows should
be treated as zero.

## Examples

``` r
library(Matrix)
mat <- rsparsematrix(3, 4, density = 0.5)
rownames(mat) <- c("A", "B", "C")
synced <- sync_sparse_rows(mat, target_rownames = c("A", "B", "C", "D"))
dim(synced)  # Should be 4 x 4
#> [1] 4 4
```
