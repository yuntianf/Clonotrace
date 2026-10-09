# long2wide

Convert Long Format Data Frame to Symmetric or Asymmetric Wide Matrix

## Usage

``` r
long2wide(long, row_names_from, col_names_from, values_from, symmetric = FALSE)
```

## Arguments

- long:

  A data frame in long format with at least three columns: row ID,
  column ID, and value.

- row_names_from:

  Name of the column in \`long\` to use as row names in the output.

- col_names_from:

  Name of the column in \`long\` to use as column names in the output.

- values_from:

  Name of the column in \`long\` to use as cell values in the output
  matrix.

- symmetric:

  Logical. If \`TRUE\`, the output will be forced to be symmetric by
  duplicating reversed entries.

## Value

A wide-format data frame or matrix where row names and column names are
as specified, and cell values come from \`values_from\`.

## Details

Converts a long-format data frame with row/column/value identifiers into
a wide-format matrix-style data frame. Optionally enforces matrix
symmetry by duplicating and merging reversed entries.

## Examples

``` r
long_df <- data.frame(from = c("A", "B"), to = c("B", "A"), value = c(1, 2))
mat <- long2wide(long_df, row_names_from = "from", col_names_from = "to", values_from = "value")
print(mat)
#>    B  A
#> A  1 NA
#> B NA  2
```
