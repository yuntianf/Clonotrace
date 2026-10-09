# long2sparse

Convert Long-Format Data Frame to Sparse Matrix

## Usage

``` r
long2sparse(
  long,
  row_names_from,
  col_names_from,
  values_from,
  unique_rows = NULL,
  unique_cols = NULL,
  symmetric = FALSE
)
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

- unique_rows:

  Optional character or factor vector defining the full set and order of
  row names.

- unique_cols:

  Optional character or factor vector defining the full set and order of
  column names.

- symmetric:

  Logical. If \`TRUE\`, creates a symmetric matrix by adding reversed
  entries (default: FALSE).

## Value

A sparse matrix of class \`dgCMatrix\` with row and column names.

## Details

Converts a long-format data frame representing pairwise relationships
into a sparse matrix (\`dgCMatrix\`), with optional symmetry enforcement
and customizable row/column ordering.

## Examples

``` r
long_df <- data.frame(from = c("A", "B"), to = c("B", "C"), value = c(1, 2))
spmat <- long2sparse(long_df, "from", "to", "value")
print(spmat)
#> 2 x 2 sparse Matrix of class "dgCMatrix"
#>   B C
#> A 1 .
#> B . 2
```
