# long2square

Convert Long Format Data Frame to a Square Matrix

## Usage

``` r
long2square(
  long,
  row_names_from,
  col_names_from,
  values_from,
  symmetric = TRUE,
  na.fill = NA,
  nodes = NULL
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

- symmetric:

  Logical. If \`TRUE\`, the output will be forced to be symmetric by
  duplicating reversed entries.

- na.fill:

  Value to fill in for missing matrix entries (default: NA).

- nodes:

  Optional. A character vector specifying the full set of row/column
  labels.

## Value

A square numeric matrix with rows and columns defined by the node set.

## Details

Converts a long-format data frame (pairwise relationships) into a square
numeric matrix. Optionally enforces symmetry and fills missing values.

## Examples

``` r
if (FALSE) long_df <- data.frame(from = c("A", "B"), 
        to = c("B", "A"), value = c(1, 2))
square_mat <- long2square(long_df, 
    row_names_from = "from", 
      col_names_from = "to", 
      values_from = "value")
#> Error: object 'long_df' not found
print(square_mat) # \dontrun{}
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'square_mat' not found
```
