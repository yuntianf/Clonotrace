# build_edges

Build Temporal Edge Table from Time-Series Data Frame

## Usage

``` r
build_edges(df, id_col, order_col, feature_cols = NULL)
```

## Arguments

- df:

  A data frame with time-ordered entries per group.

- id_col:

  Column name (string) indicating group identity (e.g., cell or clone
  ID).

- order_col:

  Column name (string) indicating time or order.

- feature_cols:

  Optional. Vector of column names to use as node features. If \`NULL\`,
  all columns except \`id_col\` and \`order_col\` are used.

## Value

A data frame where each row represents a temporal edge with start/end
times and next-step features.

## Details

For each group defined by \`id_col\`, constructs edges between
consecutive time points based on \`order_col\`. Adds forward-shifted
features (e.g., lead values) to represent transitions.

## Examples

``` r
df <- data.frame(id = c(1, 1, 2, 2), time = c(1, 2, 1, 2), val = c(5, 6, 3, 4))
edges <- build_edges(df, id_col = "id", order_col = "time")
print(edges)
#> # A tibble: 2 × 5
#> # Groups:   id [2]
#>      id Start   End   val Next_val
#>   <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1     1     1     2     5        6
#> 2     2     1     2     3        4
```
