# long_symmetry

Enforce Symmetry in a Long-Format Pairwise Data Frame

## Usage

``` r
long_symmetry(long, row_names_from, col_names_from)
```

## Arguments

- long:

  A long-format data frame with at least two ID columns and one or more
  data columns.

- row_names_from:

  Name of the column representing the source node or element.

- col_names_from:

  Name of the column representing the target node or element.

## Value

A symmetric long-format data frame with reversed entries added and
sorted.

## Details

Ensures a symmetric representation of long-format data by adding
reversed pairs and removing duplicates. Useful for converting pairwise
relationships (e.g., distances or similarities) into symmetric forms.

## Examples

``` r
long_df <- data.frame(from = c("A", "B"), to = c("B", "C"), value = c(1, 2))
sym_df <- long_symmetry(long_df, row_names_from = "from", col_names_from = "to")
print(sym_df)
#>   from to value
#> 1    A  B     1
#> 2    B  A     1
#> 3    B  C     2
#> 4    C  B     2
```
