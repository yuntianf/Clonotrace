# bin_filter

Iterative Binning Filter to Remove Sparse Bins

## Usage

``` r
bin_filter(x, col, thresh = 10, breaks = 100)
```

## Arguments

- x:

  A data frame.

- col:

  A string specifying the column to bin and filter by.

- thresh:

  Minimum number of rows required per bin (default: 10).

- breaks:

  Number of breaks to use in \`cut()\` (default: 100).

## Value

A filtered data frame with the same columns as \`x\`, excluding the
temporary bin column.

## Details

Filters rows from a data frame such that each bin of a numeric column
contains at least \`thresh\` observations. Binning is performed using
\`cut()\` over a specified number of breaks. The process repeats until
all bins meet the threshold.
