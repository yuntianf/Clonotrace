# connectivity_coord

Convert Connectivity Matrix into Edge Coordinates for Plotting

## Usage

``` r
connectivity_coord(coord, connectivity, dims = c(1, 2))
```

## Arguments

- coord:

  A numeric matrix or data frame of coordinates (rows = points, columns
  = dimensions).

- connectivity:

  A square numeric or sparse matrix representing connectivity between
  points.

- dims:

  A numeric vector of length 2 specifying which columns of \`coord\` to
  use for x and y coordinates.

## Value

A data frame with one row per edge

## Details

Transforms a connectivity matrix and a coordinate matrix into a data
frame containing edge coordinates. This is useful for plotting
graph-like structures (e.g., k-NN graphs or SNNs) using
\`ggplot2::geom_segment\`.
