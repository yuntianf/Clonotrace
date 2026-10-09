# dis_point_to_edge

Compute Distance from a Point to a Line Segment

## Usage

``` r
dis_point_to_edge(point, edge_start, edge_end)
```

## Arguments

- point:

  A numeric vector representing a point.

- edge_start:

  A numeric vector representing the start of the edge.

- edge_end:

  A numeric vector representing the end of the edge.

## Value

A numeric vector of length 2: the shortest distance and the scalar
projection (t in \[0,1\]).

## Details

Calculates the shortest distance from a point to a line segment in 2D or
3D space, and returns the scalar projection position along the segment.

## Examples

``` r
d <- dis_point_to_edge(c(0, 0), c(1, 0), c(1, 1))
print(d)  # Expect ~1.0 distance to the vertical segment
#> [1] 1 0
```
