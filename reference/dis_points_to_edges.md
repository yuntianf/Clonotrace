# dis_points_to_edges

Compute Distances from Multiple Points to Multiple Edges

## Usage

``` r
dis_points_to_edges(points, edges)
```

## Arguments

- points:

  A data frame or matrix where each row is a point (2D or 3D).

- edges:

  A list of edge definitions, where each element is a 2-row matrix: row
  1 = start, row 2 = end.

## Value

A list with two components:

- map:

  Matrix of scalar projection values (t in \[0,1\]) for each point–edge
  pair.

- dis:

  Matrix of shortest distances for each point–edge pair.

## Details

Computes the distance from each point to each edge, and also returns the
scalar projection for each edge to understand where the closest point
falls along the segment.

## Examples

``` r
pts <- data.frame(x = rnorm(3), y = rnorm(3))
edge1 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
edge2 <- matrix(c(0, 1, 1, 0), ncol = 2, byrow = TRUE)
dists <- dis_points_to_edges(pts, list(edge1, edge2))
str(dists)
#> List of 2
#>  $ map: num [1:3, 1:2] 0.4806 2.6259 0.7426 0.0183 2.4317 ...
#>  $ dis: num [1:3, 1:2] 0.513 0 0 0.16 0 ...
```
