# dismat_mst

Compute the Minimum Spanning Tree (MST) from a Distance Matrix

## Usage

``` r
dismat_mst(mat)
```

## Arguments

- mat:

  A symmetric numeric distance matrix.

## Value

A data frame with columns \`from\`, \`to\`, and \`weight\` representing
MST edges.

## Details

Converts a full square distance matrix into a graph and computes the
minimum spanning tree. Returns the MST edges as a long-format data
frame.
