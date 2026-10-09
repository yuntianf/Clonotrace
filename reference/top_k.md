# top_k

Extract Top-k Elements from a Numeric Vector

## Usage

``` r
top_k(v, k, decreasing = FALSE)
```

## Arguments

- v:

  A numeric vector.

- k:

  Integer. Number of top elements to return.

- decreasing:

  Logical. If \`TRUE\`, returns the \`k\` largest values (default:
  \`FALSE\` = smallest).

## Value

A numeric vector of length \`k\` containing the top \`k\` values.

## Details

Returns the top \`k\` elements from a numeric vector, sorted in
ascending or descending order.
