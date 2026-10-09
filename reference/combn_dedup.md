# combn_dedup

Remove Duplicate Combinations (Order-Invariant)

## Usage

``` r
combn_dedup(combn)
```

## Arguments

- combn:

  A two-column matrix where each row is a pairwise combination.

## Value

A logical vector indicating whether each row is a unique
(order-insensitive) combination.

## Details

Given a matrix of combinations (e.g., from \`combn()\`), this function
removes duplicates where row elements are identical up to ordering.

## Examples

``` r
combn_mat <- matrix(c(1,2,2,1,3,4), ncol = 2, byrow = TRUE)
unique_flags <- combn_dedup(combn_mat)
print(unique_flags)
#> [1]  TRUE FALSE  TRUE
```
