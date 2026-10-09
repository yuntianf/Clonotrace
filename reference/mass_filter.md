# mass_filter

Filter Vector Elements by Cumulative Mass Threshold

## Usage

``` r
mass_filter(mass, thresh = 0.9)
```

## Arguments

- mass:

  A numeric vector of values (e.g., weights, similarities, or
  densities).

- thresh:

  Numeric. Proportion of total cumulative mass to retain (default: 0.9).

## Value

A numeric vector of the same length as \`mass\`, with low-contributing
elements set to zero.

## Details

Retains the top elements in a numeric vector based on cumulative mass.
Elements contributing least to the total are zeroed out until the
specified proportion of total mass (\`thresh\`) is retained.

This function normalizes the input to sum to 1, ranks the elements by
size, and sets to zero those contributing below the cutoff needed to
achieve \`thresh\` cumulative mass. Useful for sparsifying matrices.

## Examples

``` r
x <- runif(10)
filtered <- mass_filter(x, thresh = 0.8)
print(filtered)
#>  [1] 0.0000000 0.6598604 0.9352709 0.0000000 0.0000000 0.5649702 0.8195636
#>  [8] 0.0000000 0.9539901 0.0000000
```
