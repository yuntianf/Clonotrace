# ceil_digit

Round Up a Number to a Specified Number of Decimal Places

## Usage

``` r
ceil_digit(x, n)
```

## Arguments

- x:

  A numeric vector to round up.

- n:

  Integer. Number of decimal places to retain.

## Value

A numeric vector rounded up to \`n\` digits.

## Details

Rounds numeric values upward to the nearest multiple at \`n\` decimal
places.

## Examples

``` r
ceil_digit(3.1415, 2)  # 3.15
#> [1] 3.15
ceil_digit(c(2.001, 2.009), 2)  # c(2.01, 2.01)
#> [1] 2.01 2.01
```
