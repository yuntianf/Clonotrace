# ridge_regression

Ridge Regression for Gene Expression Modeling

## Usage

``` r
ridge_regression(X, G, lambda = 1e-04)
```

## Arguments

- X:

  A numeric matrix of predictors (cells × features).

- G:

  A numeric response matrix or vector. Must match number of rows in
  \`X\` (cells × genes or 1).

- lambda:

  Non-negative regularization strength for ridge penalty (default:
  1e-4).

## Value

A list with components:

- coef:

  Estimated coefficient matrix (features × genes)

- rss:

  Residual sum of squares per gene

## Details

Performs ridge regression using a closed-form solution. Can be applied
to model either a single gene (vector response) or multiple genes
(matrix response) using a shared regularization penalty.
