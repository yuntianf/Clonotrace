# soft_cluster_gam_fit

Fit Soft Cluster-Specific GAMs with Ridge Regularization

## Usage

``` r
soft_cluster_gam_fit(G, t, P, df = 5, lambda = 1e-04, test = "F")
```

## Arguments

- G:

  A matrix of gene expression values (genes × cells).

- t:

  A numeric vector of pseudotime values (length = number of cells).

- P:

  A soft cluster assignment matrix (cells × clusters), where rows sum to
  1.

- df:

  Degrees of freedom for spline basis (default: 5).

- lambda:

  Ridge regularization strength (default: 1e-4).

- test:

  Statistical test to use for model comparison: \`"F"\` (default) or
  \`"LRT"\` (likelihood ratio test).

## Value

A list with components:

- stat:

  A data frame with test statistic, p-value, mean expression difference,
  and Cohen's d for each gene

- df:

  Degrees of freedom used for test

- coef:

  Fitted coefficients from the full model

- design_null:

  Design matrix for the null model

- design_full:

  Design matrix for the full model

## Details

Fits a generalized additive model (GAM) to gene expression along
pseudotime using soft cluster assignments. The model compares a shared
smooth (null model) to cluster-specific smooths (full model), using
ridge regression and F or likelihood ratio tests.

## Examples

``` r
# Example pseudotime, profiles, and expression matrix setup
set.seed(123)
t <- seq(0, 1, length.out = 100)
P <- matrix(runif(300), nrow = 100)
P <- P / rowSums(P)
G <- matrix(rnorm(500), nrow = 5)
fit <- soft_cluster_gam_fit(G, t, P, df = 3)
head(fit$stat)
#>        stat      pval    mean_diff        cohen
#> 1 0.7847855 0.6172285  0.785680855  0.708225587
#> 2 0.4980950 0.8545586  0.661185867  0.778625761
#> 3 0.7429910 0.6533914 -0.621810882 -0.700928364
#> 4 1.3263654 0.2410829  0.005949389  0.005999336
#> 5 1.5342015 0.1569317 -0.218150387 -0.235564932
```
