# clone_2_ot

Compute Optimal Transport Distance Between Two Clones

## Usage

``` r
clone_2_ot(distance, group1_mass, group2_mass)
```

## Arguments

- distance:

  A numeric cost matrix of pairwise distances (rows = group 1 cells,
  cols = group 2 cells).

- group1_mass:

  A numeric vector of mass/probabilities for group 1 (length must equal
  number of rows in \`distance\`).

- group2_mass:

  A numeric vector of mass/probabilities for group 2 (length must equal
  number of columns in \`distance\`).

## Value

A single numeric value representing the OT distance between the two
groups.

## Details

Calculates the optimal transport (OT) distance between two sets of cells
(e.g., clones) based on a cost matrix and their mass distributions. Uses
the \`transport\` package for solving the OT plan.
