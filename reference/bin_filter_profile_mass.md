# bin_filter_profile_mass

Filter Cells by Profile Mass within Pseudotime Bins

## Usage

``` r
bin_filter_profile_mass(mass, time, thresh = 5, binsize = 0.005)
```

## Arguments

- mass:

  A numeric matrix or data frame with two columns: one for the target
  profile mass and one for the control mass.

- time:

  A numeric vector of pseudotime values (must match number of rows in
  \`mass\`).

- thresh:

  Minimum total mass required in both target and control groups for a
  bin to be retained (default: 5).

- binsize:

  Width of pseudotime bins used in \`cut()\` (default: 0.005).

## Value

A logical vector indicating which rows (cells) fall into retained
pseudotime bins.

## Details

Applies a binning strategy over pseudotime and retains only bins where
both target and control profile mass exceed a given threshold. Useful
for identifying regions of sufficient information for differential
analysis.
