# clone_partition

Partition Clones Based on Containment Similarity With Centers From
Farthest Point Sampling

## Usage

``` r
clone_partition(clone_matrix, k = 10, similarity_threshold = 0)
```

## Arguments

- clone_matrix:

  A cells × clones matrix, where non-zero entries indicate clone
  membership.

- k:

  Integer. Number of clone partitions (default: 10).

- similarity_threshold:

  Numeric threshold for assigning a clone to a group based on similarity
  (default: 0).

## Value

A named list of clone ID vectors, one per partition.

## Details

Groups clones into \`k\` roughly balanced partitions based on
containment similarity. Useful for parallelizing clone-level
computations (e.g., OT) by minimizing inter-group redundancy.
