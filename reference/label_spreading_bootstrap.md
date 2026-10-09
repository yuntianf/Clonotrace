# Bootstrap-Based Stability Estimation for Label Propagation

Memory-conscious bootstrap stability estimation for label spreading

## Usage

``` r
label_spreading_bootstrap(
  adj,
  labels,
  refer = NULL,
  alpha = 0.8,
  sample_rate = 0.8,
  sample_n = 50,
  backend = c("blocked", "dense"),
  block_size = 128,
  tol = 0.001,
  max_iter = 100,
  epsilon = 0,
  refer_h5 = NULL,
  tmpdir = tempdir(),
  verbose = TRUE,
  ...
)
```

## Arguments

- adj:

  A square adjacency matrix (preferably sparse) representing the graph.

- labels:

  An integer vector of length equal to the number of nodes. Use \`NA\`
  for unlabeled entries.

- refer:

  Optional dense reference soft-label matrix (N x C). This legacy
  argument is supported only by `backend = "dense"`. If supplied while
  `backend` is omitted, the function selects the dense backend with a
  warning.

- alpha:

  Float in (0, 1). The propagation coefficient controlling the balance
  between prior and propagated labels (default: 0.9).

- sample_rate:

  Fraction of labeled nodes used in each bootstrap sample (default:
  0.8).

- sample_n:

  Number of bootstrap replicates (default: 50).

- backend:

  Bootstrap implementation. The default is `"blocked"`; use `"dense"`
  only for small inputs or legacy reproduction.

- block_size:

  Number of label classes processed in one HDF5-backed block.

- tol:

  Convergence threshold (default: 1e-3).

- max_iter:

  Maximum number of iterations for propagation (default: 100).

- epsilon:

  Small prior assigned to unlabeled entries (default: 0).

- refer_h5:

  Optional reusable HDF5 reference for the blocked backend. A
  caller-supplied file is never deleted by this function.

- tmpdir:

  Directory for temporary HDF5 files used by the blocked backend.

- verbose:

  Logical. Whether to print progress and convergence status (default:
  TRUE).

- ...:

  Additional arguments passed to
  [`label_spreading()`](https://yuntianf.github.io/Clonotrace/reference/label_spreading.md)
  when `backend = "dense"`.

## Value

A list with `prob`, the full-reference probability matrix, and
`deviance`, the per-cell mean L1 bootstrap deviation.

## Details

Repeatedly applies label spreading to subsampled label sets and returns
the full-label probability estimate together with a node-level bootstrap
deviance. The default `backend = "blocked"` processes label classes in
bounded blocks and streams intermediate matrices through HDF5. This
avoids retaining a full cells-by-label matrix for every concurrent
bootstrap.

The blocked backend still returns one dense cells-by-label probability
matrix for compatibility. Its peak memory therefore includes that final
matrix, the sparse graph, and one class block per concurrent worker.
Reduce the number of future workers when either the cell count or clone
count is large. The former dense implementation remains available with
`backend = "dense"` for small inputs and reproducibility checks.

## Examples

``` r
if (FALSE) { # \dontrun{
adj <- Matrix::rsparsematrix(100, 100, density = 0.05)
labels <- rep(NA, 100)
labels[1:10] <- sample(1:3, 10, replace = TRUE)
result <- label_spreading_bootstrap(adj, labels, block_size = 32)
str(result)
} # }
```
