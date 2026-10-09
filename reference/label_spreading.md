# label_spreading

Label Propagation via Iterative Graph-Based Spreading

## Usage

``` r
label_spreading(
  adj,
  labels,
  label_n = NULL,
  alpha = 0.9,
  max_iter = 100,
  tol = 0.001,
  epsilon = 0,
  verbose = TRUE
)
```

## Arguments

- adj:

  A square adjacency matrix (preferably sparse) representing the graph.

- labels:

  An integer vector of length equal to the number of nodes. Use \`NA\`
  for unlabeled entries.

- label_n:

  Optional. The number of label classes. If \`NULL\`, inferred as
  \`max(labels, na.rm = TRUE)\`.

- alpha:

  Float in (0, 1). The propagation coefficient controlling the balance
  between prior and propagated labels (default: 0.9).

- max_iter:

  Maximum number of iterations for propagation (default: 100).

- tol:

  Convergence threshold (default: 1e-3).

- epsilon:

  Small prior assigned to unlabeled entries (default: 0).

- verbose:

  Logical. Whether to print progress and convergence status (default:
  TRUE).

## Value

A matrix of size N x C where each row contains soft label probabilities
for a node across \`C\` classes.

## Details

Performs semi-supervised label propagation over a graph using a sparse
or dense adjacency matrix. The algorithm propagates known labels across
the graph structure, allowing soft label assignment for unlabeled nodes.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
adj <- Matrix::rsparsematrix(100, 100, density = 0.05)
labels <- rep(NA, 100)
labels[1:10] <- sample(1:3, 10, replace = TRUE)
prob_matrix <- label_spreading(adj, labels)
} # }
```
