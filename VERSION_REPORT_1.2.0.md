# Clonotrace 1.2.0 integration report

Date: 2026-09-21  
Package source:
`C:/Users/yunti/OneDrive/Documentos/Clonotrace revision/Clonotrace-revised`  
Built artifact:
`C:/Users/yunti/OneDrive/Documentos/Clonotrace revision/Clonotrace_1.2.0.tar.gz`

## Release summary

Version 1.2.0 consolidates the methodological and engineering revisions
made during the benchmark work into the main Clonotrace package. The two
principal changes are:

1.  the revised clone distance is now an audited, reusable package
    method and is the default result returned by
    [`clone_disance()`](https://yuntianf.github.io/Clonotrace/reference/clone_disance.md);
    and
2.  the memory-conscious blocked/HDF5 bootstrap is now the default
    [`label_spreading_bootstrap()`](https://yuntianf.github.io/Clonotrace/reference/label_spreading_bootstrap.md)
    backend.

The unrevised Clonotrace distance and the legacy dense bootstrap remain
available explicitly, so old results can still be reproduced and the
revised distance can be evaluated against a raw-distance ablation.

## Improvements accumulated to date

| Version | Area | Main improvement |
|----|----|----|
| 1.0.0 | Initial method | Clone-level modeling of lineage-resolved single-cell RNA data. |
| 1.1.0 | Scale and downstream analysis | Blocked label spreading, HDF5 helpers, clone diffusion pseudotime, cell-clone coembedding, profile enrichment, and pseudotime-aware differential expression. |
| 1.1.1 | Correctness and reproducibility | Validated distance conversion, robust checkpoints, exact/approximate separation, exact-OT safety limits, stable clone identifiers, explicit coembedding inputs, and clean-session tests/CI. |
| 1.1.2 | Cell-graph metric | Corrected the SNN similarity-to-cost conversion used by shortest paths and versioned the cell-graph cache. |
| 1.2.0 | Clone-distance revision and bootstrap default | Added local distance rescaling plus an optional clone-centroid component; made the revised distance and blocked bootstrap the defaults while retaining raw/legacy modes. |

### 1. Correct cell-graph distance

The cell graph is constructed as an SNN graph, whose edge weights are
shared neighbor similarities. Shortest-path routines interpret larger
edge weights as longer costs, so the similarity must first be converted
to a cost. The package now assigns the transformed graph returned by
[`igraph::set_edge_attr()`](https://r.igraph.org/reference/set_edge_attr.html)
and uses

``` text
edge cost = exp(-shared-neighbor similarity)
```

Thus, stronger local transcriptomic similarity gives a shorter traversal
cost. The graph-checkpoint metadata records this transform, preventing
an older graph with untransformed similarities from being silently
reused.

### 2. Safer exact and approximate clone distances

- Exact OT and approximate nearest-neighbor distances have separate
  checkpoint files: `clone_graph_dis_exact.rds` and
  `clone_graph_dis.rds`.
- Checkpoints carry input and parameter metadata and are written
  atomically. An incompatible cache now produces an actionable error
  rather than being silently reused.
- Exact mode has conservative cell- and clone-count guards. A large run
  must be acknowledged explicitly with `force_exact = TRUE`.
- Disconnected clone-bearing graph components and non-finite distances
  are detected explicitly.
- [`long2square()`](https://yuntianf.github.io/Clonotrace/reference/long2square.md)
  now rejects missing, duplicate, or conflicting clone pairs.
- Approximate graph-distance calculation supports clone blocking to
  limit peak memory.

These changes make exact OT usable as a reference calculation without
treating it as the routine large-dataset default.

### 3. Revised clone distance

The new exported function
[`revise_clone_distance()`](https://yuntianf.github.io/Clonotrace/reference/revise_clone_distance.md)
integrates the `revised-distance` branch into the main package. Starting
with a raw Clonotrace clone-distance matrix, it performs three steps.

#### Local rescaling

For clone `i`, let `sigma_i` be its distance to the `local_k`-th nearest
clone. The locally rescaled distance is

``` text
D_local(i,j) = D_raw(i,j) / sqrt(sigma_i * sigma_j)
```

This reduces distortions caused by unequal clone density in different
parts of the clone manifold.

#### Clone-centroid distance

The method calculates Euclidean distance between clone centroids in the
input RNA embedding. Exact observed clone labels are preferred;
otherwise the cell-by-clone probability matrix supplies weights.

#### Fate-blind blend selection

After pairwise-median normalization, the final distance is

``` text
D_revised = (1 - s) * D_local + s * D_centroid
```

If `s` is not fixed by the user, Clonotrace selects it by maximizing
Spearman agreement with clone RNA-cluster-composition differences over
multiple Leiden resolutions, then uses the median selected weight. This
selection does not use fate labels or benchmark outcomes, but it is
RNA-cluster-informed. Therefore, benchmark reports should retain raw,
local-rescaled, and full-revision ablations where relevant.

Additional integration details:

- [`clone_disance()`](https://yuntianf.github.io/Clonotrace/reference/clone_disance.md)
  now returns the revised distance in `dis` by default.
- The raw exact/approximate result is retained in `dis_clonotrace`.
- `revised = FALSE` reproduces the pre-1.2 distance behavior.
- Named and positional clone identifiers are supported.
- A named clone-to-sample vector can restrict normalization and
  automatic weight selection to within-sample pairs.
- A fixed `mean_weight` supports controlled ablations.
- Selected weights, per-resolution choices, selection curves, raw
  distance, and local scales are retained as attributes for auditing.
- Expensive raw OT/NN checkpoints are independent of the inexpensive
  revision parameters, so tuning the revision does not recompute the raw
  distance.

### 4. Blocked label-propagation bootstrap is now the default

[`label_spreading_bootstrap()`](https://yuntianf.github.io/Clonotrace/reference/label_spreading_bootstrap.md)
now dispatches to the HDF5-backed blocked implementation unless
`backend = "dense"` is requested explicitly. The blocked implementation:

- fixes the global label-to-column mapping once, so bootstrap subsets
  cannot reorder clone classes;
- processes clone classes in bounded blocks;
- streams bootstrap/reference comparisons through HDF5;
- preserves an HDF5 reference cache supplied by the caller;
- removes only internally owned temporary files; and
- guards deviance normalization when a cell has no eligible replicate.

This resolves the main avoidable source of bootstrap peak memory. It
does not make the whole workflow disk-only: the final cells-by-clones
probability matrix is still returned densely for compatibility. Worker
count and block size must therefore still be selected conservatively.

### 5. Downstream and interface improvements

- [`cell_clone_coembed()`](https://yuntianf.github.io/Clonotrace/reference/cell_clone_coembed.md)
  requires its cell-by-clone probability input explicitly, avoiding
  hidden dependence on a global variable.
- Leiden helpers avoid ambiguous S4 `dist` dispatch and preserve clone
  names.
- Clone-level diffusion pseudotime is available through
  [`clone_dpt()`](https://yuntianf.github.io/Clonotrace/reference/clone_dpt.md).
- Profile enrichment and pseudotime-aware profile
  differential-expression helpers are included.
- Documentation now distinguishes clone-profile identity from cell-type
  annotation and records the new defaults and their compatibility
  switches.
- Cross-platform R package checks and isolated-install validation are
  defined in GitHub Actions.

## Compatibility notes

### Intentional default changes

1.  `clone_disance(..., revised = TRUE)` is now the default. Use
    `revised = FALSE` for an unrevised exact-OT or approximate distance.
2.  `label_spreading_bootstrap(..., backend = "blocked")` is now the
    default. Use `backend = "dense"` only for legacy reproduction or
    small inputs.

The misspelled public function name
[`clone_disance()`](https://yuntianf.github.io/Clonotrace/reference/clone_disance.md)
is retained to avoid an API-breaking rename in this release.

### Recommended benchmark reporting

For a method-comparison figure, report at least:

- raw Clonotrace distance;
- locally rescaled distance; and
- full revised distance.

Automatic blend selection is fate-blind, but because it uses
multiresolution RNA structure it should not be described as completely
independent of RNA cluster organization.

## Validation performed for this integration

| Check | Result |
|----|----|
| Parse package R sources and base regression file | Passed |
| Build `Clonotrace_1.2.0.tar.gz` without rebuilding the vignette | Passed |
| Install the tarball into an isolated WSL R library | Passed |
| Base-R regression suite | Passed |
| Automatic internal Leiden ladder with current igraph membership type | Passed after adding explicit integer coercion |
| CARLIN real-data revised-distance application (310 clones; 47,895 pairs) | Passed; finite symmetric output, selected centroid weight 0.10 |
| Standard `R CMD check --no-manual --no-build-vignettes` | 0 errors, 2 warnings; both warnings only reflect the intentionally skipped vignette build and missing `inst/doc` output |
| Running code extracted from the vignette | Passed after external-data chunks were marked non-executing |
| `testthat` suite | Not run locally because `testthat` is not installed in the selected WSL R environment; CI is configured to install and run it |

The base-R regression suite covers the corrected SNN cost,
long-to-square validation, cache validation and atomic writes,
clone-identifier preservation, raw approximate and exact distance
workflows, revised-distance construction, the internal automatic Leiden
ladder used for centroid-weight selection, and the blocked-bootstrap
default.

A full local vignette render was not possible because Pandoc is absent
from the selected WSL environment. The GitHub Actions workflow installs
Pandoc before building and checking the package, so the final repository
check should render the vignette and remove the two local `inst/doc`
warnings.

## Deliberately not included in Clonotrace 1.2.0

- Dataset-specific benchmark scoring, diagnostic figures, and competitor
  adapters remain analysis code and belong in a separate `ClonoBench`
  package if they are to be generalized.
- Experimental disk-backed exact-OT/graph-flow prototypes remain outside
  this release until their API, numerical equivalence, and failure
  recovery are validated.
- Simulation-specific fate fields and benchmark ground-truth
  construction are not part of the Clonotrace package.

This boundary keeps Clonotrace focused on the reusable biological method
while allowing benchmark machinery to evolve independently.

## Files most relevant to the release

- `R/clone_dis.R`: corrected graph costs, validated checkpoints,
  exact/approx workflows, and revised-distance dispatch.
- `R/revised_distance.R`: local rescaling, centroid distance, and
  automatic or fixed blend selection.
- `R/label_propagation.R`: blocked-default bootstrap and legacy dense
  backend.
- `NEWS.md`: concise user-facing version history.
- `vignettes/Clonotrace.Rmd`: workflow examples and compatibility
  guidance.
- `tests/regression-base.R` and `tests/testthat/`: release regression
  tests.
