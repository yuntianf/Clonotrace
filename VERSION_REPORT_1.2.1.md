# Clonotrace 1.2.1 compatibility note

Version 1.2.1 changes one default introduced in version 1.2.0: clone-specific
local rescaling is now optional and disabled by default.

The revised distance returned by `clone_disance(revised = TRUE)` is therefore
constructed from:

1. the raw exact-OT or approximate Clonotrace clone distance, retained up to a
   single global median normalization; and
2. Euclidean distance between observed or probability-weighted clone
   centroids, with an optional fixed weight or a weight selected from
   multiresolution RNA-cluster composition.

Set `local_rescale = TRUE` in `clone_disance()` or
`revise_clone_distance()` to reproduce the version-1.2.0 local-rescaling step.
Because local rescaling can alter clone-neighbour rankings, it should be
reported as a sensitivity analysis rather than treated as a change of units.

For method comparisons, retain the raw Clonotrace distance (`dis_clonotrace`
or `revised = FALSE`) and report the centroid-blended distance as a separate
ablation. Automatic centroid-weight selection is fate-blind but is informed by
transcriptomic cluster structure; it must not be described as independent of
RNA-derived state composition.

## Disk-backed exact OT

Version 1.2.1 also includes the validated exact-OT backend used for the
multi-dataset benchmark. `build_graph_distance_store()` computes cell-graph
shortest-path distances in resumable row blocks and writes them to a row-major
binary store. `graph_clone_ot_disk()` reads bounded source/target blocks from
that store, calls the existing unregularized `clone_2_ot()` solver, balances
source-clone work across a configured `future` plan, and writes validated
per-source checkpoints.

This is a storage and execution improvement, not a new distance definition.
It exactly preserves the graph-Wasserstein objective and the input clone
probability masses. Its main advantages are lower peak RAM, resumability, and
reuse of the expensive shortest-path calculation. The full distance store
still occupies `8 * n_endpoint^2` bytes, and small jobs may run slightly more
slowly than the dense in-memory implementation. Approximate distance therefore
remains the routine default for large exploratory analyses.

## Final validation on 2026-10-08

Validation was repeated from a clean source tarball under Ubuntu 24.04 and R
4.3.3, using one build worker and single-threaded numerical libraries to keep
memory use bounded.

- All ten R source files parsed successfully and every package in `Imports`
  was available.
- `R CMD build --no-manual` completed successfully and rebuilt the vignette.
- `R CMD check --as-cran --no-manual` completed with **0 errors, 0 warnings,
  and 3 notes**. The regression suite reported **56 passes, 0 failures, and one
  CRAN-only smoke-test skip**. The notes identify a new submission/one URL that
  rejects automated checks, the bundled example-data size, and an environment
  time-verification limitation; none indicates a runtime or API failure.
- BiocCheck completed with 3 errors, 6 warnings, and 16 notes. Its blocking
  items concern Bioconductor submission policy rather than package execution:
  the 25 MB source tarball exceeds the Bioconductor size limit, fewer than 80%
  of exported man pages contain runnable examples, and the maintainer email is
  not registered on the Bioconductor support site. It also flags the largely
  non-evaluated external-data vignette, individual files over 5 MB, a concise
  DESCRIPTION, and style recommendations.

The package is therefore validated for the present GitHub release workflow.
A formal Bioconductor submission would require a separate packaging pass to
move or reduce bundled example data, expand runnable examples and vignette
evaluation, and complete the maintainer-registration requirements.

The disk-backed exact-OT addition was independently revalidated on 2026-10-09.
Serial and two-worker calculations matched the dense shortest-path OT result
to numerical tolerance, checkpoint reuse was verified, the complete package
test suite passed, the vignette compiled to HTML, and `R CMD check` completed
with status **OK**.
