# Clonotrace 1.2.1

* Accelerated approximate clone distances by resolving clone supports once
  and reusing source-clone shortest-path distances across target clones.
  `graph_clone_nn()` and `clone_disance(exact = FALSE)` accept
  `distance_cache_mb` (default 64 MiB per worker); larger source matrices
  fall back to the existing target-clone blocks. Set it to 0 to disable
  distance caching. The nearest-neighbor statistic is unchanged.
* Made clone-specific local distance rescaling optional and disabled it by
  default. Set `local_rescale = TRUE` in `revise_clone_distance()` or
  `clone_disance()` to recover the 1.2.0 behavior. The unscaled input clone
  geometry is otherwise retained up to a global median normalization.
* Added a disk-backed exact-OT workflow through
  `build_graph_distance_store()` and `graph_clone_ot_disk()`. Cell-graph
  shortest-path distances are written in resumable row blocks, and exact OT
  reads bounded cost blocks instead of retaining the complete distance matrix
  in memory. Source-clone work can be split across a configured `future` plan,
  and validated per-source checkpoints can be reused across interrupted runs.
  This preserves the exact graph-Wasserstein objective; it reduces peak RAM
  and repeated work but still requires quadratic disk space in the number of
  stored endpoint cells.
* Made the disk-backed workflow the implementation used by
  `clone_disance(exact = TRUE)`. The former dense implementation remains
  available as `exact_backend = "memory"` for small legacy reproductions.

# Clonotrace 1.2.0

* Fixed the internal multiresolution Leiden ladder to accept numeric
  membership vectors returned by current igraph releases while retaining
  integer cluster identifiers. This restores automatic centroid-weight
  selection in `revise_clone_distance()` when `clusters` is not supplied.

* Added `revise_clone_distance()`: local density rescaling of the original
  clone distance followed by an optional blend with clone-centroid distance.
  The automatic centroid weight is selected without fate or cell-type labels
  by agreement with multiresolution RNA-cluster-composition differences.
* Made the revised distance the default output of `clone_disance()`. The raw
  exact-OT or approximate distance is retained in `dis_clonotrace`; set
  `revised = FALSE` for the pre-1.2 behavior.
* Added support for observed clone identities, probability-weighted clone
  centroids, within-sample weight selection, fixed-weight ablations, named or
  positional long-form clone identifiers, and auditable selection attributes.
* Made the HDF5-backed blocked implementation the default backend of
  `label_spreading_bootstrap()`. The legacy dense backend remains available as
  `backend = "dense"`, and caller-owned HDF5 reference caches are preserved.
* Added regression tests for revised-distance construction, `clone_disance()`
  integration, blocked-bootstrap dispatch, and HDF5-cache ownership.
* Made the external-data tutorial non-executing during package builds, so clean
  checks do not depend on an unbundled `hematopoiesis.rds` file.

# Clonotrace 1.1.2

* Fixed SNN similarity-to-shortest-path cost conversion. `igraph::set_edge_attr()`
  returns a modified graph and does not mutate the input graph in place, so the
  returned graph is now assigned explicitly. Larger shared-neighbor overlap now
  correctly produces shorter traversal costs through `exp(-overlap)`.
* Versioned cell-graph checkpoints with an explicit edge-cost transform so
  checkpoints created with untransformed SNN similarities are not reused.

# Clonotrace 1.1.1

* Fixed `long2square()` conversion and added explicit validation for malformed,
  duplicate, or conflicting pairwise inputs.
* Corrected clone-distance checkpoint reuse, added parameter/input validation,
  atomic checkpoint writes, and separate exact/approximate result files.
* Added conservative safeguards before large exact clone-distance calculations.
* Removed Leiden clustering's dependence on ambiguous S4 `dist` dispatch and
  preserved clone identifiers in clustering output.
* Made `cell_clone_prob` an explicit argument to `cell_clone_coembed()`.
* Added clean-session regression tests and continuous-integration workflows.

# Clonotrace 1.1.0

* Added memory-efficient blocked label spreading (`label_spreading_blocked`,
  `label_spreading_bootstrap_blocked`) with optional HDF5 streaming.
* Added HDF5 helper utilities: `create_hdf5_matrix`, `normalize_hdf5_rows`.
* Added `clone_dpt` for clone-level diffusion pseudotime inference.
* Added `cell_clone_coembed` for clone-aware cell-cell coembedding.
* Added `soft_cluster_gam_fit` and `profile_cluster_DEG` for pseudotime-aware
  differential expression analysis.
* Added `cluster_profile_enrich` for permutation-based profile enrichment
  testing.
* Added `leiden_embedding` and `leiden_dis` for graph-based clustering.
* Various utility functions for sparse matrix manipulation, kNN graph
  construction, and optimal transport-based clone distances.

# Clonotrace 1.0.0

* Initial release.
