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
