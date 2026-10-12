# Clonotrace

Clonotrace is a package for multimodal analysis of single cell RNA-seq
data with lineage information, including lineaged-barcoded cells and
immune cells. Clonotrace generates an informative clone embedding to
show the clonal dynamics, and extend to pseudotime inference and fate
driving gene detection. ![main_figure_1_revise2
\_cut](https://github.com/user-attachments/assets/1daa8fb3-56fa-4ac9-82d3-652d6d9a7077)

## Installation

Clonotrace uses packages from both CRAN and Bioconductor. Install it
through `BiocManager`, which resolves Bioconductor dependencies such as
`rhdf5` and `bluster` in a compatible release:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

BiocManager::install(c("rhdf5", "bluster"))
remotes::install_github("yuntianf/Clonotrace", dependencies = TRUE)
```

To verify the installation in a clean R session:

``` r

library(Clonotrace)
BiocManager::valid()
```

The Python projection of the optional clone-weighted cell embedding
requires the packages in `requirements.txt`. Install them into an
environment you manage, then point `reticulate` to that environment if
it is not selected automatically:

``` shell
pip install -r requirements.txt
```

The core Clonotrace workflow does not require Python.

## Version 1.2.1 highlights

- [`clone_disance()`](https://yuntianf.github.io/Clonotrace/reference/clone_disance.md)
  returns a revised clone distance by default: the original Clonotrace
  distance, globally median-normalized, is blended with clone-centroid
  distance at a fate-blind, multiresolution RNA-structure-selected
  weight. Clone-specific local rescaling is optional and disabled by
  default in version 1.2.1 because it can change neighbour rankings; set
  `local_rescale = TRUE` only for a documented sensitivity analysis or
  to reproduce version 1.2.0. The raw distance remains available in
  `dis_clonotrace` or via `revised = FALSE`.
- [`label_spreading_bootstrap()`](https://yuntianf.github.io/Clonotrace/reference/label_spreading_bootstrap.md)
  uses the memory-conscious HDF5-backed blocked implementation by
  default; use `backend = "dense"` only for legacy reproduction or small
  inputs.
- `clone_disance(exact = TRUE)` now uses the disk-backed exact-OT
  workflow by default. Its reusable components are
  [`build_graph_distance_store()`](https://yuntianf.github.io/Clonotrace/reference/build_graph_distance_store.md)
  and
  [`graph_clone_ot_disk()`](https://yuntianf.github.io/Clonotrace/reference/graph_clone_ot_disk.md).
  The workflow computes the same unregularized graph-Wasserstein
  distance as the in-memory implementation while replacing the full
  in-memory cell-distance matrix with a resumable binary store and
  bounded reads. The tradeoff is quadratic disk space in the stored
  clone-bearing cells; approximate distance remains the recommended
  default for routine large datasets. Set `exact_backend = "memory"`
  only to reproduce the former dense implementation on a small input.

See [NEWS.md](https://yuntianf.github.io/Clonotrace/NEWS.md),
[VERSION_REPORT_1.2.0.md](https://yuntianf.github.io/Clonotrace/VERSION_REPORT_1.2.0.md),
and
[VERSION_REPORT_1.2.1.md](https://yuntianf.github.io/Clonotrace/VERSION_REPORT_1.2.1.md)
for the complete change and compatibility reports.

## Documentation

The documentation and tutorial are available at:
<https://yuntianf.github.io/Clonotrace/>

The tutorial recommends approximate clone-distance calculation for
routine analyses and documents both validated checkpoints and the
disk-backed exact-OT workflow for reference calculations. A compiled
copy is included at
[`vignettes/Clonotrace.html`](https://yuntianf.github.io/Clonotrace/vignettes/Clonotrace.md);
the GitHub Pages site can be built into `docs/` with
[`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html).

## Citation

If you use Clonotrace for published work, please cite our manuscript:

``` r
Fu Y, Mathew D, Wang M, et al. Deciphering Cell Fate and Clonal Dynamics via Integrative Single-Cell Lineage Modeling[J]. bioRxiv, 2025: 2025.09. 01.673503.
```
