# Clonotrace
<!-- badges: start -->
[![R package checks](https://github.com/yuntianf/Clonotrace/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yuntianf/Clonotrace/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Clonotrace is a package for multimodal analysis of single cell RNA-seq data with lineage information, including lineaged-barcoded cells and immune cells. Clonotrace generates an informative clone embedding to show the clonal dynamics, and extend to pseudotime inference and fate driving gene detection.
<img width="9922" height="7207" alt="main_figure_1_revise2 _cut" src="https://github.com/user-attachments/assets/1daa8fb3-56fa-4ac9-82d3-652d6d9a7077" />

## Installation

Clonotrace uses packages from both CRAN and Bioconductor. Install it through
`BiocManager`, which resolves Bioconductor dependencies such as `rhdf5` and
`bluster` in a compatible release:

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

The Python projection of the optional clone-weighted cell embedding requires
the packages in `requirements.txt`. Install them into an environment you
manage, then point `reticulate` to that environment if it is not selected
automatically:

``` shell
pip install -r requirements.txt
```

The core Clonotrace workflow does not require Python.

## Version 1.2.1 highlights

- `clone_disance()` returns a revised clone distance by default: the original
  Clonotrace distance, globally median-normalized, is blended with
  clone-centroid distance at a fate-blind, multiresolution
  RNA-structure-selected weight. Clone-specific local rescaling is optional
  and disabled by default in version 1.2.1 because it can change neighbour
  rankings; set `local_rescale = TRUE` only for a documented sensitivity
  analysis or to reproduce version 1.2.0. The raw distance remains available
  in `dis_clonotrace` or via `revised = FALSE`.
- `label_spreading_bootstrap()` uses the memory-conscious HDF5-backed blocked
  implementation by default; use `backend = "dense"` only for legacy
  reproduction or small inputs.
- `clone_disance(exact = TRUE)` now uses the disk-backed exact-OT workflow by
  default. Its reusable components are `build_graph_distance_store()` and
  `graph_clone_ot_disk()`. The workflow
  computes the same unregularized graph-Wasserstein distance as the in-memory
  implementation while replacing the full in-memory cell-distance matrix with
  a resumable binary store and bounded reads. The tradeoff is quadratic disk
  space in the stored clone-bearing cells; approximate distance remains the
  recommended default for routine large datasets. Set
  `exact_backend = "memory"` only to reproduce the former dense implementation
  on a small input.

See [NEWS.md](NEWS.md), [VERSION_REPORT_1.2.0.md](VERSION_REPORT_1.2.0.md), and
[VERSION_REPORT_1.2.1.md](VERSION_REPORT_1.2.1.md) for the complete change and
compatibility reports.

## Documentation

The documentation and tutorial are available at:
https://yuntianf.github.io/Clonotrace/

The tutorial recommends approximate clone-distance calculation for routine
analyses and documents both validated checkpoints and the disk-backed exact-OT
workflow for reference calculations. A compiled copy is included at
[`vignettes/Clonotrace.html`](vignettes/Clonotrace.html); the GitHub Pages site
can be built into `docs/` with `pkgdown::build_site()`.

## ClonoView

ClonoView is an interactive viewer for a Clonotrace run, and it comes with the package —
installing Clonotrace installs it. Two panels side by side, the clone embedding on the
left and the single-cell expression embedding on the right, with the colouring, the
profile filter and the timepoint acting on both at once.

### Try it first

A real run ships with the package, so you can see the viewer working before preparing
anything of your own:

``` r
library(Clonotrace)
clonoview_demo()          # bundles the example, serves it, opens a browser
clonoview_stop()          # when you are done
```

That takes a few seconds. The example is the NSCLC anti-PD-1 plus JAK1-inhibitor trial of
Mathew et al. — 8 patients over 5 treatment cycles, 150,150 cells, 916 clone units, 16
cell types. It is not a stored screenshot or a pre-built result: `clonoview_demo()` runs
the same `clonoview_bundle()` you would run on your own data, so if it works, the whole
path works. `clonoview_example()` returns the four objects it uses, which is the quickest
way to see the shape your own inputs need.

### Your own run

Bundle it once, then view it:

``` r
clonoview_bundle(cell_embedding   = umap,
                 cell_meta        = cell_meta,
                 clone_cluster    = clone_cluster,      # from leiden_dis(if_umap = TRUE)
                 cell_profile_prob = cell_profile_prob,
                 name             = "my_run",
                 cell_type_column = "Cell.type.annotation")
clonoview("my_run")     # opens in a browser; clonoview_stop() when done
```

Those four arguments are the four things the viewer draws, and nothing else is required.
A longitudinal study adds its structure through `clone_id_regex`, which says how to read
the patient, timepoint and receptor out of a clone identifier:

``` r
clonoview_bundle(..., clone_id_regex = "^(?P<patient>[^|]+)\\|(?P<visit>[^|]+)\\|(?P<receptor>.+)$",
                 patient_column = "patient", visit_column = "visit")
```

### The tables

Beside the two panels there is a clone table and a sample table, one row per clone and
per sample, with the cell counts and the sampling schedule. **Click any column heading to
sort by it, and again to reverse** — blanks sink to the bottom either way, so an empty
cell never displaces real data. **Click a row to show it, shift-click to add more**, and
up to eight hand-picked clones each get their own **marker shape** on the clone map, shown
beside the row so the table and the plot name the clone the same way. Shape is used rather
than colour because colour already means the cell-state profile, in that panel as in every
other; a square and a circle of the same colour still agree about the profile and disagree
about the clone.

### Notes

The viewer serves on the loopback address only. If R is running on a cluster, forward the
port from your own machine with the `ssh -L` command `clonoview()` prints, then open the
URL there.

ClonoView needs `matplotlib` in addition to the packages in `requirements.txt`, and uses
the same Python interpreter as the co-embedding step — so a `reticulate::use_python()`
call covers both. See `?clonoview_bundle` for the full set of options.

## Citation

If you use Clonotrace for published work, please cite our manuscript:

``` r
Fu Y, Mathew D, Wang M, et al. Deciphering Cell Fate and Clonal Dynamics via Integrative Single-Cell Lineage Modeling[J]. bioRxiv, 2025: 2025.09. 01.673503.
```
