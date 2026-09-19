# Clonotrace
<!-- badges: start -->
<!-- badges: end -->

Clonotrace is a package for multimodal analysis of single cell RNA-seq data with lineage information, including lineaged-barcoded cells and immune cells. Clonotrace generates an informative clone embedding to show the clonal dynamics, and extend to pseudotime inference and fate driving gene detection.
<img width="9922" height="7207" alt="main_figure_1_revise2 _cut" src="https://github.com/user-attachments/assets/1daa8fb3-56fa-4ac9-82d3-652d6d9a7077" />

## Installation

The easiest way to get Clonotrace is to install it from Github:

``` r
# install.packages("devtools")
devtools::install_github("yuntianf/Clonotrace",dependencies=TRUE)
```

If you want to visualize the clone-weighted cell embedding, you also need to install several python packages in the `requirements.txt` (https://github.com/yuntianf/Clonotrace/blob/main/requirements.txt), which can be done by:
```
pip install -r requirements.txt
```
in your python3 environment.

The python version of this package can be found in https://github.com/yizhou0201/Clonotrace_python

## Instruction
To view the documetation and tutorial of Clonotrace, you could visit:
https://yuntianf.github.io/Clonotrace/

### The revised clone distance (default)
`clone_disance()` now returns the revised Clonotrace distance by default. The graph distance between clones is locally rescaled, then combined with the distance between clone means. The weight is chosen from the data without cell-type labels or cross-validation (`?revise_clone_distance`). Pass the observed clone labels through `clone_labels = cell_clone`. `revised = FALSE` gives the earlier graph distance, and `revise_clone_distance()` applies the revision to a distance you have already computed.

## Citation

If you use Longcell for published work, please cite our manuscript:

``` r
Fu Y, Mathew D, Wang M, et al. Deciphering Cell Fate and Clonal Dynamics via Integrative Single-Cell Lineage Modeling[J]. bioRxiv, 2025: 2025.09. 01.673503.
```
