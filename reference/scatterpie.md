# scatterpie

Plot Cluster Compositions as Scatterpie Chart with Connectivity Edges

## Usage

``` r
scatterpie(
  scatter_coord,
  composition,
  connectivity = NULL,
  connectivity_thresh = 0.5,
  dims = c("umap_1", "umap_2"),
  cluster_col = "cluster",
  edge_color = "lightgrey",
  edge_alpha = 1,
  label_size = 5
)
```

## Arguments

- scatter_coord:

  A data frame of coordinates (e.g., UMAP or spatial centroids), one row
  per cluster.

- composition:

  A numeric matrix or data frame (same number of rows as
  \`scatter_coord\`) where each row gives the proportions or counts of
  categories.

- connectivity:

  Optional. A square matrix indicating edge weights (e.g., similarity or
  shared connectivity) between clusters.

- connectivity_thresh:

  Numeric threshold above which edges are shown (default: 0.5).

- dims:

  A character vector of length 2 specifying the coordinate column names
  in \`scatter_coord\` (default: \`c("umap_1", "umap_2")\`).

- cluster_col:

  Character. The column in \`scatter_coord\` that contains cluster or
  group labels for annotation (default: \`"cluster"\`).

- edge_color:

  Color for the connectivity edges (default: \`"lightgrey"\`).

- edge_alpha:

  Alpha transparency for the edges (default: 1).

- label_size:

  Numeric size for the cluster label text (default: 5).

## Value

A list of \`ggplot2\` layers (can be added to a base ggplot using \`+\`)
including:

- geom_segment:

  (optional) edges between connected clusters

- geom_scatterpie:

  pie chart layer showing cluster composition

- geom_text_repel:

  cluster labels

- theme_classic:

  plot styling

## Details

Creates a \`ggplot2\` layer list showing pie charts at spatial or
embedding coordinates, where each pie represents the relative
composition of categories (e.g., cell types, fates) per cluster or
region. Optionally overlays connectivity edges between clusters based on
a provided adjacency matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scatterpie)
coords <- data.frame(umap_1 = rnorm(5), umap_2 = rnorm(5), cluster = letters[1:5])
comp <- matrix(runif(5 * 3), nrow = 5)
colnames(comp) <- c("A", "B", "C")
comp <- comp / rowSums(comp)
scatterpie(coords, comp)
} # }
```
