# dimplot

Plot 2D Embedding with Annotations, Connectivity Edges, and Optional
Labels

## Usage

``` r
dimplot(
  embedding,
  annot,
  color_by,
  alpha_by = NULL,
  connectivity = NULL,
  label = TRUE,
  dims = c(1, 2),
  connectivity_thresh = 0.1,
  label_size = 5,
  label_type = "text",
  label_color = "black",
  box.padding = 0.25,
  point.padding = 1e-06,
  raster_thresh = 10000,
  ...
)
```

## Arguments

- embedding:

  A numeric matrix or data frame of embeddings (rows = cells, columns =
  dimensions).

- annot:

  A data frame of annotations with rownames matching \`embedding\`.

- color_by:

  The name of the column in \`annot\` used to color points.

- alpha_by:

  Optional. The name of the column in \`annot\` used to set point
  transparency.

- connectivity:

  Optional. A square matrix indicating connectivity between groups
  (e.g., clusters).

- label:

  Logical. Whether to add text or label annotations at group centers.

- dims:

  Integer vector of length 2 specifying which dimensions to plot
  (default: c(1, 2)).

- connectivity_thresh:

  Threshold for filtering weak edges in the connectivity matrix
  (default: 0.1).

- label_size:

  Numeric. Size of label text (default: 5).

- label_type:

  Character. One of \`"text"\` or \`"label"\`, determines the type of
  label used (default: "text").

- label_color:

  Color used for label text (default: "black").

- box.padding:

  Padding around label boxes (default: 0.25).

- point.padding:

  Padding between labels and points (default: 1e-6).

- raster_thresh:

  Integer. Threshold above which rasterization is applied to speed up
  rendering (default: 10000).

- ...:

  Additional arguments passed to \`geom_point()\` for the main data
  points.

## Value

A \`ggplot\` object displaying the embedding with color and optional
annotations, connectivity, and labels.

## Details

This function generates a 2D scatter plot of a given embedding (e.g.,
UMAP or PCA) with points colored by annotation, optionally alpha-scaled,
labeled, and overlaid with edges based on a connectivity matrix. It is
designed for use in single-cell or embedding-based visualizations with
optional rasterization support for large datasets.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with dummy UMAP and cluster annotation
umap <- matrix(rnorm(200), ncol = 2)
rownames(umap) <- paste0("Cell", 1:100)
annot <- data.frame(cluster = sample(letters[1:4], 100, TRUE))
rownames(annot) <- rownames(umap)
dimplot(umap, annot, color_by = "cluster")
} # }
```
