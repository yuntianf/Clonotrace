


#' @title connectivity_coord
#' @description Convert Connectivity Matrix into Edge Coordinates for Plotting
#' @details
#'
#' Transforms a connectivity matrix and a coordinate matrix into a data frame containing edge coordinates.
#' This is useful for plotting graph-like structures (e.g., k-NN graphs or SNNs) using `ggplot2::geom_segment`.
#'
#' @param coord A numeric matrix or data frame of coordinates (rows = points, columns = dimensions).
#' @param connectivity A square numeric or sparse matrix representing connectivity between points.
#' @param dims A numeric vector of length 2 specifying which columns of `coord` to use for x and y coordinates.
#'
#' @return A data frame with one row per edge
#'
#' @importFrom dplyr left_join
#' @importFrom Matrix sparseMatrix
#'
connectivity_coord = function(coord,connectivity,dims = c(1,2)){
  if(nrow(coord) != nrow(connectivity)){
    stop("The number of points in the coord and the connectivity doesn't match!")
  }
  if(length(dims) != 2){
    stop("This function is for plot and only 2 dimnesions are allowed.")
  }
  coord = as.data.frame(coord[,dims,drop = FALSE])
  coord$id = 1:nrow(coord)

  diag(connectivity) = 0
  connectivity = as(connectivity, "sparseMatrix")
  connectivity = as.data.frame(summary(connectivity))
  connectivity = left_join(connectivity,coord,by = c("i" = "id"))
  connectivity = left_join(connectivity,coord,by = c("j" = "id"))
  colnames(connectivity)[(ncol(connectivity)-3):ncol(connectivity)] = c("i_x","i_y","j_x","j_y")

  return(connectivity)
}


#' @title dimplot
#' @description Plot 2D Embedding with Annotations, Connectivity Edges, and Optional Labels
#' @details
#'
#' This function generates a 2D scatter plot of a given embedding (e.g., UMAP or PCA)
#' with points colored by annotation, optionally alpha-scaled, labeled, and overlaid
#' with edges based on a connectivity matrix. It is designed for use in single-cell
#' or embedding-based visualizations with optional rasterization support for large datasets.
#'
#' @param embedding A numeric matrix or data frame of embeddings (rows = cells, columns = dimensions).
#' @param annot A data frame of annotations with rownames matching `embedding`.
#' @param color_by The name of the column in `annot` used to color points.
#' @param alpha_by Optional. The name of the column in `annot` used to set point transparency.
#' @param connectivity Optional. A square matrix indicating connectivity between groups (e.g., clusters).
#' @param label Logical. Whether to add text or label annotations at group centers.
#' @param dims Integer vector of length 2 specifying which dimensions to plot (default: c(1, 2)).
#' @param connectivity_thresh Threshold for filtering weak edges in the connectivity matrix (default: 0.1).
#' @param label_size Numeric. Size of label text (default: 5).
#' @param label_type Character. One of `"text"` or `"label"`, determines the type of label used (default: "text").
#' @param label_color Color used for label text (default: "black").
#' @param box.padding Padding around label boxes (default: 0.25).
#' @param point.padding Padding between labels and points (default: 1e-6).
#' @param raster_thresh Integer. Threshold above which rasterization is applied to speed up rendering (default: 10000).
#' @param ... Additional arguments passed to `geom_point()` for the main data points.
#'
#' @return A `ggplot` object displaying the embedding with color and optional annotations, connectivity, and labels.
#'
#' @import ggplot2
#' @importFrom dplyr group_by_at summarise filter_at across filter
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom ggrastr rasterize
#'
#' @examples
#' \dontrun{
#' # Example with dummy UMAP and cluster annotation
#' umap <- matrix(rnorm(200), ncol = 2)
#' rownames(umap) <- paste0("Cell", 1:100)
#' annot <- data.frame(cluster = sample(letters[1:4], 100, TRUE))
#' rownames(annot) <- rownames(umap)
#' dimplot(umap, annot, color_by = "cluster")
#' }
#'
#' @export

dimplot = function(embedding, annot, color_by, alpha_by = NULL,connectivity = NULL, label = TRUE,
                   dims = c(1, 2),connectivity_thresh = 0.1,label_size = 5,label_type = "text",label_color = "black",
                   box.padding = 0.25,point.padding = 1e-6,
                   raster_thresh = 10000,...) {
  embedding = as.data.frame(embedding[, dims, drop = FALSE])
  annot = as.data.frame(annot[,c(color_by,alpha_by),drop = FALSE])
  embedding = cbind(embedding, annot[rownames(embedding), c(color_by,alpha_by),drop = FALSE])

  aes_base <- aes(
    x = !!sym(colnames(embedding)[1]),
    y = !!sym(colnames(embedding)[2]),
    color = !!sym(color_by)
  )
  aes_full <- if (!is.null(alpha_by)) {
    modifyList(aes_base, aes(alpha = !!sym(alpha_by)))
  } else {
    aes_base
  }


  out = ggplot() +
    geom_point(data = embedding,
               mapping = aes(x = !!sym(colnames(embedding)[1]), y = !!sym(colnames(embedding)[2])),
               color = "lightgrey",size = 0.1,alpha = 0.5)+
    geom_point(data = embedding %>% filter_at(color_by,~!is.na(.)),
               mapping = aes_full, ...) +
    theme_classic() +
    theme(text = element_text(size = 15))#+guides(colour = guide_legend(override.aes = list(size=5)))

  if(nrow(embedding) > raster_thresh){
    out = rasterize(out, layers='Point', dpi=300)
  }

  if(label | !is.null(connectivity)){
    center_coord = embedding %>%
      group_by_at(color_by) %>%
      summarise(across(c(colnames(embedding)[1:2]), ~ median(.x, na.rm = TRUE)), count = n())
    center_coord = na.omit(center_coord)

    if (!is.null(connectivity)) {
      edge_coord = connectivity_coord(center_coord, connectivity,dims = c(2,3))
      edge_coord = edge_coord %>% filter(x >= connectivity_thresh)
      out = out +
        geom_segment(data = edge_coord,
                     aes(x = i_x, y = i_y, xend = j_x, yend = j_y, linewidth = x), color = "Honeydew",alpha = 0.75)+
        geom_point(data = center_coord,
                   aes(x = !!sym(colnames(embedding)[1]), y = !!sym(colnames(embedding)[2]), size = log(count)))
    }

    if (label) {
      label_type = match.arg(label_type,choices = c("text","label"))

      if(label_type == "text"){
        out = out +
          ggrepel::geom_text_repel(data = center_coord,
                                   mapping = aes(x = !!sym(colnames(embedding)[1]), y = !!sym(colnames(embedding)[2]), label = !!sym(color_by)),
                                   size = label_size,point.padding = point.padding,box.padding = box.padding,
                                   color = label_color)
      }
      else{
        out = out +
          ggrepel::geom_label_repel(data = center_coord,
                                   mapping = aes(x = !!sym(colnames(embedding)[1]), y = !!sym(colnames(embedding)[2]), label = !!sym(color_by)),
                                   size = label_size,point.padding = point.padding,box.padding = box.padding,
                                   color = label_color)
      }
    }
  }

  return(out)
}


#' @title scatterpie
#' @description Plot Cluster Compositions as Scatterpie Chart with Connectivity Edges
#' @details
#'
#' Creates a `ggplot2` layer list showing pie charts at spatial or embedding coordinates, where each pie represents
#' the relative composition of categories (e.g., cell types, fates) per cluster or region. Optionally overlays connectivity
#' edges between clusters based on a provided adjacency matrix.
#'
#' @param scatter_coord A data frame of coordinates (e.g., UMAP or spatial centroids), one row per cluster.
#' @param composition A numeric matrix or data frame (same number of rows as `scatter_coord`) where each row gives the proportions or counts of categories.
#' @param connectivity Optional. A square matrix indicating edge weights (e.g., similarity or shared connectivity) between clusters.
#' @param connectivity_thresh Numeric threshold above which edges are shown (default: 0.5).
#' @param dims A character vector of length 2 specifying the coordinate column names in `scatter_coord` (default: `c("umap_1", "umap_2")`).
#' @param cluster_col Character. The column in `scatter_coord` that contains cluster or group labels for annotation (default: `"cluster"`).
#' @param edge_color Color for the connectivity edges (default: `"lightgrey"`).
#' @param edge_alpha Alpha transparency for the edges (default: 1).
#' @param label_size Numeric size for the cluster label text (default: 5).
#'
#' @return A list of `ggplot2` layers (can be added to a base ggplot using `+`) including:
#' \describe{
#'   \item{geom_segment}{(optional) edges between connected clusters}
#'   \item{geom_scatterpie}{pie chart layer showing cluster composition}
#'   \item{geom_text_repel}{cluster labels}
#'   \item{theme_classic}{plot styling}
#' }
#'
#' @importFrom ggplot2 geom_segment aes theme_classic
#' @importFrom ggrepel geom_text_repel
#' @importFrom scatterpie geom_scatterpie
#' @importFrom dplyr filter
#'
#' @examples
#' \dontrun{
#' library(scatterpie)
#' coords <- data.frame(umap_1 = rnorm(5), umap_2 = rnorm(5), cluster = letters[1:5])
#' comp <- matrix(runif(5 * 3), nrow = 5)
#' colnames(comp) <- c("A", "B", "C")
#' comp <- comp / rowSums(comp)
#' scatterpie(coords, comp)
#' }
#'
#' @export
scatterpie <- function(scatter_coord, composition,
                       connectivity = NULL, connectivity_thresh = 0.5,
                       dims = c("umap_1", "umap_2"), cluster_col = "cluster",
                       edge_color= "lightgrey",edge_alpha = 1,
                       label_size = 5) {
  if (nrow(scatter_coord) != nrow(composition)) {
    stop("The number of elements of the coordinates and the compositions should match!")
  }

  layers <- list()

  if (!is.null(connectivity)) {
    edge_coord <- connectivity_coord(scatter_coord, connectivity, dims = c(2, 3))
    edge_coord <- edge_coord %>% filter(x >= connectivity_thresh)
    layers[[length(layers) + 1]] <- geom_segment(
      data = edge_coord,
      aes(x = i_x, y = i_y, xend = j_x, yend = j_y, linewidth = x),
      color = edge_color,alpha = edge_alpha
    )
  }

  coord_comp <- as.data.frame(cbind(scatter_coord, composition))
  coord_comp$count <- rowSums(composition)

  layers[[length(layers) + 1]] <- geom_scatterpie(
    data = coord_comp,
    aes(x = !!sym(dims[1]), y = !!sym(dims[2]), r = log(count) / 12),
    cols = colnames(composition)
  )

  layers[[length(layers) + 1]] <- ggrepel::geom_text_repel(
    data = scatter_coord,
    aes(x = !!sym(dims[1]), y = !!sym(dims[2]), label = !!sym(cluster_col)),
    size = label_size
  )

  layers[[length(layers) + 1]] <- theme_classic()

  return(layers)
}

