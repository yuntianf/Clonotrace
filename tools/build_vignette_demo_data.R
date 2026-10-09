source_path <- Sys.getenv("CLONOTRACE_VIGNETTE_SOURCE", unset = "")
output_path <- Sys.getenv(
  "CLONOTRACE_VIGNETTE_OUTPUT",
  unset = file.path("inst", "extdata", "vignette_demo_data.rds")
)

if (!nzchar(source_path) || !file.exists(source_path)) {
  stop("Set CLONOTRACE_VIGNETTE_SOURCE to the coSpar Seurat RDS.")
}
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("SeuratObject is required to build the vignette data.")
}

object <- readRDS(source_path)
pca <- SeuratObject::Embeddings(object, reduction = "pca")
pca <- pca[, seq_len(min(15L, ncol(pca))), drop = FALSE]
umap <- SeuratObject::Embeddings(object, reduction = "umap")
metadata <- object[[]]

required <- c(
  "Time.point", "assigned_lineage", "Cell.type.annotation", "paga"
)
missing <- setdiff(required, colnames(metadata))
if (length(missing)) {
  stop("The source object is missing metadata: ", paste(missing, collapse = ", "))
}

cell_meta <- metadata[, required, drop = FALSE]
cell_meta$clone <- paste(cell_meta$assigned_lineage, cell_meta$Time.point, sep = "-")
cell_meta$cluster <- factor(cell_meta$paga)
cell_meta <- cell_meta[, c("Time.point", "Cell.type.annotation", "clone", "cluster")]

common <- Reduce(intersect, list(rownames(pca), rownames(umap), rownames(cell_meta)))
pca <- pca[common, , drop = FALSE]
umap <- as.data.frame(umap[common, , drop = FALSE])
cell_meta <- cell_meta[common, , drop = FALSE]

demo <- list(
  pca = pca,
  umap = umap,
  cell_meta = cell_meta,
  source = list(
    dataset = "Weinreb et al. lineage-traced hematopoiesis",
    doi = "10.1126/science.aaw3381",
    cells = length(common),
    pca_dimensions = ncol(pca)
  )
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(demo, output_path, compress = "xz")
message("Wrote ", normalizePath(output_path), " (", file.info(output_path)$size, " bytes)")
