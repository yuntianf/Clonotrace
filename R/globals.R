# Definitions for functions resolved dynamically (for example by reticulate)
# and data-mask column names used by dplyr/ggplot2. Keeping these declarations
# together allows `R CMD check` to distinguish intentional non-standard
# evaluation from accidental undeclared symbols.
utils::globalVariables(c(
  "%^%", "Next_Time", "bin", "cluster", "control", "count", "dis", "dist",
  "everything", "from", "i", "i_x", "i_y", "j", "j_x", "j_y", "jaccard",
  "mds_from_knn_py", "name", "node1", "node2", "pcali", "ratio", "target",
  "to", "umap_from_knn_py", "weight", "x"
))
