## ----setup, include=FALSE-----------------------------------------------------
options(future.globals.maxSize = 4 * 1024^3)

demo_data_path <- Sys.getenv(
  "CLONOTRACE_DEMO_DATA",
  unset = "../../hematopoiesis.rds"
)
demo_data_available <- file.exists(demo_data_path)

knitr::opts_chunk$set(
  eval = demo_data_available,
  message = FALSE,
  warning = FALSE,
  fig.align = "center"
)

cached_demo_figure <- function(filename) {
  path <- system.file("extdata", "demo_files", filename, package = "Clonotrace")
  if (!nzchar(path)) stop("Bundled vignette figure not found: ", filename)
  path
}

## -----------------------------------------------------------------------------
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(Matrix)
# library(future)
# library(future.apply)

## -----------------------------------------------------------------------------
# library(Clonotrace)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("devtools")
# devtools::install_github("yuntianf/Clonotrace")

## -----------------------------------------------------------------------------
# seurat_object = readRDS(demo_data_path)

## -----------------------------------------------------------------------------
# pca = seurat_object@reductions$pca@cell.embeddings
# umap = as.data.frame(seurat_object@reductions$umap@cell.embeddings)
# cell_meta = seurat_object@meta.data

## -----------------------------------------------------------------------------
# table(cell_meta$Time.point)

## -----------------------------------------------------------------------------
# cell_type_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "Cell.type.annotation",size = 0.1,label_size = 5)+
# NoLegend()+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank())+
# xlab("UMAP 1")+ylab("UMAP 2")
# cell_type_umap

## ----cached-cell-type-umap, eval=TRUE, echo=FALSE-----------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_10_0.png"))

## -----------------------------------------------------------------------------
# cell_cluster_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "cluster",size = 0.1,label_size = 5)+
# NoLegend()+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank())+
# xlab("UMAP 1")+ylab("UMAP 2")
# cell_cluster_umap

## ----cached-cell-cluster-umap, eval=TRUE, echo=FALSE--------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_12_0.png"))

## -----------------------------------------------------------------------------
# cell_knn = embedding2knn(embedding = as.matrix(pca),k = 30,mode = "connectivity",if_self = FALSE)
# cell_knn = compute_transition(cell_knn)

## -----------------------------------------------------------------------------
# clone_size = cell_meta %>% group_by(clone) %>% summarise(count = n())
# expanded_clones = clone_size %>% filter(count >=10)

## -----------------------------------------------------------------------------
# clone_size_hist = ggplot(clone_size)+
# geom_histogram(aes(x = sqrt(count)),binwidth = 1,color = "black",fill = "steelblue")+
# geom_vline(xintercept = sqrt(10),linetype = "dashed",color = "coral",linewidth = 1)+
# theme_classic()+
# xlab("sqrt(clone size)")+
# ylab("clone count")
# clone_size_hist

## ----cached-clone-size-histogram, eval=TRUE, echo=FALSE-----------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_18_0.png"))

## -----------------------------------------------------------------------------
# cell_clone = data.frame(cell = rownames(cell_meta),clone = cell_meta[,"clone"]) %>%
#             mutate(clone = if_else(clone %in% expanded_clones$clone, clone, NA_character_))

## -----------------------------------------------------------------------------
# start = Sys.time()
# plan(multisession,workers = 2)
# clone_labels = label_spreading_bootstrap(adj = cell_knn,labels = cell_clone$clone,
#                                          alpha = 0.6, sample_rate = 0.8,
#                                          sample_n = 48, backend = "blocked",
#                                          block_size = 128)
# end = Sys.time()
# print(end-start)

## -----------------------------------------------------------------------------
# cell_clone_prob_raw = clone_labels[[1]]
# deviance = clone_labels[[2]]

## -----------------------------------------------------------------------------
# rownames(cell_clone_prob_raw) = rownames(cell_meta)

## -----------------------------------------------------------------------------
# deviance_hist = ggplot()+
# geom_histogram(aes(x = deviance),fill = "coral",color = "black",binwidth = 0.02)+
# geom_vline(xintercept = 0.3,linetype = "dashed")+
# theme_classic()+
# ylab("cell count")
# deviance_hist

## ----cached-deviance-histogram, eval=TRUE, echo=FALSE-------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_24_0.png"))

## -----------------------------------------------------------------------------
# cell_meta$deviance = deviance

## -----------------------------------------------------------------------------
# deviance_umap = dimplot(umap,annot = cell_meta %>% filter(deviance > 0.3),color_by = "deviance",label = FALSE,size = 0.1)+
# scale_color_viridis_c(name = "deviance",option = "plasma")+
# theme_classic()+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank())+
# ggtitle("Cells with high deviance")
# deviance_umap

## ----cached-deviance-umap, eval=TRUE, echo=FALSE------------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_27_0.png"))

## -----------------------------------------------------------------------------
# rownames(cell_clone_prob_raw) = rownames(cell_meta)

## -----------------------------------------------------------------------------
# cell_clone_prob = cell_clone_prob_raw[deviance < 0.3,]
# cell_clone_prob = cell_clone_prob/rowSums(cell_clone_prob)

## -----------------------------------------------------------------------------
# cell_clone_prob = mat_sparsify(mat = cell_clone_prob,row_mass = 0.9,col_mass = 0.9)
# cell_clone_prob = cell_clone_prob/rowSums(cell_clone_prob)
# cell_clone_prob = Matrix(cell_clone_prob, sparse = TRUE)

## -----------------------------------------------------------------------------
# colnames(cell_clone_prob) = names(table(cell_clone$clone))

## -----------------------------------------------------------------------------
# example_clone = "Lineage-1193-6"

## -----------------------------------------------------------------------------
# clone_example_umap_raw = dimplot(umap,annot = cell_meta %>% filter(clone == example_clone) %>% mutate(group = "1"),
#                                  color_by = "group",label = FALSE,size = 2)+
# theme_classic()+
# NoLegend()+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank())
# 
# clone_example_umap_raw

## ----cached-raw-clone-example, eval=TRUE, echo=FALSE--------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_35_0.png"))

## -----------------------------------------------------------------------------
# cell_meta$prob = NA
# cell_meta[rownames(cell_clone_prob),]$prob = cell_clone_prob[,example_clone]

## -----------------------------------------------------------------------------
# umap$prob = cell_meta$prob

## -----------------------------------------------------------------------------
# clone_example_umap_smooth = ggplot(mapping = aes(x= umap_1,y = umap_2))+
# geom_point(data = umap,color = "lightgrey",size = 0.1)+
# geom_point(data = umap %>% filter(prob > 0.025) %>% arrange(prob),
#            aes(color = prob),size = 2,alpha = 0.5)+
# scale_color_viridis_c(name = "prob",option = "turbo")+
# theme_classic()+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank())
# clone_example_umap_smooth = ggrastr::rasterize(clone_example_umap_smooth,layers='Point', dpi=300)
# clone_example_umap_smooth

## ----cached-smoothed-clone-example, eval=TRUE, echo=FALSE---------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_38_0.png"))

## ----eval = FALSE-------------------------------------------------------------
# # not run
# plan(multisession, workers = 2)
# ### exact OT distance
# clone_dis = clone_disance(
#   as.matrix(pca), cell_clone_prob, outpath = "./",
#   exact = TRUE, exact_backend = "disk", force_exact = TRUE,
#   revised = FALSE, cores = 2, distance_block_size = 256,
#   target_cell_limit = 5000,
#   checkpoint_dir = "./exact_ot_checkpoints"
# )
# ### nearest neighbor distance approximation
# clone_dis = clone_disance(
#   as.matrix(pca), cell_clone_prob, outpath = "./",
#   exact = FALSE, revised = FALSE
# )

## -----------------------------------------------------------------------------
# file_path <- system.file("extdata", "clone_graph_dis.tsv", package = "Clonotrace")
# clone_dis = read.table(file_path,sep = "\t")
# 
# clone_dis = long2square(long = as.data.frame(clone_dis),row_names_from = "group1",
#                         col_names_from = "group2",values_from = "dis",symmetric = T)
# diag(clone_dis) = 0
# rownames(clone_dis) = colnames(clone_dis) = colnames(cell_clone_prob)

## -----------------------------------------------------------------------------
# set.seed(1230)
# clone_cluster = leiden_dis(dismat = clone_dis,k = 20,resolution = 0.5,if_umap = TRUE)
# clone_cluster_center = clone_cluster %>% group_by(cluster) %>% summarise_all(median)

## -----------------------------------------------------------------------------
# clone_cluster$mass = colSums(cell_clone_prob)

## -----------------------------------------------------------------------------
# clone_leiden_umap = ggplot(mapping = aes(x = umap1,y = umap2))+
# geom_point(data = clone_cluster,aes(fill = cluster,size = log2(mass)),pch=21,color = "black",alpha = 0.75)+
# ggrepel::geom_label_repel(data = clone_cluster_center,aes(label = cluster),size = 8)+
# theme_bw()+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank())+
# guides(fill = guide_legend(override.aes = list(size=5)))+
# guides(fill = "none")+
# theme(legend.position = "bottom")+
# xlab("UMAP 1")+ylab("UMAP 2")
# clone_leiden_umap

## ----cached-clone-umap, eval=TRUE, echo=FALSE---------------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_50_0.png"))

## -----------------------------------------------------------------------------
# clone_leiden_dm <- diffusionMap::diffuse(clone_dis,maxdim = 30)
# clone_leiden_dm_coord = as.data.frame(clone_leiden_dm$X)
# colnames(clone_leiden_dm_coord) = paste("dm",1:ncol(clone_leiden_dm_coord),sep = "")
# rownames(clone_leiden_dm_coord) = rownames(clone_cluster)

## -----------------------------------------------------------------------------
# dm_cluster_center = clone_leiden_dm_coord %>% group_by(clone_cluster$cluster) %>% summarise_all(median)
# colnames(dm_cluster_center)[1] = "cluster"

## -----------------------------------------------------------------------------
# clone_leiden_diffmap = ggplot(mapping = aes(x = dm1,y = dm2))+
# geom_point(data = clone_leiden_dm_coord,aes(fill = clone_cluster$cluster,size = log2(clone_cluster$mass)),
#            pch=21,color = "black",alpha = 0.75)+
# ggrepel::geom_label_repel(data = dm_cluster_center,aes(label = cluster),size = 8)+
# theme_bw()+
# theme(axis.text = element_blank(),    # Remove axis numbers (labels)
#       axis.ticks = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank())+
# guides(fill = "none",size = guide_legend(title = "log2(mass)"))+
# theme(legend.position = "bottom")+
# xlab("dm 1")+ylab("dm 2")
# clone_leiden_diffmap

## ----cached-clone-diffusion-map, eval=TRUE, echo=FALSE------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_54_0.png"))

## -----------------------------------------------------------------------------
# clone_cluster = clone_cluster %>% mutate(clone = rownames(clone_cluster)) %>%
#                                   mutate(lineage = substr(clone,start = 1,stop = nchar(clone)-2),
#                                   time = substr(clone,start = nchar(clone),stop = nchar(clone))) %>%
#                                   mutate(time = paste("day",time,sep = ""))

## -----------------------------------------------------------------------------
# clone_link = long2wide(clone_cluster %>% dplyr::select(lineage,time,clone),
#                        row_names_from = "lineage",
#                        col_names_from = "time",
#                        values_from = "clone",
#                        symmetric = FALSE)
# clone_link = na.omit(clone_link)

## -----------------------------------------------------------------------------
# colnames(clone_link) = c("start","end")

## -----------------------------------------------------------------------------
# clone_link_coord = cbind(clone_link,
#                          clone_leiden_dm_coord[clone_link$start,c("dm1","dm2")],
#                          clone_leiden_dm_coord[clone_link$end,c("dm1","dm2")],
#                          clone_cluster[clone_link$start,"cluster"],
#                          clone_cluster[clone_link$end,"cluster"])
# colnames(clone_link_coord) = c("start","end","xdm1","xdm2","ydm1","ydm2","x_cluster","y_cluster")

## -----------------------------------------------------------------------------
# sub_annot = clone_cluster %>% filter(lineage %in% rownames(clone_link_coord))
# sub_dm_coord = clone_leiden_dm_coord[rownames(sub_annot),]

## -----------------------------------------------------------------------------
# clone_leiden_diffmap_arrow = ggplot()+
# geom_point(data = clone_leiden_dm_coord,aes(x = dm1,y = dm2),size = 0.2,color = "lightgrey")+
# geom_point(data = sub_dm_coord,aes(x = dm1,y = dm2,fill = sub_annot$time,size = log2(sub_annot$mass)),shape = 21)+
# geom_segment(data = clone_link_coord,aes(x = xdm1,y = xdm2,xend = ydm1,yend = ydm2),
#              arrow = arrow(length = unit(0.03, "npc")),linewidth = 0.25,alpha = 0.75)+
# scale_fill_viridis_d(name = "time")+
# theme_bw()+
# theme(axis.text = element_blank(),    # Remove axis numbers (labels)
#       axis.ticks = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank())+
# guides(fill = guide_legend(override.aes = list(size=5)),
#        size = guide_legend(title = "log2(mass)"))
# clone_leiden_diffmap_arrow

## ----cached-lineage-arrows, eval=TRUE, echo=FALSE-----------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_61_0.png"))

## -----------------------------------------------------------------------------
# clone_profile = clone_cluster %>% dplyr::select(clone,cluster) %>% mutate(flag = 1)
# clone_profile = long2sparse(long = clone_profile,row_names_from = "clone",col_names_from = "cluster",values_from = "flag")

## -----------------------------------------------------------------------------
# cell_profile_prob = cell_clone_prob %*% clone_profile[colnames(cell_clone_prob),]
# rownames(cell_profile_prob) = rownames(cell_clone_prob)

## -----------------------------------------------------------------------------
# cell_meta$profile = NA
# cell_meta[rownames(cell_profile_prob),]$profile = apply(cell_profile_prob,1,function(x){
#     if(max(x) > 0.5){
#         return(which(x == max(x))[1])
#     }
#     else{
#         return(NA)
#     }
# })
# cell_meta = cell_meta %>% mutate(profile = as.factor(profile))

## -----------------------------------------------------------------------------
# cell_meta$profile_prob = NA
# cell_meta[rownames(cell_profile_prob),]$profile_prob = apply(cell_profile_prob,1,max)

## -----------------------------------------------------------------------------
# cell_profile_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "profile",
#                             size = 0.1,alpha_by = "profile_prob",label = TRUE)+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank())+
# xlab("UMAP 1")+ylab("UMAP 2")+
# guides(alpha = "none",color = guide_legend(override.aes = list(size=5)))
# cell_profile_umap

## ----cached-profile-projection, eval=TRUE, echo=FALSE-------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_70_0.png"))

## ----eval = FALSE-------------------------------------------------------------
# clone_mds = MASS::isoMDS(as.matrix(clone_dis), k=30)
# clone_embedding = clone_mds$points
# colnames(clone_embedding) = paste("mds",1:ncol(clone_embedding),sep = "_")
# rownames(clone_embedding) = rownames(clone_cluster)

## -----------------------------------------------------------------------------
# file_path <- system.file("extdata", "clone_mds.tsv", package = "Clonotrace")
# clone_embedding = read.table(file_path,sep = "\t",header = TRUE)

## -----------------------------------------------------------------------------
# clone_t = clone_dpt(clone_embedding = clone_embedding,cell_meta = cell_meta,
#                     clone_col = "clone",cluster_col = "cluster",start_cluster = "0")

## -----------------------------------------------------------------------------
# clone_cluster$dpt = clone_t

## -----------------------------------------------------------------------------
# clone_dpt_diffmap = ggplot(mapping = aes(x = dm1,y = dm2))+
# geom_point(data = clone_leiden_dm_coord,aes(fill = clone_cluster$dpt,size = log2(clone_cluster$mass)),
#            pch=21,color = "black",alpha = 0.75)+
# theme_bw()+
# theme(axis.text = element_blank(),    # Remove axis numbers (labels)
#       axis.ticks = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank())+
# guides(fill = guide_legend(override.aes = list(size=5)))+
# scale_fill_viridis_c(name = "clone t")+
# guides(size = "none")+
# theme(legend.position = "bottom")+
# xlab("dm 1")+ylab("dm 2")
# clone_dpt_diffmap

## ----cached-clone-pseudotime, eval=TRUE, echo=FALSE---------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_80_0.png"))

## -----------------------------------------------------------------------------
# cell_meta$cell_t = NA
# cell_meta[rownames(cell_clone_prob),]$cell_t = cell_clone_prob[,clone_cluster$clone] %*% clone_cluster$dpt

## -----------------------------------------------------------------------------
# cell_dpt_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "cell_t",
#                             size = 0.1,alpha_by = NULL,label = FALSE)+
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank())+
# xlab("UMAP 1")+ylab("UMAP 2")+
# scale_color_viridis_c(name = "cell t")
# cell_dpt_umap

## ----cached-cell-pseudotime, eval=TRUE, echo=FALSE----------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_83_0.png"))

## -----------------------------------------------------------------------------
# cell_k = 30
# clone_k = 15
# cell_embedding = pca[rownames(cell_clone_prob),]

## -----------------------------------------------------------------------------
# plan(multisession, workers = 2)
# coembed_dis = cell_clone_coembed(
#   cell_embedding, clone_embedding, cell_clone_prob,
#   cell_k = cell_k, clone_k = clone_k
# )

## -----------------------------------------------------------------------------
# # Optionally set RETICULATE_PYTHON before starting R.
# library(reticulate)
# python_path <- Sys.getenv("RETICULATE_PYTHON")
# if (nzchar(python_path)) use_python(python_path, required = TRUE)

## -----------------------------------------------------------------------------
# script_path <- system.file("py", "embedding_from_kNN.py", package = "Clonotrace")
# source_python(script_path)

## -----------------------------------------------------------------------------
# cell_coembed_umap = umap_from_knn(coembed_dis,n_neighbors = 10,seed = 512)

## -----------------------------------------------------------------------------
# coembed_umap_cluster = dimplot(embedding = cell_coembed_umap,annot = cell_meta,alpha_by = NULL,label_size = 8,
#                           color_by = "cluster",size = 0.1,alpha = 1)+
# theme(axis.text = element_blank(),    # Remove axis numbers (labels)
#       axis.ticks = element_blank())+
# xlab("UMAP 1")+ylab("UMAP 2")+
# NoLegend()
# coembed_umap_cluster

## ----cached-coembedding-clusters, eval=TRUE, echo=FALSE-----------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_92_0.png"))

## -----------------------------------------------------------------------------
# coembed_umap_profile = dimplot(embedding = cell_coembed_umap,annot = cell_meta,alpha_by = NULL,label_size = 8,
#                           color_by = "profile",size = 0.1,alpha = 1)+
# theme(axis.text = element_blank(),    # Remove axis numbers (labels)
#       axis.ticks = element_blank())+
# xlab("UMAP 1")+ylab("UMAP 2")+
# NoLegend()
# coembed_umap_profile

## ----cached-coembedding-profiles, eval=TRUE, echo=FALSE-----------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_94_0.png"))

## -----------------------------------------------------------------------------
# enrich = cluster_profile_enrich(cell_profile_prob,cell_meta[rownames(cell_profile_prob),"cluster"],permute_n  = 1000)

## -----------------------------------------------------------------------------
# mass_long = wide2long(enrich[[1]])
# mass_long = mass_long %>% mutate(i = rownames(enrich[[1]])[i],
#                                  j = colnames(enrich[[1]])[j])
# colnames(mass_long) = c("cluster","profile","mass")
# mass_long = mass_long %>% group_by(cluster) %>% mutate(freq = mass/sum(mass))
# mass_long$pvalue = c(t(enrich[[2]]))

## -----------------------------------------------------------------------------
# mass_long$signif = "ns"
# mass_long$signif[mass_long$pvalue< 0.05] = "*"
# mass_long$signif[mass_long$pvalue< 0.01] = "**"
# mass_long$signif[mass_long$pvalue< 0.001] = "***"

## -----------------------------------------------------------------------------
# cell_type_profile_enrich_heatmap = ggplot(mass_long,aes(x = cluster,y = profile))+
# geom_tile(aes(fill = freq),color = "black")+
# scale_fill_viridis_c(option = "turbo")+
# geom_text(data = mass_long %>% filter(signif != "ns"),aes(label = signif),size = 10)+
# theme_classic()+
# xlab("cell cluster")+
# ylab("clone profile")
# cell_type_profile_enrich_heatmap

## ----cached-profile-enrichment, eval=TRUE, echo=FALSE-------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_102_0.png"))

## -----------------------------------------------------------------------------
# exprs = seurat_object@assays$RNA$data

## -----------------------------------------------------------------------------
# start = Sys.time()
# DEG_1_4 = profile_cluster_DEG(profile = "1",cluster = "4",exprs = exprs,cell_meta = cell_meta,cell_profile_prob = cell_profile_prob)
# end = Sys.time()
# print(end - start)

## -----------------------------------------------------------------------------
# DEG_1_4_sig = DEG_1_4$stat %>% filter(padj < 0.05)
# nrow(DEG_1_4_sig)

## -----------------------------------------------------------------------------
# DEG_1_4_sig %>% top_n(abs(cohen),n = 10)

## -----------------------------------------------------------------------------
# test_object = DEG_1_4

## -----------------------------------------------------------------------------
# null_fit = test_object$design_null %*% test_object$coef[3:4,]
# full_fit = test_object$design_null %*% test_object$coef[1:2,]
# rownames(null_fit) = rownames(full_fit) = test_object$cell

## -----------------------------------------------------------------------------
# example_gene = "Dmkn"

## -----------------------------------------------------------------------------
# cell_meta$exprs = exprs[example_gene,rownames(cell_meta)]

## -----------------------------------------------------------------------------
# cell_meta$null_exprs = NA
# cell_meta[test_object$cell,]$null_exprs = null_fit[,example_gene]
# cell_meta$full_exprs = NA
# cell_meta[test_object$cell,]$full_exprs = full_fit[,example_gene]

## -----------------------------------------------------------------------------
# gene_fit_cluster = ggplot(cell_meta[test_object$cell,])+
#         geom_point(aes(x = cell_t,y = exprs,fill = cell_profile_prob[test_object$cell,"1"]),size = 2,shape = 21,color = "black")+
#         geom_line(aes(x = cell_t,y = null_exprs),color = "steelblue",linewidth = 1.5)+
#         geom_line(aes(x = cell_t,y = full_exprs),color = "coral",linewidth = 1.5)+
#         scale_fill_gradient2(name = "profile 1",low = "steelblue",mid = "whitesmoke",high = "coral",midpoint = 0.5,limit = c(0,1))+
#         theme_classic()+
#         xlab("pseudotime")+ylab(paste(example_gene,"exprs",sep = " "))
# gene_fit_cluster

## ----cached-gene-fit, eval=TRUE, echo=FALSE-----------------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_117_0.png"))

## -----------------------------------------------------------------------------
# feature_scatter = dimplot(embedding = umap,annot = t(exprs),color_by = example_gene,label  = FALSE,size = 0.1)+
# scale_color_viridis_c(name = "",option = "plasma")+
# theme(axis.text = element_blank(),    # Remove axis numbers (labels)
#       axis.ticks = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "right")+
# ggtitle(example_gene)+
# xlab("UMAP 1")+ylab("UMAP 2")
# feature_scatter

## ----cached-gene-umap, eval=TRUE, echo=FALSE----------------------------------
if (!demo_data_available) knitr::include_graphics(cached_demo_figure("demo_119_0.png"))

