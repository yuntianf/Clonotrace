```R
options(future.globals.maxSize = 4 * 1024^3)
```


```R
library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)
```

    Loading required package: SeuratObject
    
    Loading required package: sp
    
    
    Attaching package: ‘SeuratObject’
    
    
    The following object is masked from ‘package:base’:
    
        intersect
    
    
    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    



```R
library(Clonotrace)
```

### data

We first read in data. The demo data is from a public in vitro longitudinal experiment to monitor the hematopoiesis (https://www.science.org/doi/10.1126/science.aaw3381). You can download the demo data from: https://upenn.box.com/s/9bxv50lueelrrf5zv4ked20dlkalrtir


```R
seurat_object = readRDS("../../hematopoiesis.rds")
```


```R
pca = seurat_object@reductions$pca@cell.embeddings
umap = as.data.frame(seurat_object@reductions$umap@cell.embeddings)
cell_meta = seurat_object@meta.data
```

Within this data, there are 34782 cells are sampled at 3 time points.


```R
table(cell_meta$Time.point)
```


    
        2     4     6 
     3407 11101 20274 


This data has a simple cell type compistion with only 3 cell types:


```R
cell_type_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "Cell.type.annotation",size = 0.1,label_size = 5)+
NoLegend()+
theme(axis.text = element_blank(),  
      axis.ticks = element_blank())+
xlab("UMAP 1")+ylab("UMAP 2")
cell_type_umap
```


    
![png](demo_files/demo_10_0.png)
    


For a higher resolustion in the downstream analysis, we further cluster cells into 9 cell clusters.


```R
cell_cluster_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "cluster",size = 0.1,label_size = 5)+
NoLegend()+
theme(axis.text = element_blank(),  
      axis.ticks = element_blank())+
xlab("UMAP 1")+ylab("UMAP 2")
cell_cluster_umap
```


    
![png](demo_files/demo_12_0.png)
    


### clone label spreading

We first build a cell kNN graph based within the PCA space:


```R
cell_knn = embedding2knn(embedding = as.matrix(pca),k = 30,mode = "connectivity",if_self = FALSE)
cell_knn = compute_transition(cell_knn)
```

Here we define a clone as a population of cells which share the same lineage barcode at the same time point. Based on this definition there are 8108 clones in total. 802 out of them expanded over 10 cells.


```R
clone_size = cell_meta %>% group_by(clone) %>% summarise(count = n())
expanded_clones = clone_size %>% filter(count >=10)
```


```R
clone_size_hist = ggplot(clone_size)+
geom_histogram(aes(x = sqrt(count)),binwidth = 1,color = "black",fill = "steelblue")+
geom_vline(xintercept = sqrt(10),linetype = "dashed",color = "coral",linewidth = 1)+
theme_classic()+
xlab("sqrt(clone size)")+
ylab("clone count")
clone_size_hist
```


    
![png](demo_files/demo_18_0.png)
    


Then we could spread the clone label in the cell graph to smooth the cell density in those expanding clones. There is a bootstrap during this process to filter out propogation with high deviance, usually in regions with high density of lowly expanded or unexpanded clones. So it's better to run this in parallel.


```R
cell_clone = data.frame(cell = rownames(cell_meta),clone = cell_meta[,"clone"]) %>% 
            mutate(clone = if_else(clone %in% expand_clones$clone, clone, NA_character_))
```


```R
start = Sys.time()
plan(multisession,workers = 8)
clone_labels = label_spreading_bootstrap(adj = cell_knn,labels = cell_clone$clone,
                                              alpha = 0.6, sample_rate = 0.8,sample_n = 48)
end = Sys.time()
print(end-start)
```

    Time difference of 6.385496 mins



```R
cell_clone_prob_raw = clone_labels[[1]]
deviance = clone_labels[[2]]
```


```R
rownames(cell_clone_prob_raw) = rownames(cell_meta)
```


```R
deviance_hist = ggplot()+
geom_histogram(aes(x = deviance),fill = "coral",color = "black",binwidth = 0.02)+
geom_vline(xintercept = 0.3,linetype = "dashed")+
theme_classic()+
ylab("cell count")
deviance_hist
```


    
![png](demo_files/demo_24_0.png)
    


We can visualize the deviance of clone label spreading in each cell via cell UMAP, and here we use the 0.3 as the upper bound threshold for deviance. Cells with high deviance are enriched in undifferentiated cells, which due to they haven't proliferated enough to expand their clones.


```R
cell_meta$deviance = deviance
```


```R
deviance_umap = dimplot(umap,annot = cell_meta %>% filter(deviance > 0.3),color_by = "deviance",label = FALSE,size = 0.1)+
scale_color_viridis_c(name = "deviance",option = "plasma")+
theme_classic()+
theme(axis.text = element_blank(),  
      axis.ticks = element_blank())+
ggtitle("Cells with high deviance")
deviance_umap
```


    
![png](demo_files/demo_27_0.png)
    


Then we filter out cells with high deviance, and sparsify the cell-to-clone probablity by keeping top 90% mass


```R
rownames(cell_clone_prob) = rownames(cell_meta)
```

      [[ suppressing 802 column names ‘Lineage-100-4’, ‘Lineage-1042-4’, ‘Lineage-1042-6’ ... ]]
    



    6 x 802 sparse Matrix of class "dgCMatrix"
                                                                                  
    [1,] . .         . .          . . .           . . . . . . . 0.00531143 . . . .
    [2,] . .         . .          . . .           . . . . . . . .          . . . .
    [3,] . .         . .          . . .           . . . . . . . .          . . . .
    [4,] . .         . .          . . .           . . . . . . . .          . . . .
    [5,] . .         . 0.00476413 . . 0.004796612 . . . . . . . .          . . . .
    [6,] . 0.0036061 . .          . . .           . . . . . . . 0.58765214 . . . .
                                                                                  
    [1,] . . . .           . 0.008138689 .          . . .           . . .         
    [2,] . . . .           . .           .          . . .           . . .         
    [3,] . . . .           . .           .          . . .           . . .         
    [4,] . . . .           . .           .          . . .           . . 0.00755937
    [5,] . . . 0.004912297 . .           0.02908699 . . .           . . .         
    [6,] . . . .           . .           .          . . 0.008031904 . . .         
                                                                               
    [1,] . .           . . .           . . . . .           . . . . . .         
    [2,] . .           . . .           . . . . .           . . . . . .         
    [3,] . .           . . .           . . . . .           . . . . . .         
    [4,] . .           . . .           . . . . .           . . . . . .         
    [5,] . .           . . 0.004988477 . . . . 0.006318303 . . . . . .         
    [6,] . 0.004231581 . . .           . . . . .           . . . . . 0.01170913
                                                                             
    [1,] .           . .           . . . .          . .           0.006685797
    [2,] .           . .           . . . 0.01207639 . .           .          
    [3,] .           . .           . . . .          . .           .          
    [4,] .           . .           . . . .          . .           .          
    [5,] 0.005420687 . 0.003751568 . . . .          . .           .          
    [6,] 0.006794085 . .           . . . .          . 0.007406298 .          
                                                                                  
    [1,] 0.007820217 . .           .          . . . . . . . . . . .          . . .
    [2,] .           . .           .          . . . . . . . . . . .          . . .
    [3,] .           . .           .          . . . . . . . . . . 0.03908448 . . .
    [4,] .           . .           0.01761085 . . . . . . . . . . .          . . .
    [5,] .           . 0.007496687 .          . . . . . . . . . . .          . . .
    [6,] .           . .           .          . . . . . . . . . . .          . . .
                                                                            
    [1,] .          . . . . 0.00770244 .           . . . . . . .           .
    [2,] 0.01065861 . . . . .          .           . . . . . . 0.019026429 .
    [3,] .          . . . . .          .           . . . . . . .           .
    [4,] .          . . . . .          0.006045595 . . . . . . .           .
    [5,] .          . . . . .          .           . . . . . . 0.006644126 .
    [6,] .          . . . . .          .           . . . . . . .           .
                                                                                  
    [1,] .          . . . .           .          . . .         . . . . . . . . . .
    [2,] .          . . . 0.004488529 .          . . .         . . . . . . . . . .
    [3,] .          . . . .           .          . . .         . . . . . . . . . .
    [4,] .          . . . .           .          . . 0.0154877 . . . . . . . . . .
    [5,] 0.01381945 . . . .           0.01197759 . . .         . . . . . . . . . .
    [6,] .          . . . .           .          . . .         . . . . . . . . . .
                                                                                
    [1,] . . . . . . . .          . . .          . .           .           . . .
    [2,] . . . . . . . .          . . .          . 0.004286043 .           . . .
    [3,] . . . . . . . .          . . .          . .           .           . . .
    [4,] . . . . . . . .          . . .          . .           .           . . .
    [5,] . . . . . . . 0.01642447 . . 0.01200701 . .           0.006009263 . . .
    [6,] . . . . . . . .          . . .          . .           .           . . .
                                                                                   
    [1,] .           . . . . . .           . . .          .           . . . . . . .
    [2,] .           . . . . . .           . . .          .           . . . . . . .
    [3,] .           . . . . . .           . . .          .           . . . . . . .
    [4,] .           . . . . . .           . . .          .           . . . . . . .
    [5,] 0.005034679 . . . . . 0.006497433 . . 0.01225052 0.008139744 . . . . . . .
    [6,] .           . . . . . .           . . .          .           . . . . . . .
                                                                              
    [1,] . . . .          . . . . . .          .          . . . .          . .
    [2,] . . . 0.01790191 . . . . . .          .          . . . .          . .
    [3,] . . . .          . . . . . .          .          . . . .          . .
    [4,] . . . .          . . . . . .          .          . . . .          . .
    [5,] . . . .          . . . . . 0.01562914 0.01593801 . . . .          . .
    [6,] . . . .          . . . . . .          .          . . . 0.01392533 . .
                                                                             
    [1,] .          . . .          .           . . . .          . . . . . . .
    [2,] .          . . .          .           . . . .          . . . . . . .
    [3,] .          . . .          .           . . . .          . . . . . . .
    [4,] .          . . .          .           . . . 0.01476958 . . . . . . .
    [5,] 0.01254928 . . 0.01202797 0.007228893 . . . .          . . . . . . .
    [6,] .          . . .          .           . . . .          . . . . . . .
                                                                                  
    [1,] .          .          . . . . .          . . . . 0.004697769 .          .
    [2,] .          .          . . . . .          . . . . .           .          .
    [3,] .          0.04163761 . . . . .          . . . . .           .          .
    [4,] .          .          . . . . .          . . . . .           .          .
    [5,] 0.01350831 .          . . . . 0.01178876 . . . . 0.004671211 0.02783541 .
    [6,] .          .          . . . . .          . . . . .           .          .
                                                                                 
    [1,] . . .           .         . . . .          .           . . . .          
    [2,] . . .           .         . . . .          .           . . . .          
    [3,] . . .           0.0265662 . . . .          .           . . . .          
    [4,] . . .           .         . . . .          .           . . . .          
    [5,] . . 0.007127276 .         . . . 0.00833106 0.004019925 . . . 0.003711077
    [6,] . . .           .         . . . .          .           . . . .          
                                                                                   
    [1,] .          . .           . . . 0.01738502 . . . . . . . . .          . . .
    [2,] .          . .           . . . .          . . . . . . . . .          . . .
    [3,] .          . .           . . . .          . . . . . . . . .          . . .
    [4,] .          . .           . . . .          . . . . . . . . .          . . .
    [5,] 0.01283954 . .           . . . .          . . . . . . . . 0.01674855 . . .
    [6,] .          . 0.004744235 . . . .          . . . . . . . . .          . . .
                                                                                 
    [1,] .           .           . .           . .           . . . . . .         
    [2,] .           .           . .           . .           . . . . . .         
    [3,] .           .           . 0.004843621 . .           . . . . . .         
    [4,] .           .           . .           . .           . . . . . .         
    [5,] 0.009825711 0.007512789 . .           . 0.004162588 . . . . . 0.01263067
    [6,] .           .           . .           . .           . . . . . .         
                                                                          
    [1,] .           0.01006134 . .           .           . .          . .
    [2,] .           .          . 0.008418974 0.003872038 . 0.02328671 . .
    [3,] .           .          . .           .           . .          . .
    [4,] .           .          . .           .           . 0.01586028 . .
    [5,] 0.008257124 .          . .           .           . .          . .
    [6,] .           .          . .           .           . .          . .
                                                                                 
    [1,] .           .           . . 0.005635643 . . . .           . .          .
    [2,] 0.008700377 0.007854785 . . .           . . . .           . .          .
    [3,] .           .           . . .           . . . .           . 0.07413213 .
    [4,] .           .           . . .           . . . .           . .          .
    [5,] 0.006230622 .           . . .           . . . 0.009704865 . .          .
    [6,] .           .           . . .           . . . 0.004697517 . .          .
                                                                                
    [1,] 0.007300794 . . . . .          . 0.004484434 . . . . . . .          . .
    [2,] .           . . . . .          . .           . . . . . . .          . .
    [3,] .           . . . . 0.00439089 . .           . . . . . . .          . .
    [4,] .           . . . . .          . .           . . . . . . 0.01261964 . .
    [5,] .           . . . . .          . .           . . . . . . .          . .
    [6,] .           . . . . .          . 0.010150794 . . . . . . .          . .
                                                                                  
    [1,] .           . .           . . 0.009052554 . . .           . . . . . . . .
    [2,] .           . .           . . .           . . 0.004859006 . . . . . . . .
    [3,] .           . .           . . .           . . .           . . . . . . . .
    [4,] .           . .           . . .           . . .           . . . . . . . .
    [5,] 0.006626202 . 0.005760232 . . .           . . .           . . . . . . . .
    [6,] .           . .           . . .           . . .           . . . . . . . .
                                                                                  
    [1,] .          . . .           . . . . .           .          .           . .
    [2,] .          . . .           . . . . .           .          .           . .
    [3,] .          . . 0.004799538 . . . . .           .          .           . .
    [4,] .          . . .           . . . . .           .          .           . .
    [5,] 0.02008033 . . .           . . . . 0.006032500 0.01654425 .           . .
    [6,] .          . . .           . . . . 0.003962506 .          0.009633478 . .
                                                                           
    [1,] . . . . . 0.008453072 . . 0.010149159 . . . .           .         
    [2,] . . . . . .           . . .           . . . 0.004282881 .         
    [3,] . . . . . .           . . .           . . . .           0.01039887
    [4,] . . . . . .           . . .           . . . .           .         
    [5,] . . . . . .           . . .           . . . .           .         
    [6,] . . . . . .           . . 0.007447468 . . . .           .         
                                                                              
    [1,] .           .           .          .          . . . . . . .          
    [2,] .           .           .          .          . . . . . . .          
    [3,] .           .           0.03301905 0.01197947 . . . . . . 0.004073551
    [4,] .           .           .          .          . . . . . . .          
    [5,] 0.008662538 0.007303984 .          .          . . . . . . .          
    [6,] .           .           .          .          . . . . . . 0.003926391
                                                                                   
    [1,] .          . . . .           . . . .           .          . . . .         
    [2,] .          . . . 0.004311561 . . . 0.005945972 0.01396755 . . . .         
    [3,] 0.04970113 . . . .           . . . .           .          . . . .         
    [4,] .          . . . 0.003645155 . . . .           .          . . . .         
    [5,] .          . . . .           . . . .           .          . . . 0.01660348
    [6,] .          . . . .           . . . .           .          . . . .         
                                                                            
    [1,] . . . . . . . . . .           . . . . . .          .          . . .
    [2,] . . . . . . . . . 0.009461547 . . . . . .          0.01654894 . . .
    [3,] . . . . . . . . . .           . . . . . .          .          . . .
    [4,] . . . . . . . . . .           . . . . . .          .          . . .
    [5,] . . . . . . . . . .           . . . . . 0.01770614 .          . . .
    [6,] . . . . . . . . . 0.007478691 . . . . . .          .          . . .
                                                                                  
    [1,] .          . . . . . . .           .           . . . .          . . . . .
    [2,] .          . . . . . . .           .           . . . .          . . . . .
    [3,] 0.01387176 . . . . . . 0.008256035 .           . . . .          . . . . .
    [4,] .          . . . . . . .           .           . . . .          . . . . .
    [5,] .          . . . . . . .           0.005338825 . . . 0.01003074 . . . . .
    [6,] .          . . . . . . .           .           . . . .          . . . . .
                                                                               
    [1,] .          .          . . . .          . .           . . . .          
    [2,] .          .          . . . .          . .           . . . 0.006320271
    [3,] .          .          . . . 0.06470884 . .           . . . .          
    [4,] .          .          . . . .          . .           . . . .          
    [5,] 0.01498045 .          . . . .          . 0.005609292 . . . .          
    [6,] .          0.01006873 . . . .          . .           . . . .          
                                                                             
    [1,] .          .          .          .           . .          .         
    [2,] .          0.01949232 .          0.004978577 . .          .         
    [3,] 0.01826152 .          .          .           . 0.03315438 .         
    [4,] .          .          .          .           . .          .         
    [5,] .          .          0.01007581 .           . .          0.02198379
    [6,] .          .          .          .           . .          .         
                                                                                  
    [1,] .          .           . .           . . . . . 0.01087907 . 0.005693848 .
    [2,] .          .           . .           . . . . . .          . .           .
    [3,] 0.04269007 .           . .           . . . . . .          . .           .
    [4,] .          .           . .           . . . . . .          . .           .
    [5,] .          0.006074536 . 0.009062067 . . . . . .          . .           .
    [6,] .          .           . 0.023330702 . . . . . .          . .           .
                                                                            
    [1,] .          . . .           . . 0.006867066 . .          .          
    [2,] .          . . .           . . .           . .          0.008242525
    [3,] 0.02250997 . . .           . . .           . .          .          
    [4,] .          . . .           . . .           . .          .          
    [5,] .          . . 0.007671263 . . .           . 0.02362616 .          
    [6,] .          . . .           . . .           . .          .          
                                                                             
    [1,] .           .          . .          . . . . . . . . . . .          .
    [2,] .           .          . 0.02811358 . . . . . . . . . . .          .
    [3,] 0.005839535 0.01975016 . .          . . . . . . . . . . .          .
    [4,] .           .          . .          . . . . . . . . . . .          .
    [5,] .           .          . .          . . . . . . . . . . .          .
    [6,] 0.006320266 .          . .          . . . . . . . . . . 0.01404636 .
                                                                                  
    [1,] .         . . .           . . .           . .           .           . . .
    [2,] .         . . 0.008693621 . . .           . .           0.005874137 . . .
    [3,] 0.0142553 . . .           . . 0.007282505 . .           .           . . .
    [4,] .         . . .           . . .           . .           .           . . .
    [5,] .         . . .           . . .           . .           0.004083009 . . .
    [6,] .         . . .           . . .           . 0.008829159 .           . . .
                                                                          
    [1,] . .          .           . . . .           .          0.006068631
    [2,] . .          .           . . . 0.536427760 .          .          
    [3,] . .          .           . . . .           .          .          
    [4,] . .          0.008635409 . . . 0.007007594 0.02621428 .          
    [5,] . 0.01670763 .           . . . .           .          .          
    [6,] . .          .           . . . .           .          .          
                                                                                  
    [1,] .           . . . .           .           . . . .           . .          
    [2,] .           . . . .           .           . . . .           . .          
    [3,] .           . . . .           .           . . . .           . .          
    [4,] .           . . . .           .           . . . .           . .          
    [5,] 0.004333847 . . . .           .           . . . .           . .          
    [6,] .           . . . 0.005734107 0.007855125 . . . 0.008403258 . 0.008617288
                                                                               
    [1,] .           . . . . . . .           . . . . . .         . . .         
    [2,] 0.009459307 . . . . . . .           . . . . . .         . . .         
    [3,] 0.009618967 . . . . . . 0.006098846 . . . . . .         . . 0.02219138
    [4,] 0.004367379 . . . . . . .           . . . . . .         . . .         
    [5,] .           . . . . . . .           . . . . . 0.0156496 . . .         
    [6,] .           . . . . . . .           . . . . . .         . . .         
                                                                           
    [1,] .           .           . .           . . .           . .         
    [2,] .           .           . .           . . .           . .         
    [3,] .           0.037537577 . .           . . .           . .         
    [4,] .           0.009566316 . .           . . .           . .         
    [5,] .           .           . .           . . .           . .         
    [6,] 0.008563702 .           . 0.008744105 . . 0.007512619 . 0.01519996
                                                                           
    [1,] 0.0152339 . .           .           . .          . .           . .
    [2,] .         . .           .           . .          . .           . .
    [3,] .         . .           .           . .          . 0.004946337 . .
    [4,] .         . .           0.005147847 . .          . .           . .
    [5,] .         . 0.004423866 .           . 0.01285646 . .           . .
    [6,] 0.0111323 . .           .           . .          . .           . .
                                                                                
    [1,] .           . . .          . . . . . . . . . . .          . 0.007340281
    [2,] 0.005972213 . . .          . . . . . . . . . . .          . .          
    [3,] .           . . .          . . . . . . . . . . 0.09689098 . .          
    [4,] .           . . .          . . . . . . . . . . .          . .          
    [5,] .           . . .          . . . . . . . . . . .          . .          
    [6,] .           . . 0.01013857 . . . . . . . . . . .          . .          
                                                                                
    [1,] .           .           . . . . . . . . . . .          .          . . .
    [2,] .           0.005196752 . . . . . . . . . . .          .          . . .
    [3,] .           .           . . . . . . . . . . .          .          . . .
    [4,] .           .           . . . . . . . . . . .          .          . . .
    [5,] 0.003744214 .           . . . . . . . . . . 0.01910867 0.01304968 . . .
    [6,] .           .           . . . . . . . . . . .          .          . . .
                                                                              
    [1,] .           . . . .           .           . . . . . . .           . .
    [2,] .           . . . .           0.006517798 . . . . . . .           . .
    [3,] .           . . . .           .           . . . . . . .           . .
    [4,] .           . . . .           .           . . . . . . .           . .
    [5,] 0.004470203 . . . .           .           . . . . . . 0.008451937 . .
    [6,] .           . . . 0.003724789 .           . . . . . . .           . .
                                                                               
    [1,] .          .          . .          .          . . . .          . . . .
    [2,] .          .          . .          .          . . . .          . . . .
    [3,] .          0.00821102 . .          .          . . . .          . . . .
    [4,] .          .          . .          .          . . . .          . . . .
    [5,] 0.01925678 .          . 0.03017141 0.01009571 . . . .          . . . .
    [6,] .          .          . .          .          . . . 0.01048881 . . . .
                                                                                
    [1,] .           . . .          . . .          . .          . .          . .
    [2,] 0.003720756 . . .          . . .          . .          . 0.02157850 . .
    [3,] .           . . .          . . .          . .          . .          . .
    [4,] .           . . 0.00613788 . . .          . .          . .          . .
    [5,] .           . . .          . . .          . 0.01083527 . 0.00499292 . .
    [6,] .           . . .          . . 0.01292093 . .          . .          . .
                                                                               
    [1,] .           . 0.01067208 . . .          .           . . .          . .
    [2,] .           . .          . . .          0.005193962 . . .          . .
    [3,] .           . .          . . .          .           . . .          . .
    [4,] .           . .          . . .          .           . . .          . .
    [5,] .           . .          . . 0.01399901 .           . . 0.02639005 . .
    [6,] 0.009388794 . .          . . .          .           . . .          . .
                                                                             
    [1,] .           . . .           . . . . . . .           .          . . .
    [2,] .           . . .           . . . . . . .           .          . . .
    [3,] .           . . 0.008619635 . . . . . . .           .          . . .
    [4,] 0.005160685 . . .           . . . . . . .           .          . . .
    [5,] .           . . .           . . . . . . .           0.01270369 . . .
    [6,] .           . . .           . . . . . . 0.006894525 .          . . .
                                                                                   
    [1,] .           .          .           .           . .           .         . .
    [2,] .           .          .           .           . 0.009277199 0.0158077 . .
    [3,] .           0.01142365 .           .           . .           .         . .
    [4,] .           .          .           0.009287344 . .           .         . .
    [5,] 0.006549919 .          .           .           . .           .         . .
    [6,] .           .          0.007480508 .           . .           .         . .
                                                                         
    [1,] . .           .           .          . . . . . . .           . .
    [2,] . .           .           .          . . . . . . .           . .
    [3,] . .           .           .          . . . . . . 0.008064689 . .
    [4,] . 0.004889377 0.008989139 0.01467654 . . . . . . .           . .
    [5,] . .           0.007050021 .          . . . . . . .           . .
    [6,] . .           .           .          . . . . . . .           . .
                                                                                 
    [1,] .           . . . . . . . .           . 0.01301946 . . . .           . .
    [2,] .           . . . . . . . .           . .          . . . .           . .
    [3,] .           . . . . . . . .           . .          . . . .           . .
    [4,] 0.004895602 . . . . . . . 0.007407242 . .          . . . .           . .
    [5,] .           . . . . . . . .           . .          . . . .           . .
    [6,] .           . . . . . . . .           . .          . . . 0.009963157 . .
                                                                                   
    [1,] .           . . .           . . .          . . . . . .           . . . . .
    [2,] .           . . .           . . .          . . . . . .           . . . . .
    [3,] .           . . 0.004672837 . . 0.07022464 . . . . . .           . . . . .
    [4,] .           . . .           . . .          . . . . . .           . . . . .
    [5,] 0.004618227 . . .           . . .          . . . . . 0.006221447 . . . . .
    [6,] .           . . .           . . .          . . . . . .           . . . . .
                                                                                  
    [1,] . . . . . . . . 0.005132438 . .          .          . .           . . . .
    [2,] . . . . . . . . .           . 0.01688715 .          . .           . . . .
    [3,] . . . . . . . . 0.004537015 . .          0.01270395 . .           . . . .
    [4,] . . . . . . . . .           . .          .          . .           . . . .
    [5,] . . . . . . . . .           . .          .          . .           . . . .
    [6,] . . . . . . . . .           . .          .          . 0.006906538 . . . .
                                                                              
    [1,] . . . . . . .           .          . . . .          . .           . .
    [2,] . . . . . . .           .          . . . .          . .           . .
    [3,] . . . . . . .           .          . . . 0.05124961 . .           . .
    [4,] . . . . . . 0.004524288 .          . . . .          . .           . .
    [5,] . . . . . . .           .          . . . .          . .           . .
    [6,] . . . . . . .           0.00566217 . . . .          . 0.004999911 . .
                                                                        
    [1,] .          . .         . . .         .           .          . .
    [2,] .          . .         . . .         0.011383620 .          . .
    [3,] .          . .         . . 0.0212285 .           .          . .
    [4,] .          . 0.7501666 . . .         .           .          . .
    [5,] 0.01447602 . .         . . .         0.006928652 0.02381266 . .
    [6,] .          . .         . . .         .           .          . .
                                                                           
    [1,] .           . . . 0.004607344 .          0.7965106 .           . .
    [2,] .           . . . .           .          .         .           . .
    [3,] .           . . . .           .          .         .           . .
    [4,] .           . . . .           .          .         .           . .
    [5,] .           . . . 0.009767009 0.01072186 .         0.006172673 . .
    [6,] 0.006842767 . . . 0.005836133 .          .         .           . .
                                                                                   
    [1,] .          . . . . .          .           .           .          . . . . .
    [2,] .          . . . . 0.02999197 .           0.004815177 .          . . . . .
    [3,] .          . . . . .          .           .           .          . . . . .
    [4,] .          . . . . .          .           .           .          . . . . .
    [5,] .          . . . . .          0.005683789 0.010187096 .          . . . . .
    [6,] 0.00700588 . . . . .          .           .           0.01296426 . . . . .
                                                                          
    [1,] .          . . .          . . 0.005096926 . . . . . .           .
    [2,] 0.01075147 . . .          . . .           . . . . . .           .
    [3,] .          . . .          . . .           . . . . . .           .
    [4,] .          . . .          . . .           . . . . . .           .
    [5,] .          . . 0.01010212 . . .           . . . . . .           .
    [6,] .          . . .          . . 0.010909978 . . . . . 0.003597401 .
                                                                          
    [1,] .           . .           . . . .           . . . . .           .
    [2,] .           . 0.006526198 . . . .           . . . . .           .
    [3,] .           . .           . . . 0.008523518 . . . . .           .
    [4,] .           . .           . . . .           . . . . 0.004646291 .
    [5,] .           . .           . . . .           . . . . .           .
    [6,] 0.005229175 . .           . . . .           . . . . .           .
                                                                                 
    [1,] .          . .          . . . . .           .           .          . . .
    [2,] .          . .          . . . . .           0.005313191 .          . . .
    [3,] .          . .          . . . . .           .           .          . . .
    [4,] .          . 0.01890981 . . . . .           .           0.00577222 . . .
    [5,] 0.01999313 . .          . . . . .           .           .          . . .
    [6,] .          . .          . . . . 0.008423066 .           .          . . .
                                                                                   
    [1,] .          .         . . . . . . . . . .          . . . .          . . . .
    [2,] 0.03351598 .         . . . . . . . . . .          . . . .          . . . .
    [3,] .          0.0151425 . . . . . . . . . .          . . . 0.01131623 . . . .
    [4,] .          .         . . . . . . . . . .          . . . .          . . . .
    [5,] .          .         . . . . . . . . . .          . . . .          . . . .
    [6,] .          .         . . . . . . . . . 0.01086732 . . . .          . . . .
                             
    [1,] . . . . . .         
    [2,] . . . . . .         
    [3,] . . . . . 0.03159152
    [4,] . . . . . .         
    [5,] . . . . . .         
    [6,] . . . . . .         



```R
cell_clone_prob = cell_clone_prob_raw[deviance < 0.3,]
cell_clone_prob = cell_clone_prob/rowSums(cell_clone_prob)
```


```R
cell_clone_prob = mat_sparsify(mat = cell_clone_prob,row_mass = 0.9,col_mass = 0.9)
cell_clone_prob = cell_clone_prob/rowSums(cell_clone_prob)
cell_clone_prob = Matrix(cell_clone_prob, sparse = TRUE)
```


```R
colnames(cell_clone_prob) = names(table(cell_clone$clone))
```

We can use an example to show how the cell clone assignment is smoothed during this process:


```R
example_clone = "Lineage-1193-6"
```


```R
clone_example_umap_raw = dimplot(umap,annot = cell_meta %>% filter(clone == example_clone) %>% mutate(group = "1"),
                                 color_by = "group",label = FALSE,size = 2)+
theme_classic()+
NoLegend()+
theme(axis.text = element_blank(),    
      axis.ticks = element_blank())

clone_example_umap_raw
```


    
![png](demo_files/demo_35_0.png)
    



```R
cell_meta$prob = NA
cell_meta[rownames(cell_clone_prob),]$prob = cell_clone_prob[,example_clone]
```


```R
umap$prob = cell_meta$prob
```


```R
clone_example_umap_smooth = ggplot(mapping = aes(x= umap_1,y = umap_2))+
geom_point(data = umap,color = "lightgrey",size = 0.1)+
geom_point(data = umap %>% filter(prob > 0.025) %>% arrange(prob),
           aes(color = prob),size = 2,alpha = 0.5)+
scale_color_viridis_c(name = "prob",option = "turbo")+
theme_classic()+
theme(axis.text = element_blank(),   
      axis.ticks = element_blank())
clone_example_umap_smooth = ggrastr::rasterize(clone_example_umap_smooth,layers='Point', dpi=300)
clone_example_umap_smooth
```


    
![png](demo_files/demo_38_0.png)
    


### Clonewise transcriptomic distance

After we smooth the cell density for each expanded clone over the cell transcriptomic embedding, we can measure their cell distribution difference via optimal transportation distance. 

We provide two methods for this, the first one calculate the exact optimal transportation distance, and the second one is a faster approximation, but less accurate. You can switch the modes by the `exact` parameter in the `clone_disance` function.
### not run
plan(multisession, workers = 8)
### exact OT distance
clone_dis = clone_disance(as.matrix(pca),cell_clone_prob,outpath = "./",exact = TRUE, cores = 8)
### nearest neighbor distance approximation
clone_dis = clone_disance(as.matrix(pca),cell_clone_prob,outpath = "./",exact = FALSE)
This step is time and memory consuming, which can take hours even with the approximation. We strongly recommend to put this step in a server or cluster and run it with multiple cores. (As a reference, the approximation method takes 40 minutes to run on the HPC paralleled with 4 cores using 8 Gb RAM).

Here for time efficiency, we can use the pre-computed clone-wise distance for downstream applicaitons.


```R
file_path <- system.file("extdata", "clone_graph_dis.rds", package = "Clonotrace")
clone_dis = readRDS(file_path)
```

### clone embedding visualization

We first cluster clones using leiden algorithm and visualize the clone community in the UMAP. 


```R
set.seed(1230)
clone_cluster = leiden_dis(dismat = clone_dis,k = 20,resolution = 0.5,if_umap = TRUE)
clone_cluster_center = clone_cluster %>% group_by(cluster) %>% summarise_all(median)
```


```R
clone_cluster$mass = colSums(cell_clone_prob)
```


```R
clone_leiden_umap = ggplot(mapping = aes(x = umap1,y = umap2))+
geom_point(data = clone_cluster,aes(fill = cluster,size = log2(mass)),pch=21,color = "black",alpha = 0.75)+
ggrepel::geom_label_repel(data = clone_cluster_center,aes(label = cluster),size = 8)+
theme_bw()+
theme(axis.text = element_blank(),   
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
guides(fill = guide_legend(override.aes = list(size=5)))+
guides(fill = "none")+
theme(legend.position = "bottom")+
xlab("UMAP 1")+ylab("UMAP 2")
clone_leiden_umap
```


    
![png](demo_files/demo_50_0.png)
    


We can also visualize the clone embedding in other low dimension reduction like diffusion map:


```R
clone_leiden_dm <- diffusionMap::diffuse(clone_dis,maxdim = 30)
clone_leiden_dm_coord = as.data.frame(clone_leiden_dm$X)
colnames(clone_leiden_dm_coord) = paste("dm",1:ncol(clone_leiden_dm_coord),sep = "")
rownames(clone_leiden_dm_coord) = rownames(clone_cluster)
```

    Performing eigendecomposition
    Computing Diffusion Coordinates
    Used default value: 11 dimensions
    Elapsed time: 1.068 seconds



```R
dm_cluster_center = clone_leiden_dm_coord %>% group_by(clone_cluster$cluster) %>% summarise_all(median)
colnames(dm_cluster_center)[1] = "cluster"
```


```R
clone_leiden_diffmap = ggplot(mapping = aes(x = dm1,y = dm2))+
geom_point(data = clone_leiden_dm_coord,aes(fill = clone_cluster$cluster,size = log2(clone_cluster$mass)),
           pch=21,color = "black",alpha = 0.75)+
ggrepel::geom_label_repel(data = dm_cluster_center,aes(label = cluster),size = 8)+
theme_bw()+
theme(axis.text = element_blank(),    # Remove axis numbers (labels)
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
guides(fill = "none",size = guide_legend(title = "log2(mass)"))+
theme(legend.position = "bottom")+
xlab("dm 1")+ylab("dm 2")
clone_leiden_diffmap
```


    
![png](demo_files/demo_54_0.png)
    


Diffusion map is good to layout continous movement of clones over time.


```R
clone_cluster = clone_cluster %>% mutate(clone = rownames(clone_cluster)) %>% 
                                  mutate(lineage = substr(clone,start = 1,stop = nchar(clone)-2),
                                  time = substr(clone,start = nchar(clone),stop = nchar(clone))) %>% 
                                  mutate(time = paste("day",time,sep = ""))
```


```R
clone_link = long2wide(clone_cluster %>% dplyr::select(lineage,time,clone),
                       row_names_from = "lineage",
                       col_names_from = "time",
                       values_from = "clone",
                       symmetric = FALSE)
clone_link = na.omit(clone_link)
```


```R
colnames(clone_link) = c("start","end")
```


```R
clone_link_coord = cbind(clone_link,
                         clone_leiden_dm_coord[clone_link$start,c("dm1","dm2")],
                         clone_leiden_dm_coord[clone_link$end,c("dm1","dm2")],
                         clone_cluster[clone_link$start,"cluster"],
                         clone_cluster[clone_link$end,"cluster"])
colnames(clone_link_coord) = c("start","end","xdm1","xdm2","ydm1","ydm2","x_cluster","y_cluster")
```


```R
sub_annot = clone_cluster %>% filter(lineage %in% rownames(clone_link_coord))
sub_dm_coord = clone_leiden_dm_coord[rownames(sub_annot),]
```


```R
clone_leiden_diffmap_arrow = ggplot()+
geom_point(data = clone_leiden_dm_coord,aes(x = dm1,y = dm2),size = 0.2,color = "lightgrey")+
geom_point(data = sub_dm_coord,aes(x = dm1,y = dm2,fill = sub_annot$time,size = log2(sub_annot$mass)),shape = 21)+
geom_segment(data = clone_link_coord,aes(x = xdm1,y = xdm2,xend = ydm1,yend = ydm2),
             arrow = arrow(length = unit(0.03, "npc")),linewidth = 0.25,alpha = 0.75)+
scale_fill_viridis_d(name = "time")+
theme_bw()+
theme(axis.text = element_blank(),    # Remove axis numbers (labels)
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
guides(fill = guide_legend(override.aes = list(size=5)),
       size = guide_legend(title = "log2(mass)"))
clone_leiden_diffmap_arrow
```


    
![png](demo_files/demo_61_0.png)
    


### profile projection in cell embedding

After we cluster the clones into profiles, with the cell clone probabilty, we can also map cells to profiles:


```R
clone_profile = clone_cluster %>% dplyr::select(clone,cluster) %>% mutate(flag = 1)
clone_profile = long2sparse(long = clone_profile,row_names_from = "clone",col_names_from = "cluster",values_from = "flag")
```


```R
cell_profile_prob = cell_clone_prob %*% clone_profile[colnames(cell_clone_prob),]
rownames(cell_profile_prob) = rownames(cell_clone_prob)
```

For cells with profile probability higher than 0.5, we assign it with a profile label:


```R
cell_meta$profile = NA
cell_meta[rownames(cell_profile_prob),]$profile = apply(cell_profile_prob,1,function(x){
    if(max(x) > 0.5){
        return(which(x == max(x))[1])
    }
    else{
        return(NA)
    }
})
cell_meta = cell_meta %>% mutate(profile = as.factor(profile))
```


```R
cell_meta$profile_prob = NA
cell_meta[rownames(cell_profile_prob),]$profile_prob = apply(cell_profile_prob,1,max)
```

Then we can visualize how profile flow over the cell embedding:


```R
cell_profile_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "profile",
                            size = 0.1,alpha_by = "profile_prob",label = TRUE)+
theme(axis.text = element_blank(), 
      axis.ticks = element_blank())+
xlab("UMAP 1")+ylab("UMAP 2")+
guides(alpha = "none",color = guide_legend(override.aes = list(size=5)))
cell_profile_umap
```


    
![png](demo_files/demo_70_0.png)
    


Obvisouly, profile 5 and 6 monopotently differentiate to monocytes and profile 3 and 7 monopotently differentiate to netriphil, while profile 1 and 4 denote the bipotent differenitation.

### Clone level pseudotime

We first generate the clone level pseudotime by diffusion pseudotime in the clone embedding. Here we use the clone distance projection in the MDS space as the clone embedding.


```R
clone_mds = MASS::isoMDS(as.matrix(clone_dis), k=30)
clone_embedding = clone_mds$points
colnames(clone_embedding) = paste("mds",1:ncol(clone_embedding),sep = "_")
```

    initial  value 6.673968 
    iter   5 value 3.757695
    iter  10 value 2.922664
    iter  15 value 2.544110
    iter  20 value 2.367482
    iter  25 value 2.291132
    iter  30 value 2.251516
    iter  35 value 2.225022
    iter  40 value 2.207040
    final  value 2.194391 
    converged



```R
rownames(clone_embedding) = rownames(clone_cluster)
```

Here we choose the cell cluster 0 as the starting cluster for diffusion pseudotime. The `clone_dpt` function would automatically choose the clone which has the highest enrichment in this starting cluster to be the root of diffusion pseudotime.


```R
clone_t = clone_dpt(clone_embedding = clone_embedding,cell_meta = cell_meta,
                    clone_col = "clone",cluster_col = "cluster",start_cluster = "0")
```


```R
clone_cluster$dpt = clone_t
```

We can then visualize the clone level pseudotime in the diffusion map:


```R
clone_dpt_diffmap = ggplot(mapping = aes(x = dm1,y = dm2))+
geom_point(data = clone_leiden_dm_coord,aes(fill = clone_cluster$dpt,size = log2(clone_cluster$mass)),
           pch=21,color = "black",alpha = 0.75)+
theme_bw()+
theme(axis.text = element_blank(),    # Remove axis numbers (labels)
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
guides(fill = guide_legend(override.aes = list(size=5)))+
scale_fill_viridis_c(name = "clone t")+
guides(size = "none")+
theme(legend.position = "bottom")+
xlab("dm 1")+ylab("dm 2")
clone_dpt_diffmap
```


    
![png](demo_files/demo_80_0.png)
    


Then we can smooth the clone level pseudotime over the cell embedding by multiplying the clone pseudotime to cell clone probability. The pseudotime inference is necessary for the fate driving gene detection in the downstream analysis. 


```R
cell_meta$cell_t = NA
cell_meta[rownames(cell_clone_prob),]$cell_t = cell_clone_prob[,clone_cluster$clone] %*% clone_cluster$dpt
```


```R
cell_dpt_umap = dimplot(embedding = umap,annot = cell_meta,color_by = "cell_t",
                            size = 0.1,alpha_by = NULL,label = FALSE)+
theme(axis.text = element_blank(), 
      axis.ticks = element_blank())+
xlab("UMAP 1")+ylab("UMAP 2")+
scale_color_viridis_c(name = "cell t")
cell_dpt_umap
```


    
![png](demo_files/demo_83_0.png)
    


### clone weighted cell embedding (optional)

By integrating information from both the clone embedding and the cell embedding, we developed a new visualization method called clone-weighted cell embedding. In this representation, even if two cells are close in the transcriptomic (cell) embedding, they will be pulled apart in the clone-weighted embedding if their corresponding clones are distant in clone space. This embedding can better illustrate the clone level heterogeneity.


```R
cell_k = 30
clone_k = 15
cell_embedding = pca[rownames(cell_clone_prob),]
```


```R
plan(multisession, workers = 8)
coembed_dis = cell_clone_coembed(cell_embedding,clone_embedding)
```

This function `cell_clone_coembed` would return a sparse distance matrix, which denotes a weighted kNN network between cells. We can project this kNN network to UMAP for visulization. This step requires a python script.


```R
# set you python path
use_python("~/softwares/anaconda3/envs/py3.8/bin/python")
```


```R
script_path <- system.file("py", "embedding_from_kNN.py", package = "Clonotrace")
source_python("../py/embedding_from_kNN.py")
```


```R
cell_coembed_umap = umap_from_knn(coembed_dis,n_neighbors = 10,seed = 512)
```


```R
coembed_umap_cluster = dimplot(embedding = cell_coembed_umap,annot = cell_meta,alpha_by = NULL,label_size = 8,
                          color_by = "cluster",size = 0.1,alpha = 1)+
theme(axis.text = element_blank(),    # Remove axis numbers (labels)
      axis.ticks = element_blank())+
xlab("UMAP 1")+ylab("UMAP 2")+
NoLegend()
coembed_umap_cluster
```


    
![png](demo_files/demo_92_0.png)
    


Obviously, there appear two boundaries within the two branchs for monocytes (cluster 4, 3, 5) and neutriphils (cluster 1, 6, 7), which denote the clone level heterogenity. Once we colored cells by the profile identity, we can see the two boundaries seperate the bipoent and unipotent profiles.


```R
coembed_umap_profile = dimplot(embedding = cell_coembed_umap,annot = cell_meta,alpha_by = NULL,label_size = 8,
                          color_by = "profile",size = 0.1,alpha = 1)+
theme(axis.text = element_blank(),    # Remove axis numbers (labels)
      axis.ticks = element_blank())+
xlab("UMAP 1")+ylab("UMAP 2")+
NoLegend()
coembed_umap_profile
```


    
![png](demo_files/demo_94_0.png)
    


### profile enrichment in each cluster

In the dowsntream analysis we will focus on the early bipotent profiles (profile) and detect early fate-associated genes. To achieve this goal we first do an enrichment test between profiles and clusters to show in which cluster profile 1 has significant high density. And then within the target cluster, we compare the gene expression between profile 1 and other profiles to identify early fate-associated genes.

We first do the cell cluster profile enrichment using permutation test:


```R
enrich = cluster_profile_enrich(cell_profile_prob,cell_meta[rownames(cell_profile_prob),"cluster"],permute_n  = 1000)
```

We can then visualize the enrichment result using a heatmap:


```R
mass_long = wide2long(enrich[[1]])
mass_long = mass_long %>% mutate(i = rownames(enrich[[1]])[i],
                                 j = colnames(enrich[[1]])[j])
colnames(mass_long) = c("cluster","profile","mass")
mass_long = mass_long %>% group_by(cluster) %>% mutate(freq = mass/sum(mass))
mass_long$pvalue = c(t(enrich[[2]]))
```


```R
mass_long$signif = "ns"
mass_long$signif[mass_long$pvalue< 0.05] = "*"
mass_long$signif[mass_long$pvalue< 0.01] = "**"
mass_long$signif[mass_long$pvalue< 0.001] = "***"
```


```R
cell_type_profile_enrich_heatmap = ggplot(mass_long,aes(x = cluster,y = profile))+
geom_tile(aes(fill = freq),color = "black")+
scale_fill_viridis_c(option = "turbo")+
geom_text(data = mass_long %>% filter(signif != "ns"),aes(label = signif),size = 10)+
theme_classic()+
xlab("cell cluster")+
ylab("clone profile")
cell_type_profile_enrich_heatmap
```


    
![png](demo_files/demo_102_0.png)
    


From the heatmap we can see that our target profile 1 is enriched in 4 cell clusters, including 1, 3, 4 and 8, within which cluster 1 and 4 are at early differentiation stage. Here, as an example, we will identify genes differentially expressed in profile 1 in cluter 4.

### profile DEG within clusters


```R
exprs = seurat_object@assays$RNA$data
```


```R
start = Sys.time()
DEG_1_4 = profile_cluster_DEG(profile = "1",cluster = "4",exprs = exprs,cell_meta = cell_meta,cell_profile_prob = cell_profile_prob)
end = Sys.time()
print(end - start)
```

    Time difference of 7.764702 mins



```R
DEG_1_4_sig = DEG_1_4$stat %>% filter(padj < 0.05)
nrow(DEG_1_4_sig)
```


2165


In total we identified 2165 genes which are differentially expressed in profile 1 in cluster 4.

Here we can take a look at the top 10 genes with the largest effect size.


```R
DEG_1_4_sig %>% top_n(abs(cohen),n = 10)
```


<table class="dataframe">
<caption>A data.frame: 10 × 6</caption>
<thead>
	<tr><th></th><th scope=col>stat</th><th scope=col>pval</th><th scope=col>mean_diff</th><th scope=col>cohen</th><th scope=col>pcali</th><th scope=col>padj</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Dmkn</th><td>386.7920</td><td>9.632470e-152</td><td> 2.882819</td><td> 1.2346078</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Elane</th><td>326.3879</td><td>3.829860e-130</td><td> 3.278659</td><td> 1.1534887</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Fth1</th><td>209.3222</td><td> 2.394921e-86</td><td>-6.708008</td><td>-0.9494978</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>H2afy</th><td>292.5566</td><td>9.316957e-118</td><td>-6.653452</td><td>-1.1008000</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Lcn2</th><td>360.6077</td><td>1.851466e-142</td><td> 1.686889</td><td> 1.0452569</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Ly6a</th><td>268.3549</td><td>9.164804e-109</td><td>-2.698701</td><td>-1.0383824</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Olfm1</th><td>256.0305</td><td>3.846252e-104</td><td>-1.223922</td><td>-0.9836927</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Prtn3</th><td>261.5852</td><td>3.146175e-106</td><td> 2.319437</td><td> 1.0212434</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>S100a8</th><td>331.2332</td><td>6.707127e-132</td><td> 9.909431</td><td> 0.9617070</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Wfdc21</th><td>589.5540</td><td>7.449556e-220</td><td> 2.966061</td><td> 1.3873982</td><td>0</td><td>0</td></tr>
</tbody>
</table>



We can take `Dmkn` as an example to show how the linear regression model fits for the profile 1 and other profiles in cluster 4 along pseudotime:


```R
test_object = DEG_1_4
```


```R
null_fit = test_object$design_null %*% test_object$coef[3:4,]
full_fit = test_object$design_null %*% test_object$coef[1:2,]
rownames(null_fit) = rownames(full_fit) = test_object$cell
```


```R
example_gene = "Dmkn"
```


```R
cell_meta$exprs = exprs[example_gene,rownames(cell_meta)]
```


```R
cell_meta$null_exprs = NA
cell_meta[test_object$cell,]$null_exprs = null_fit[,example_gene]
cell_meta$full_exprs = NA
cell_meta[test_object$cell,]$full_exprs = full_fit[,example_gene]
```


```R
gene_fit_cluster = ggplot(cell_meta[test_object$cell,])+
        geom_point(aes(x = cell_t,y = exprs,fill = cell_profile_prob[test_object$cell,"1"]),size = 2,shape = 21,color = "black")+
        geom_line(aes(x = cell_t,y = null_exprs),color = "steelblue",linewidth = 1.5)+
        geom_line(aes(x = cell_t,y = full_exprs),color = "coral",linewidth = 1.5)+
        scale_fill_gradient2(name = "profile 1",low = "steelblue",mid = "whitesmoke",high = "coral",midpoint = 0.5,limit = c(0,1))+
        theme_classic()+
        xlab("pseudotime")+ylab(paste(example_gene,"exprs",sep = " "))
gene_fit_cluster
```


    
![png](demo_files/demo_117_0.png)
    


We can further visualize the expression for `Dmkn` over the cell UMAP:


```R
feature_scatter = dimplot(embedding = umap,annot = t(exprs),color_by = example_gene,label  = FALSE,size = 0.1)+
scale_color_viridis_c(name = "",option = "plasma")+
theme(axis.text = element_blank(),    # Remove axis numbers (labels)
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom")+
ggtitle(example_gene)+
xlab("UMAP 1")+ylab("UMAP 2")
feature_scatter
```


    
![png](demo_files/demo_119_0.png)
    



```R

```
