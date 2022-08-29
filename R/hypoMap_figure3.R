
##########
### Load & Prepare
##########

results_path_figure3 = "figure_outputs/figure_3/"
system(paste0("mkdir -p ",results_path_figure3))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
library(scales)
source("R/utility_functions.R")
source("R/plot_functions.R")

# where to find large data objects (need to change to local dir !)
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/"

# load seurat objects via large_data_path
load_required_files(large_data_path = large_data_path)
hypoMap_v2_seurat@meta.data$Author_Class_Curated[hypoMap_v2_seurat@meta.data$Author_Class_Curated=="Differentiating"] = "Dividing"

# load colors 
load_colors()
getLongPalette = colorRampPalette(long_palette_strong)
getOkabeItoPalette = colorRampPalette(short_palette)


## plotting
load_plot_params()



##########
### Figure 3 sankey plot kim
##########


# check Kim et all anoatations vs my clustering:
subset_seurat_kim10x = subset(hypoMap_v2_seurat,subset = Dataset == "Kim10x")
# DimPlot(subset_seurat_kim10x,group.by = "C25",label=TRUE,label.size = 3)+NoLegend()

# create cluster overview using mapscvi:
overview_clustering_kim = mapscvi::compare_clustering(query_seura_object = subset_seurat_kim10x,
                                                      clustering_1 = "Author_CellType",
                                                      clustering_2 = "C286_named",
                                                      min_cells = 5,
                                                      min_pct = 0.05,
                                                      return_data=TRUE)
# sankey
# i focus on the vmh part of the clusters
clustering_2_filter = unique(subset_seurat_kim10x@meta.data$C286_named[subset_seurat_kim10x@meta.data$C25 == "C25-3"]) # neurons only
#clustering_2_filter = unique(subset_seurat@meta.data$leiden_clusters_27[subset_seurat@meta.data$leiden_clusters_0.01 == "1"]) # neurons only
clustering_1_filter = NULL
sankey_clusters_kim = mapscvi::plot_sankey_comparison(overview_clustering_kim,
                                                      clustering_1_filter = clustering_1_filter,
                                                      clustering_2_filter = clustering_2_filter,
                                                      text_size=20,
                                                      use_and = FALSE,
                                                      light_factor = 0.45)
sankey_clusters_kim

# # save sankeys
library(webshot)
# https://stackoverflow.com/questions/65158327/how-can-i-save-a-networkd3sankeynetwork-into-a-static-image-automatically-vi
networkD3::saveNetwork(sankey_clusters_kim, paste0(results_path_figure3,"sankey_clusters_kim.html"))
# convert it
# need: webshot::install_phantomjs()
webshot::webshot(paste0(results_path_figure3,"sankey_clusters_kim.html"),file=paste0(results_path_figure3,"sankey_clusters_kim.png"), vwidth = 1000, vheight = 900)
webshot::webshot(paste0(results_path_figure3,"sankey_clusters_kim.html"),file=paste0(results_path_figure3,"sankey_clusters_kim.pdf"), vwidth = 1000, vheight = 900)

# source data
overview_clustering_kim_source = overview_clustering_kim[,1:7] %>% dplyr::filter(clustering_1 %in% clustering_1_filter | clustering_2 %in% clustering_2_filter)
data.table::fwrite(overview_clustering_kim_source,paste0(results_path_figure3,"source_figure_3_e_clustering_kim.txt"),sep="\t")

##########
### Figure 3 sankey plot chen
##########


# check Kim et all anoatations vs my clustering:
subset_seurat_chen = subset(hypoMap_v2_seurat,subset = Dataset == "ChenDropseq")
# DimPlot(subset_seurat_chen,group.by = "C25",label=TRUE,label.size = 3)+NoLegend()

# create cluster overview using mapscvi:
overview_clustering_chen = mapscvi::compare_clustering(query_seura_object = subset_seurat_chen,
                                                       clustering_1 = "Author_CellType",
                                                       clustering_2 = "C286_named",
                                                       min_cells = 5,
                                                       min_pct = 0.05,
                                                       return_data=TRUE)
# sankey
# i focus on the vmh part of the clusters
clustering_2_filter = unique(subset_seurat_chen@meta.data$C286_named[subset_seurat_chen@meta.data$C25 == "C25-3"]) # neurons only
#clustering_2_filter = unique(subset_seurat@meta.data$leiden_clusters_27[subset_seurat@meta.data$leiden_clusters_0.01 == "1"]) # neurons only
clustering_1_filter = NULL
sankey_clusters_chen = mapscvi::plot_sankey_comparison(overview_clustering_chen,
                                                       clustering_1_filter = clustering_1_filter,
                                                       clustering_2_filter = clustering_2_filter,
                                                       text_size=20,
                                                       use_and = FALSE,
                                                       light_factor = 0.45)
sankey_clusters_chen

# # save sankeys
library(webshot)
# https://stackoverflow.com/questions/65158327/how-can-i-save-a-networkd3sankeynetwork-into-a-static-image-automatically-vi
networkD3::saveNetwork(sankey_clusters_chen, paste0(results_path_figure3,"sankey_clusters_chen.html"))
# convert it
# need: webshot::install_phantomjs()
webshot::webshot(paste0(results_path_figure3,"sankey_clusters_chen.html"),file=paste0(results_path_figure3,"sankey_clusters_chen.png"), vwidth = 1000, vheight = 900)
webshot::webshot(paste0(results_path_figure3,"sankey_clusters_chen.html"),file=paste0(results_path_figure3,"sankey_clusters_chen.pdf"), vwidth = 1000, vheight = 900)

# source data
overview_clustering_chen_source = overview_clustering_chen[,1:7] %>% dplyr::filter(clustering_1 %in% clustering_1_filter | clustering_2 %in% clustering_2_filter)
data.table::fwrite(overview_clustering_chen_source,paste0(results_path_figure3,"source_figure_3_d_clustering_chen.txt"),sep="\t")

##########
### Figure 3 orientation umap:
##########

# # orientation umap:
vmh_subset = subset(hypoMap_v2_seurat,subset = C25 %in% "C25-3" & umapscvi_1 > -5 & umapscvi_1 < 1.5 & umapscvi_2 > -1.5 & umapscvi_2 < 7 )
#chen_vmh_subset = subset(hypoMap_v2_seurat,subset = C25 %in% "C25-2" & Dataset %in% c("ChenDropseq") & umapscvi_1 < -5 & umapscvi_2 > -3 & umapscvi_2 < 4.5 )
vmh_subset@meta.data$chen_celltypes = NA
vmh_subset@meta.data$chen_celltypes[vmh_subset@meta.data$Dataset %in% c("ChenDropseq")] = vmh_subset@meta.data$Author_CellType[vmh_subset@meta.data$Dataset %in% c("ChenDropseq")]
away= names(table(vmh_subset@meta.data$chen_celltypes)[table(vmh_subset@meta.data$chen_celltypes) < 5])
vmh_subset@meta.data$chen_celltypes[vmh_subset@meta.data$chen_celltypes %in% away] =NA
chen_vmh_subset_plot = DimPlot(vmh_subset,group.by = "chen_celltypes",reduction = paste0("umap_","scvi"),na.value = bg_col,cols=getOkabeItoPalette(7),order=TRUE,
                               label = TRUE,label.size = 5,raster = TRUE,pt.size = 3,raster.dpi = c(1536,1536))+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))
chen_vmh_subset_plot

ggsave(filename = paste0(results_path_figure3,"chen_C25_vmh_umap.png"),
       plot = chen_vmh_subset_plot, "png",dpi=400,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure3,"chen_C25_vmh_umap.pdf"),
       plot = chen_vmh_subset_plot, "pdf",dpi=400,width=200,height = 200,units="mm")

# source data UMAP
source_data_umap_c= chen_vmh_subset_plot$data %>% dplyr::mutate(Cell_ID = rownames( chen_vmh_subset_plot$data)) 
data.table::fwrite(source_data_umap_c,paste0(results_path_figure3,"source_figure_3_b_umap_chen.txt"),sep="\t")

# Kim10x
#kim_vmh_subset = subset(hypoMap_v2_seurat,subset = C25 %in% "C25-2" & Dataset %in% c("Kim10x") & umapscvi_1 < -5 & umapscvi_2 > -3 & umapscvi_2 < 4.5 )
vmh_subset@meta.data$kim_celltypes = NA
vmh_subset@meta.data$kim_celltypes[vmh_subset@meta.data$Dataset %in% c("Kim10x")] = vmh_subset@meta.data$Author_CellType[vmh_subset@meta.data$Dataset %in% c("Kim10x")]
away= names(table(vmh_subset@meta.data$kim_celltypes)[table(vmh_subset@meta.data$kim_celltypes) < 5])
vmh_subset@meta.data$kim_celltypes[vmh_subset@meta.data$kim_celltypes %in% away] =NA
kim_vmh_subset_plot = DimPlot(vmh_subset,group.by = "kim_celltypes",reduction = paste0("umap_","scvi"),na.value = bg_col,cols=getOkabeItoPalette(37),order=TRUE,
                              label = TRUE,label.size = 4.5,repel = TRUE,raster = TRUE,pt.size = 3,raster.dpi = c(1536,1536))+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))
kim_vmh_subset_plot

ggsave(filename = paste0(results_path_figure3,"kim10x_C25_vmh_umap.png"),
       plot = kim_vmh_subset_plot, "png",dpi=400,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure3,"kim10x_C25_vmh_umap.pdf"),
       plot = kim_vmh_subset_plot, "pdf",dpi=400,width=200,height = 200,units="mm")

# source data UMAP
source_data_umap_d= kim_vmh_subset_plot$data %>% dplyr::mutate(Cell_ID = rownames( kim_vmh_subset_plot$data)) 
data.table::fwrite(source_data_umap_d,paste0(results_path_figure3,"source_figure_3_c_umap_kim.txt"),sep="\t")

##########
### Figure 3 overview plots
##########

cellsh = hypoMap_v2_seurat@meta.data$Cell_ID[hypoMap_v2_seurat@meta.data$C25=="C25-3" & 
                                               hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,1] > -5 & 
                                               hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,1] < 1.5 & 
                                               hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,2] > -1.5 & 
                                               hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,2] < 7]
hypoMap_annotated_C25_plot = DimPlot(hypoMap_v2_seurat,group.by = "C25",reduction = paste0("umap_","scvi"),
                                     cells.highlight = cellsh,sizes.highlight = 0.1,cols.highlight = "#D55E00",na.value = bg_col,
                                     label = F,label.size = 4,raster = TRUE,pt.size = seurat_pt_size,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
hypoMap_annotated_C25_plot

ggsave(filename = paste0(results_path_figure3,"highlight_C25_vmh_umap.png"),
       plot = hypoMap_annotated_C25_plot, "png",dpi=400,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure3,"highlight_C25_vmh_umap.pdf"),
       plot = hypoMap_annotated_C25_plot, "pdf",dpi=400,width=200,height = 200,units="mm")


# source data UMAPs
source_data_umap_b = hypoMap_annotated_C25_plot$data %>% dplyr::mutate(Cell_ID = rownames( hypoMap_annotated_C25_plot$data)) 
data.table::fwrite(source_data_umap_b,paste0(results_path_figure3,"source_figure_3_a_umap_vmh.txt"),sep="\t")




