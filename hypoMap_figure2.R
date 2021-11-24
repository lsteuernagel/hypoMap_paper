
##########
### Load & Prepare
##########

require(tidyverse)
require(data.table)
require(Seurat)
require(ggplot2)
require(mapscvi)

### Set parameters
global_seed = 123467# seed
map_name = "hypothalamus_neurons_reference"
map_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_objects/"

results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_2/"
system(paste0("mkdir -p ",results_path))

### Load map
map_seurat_path = paste0(map_path,map_name,".rds")
neuron_map_seurat = readRDS(map_seurat_path)
#neuron_map_seurat@meta.data = cbind(neuron_map_seurat@meta.data ,neuron_map_seurat@misc$other_metdata)

## plotting
rasterize_point_size = 2.2
rasterize_pixels = 2048

##########
### tree new
##########


# http://yulab-smu.top/treedata-book/chapter2.html

# using the function in fisualization_functions.R

source("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scHarmonize/harmonization/visualization_functions.R")

# make data for first heatmap with percentages per dataset
heatmap_data = neuron_map_seurat@meta.data %>% dplyr::select(Cell_ID,Dataset,K169_pruned) %>% dplyr::group_by(K169_pruned,Dataset) %>% #dplyr::filter(predicted_Campbell!="NA") 
  dplyr::count(name = "presence")  %>% dplyr::group_by(K169_pruned) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label")) 
  tidyr::spread(key = Dataset,value=presence) 
heatmap_matrix = as.matrix(heatmap_data[,2:ncol(heatmap_data)])
rownames(heatmap_matrix) = heatmap_data$K169_pruned
heatmap_matrix[is.na(heatmap_matrix)] = 0

# make data for second heatmap with regions
heatmap_data2 = neuron_map_seurat@meta.data %>% dplyr::select(Cell_ID,suggested_region_curated,K169_pruned) %>% 
  dplyr::distinct(K169_pruned,suggested_region_curated,.keep_all=TRUE)
heatmap_data2$suggested_region_curated[heatmap_data2$suggested_region_curated=="NA"]=NA
heatmap_data2$suggested_region_curated[heatmap_data2$suggested_region_curated=="Paraventricular hypothalamic nucleus descending division"]="Paraventricular hypothalamic nucleus"
heatmap_matrix2 = as.matrix(heatmap_data2[,"suggested_region_curated"])
rownames(heatmap_matrix2) = heatmap_data2$K169_pruned
colnames(heatmap_matrix2) = "suggested_region_curated"

require(RColorBrewer)
## expand palette size
colourCount <- length(unique(heatmap_data2$suggested_region_curated)) # number of levels
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
# plot tree with heatmaps
circular_tree_heat = plot_cluster_tree(edgelist = neuron_map_seurat@misc$mrtree_edgelist,
                  heatmap_matrix=heatmap_matrix,
                  heatmap_matrix2 = heatmap_matrix2,
                  leaf_level=6,metadata=neuron_map_seurat@meta.data,
                  label_size = 2, show_genes = TRUE, legend_title_1 = "Pct", legend_title_2 = "Region",
                  matrix_offset = 0.1, matrix_width =0.4,matrix_width_2 = 0.1,heatmap_colnames = TRUE,
                  manual_off_second = 2,legend_text_size = 8,heatmap_text_size = 2,
                  heatmap_colors =c("grey90","darkred")) +
  scale_fill_brewer(palette = "Paired",na.value = "grey90")
  # scale_fill_gradient2(low="#eda09a",mid = "white",high = "#9aa9ed")
  #scale_fill_gradient2(low="#cc2118",mid = "white",high = "#302ac9",na.value = "#878787")
circular_tree_heat

#store:
ggsave(filename = paste0(results_path,"circular_tree_heat_label2.png"),
       plot = circular_tree_heat, "png",dpi=600,width=400,height = 400,units="mm")

ggsave(filename = paste0(results_path,"circular_tree_heat_label2.pdf"),
       plot = circular_tree_heat, "pdf",dpi=600,width=400,height = 400,units="mm")

##########
### campbell anno
##########
neuron_map_seurat@meta.data$new_name_col = neuron_map_seurat@meta.data$K31_named
neuron_map_seurat@meta.data$new_name_col = gsub("Slc32a1\\.|Slc17a6\\.","",neuron_map_seurat@meta.data$new_name_col)

# prepare column
campbell_names= names(table(neuron_map_seurat@meta.data$Author_CellType[neuron_map_seurat@meta.data$Dataset=="Campbell"]))[table(neuron_map_seurat@meta.data$Author_CellType[neuron_map_seurat@meta.data$Dataset=="Campbell"]) > 5]
campbell_names = campbell_names[!grepl("NA_",campbell_names)]

neuron_map_seurat@meta.data$campbell_anno_col=NA
neuron_map_seurat@meta.data$campbell_anno_col[neuron_map_seurat@meta.data$Dataset=="Campbell"] = neuron_map_seurat@meta.data$Author_CellType[neuron_map_seurat@meta.data$Dataset=="Campbell"]
neuron_map_seurat@meta.data$campbell_anno_col[! neuron_map_seurat@meta.data$campbell_anno_col %in% campbell_names] =NA
neuron_map_seurat@meta.data$campbell_anno_col = gsub("_Neurons[0-9]","",neuron_map_seurat@meta.data$campbell_anno_col)

# plot
campbell_anno_plot=DimPlot(neuron_map_seurat,group.by = "campbell_anno_col",reduction = paste0("umap_","scvi"),label = TRUE,label.size = 5,repel = TRUE,order = TRUE,na.value = "lightgrey")+
  NoLegend()+NoAxes()+scale_color_discrete(na.value="lightgrey")+ggtitle("Campbell ARH celltypes")

# change order
plot_data = campbell_anno_plot$data
campbell_anno_plot$data = plot_data[order(plot_data$campbell_anno_col,na.last = FALSE),]
#rasterize
campbell_anno_plot = rasterize_ggplot(campbell_anno_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
campbell_anno_plot

ggsave(filename = paste0(results_path,"campbell_annotations.png"),
       plot = campbell_anno_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path,"campbell_annotations.pdf"),
       plot = campbell_anno_plot, "pdf",dpi=600,width=300,height = 300,units="mm")

##########
### VIP inlay plot
##########

neuron_map_seurat@meta.data$new_name_col_169 = neuron_map_seurat@meta.data$K169_named
#neuron_map_seurat@meta.data$new_name_col_169 = gsub("Slc32a1\\.|Slc17a6\\.","",neuron_map_seurat@meta.data$new_name_col_169)

# plot
Idents(neuron_map_seurat) <- "K14_pruned"
neuron_map_seurat_vip_subset = subset(neuron_map_seurat,subset = K14_pruned == "K14-4")
neuron_map_seurat_vip_subset = subset(neuron_map_seurat_vip_subset,subset = umapscvi_1 > 1 & umapscvi_2 < -0.5)
vip_small_plot=DimPlot(neuron_map_seurat_vip_subset,group.by = "new_name_col_169",reduction = paste0("umap_","scvi"),label = TRUE,label.size = 6,repel = TRUE,order = TRUE,na.value = "lightgrey")+
  NoLegend()+NoAxes()+scale_color_discrete(na.value="lightgrey")+ggtitle("Vip neurons reference map")
vip_small_plot = rasterize_ggplot(vip_small_plot,pixel_raster = 1536,pointsize = 1.8)
vip_small_plot

ggsave(filename = paste0(results_path,"vip_small_plot.png"),
       plot = vip_small_plot, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"vip_small_plot.pdf"),
       plot = vip_small_plot, "pdf",dpi=450,width=200,height = 200,units="mm")

##########
### Campbell VIP sankey plot -- probably not !
##########

# library(mapscvi)
# 
# original_clustering = compare_clustering(neuron_map_seurat,"K169_named","Author_CellType",min_cells = 5,min_pct = 0.025,return_data=TRUE)
# campbell_clusters = unique(neuron_map_seurat@meta.data$Author_CellType[neuron_map_seurat@meta.data$Dataset=="Campbell"])
# original_clustering =  original_clustering[original_clustering$clustering_2 %in% campbell_clusters,]
# 
# clustering_2_filter = c("Rgs16/Vip_Neurons2")
# plot_sankey_comparison(original_clustering,clustering_1_filter = NULL,clustering_2_filter = clustering_2_filter,text_size=15)


##########
### Romanov 
##########

## see the romanov_scvi RMD in /scHarmonize/scripts/mapscvi or the mapscvi vignette
query_romanov_neurons = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypothalamus/romanov/mapped_data/query_romanov_neurons.rds")

# need mapscvi function!
# make plot
romanov_mapped_plot = mapscvi::plot_query_labels(query_seura_object=query_romanov_neurons,reference_seurat=neuron_map_seurat,label_col="K31_named",
                                                 label_col_query = "predicted_K31_named",overlay = TRUE,bg_col = "lightgrey",
                                                 query_pt_size = 0.4,labelonplot = TRUE,label.size=5)+
  ggtitle("Romanov et al. mapped on HypoMap")

romanov_mapped_plot = rasterize_ggplot(romanov_mapped_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
romanov_mapped_plot

ggsave(filename = paste0(results_path,"romanov_neurons_mapped.png"),
       plot = romanov_mapped_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path,"romanov_neurons_mapped.pdf"),
       plot = romanov_mapped_plot, "pdf",dpi=600,width=300,height = 300,units="mm")


