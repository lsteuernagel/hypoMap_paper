
##########
### Load & Prepare
##########

results_path_figure1 = "figure_outputs/figure_1/"
system(paste0("mkdir -p ",results_path_figure1))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
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
### Figure 1a
##########

hypoMap_annotated_classes_plot = DimPlot(hypoMap_v2_seurat,group.by = "Author_Class_Curated",reduction = paste0("umap_","scvi"),cols=getOkabeItoPalette(14),
                                         label = TRUE,label.size = 7.5,raster = TRUE,pt.size = seurat_pt_size,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# hypoMap_annotated_classes_plot = DimPlot(hypoMap_v2_seurat,group.by = "Author_Class_Curated",reduction = paste0("umap_","scvi"),
#                                          label = TRUE,label.size = 6,raster = F,pt.size = 0.2,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
#   theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# neurons_annotated_plot = rasterize_ggplot(hypoMap_annotated_classes_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
hypoMap_annotated_classes_plot

ggsave(filename = paste0(results_path_figure1,"hypoMap_annotated_celltypes.png"),
       plot = hypoMap_annotated_classes_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_figure1,"hypoMap_annotated_celltypes.pdf"),
       plot = hypoMap_annotated_classes_plot, "pdf",dpi=600,width=300,height = 300,units="mm")

# source files
figure1_sourcedata = hypoMap_annotated_classes_plot$data
figure1_sourcedata$Cell_ID = rownames(figure1_sourcedata)

##########
### Figure 1b:neuron subset UMAP
##########

# factor:
neuron_clusters = scUtils::find_children("C2-1",hypoMap_v2_seurat@misc$clustering_edgelist)
hypoMap_v2_seurat@meta.data$C66_named_neurons = hypoMap_v2_seurat@meta.data$C66_named
hypoMap_v2_seurat@meta.data$C66_named_neurons[! hypoMap_v2_seurat@meta.data$C66 %in% neuron_clusters] = NA
hypoMap_v2_seurat@meta.data$C66_named_neurons = factor(hypoMap_v2_seurat@meta.data$C66_named_neurons,levels = unique(hypoMap_v2_seurat@meta.data$C66_named_neurons))
# non-neurons in grey (just for paper plot)
hypoMap_v2_seurat@meta.data$C66_named_neurons[hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,1] < -6.5 | 
                                                (hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,1] > 2) & (hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,2] > 5) | 
                                                hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings[,2]  > 10 ] =NA

neuron_umap_plot_c66 = DimPlot(hypoMap_v2_seurat,group.by = "C66_named_neurons",reduction = paste0("umap_scvi"),cols=getOkabeItoPalette(51),order = FALSE,
                               label = TRUE,label.size = 4.5,raster = TRUE,pt.size = seurat_pt_size,repel = TRUE,na.value = bg_col,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))+ggtitle("Neurons C66 on hypothalamus reference map")
neuron_umap_plot_c66


ggsave(filename = paste0(results_path_figure1,"neuron_umap_C66.png"),
       plot = neuron_umap_plot_c66, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_figure1,"neuron_umap_C66.pdf"),
       plot = neuron_umap_plot_c66, "pdf",dpi=600,width=300,height = 300,units="mm")

# source files
figure1_sourcedata =  dplyr::left_join(figure1_sourcedata,neuron_umap_plot_c66$data %>% dplyr::mutate(Cell_ID = rownames(neuron_umap_plot_c66$data)) %>% dplyr::select(Cell_ID, C66_named_neurons),by="Cell_ID")

##########
### Figure 1c
##########

## make feature plot with multiple genes:
# c("Slc32a1","Slc17a6","Slc18a2","Trh","Vip","Crh","Npy","Cck")
#"Dbh","Hdc","Tph2","Th","Chat")
p <- FeaturePlot(hypoMap_v2_seurat,features = c("Slc32a1","Slc17a6","Th","Hdc","Sim1","Nr5a1","Tbx3","Rgs16"), combine = FALSE,raster = TRUE,order=TRUE,
                 cols = cols_for_feature_plot,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
  #p[[i]] <- scUtils::rasterize_ggplot(p[[i]],pixel_raster = 2048,pointsize = 1.8)
}
combined_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 2)
combined_feature_plots

# save
ggsave(filename = paste0(results_path_figure1,"gene_example_umaps.png"),
       plot = combined_feature_plots, "png",dpi=400,width=200,height = 400,units="mm")
ggsave(filename = paste0(results_path_figure1,"gene_example_umaps.pdf"),
       plot = combined_feature_plots, "pdf",dpi=400,width=200,height =400,units="mm")

# source files
all_gene_expression = Seurat::FetchData(hypoMap_v2_seurat,vars = c("Slc32a1","Slc17a6","Th","Hdc","Sim1","Nr5a1","Tbx3","Rgs16"))
all_gene_expression$Cell_ID = rownames(all_gene_expression)
figure1_sourcedata =  dplyr::left_join(figure1_sourcedata,all_gene_expression,by="Cell_ID")


##########
### Save source data
##########

data.table::fwrite(figure1_sourcedata,paste0(results_path_figure1,"source_figure1_umap_data.txt"),sep="\t")
