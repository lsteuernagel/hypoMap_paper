
##########
### Load & Prepare
##########

results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_1/"
system(paste0("mkdir -p ",results_path))

# load everything required
source("scripts/paper_figures_new/load_data.R")

# source
# TODO: needed ?
source("utils.R")

### load comparison data
neurons_metrics = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/documentation/hypothalamusMapNeurons_v4_comparison_457fc60c3c4f1911bcbc6c5d46127037.txt",data.table = F)

## plotting
rasterize_point_size = 2.2
rasterize_pixels = 2048

##########
### Figure 1b
##########

#text_size = 80
text_size = 30

# add information
neurons_metrics$assay=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[a-zA-Z]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
neurons_metrics$ndim=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")%>% as.numeric()
neurons_metrics$features=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..*\\.features") %>% stringr::str_replace(pattern = "\\.\\..*\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.features",replacement = "")
neurons_metrics$method = stringr::str_extract(neurons_metrics$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
neurons_metrics$features_ngenes = stringr::str_extract(neurons_metrics$reduction,pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
# filter
neurons_metrics  =  neurons_metrics %>% dplyr::filter(ndim %in% c(30,50,60,80,90,NA))

#plot
metrics_plot = ggplot2::ggplot(neurons_metrics,aes(x=mixing_score,y=purity_score,color=method))+geom_point(size=1.5)+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))+
  guides(color=guide_legend(ncol=1,override.aes = list(size=5)))+
  ggtitle("Mean mixing vs mean purity")+xlim(0,100)+ylim(0,100)#+ theme_bw()
metrics_plot

#save
ggsave(filename = paste0(results_path,"neurons_metrics_scatter.png"),
       plot = metrics_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path,"neurons_metrics_scatter.pdf"),
       plot = metrics_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

##########
### Figure 1c
##########

text_size = 20

#DimPlot(full_map_seurat,group.by = "Curated_Class",reduction = paste0("umap_","scvi"),label = TRUE,label.size = 6,raster = F,pt.size = 0.2)+NoLegend()+NoAxes() # ,cols = "polychrome"
full_celltype_plot = DimPlot(full_map_seurat,group.by = "Curated_Class",reduction = paste0("umap_","scvi"),label = TRUE,label.size = 6,raster = F,pt.size = 0.2)+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))+ggtitle("Celltypes on hypothalamus reference map")
full_celltype_plot = rasterize_ggplot(full_celltype_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
full_celltype_plot

#save
ggsave(filename = paste0(results_path,"full_map_celltypes.png"),
       plot = full_celltype_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path,"full_map_celltypes.pdf"),
       plot = full_celltype_plot, "pdf",dpi=600,width=300,height = 300,units="mm")

##########
### Figure 1d
##########

neuron_map_seurat@meta.data$new_name_col = neuron_map_seurat@meta.data$K31_named
#DimPlot(neuron_map_seurat,group.by = "new_name_col",reduction = paste0("umap_","scvi"),label = TRUE,label.size = 5,repel = TRUE)+NoLegend() # ,cols = "polychrome"

text_size = 20
neurons_annotated_plot = DimPlot(neuron_map_seurat,group.by = "new_name_col",reduction = paste0("umap_","scvi"),label = TRUE,label.size = 6,raster = F,pt.size = 0.2)+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))+ggtitle("Clusters on hypothalamus neuron reference map")
neurons_annotated_plot = rasterize_ggplot(neurons_annotated_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
neurons_annotated_plot

ggsave(filename = paste0(results_path,map_name,"_annotated_clusters.png"),
       plot = neurons_annotated_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path,map_name,"_annotated_clusters.pdf"),
       plot = neurons_annotated_plot, "pdf",dpi=600,width=300,height = 300,units="mm")



