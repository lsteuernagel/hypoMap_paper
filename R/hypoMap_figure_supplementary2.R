##########
### Load & Prepare
##########

#set path
results_path_supplementary_figure2 = "figure_outputs/figure_supplementary_2/"
system(paste0("mkdir -p ",results_path_supplementary_figure2))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
source("R/utility_functions.R")
source("R/plot_functions.R")

# where to find large data objects (need to change to local dir !)
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

# load seurat objects via large_data_path
load_required_files(large_data_path = large_data_path)

# plotting
rasterize_point_size = 2.2
rasterize_pixels = 2048
text_size = 20
bg_col = "grey90"
cols_for_feature_plot = c(bg_col,"#0b3ebd") # "#0b3ebd"

##########
### Supplemental Figure 2: metric scatter for full map + Dataset UMAPS
##########

### load comparison data full
full_metrics = data.table::fread(paste0("data_inputs/hypothalamusMapFull_v4_comparison_8af8a1cd950067bb6859cbc1225c818d.txt"),data.table = F)

#text_size = 80
x_pixels = 3000

# add information
full_metrics$assay=full_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[a-zA-Z]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
full_metrics$ndim=full_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")%>% as.numeric()
full_metrics$features=full_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..*\\.features") %>% stringr::str_replace(pattern = "\\.\\..*\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.features",replacement = "")
full_metrics$method = stringr::str_extract(full_metrics$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
full_metrics$features_ngenes = stringr::str_extract(full_metrics$reduction,pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
# filter
full_metrics  =  full_metrics %>% dplyr::filter(ndim %in% c(30,50,60,80,90,NA))
data.table::fwrite(full_metrics,paste0(results_path_supplementary_figure2,"full_metrics_curated.csv"))

#plot
full_metrics_plot = ggplot2::ggplot(full_metrics,aes(x=mixing_score,y=purity_score,color=method))+geom_point(size=1.5)+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))+
  guides(color=guide_legend(ncol=1,override.aes = list(size=5)))+
  xlab("Mixing score")+ylab("Purity score")+
  ggtitle("Mean mixing vs mean purity")+xlim(0,100)+ylim(0,100)#+ theme_bw()
full_metrics_plot

#save
ggsave(filename = paste0(results_path_supplementary_figure2,"full_metrics_scatter.png"),
       plot = full_metrics_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path_supplementary_figure2,"full_metrics_scatter.pdf"),
       plot = full_metrics_plot, "pdf",dpi=600,width=350,height = 300,units="mm")


##########
### Supplemental Figure 2: Dataset UMAPS 
##########

# by Dataset
full_datasets_plot = DimPlot(full_map_seurat,group.by = "Dataset",reduction = paste0("umap_","scvi"),label = F,raster = F,pt.size = 0.2,shuffle = TRUE)+NoAxes()+
  theme(text = element_text(size=text_size))+ggtitle("Datasets on hypothalamus reference map")+guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
full_datasets_plot = rasterize_ggplot(full_datasets_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
full_datasets_plot

#save
ggsave(filename = paste0(results_path_supplementary_figure2,"full_map_datasets.png"),
       plot = full_datasets_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path_supplementary_figure2,"full_map_datasets.pdf"),
       plot = full_datasets_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

# by BatchID
full_BatchID_plot = DimPlot(full_map_seurat,group.by = "Batch_ID",reduction = paste0("umap_","scvi"),label = F,raster = F,pt.size = 0.2,shuffle = TRUE)+NoAxes()+
  theme(text = element_text(size=text_size))+ggtitle("Batches on hypothalamus reference map")+guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
full_BatchID_plot = rasterize_ggplot(full_BatchID_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
full_BatchID_plot

# # save
# ggsave(filename = paste0(results_path_supplementary_figure2,"full_map_BatchID.png"),
#        plot = full_BatchID_plot, "png",dpi=600,width=350,height = 300,units="mm")
# ggsave(filename = paste0(results_path_supplementary_figure2,"full_map_BatchID.pdf"),
#        plot = full_BatchID_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

# neuronmap:
# by Dataset
neurons_datasets_plot = DimPlot(neuron_map_seurat,group.by = "Dataset",reduction = paste0("umap_","scvi"),label = F,raster = F,pt.size = 0.2,shuffle = TRUE)+NoAxes()+
  theme(text = element_text(size=text_size))+ggtitle("Datasets on hypothalamus neuron reference map")+guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
neurons_datasets_plot = rasterize_ggplot(neurons_datasets_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
neurons_datasets_plot

# save
ggsave(filename = paste0(results_path_supplementary_figure2,"neuronmap_datasets.png"),
       plot = neurons_datasets_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path_supplementary_figure2,"neuronmap_datasets.pdf"),
       plot = neurons_datasets_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

# by BatchID
neurons_batch_plot = DimPlot(neuron_map_seurat,group.by = "Batch_ID",reduction = paste0("umap_","scvi"),label = F,raster = F,pt.size = 0.2,shuffle = TRUE)+NoAxes()+
  theme(text = element_text(size=text_size))+ggtitle("Batches on hypothalamus neuron reference map")+guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
neurons_batch_plot = rasterize_ggplot(neurons_batch_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
neurons_batch_plot

# # save
# ggsave(filename = paste0(results_path_supplementary_figure2,map_name,"_batchid.png"),
#        plot = neurons_batch_plot, "png",dpi=600,width=350,height = 300,units="mm")
# ggsave(filename = paste0(results_path_supplementary_figure2,map_name,"_batchid.pdf"),
#        plot = neurons_batch_plot, "pdf",dpi=600,width=350,height = 300,units="mm")
