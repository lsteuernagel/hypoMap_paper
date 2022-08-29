##########
### Load & Prepare
##########

#set path
results_path_extended_figure2 = "figure_outputs/figure_extended_2/"
system(paste0("mkdir -p ",results_path_extended_figure2))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
source("R/utility_functions.R")
source("R/plot_functions.R")

# load seurat objects via large_data_path
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/"
load_required_files(large_data_path = large_data_path)

# load colors 
load_colors()
getLongPalette = colorRampPalette(long_palette_strong)
getOkabeItoPalette = colorRampPalette(short_palette)

## plotting
load_plot_params()

##########
### a: scvi metrics
##########

complete_evaluation_result = data.table::fread("data_inputs/complete_evaluation_results.txt",data.table = F)
# add information
complete_evaluation_result$method = stringr::str_extract(complete_evaluation_result$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
# pca:
complete_evaluation_result$ndim[complete_evaluation_result$method=="PCA"]=complete_evaluation_result$reduction[complete_evaluation_result$method=="PCA"] %>% stringr::str_extract(pattern="\\.[0-9]+\\.") %>% stringr::str_replace(pattern = "\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.",replacement = "")
complete_evaluation_result$features_ngenes[complete_evaluation_result$method=="PCA"]=complete_evaluation_result$reduction[complete_evaluation_result$method=="PCA"] %>% stringr::str_extract(pattern="[0-9]+\\.[0-9]+") %>% stringr::str_replace(pattern = "[0-9]+\\.",replacement = "") #%>% stringr::str_replace(pattern = "\\.",replacement = "")
# scvi
complete_evaluation_result$ndim[complete_evaluation_result$method=="scVI"]=complete_evaluation_result$reduction[complete_evaluation_result$method=="scVI"] %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
complete_evaluation_result$features_ngenes[complete_evaluation_result$method=="scVI"] = stringr::str_extract(complete_evaluation_result$reduction[complete_evaluation_result$method=="scVI"],pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
complete_evaluation_result$scvi_params[complete_evaluation_result$method=="scVI"] = complete_evaluation_result$reduction[complete_evaluation_result$method=="scVI"] %>%
  stringr::str_extract(pattern="scVI_[0-9]+_[0-9]+_0\\.[0-9]+_[0-9]+_[0-9]+") %>% stringr::str_replace(pattern = "scVI_[0-9]+_",replacement = "") #%>% as.numeric()
complete_evaluation_result$scvi_params[complete_evaluation_result$method=="PCA"] = "PCA"
complete_evaluation_result$nlayers = complete_evaluation_result$scvi_params %>% stringr::str_extract(pattern = "_[0-9]{1,}_")%>% stringr::str_replace(pattern = "_",replacement = "") %>% stringr::str_replace(pattern = "_",replacement = "")
complete_evaluation_result$nlayers = complete_evaluation_result$scvi_params %>% stringr::str_extract(pattern = "_[0-9]{1,}_")%>% stringr::str_replace(pattern = "_",replacement = "") %>% stringr::str_replace(pattern = "_",replacement = "")
complete_evaluation_result$max_epochs  = complete_evaluation_result$scvi_params %>% stringr::str_extract(pattern = "[0-9]+_0\\.")%>% stringr::str_replace(pattern = "_0\\.",replacement = "") #%>% stringr::str_replace(pattern = "_",replacement = "")
complete_evaluation_result$dropout_rate = complete_evaluation_result$scvi_params %>% stringr::str_extract(pattern = "_0\\.[0-9]+")%>% stringr::str_replace(pattern = "_",replacement = "") #%>% stringr::str_replace(pattern = "_",replacement = "")

data.table::fwrite(complete_evaluation_result,"data_inputs/complete_evaluation_results_updated.txt",sep="\t")

require(ggplot2)
hypoMap_v2_eval_plot= ggplot2::ggplot(complete_evaluation_result,aes(mixing_score,purity_score,color=method,size=asw_norm))+
  geom_point()+
  scale_color_manual(values=getOkabeItoPalette(2))+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))+
  guides(color=guide_legend(ncol=1,override.aes = list(size=5)))+
  ggtitle("Mean mixing vs mean purity")#+xlim(0,100)+ylim(0,100)#+ theme_bw()
hypoMap_v2_eval_plot

#save
ggsave(filename = paste0(results_path_extended_figure2,"v2_metrics_scatter.png"),
       plot = hypoMap_v2_eval_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure2,"v2_metrics_scatter.pdf"),
       plot = hypoMap_v2_eval_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

# source: 
source_ext_figure2_a =hypoMap_v2_eval_plot$data
data.table::fwrite(source_ext_figure2_a,paste0(results_path_extended_figure2,"source_ext_figure2_a_metrics.txt"),sep="\t")


##########
### b: scvi full vs neurons metrics
##########

full_vs_neurons_evaluation_result = data.table::fread("data_inputs/full_vs_neurons_evaluation_results.txt",data.table = F)

full_neuron_eval_plot= ggplot2::ggplot(full_vs_neurons_evaluation_result,aes(mixing_score,purity_score,color=source,size=asw_norm))+
  geom_point()+
  scale_color_manual(values=getOkabeItoPalette(5)[c(2,4)])+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))+
  guides(color=guide_legend(ncol=1,override.aes = list(size=5)))+
  ggtitle("Mean mixing vs mean purity")#+xlim(0,100)+ylim(0,100)#+ theme_bw()
full_neuron_eval_plot

#save
ggsave(filename = paste0(results_path_extended_figure2,"full_neuron_metrics_scatter.png"),
       plot = full_neuron_eval_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure2,"full_neuron_metrics_scatter.pdf"),
       plot = full_neuron_eval_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

# source: 
source_ext_figure2_b =full_neuron_eval_plot$data
data.table::fwrite(source_ext_figure2_b,paste0(results_path_extended_figure2,"source_ext_figure2_b_metrics_neuron.txt"),sep="\t")


##########
### c: boxplots scvi metrics
##########

## PROBABLY don't include this !!!!

plot_df = complete_evaluation_result[complete_evaluation_result$method=="scVI" & complete_evaluation_result$ndim %in% c("110","65","85"),]
plot_df$ndim= factor(plot_df$ndim,levels = c("65","85","110"))
plot_df$max_epochs= factor(plot_df$max_epochs,levels = c("50","100","150","200","300","400"))
plot_df$dropout_rate= factor(plot_df$dropout_rate,levels = c("0.01","0.05","0.1"))
plot_df$features_ngenes= factor(plot_df$features_ngenes,levels = as.character(sort(as.numeric(unique(plot_df$features_ngenes)))))

## plot epochs
p_asw_epochs = ggplot(plot_df, aes(x=max_epochs, y=asw_norm)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers)) +
  scale_fill_manual(values=getOkabeItoPalette(4))+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))
p_asw_epochs

p_mixing_epochs = ggplot(plot_df, aes(x=max_epochs, y=mixing_score)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(4))+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))
p_mixing_epochs

p_purity_epochs = ggplot(plot_df, aes(x=max_epochs, y=purity_score)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(4))+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))
p_purity_epochs

# TODO save if needed

#save
ggsave(filename = paste0(results_path_extended_figure2,"boxplot_asw_epochs.png"),
       plot = p_asw_epochs, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure2,"boxplot_asw_epochs.pdf"),
       plot = p_asw_epochs, "pdf",dpi=400,width=300,height = 200,units="mm")

ggsave(filename = paste0(results_path_extended_figure2,"boxplot_mixing_epochs.png"),
       plot = p_mixing_epochs, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure2,"boxplot_mixing_epochs.pdf"),
       plot = p_mixing_epochs, "pdf",dpi=400,width=300,height = 200,units="mm")

ggsave(filename = paste0(results_path_extended_figure2,"boxplot_purity_epochs.png"),
       plot = p_purity_epochs, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure2,"boxplot_purity_epochs.pdf"),
       plot = p_purity_epochs, "pdf",dpi=400,width=300,height = 200,units="mm")

# source: 
source_ext_figure2_c =plot_df
data.table::fwrite(source_ext_figure2_c,paste0(results_path_extended_figure2,"source_ext_figure2_c_scvi_tuning.txt"),sep="\t")

##########
### d: datasets on UMAP
##########

hypoMap_dataset_plot = DimPlot(hypoMap_v2_seurat,group.by = "Dataset",reduction = paste0("umap_","scvi"),cols=getOkabeItoPalette(18),shuffle = TRUE,
                                         label = F,label.size = 7,raster = TRUE,pt.size = seurat_pt_size,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+#NoLegend()+
  theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# hypoMap_annotated_classes_plot = DimPlot(hypoMap_v2_seurat,group.by = "Author_Class_Curated",reduction = paste0("umap_","scvi"),
#                                          label = TRUE,label.size = 6,raster = F,pt.size = 0.2,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
#   theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# neurons_annotated_plot = rasterize_ggplot(hypoMap_annotated_classes_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
hypoMap_dataset_plot

ggsave(filename = paste0(results_path_extended_figure2,"hypoMap_datasets.png"),
       plot = hypoMap_dataset_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure2,"hypoMap_datasets.pdf"),
       plot = hypoMap_dataset_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

##########
### e: purity cell types on UMAP
##########

# load
all_detected_celltypes =jsonlite::read_json("data_inputs/detected_celltypes.json")
all_detected_celltypes=lapply(all_detected_celltypes,unlist)

# add to seurat
library(purrr)
celltype_df <- purrr::map_df(all_detected_celltypes, ~as.data.frame(.x), .id="celltype")
colnames(celltype_df) = c("detected_celltype","Cell_ID")

temp_meta = dplyr::left_join(hypoMap_v2_seurat@meta.data,celltype_df,by=c("Cell_ID")) %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID
hypoMap_v2_seurat@meta.data=temp_meta

#plot
hypoMap_celltype_plot = DimPlot(hypoMap_v2_seurat,group.by = "detected_celltype",reduction = paste0("umap_","scvi"),cols=getOkabeItoPalette(34),order = TRUE,na.value = bg_col,
                               label = TRUE,label.size = 6,raster = TRUE,pt.size = seurat_pt_size,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
  theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# hypoMap_annotated_classes_plot = DimPlot(hypoMap_v2_seurat,group.by = "Author_Class_Curated",reduction = paste0("umap_","scvi"),
#                                          label = TRUE,label.size = 6,raster = F,pt.size = 0.2,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
#   theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# neurons_annotated_plot = rasterize_ggplot(hypoMap_annotated_classes_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
hypoMap_celltype_plot

ggsave(filename = paste0(results_path_extended_figure2,"hypoMap_celltypes_purity.png"),
       plot = hypoMap_celltype_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure2,"hypoMap_celltypes_purity.pdf"),
       plot = hypoMap_celltype_plot, "pdf",dpi=600,width=300,height = 300,units="mm")

# source: 
source_ext_figure2_d =hypoMap_dataset_plot$data
source_ext_figure2_e =hypoMap_celltype_plot$data
source_ext_figure2_e = bind_cols(source_ext_figure2_e,Dataset = source_ext_figure2_d$Dataset)
data.table::fwrite(source_ext_figure2_e,paste0(results_path_extended_figure2,"source_ext_figure2_de_umap_data.txt"),sep="\t")




