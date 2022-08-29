##########
### Load & Prepare
##########

#set path
results_path_extended_figure1 = "figure_outputs/figure_extended_1/"
system(paste0("mkdir -p ",results_path_extended_figure1))

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
### Figure 1a
##########

# TODO: update and add workflow diagram

##########
### Figure 1b
##########

# load comparison data
neurons_metrics = data.table::fread("data_inputs/hypothalamusMapNeurons_v4_comparison_457fc60c3c4f1911bcbc6c5d46127037.txt",data.table = F)


# add information
neurons_metrics$assay=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[a-zA-Z]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
neurons_metrics$ndim=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")%>% as.numeric()
neurons_metrics$features=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..*\\.features") %>% stringr::str_replace(pattern = "\\.\\..*\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.features",replacement = "")
neurons_metrics$method = stringr::str_extract(neurons_metrics$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
neurons_metrics$features_ngenes = stringr::str_extract(neurons_metrics$reduction,pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
# filter
neurons_metrics  =  neurons_metrics %>% dplyr::filter(ndim %in% c(30,50,60,80,90,NA))
data.table::fwrite(neurons_metrics,paste0(results_path_extended_figure1,"neurons_metrics_curated.csv"))

#plot
metrics_plot = ggplot2::ggplot(neurons_metrics,aes(x=mixing_score,y=purity_score,color=method))+geom_point(size=1.5)+
  theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1))+
  scale_color_manual(values=getOkabeItoPalette(6))+#(getOkabeItoPalette(6))+
  guides(color=guide_legend(ncol=1,override.aes = list(size=5)))+
  ggtitle("Mean mixing vs mean purity")+xlim(0,100)+ylim(0,100)#+ theme_bw()
metrics_plot

#save
ggsave(filename = paste0(results_path_extended_figure1,"neurons_metrics_scatter.png"),
       plot = metrics_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure1,"neurons_metrics_scatter.pdf"),
       plot = metrics_plot, "pdf",dpi=600,width=350,height = 300,units="mm")

# source: 
source_ext_figure1_b =metrics_plot$data
data.table::fwrite(source_ext_figure1_b,paste0(results_path_extended_figure1,"source_ext_figure1_b_metrics.txt"),sep="\t")

##########
### Supplemental Figure 1c: UMAPs of best methods
##########

# load subset  version:
large_data_path_v1 = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

# load seurat objects via large_data_path
load_required_files(large_data_path = large_data_path_v1,filenames = c("neuron_map_seurat"="hypothalamus_neurons_map.rds"))

## make plots
# rasterize_point_size = 1.5
# rasterize_pixels = 1536
# 

# filter
top_from_each_method=neurons_metrics %>% dplyr::filter(ndim %in% c("30","50","60","80","90",NA)) %>% dplyr::mutate(rank_sum = purity_score+mixing_score,rank_diff = abs(purity_score-mixing_score)) %>%
  dplyr::mutate(rank_col =rank_sum-rank_diff) %>% dplyr::arrange(desc(rank_col)) %>% dplyr::group_by(method) %>% dplyr::filter(silhouette_score_norm>50 | silhouette_score_norm ==max(silhouette_score_norm)) %>% dplyr::top_n(n = 1)
names_to_check = top_from_each_method$reduction
names(names_to_check) = paste0("best_",top_from_each_method$method)

## add 2 more examples:
names_to_check = c(names_to_check,"scvi_cov"="scVI_1_300_0.025_3_256_gene_zinb_cov2..scVI..90..features_RNA.log.vst.split_Batch_ID.features.750_73e0c7a2fe2b9ab7fe2f1b1afa2980d6")
names_to_check = c(names_to_check,"harmony_lowdim"="harmony_8_Batch_ID_2_0.1_130_PCA..RNA..30..RNA.log.vst.split_Batch_ID.features.750_9625657cf22d373e46bdde89031e6f00")
names_to_check

## CAREFUL THIS RUNS SOME TIME TO CALC ALL UMAPS!

# add reductions:
integration_res_path = paste0(large_data_path_v1,"best_reductions_per_method/")
for( i in 1:length(names_to_check)){
  neuron_map_seurat = add_reduction_seurat(neuron_map_seurat,integration_name=names_to_check[i],new_name=names(names_to_check)[i],
                                           integration_res_path,max_dim=200,calc_umap=TRUE,k_param_umap=30,overwrite =F,overwrite2=F)
}

# umap by Dataset or Batch_ID
p <-list()
for(i in 1:length(names_to_check)) {
  p[[i]] <- DimPlot(neuron_map_seurat,group.by = "Dataset",reduction = paste0("umap_",names(names_to_check)[i]),label=FALSE,shuffle = TRUE,cols = getOkabeItoPalette(13),
                    raster = TRUE,pt.size = 0.75,repel = TRUE,raster.dpi = c(1536,1536))+
    NoAxes() + NoLegend()+ggtitle(names(names_to_check)[i])
  # p[[i]] <- rasterize_ggplot( p[[i]],pixel_raster = rasterize_pixels, pointsize =  rasterize_point_size)
}
compareMthods_cp_batch = cowplot::plot_grid(plotlist = p,ncol=4)
compareMthods_cp_batch

#save
ggsave(filename = paste0(results_path_extended_figure1,"umap_dataset_methods.png"),
       plot = compareMthods_cp_batch, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure1,"umap_dataset_methods.pdf"),
       plot = compareMthods_cp_batch, "pdf",dpi=400,width=300,height =200,units="mm")

# source: 
for(i in 1:length(p)){
  plot_data = p[[i]]$data
  name = names(names_to_check)[[i]]
  data.table::fwrite(plot_data,paste0(results_path_extended_figure1,"source_ext_figure1_c_",name,".txt"),sep="\t")
}
