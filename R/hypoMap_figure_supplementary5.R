##########
### Load & Prepare
##########

#set path
results_path_supplementary_figure5 = "figure_outputs/figure_supplementary_5/"
system(paste0("mkdir -p ",results_path_supplementary_figure5))

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

## plotting
rasterize_point_size = 1.5
rasterize_pixels = 1536
bg_col = "grey90"
cols_for_feature_plot = c(bg_col,"#0b3ebd") # "#0b3ebd"

##########
### Figure 7: marker gene plots sn and sc ? side-by-side ?
##########

key_bactrap_genes = c("Agrp","Pomc","Lepr","Cartpt","Pnoc")# "Glp1r",

# make plots from sc data
p_sc <- FeaturePlot(neuron_map_seurat,features = key_bactrap_genes,reduction =paste0("umap_scvi"),cols = cols_for_feature_plot, combine = FALSE,order = FALSE)
for(i in 1:length(p_sc)){
  p_sc[[i]] <- p_sc[[i]] + NoLegend() + NoAxes() #+ scale_color_gradientn(colours = colorvec) 
  p_sc[[i]] <- rasterize_ggplot( p_sc[[i]],pixel_raster = rasterize_pixels, pointsize =  rasterize_point_size)
} #+scale_color_gradient(low="lightgrey",high=max_color)} # remove axes and legends
length(p_sc)
cp_bacTRAPgenes = cowplot::plot_grid(plotlist = p_sc,ncol = 2)
#cp_bacTRAPgenes

# make plots from sn data
p_sn <- FeaturePlot(query_snseq_neurons,features = key_bactrap_genes,reduction =paste0("umap_scvi"),cols = cols_for_feature_plot, combine = FALSE,order = FALSE)
for(i in 1:length(p_sn)){
  p_sn[[i]] <- p_sn[[i]] + NoLegend() + NoAxes() #+ scale_color_gradientn(colours = colorvec) 
  p_sn[[i]] <- rasterize_ggplot( p_sn[[i]],pixel_raster = rasterize_pixels, pointsize =  rasterize_point_size)
}#+scale_color_gradient(low="lightgrey",high=max_color)} # remove axes and legends
length(p_sn)
cp_bacTRAPgenes_sn = cowplot::plot_grid(plotlist = p_sn,ncol = 2)

# combine side by side
new_list = list()
for(i in 1:length(p_sc)){
  new_list[[(i*2-1)]] = p_sc[[i]]
  new_list[[(i*2)]] = p_sn[[i]]
}
combined_feature_plots = cowplot::plot_grid(plotlist = new_list,ncol = 2)
combined_feature_plots

# save
ggsave(filename = paste0(results_path_supplementary_figure5,"genes_bacTRAP_Umaps.png"),
       plot = combined_feature_plots, "png",dpi=400,width=200,height = 300,units="mm")
ggsave(filename = paste0(results_path_supplementary_figure5,"genes_bacTRAP_Umaps.pdf"),
       plot = combined_feature_plots, "pdf",dpi=400,width=200,height =300,units="mm")


