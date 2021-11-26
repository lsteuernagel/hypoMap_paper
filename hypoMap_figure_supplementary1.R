##########
### Load & Prepare
##########

#set path
results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_supplementary_1/"
system(paste0("mkdir -p ",results_path))

# load everything required
source("load_data.R")

# path with output files
data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_input/"

## plotting
rasterize_point_size = 1.1
rasterize_pixels = 1024

##########
### Supplemental Figure 1:  aucell mapped celltype cells mini umaps -- OPTIONAL
##########

# get mapped celltypes
mapped_celltypes = readRDS(paste0(data_path,"mapped_celltypes_neuronMap.rds"))

# plot
p=list()
for(i in 1:length(mapped_celltypes)) {
  p[[i]] <- DimPlot(neuron_map_seurat,reduction = "umap_scvi",cells.highlight = mapped_celltypes[[i]], sizes.highlight = 0.05,pt.size = 0.05,raster = F) + NoLegend() + NoAxes() + ggtitle(names(mapped_celltypes)[i])
  p[[i]] = rasterize_ggplot(p[[i]],pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
}
length(p)
cp_allcelltypes = cowplot::plot_grid(plotlist = p,ncol = 4)
cp_allcelltypes

ggsave(filename = paste0(results_path,"celltypes_miniUmaps.png"),
       plot = cp_allcelltypes, "png",dpi=300,width=200,height = 300,units="mm")
ggsave(filename = paste0(results_path,"celltypes_miniUmaps.pdf"),
       plot = cp_allcelltypes, "pdf",dpi=300,width=200,height =300,units="mm")

