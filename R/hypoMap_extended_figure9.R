##########
### Load & Prepare
##########

results_path_extended_figure9 = "figure_outputs/figure_extended_9/"
system(paste0("mkdir -p ",results_path_extended_figure9))


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
### Overview plot
##########

hypoMap_v2_seurat@meta.data$dummy=NA
bg_plot = DimPlot(hypoMap_v2_seurat,group.by = "dummy",reduction = paste0("umap_","scvi"),cols=getOkabeItoPalette(14),na.value = bg_col,
                  label = F,label.size = 7,raster = TRUE,pt.size = seurat_pt_size,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoLegend()+ggtitle("")+NoAxes()

ggsave(filename = paste0(results_path_extended_figure9,"hypoMap_grey_background.png"),
       plot = bg_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure9,"hypoMap_grey_background.pdf"),
       plot = bg_plot, "pdf",dpi=600,width=300,height = 300,units="mm")

##########
### Plots of bacTRAP gene combinations: Pnoc
##########

pnoc_gene_combinations = list(c("Pnoc","Sst"),c("Pnoc","Sst","Unc13c"),c("Pnoc","Sst","Nts"),c("Pnoc","Crabp1"),c("Pnoc","Crabp1","Htr3b"),c("Pnoc","Crabp1","Tmem215"))
all_pnoc_comb_plots = list()
for( i in 1:length(pnoc_gene_combinations)){
  current_comb = pnoc_gene_combinations[[i]]
  current_comb_name = paste0(current_comb,collapse = "_")
  hypoMap_v2_seurat@meta.data$tmp_score = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = current_comb)
  umap_plot =FeaturePlot(hypoMap_v2_seurat,features = "tmp_score",order=TRUE,cols = c(bg_col,"#c96410"),
                         raster = F)+NoAxes()+
    ggtitle(paste0(current_comb_name))#
  all_pnoc_comb_plots[[current_comb_name]] = umap_plot
}

names(all_pnoc_comb_plots)

# "Pnoc_Sst"
current_comb_name = "Pnoc_Sst"
current_plot = all_pnoc_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")


# "Pnoc_Sst_Unc13c"
current_comb_name = "Pnoc_Sst_Unc13c"
current_plot = all_pnoc_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")

# "Pnoc_Sst_Nts"
current_comb_name = "Pnoc_Sst_Nts"
current_plot = all_pnoc_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")

# "Pnoc_Crabp1"
current_comb_name = "Pnoc_Crabp1"
current_plot = all_pnoc_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")

# "Pnoc_Crabp1_Htr3b"
current_comb_name = "Pnoc_Crabp1_Htr3b"
current_plot = all_pnoc_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")

# "Pnoc_Crabp1_Tmem215"
current_comb_name = "Pnoc_Crabp1_Tmem215"
current_plot = all_pnoc_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure9,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")

##########
### source data
##########

#source
source_ext_figure9_bg_plot =bg_plot$data
# get all scores:
pnoc_gene_combinations = list(c("Pnoc","Sst"),c("Pnoc","Sst","Unc13c"),c("Pnoc","Sst","Nts"),c("Pnoc","Crabp1"),c("Pnoc","Crabp1","Htr3b"),c("Pnoc","Crabp1","Tmem215"))

for( i in 1:length(pnoc_gene_combinations)){
  current_comb = pnoc_gene_combinations[[i]]
  current_comb_name = paste0(current_comb,collapse = "_")
  source_ext_figure9_bg_plot$tmp_score = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = current_comb)
  colnames(source_ext_figure9_bg_plot)[colnames(source_ext_figure9_bg_plot) == "tmp_score"] = current_comb_name
}

data.table::fwrite(source_ext_figure9_bg_plot ,paste0(results_path_extended_figure9,"source_ext_figure9.txt"),sep="\t")

