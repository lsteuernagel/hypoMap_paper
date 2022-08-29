##########
### Load & Prepare
##########

results_path_extended_figure8 = "figure_outputs/figure_extended_8/"
system(paste0("mkdir -p ",results_path_extended_figure8))


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

ggsave(filename = paste0(results_path_extended_figure8,"hypoMap_grey_background.png"),
       plot = bg_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,"hypoMap_grey_background.pdf"),
       plot = bg_plot, "pdf",dpi=600,width=300,height = 300,units="mm")

##########
### Plots of bacTRAP gene combinations: Glp1r
##########

glp1r_gene_combinations = list(c("Glp1r","Ghrh"),c("Glp1r","Pomc"),c("Glp1r","Pomc","Anxa2"),c("Glp1r","Sst"),c("Glp1r","Sst","Unc13c"),
                               c("Glp1r","Tbx19","Anxa2"),c("Glp1r","Nkx2-4","Trh"),c("Glp1r","Oxt"))
all_glp1r_comb_plots = list()
for( i in 1:length(glp1r_gene_combinations)){
  current_comb = glp1r_gene_combinations[[i]]
  current_comb_name = paste0(current_comb,collapse = "_")
  hypoMap_v2_seurat@meta.data$tmp_score = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = current_comb)
  umap_plot =FeaturePlot(hypoMap_v2_seurat,features = "tmp_score",order=TRUE,cols = c(bg_col,"#c96410"),
                         raster = F)+NoAxes()+
    ggtitle(paste0(current_comb_name))
  all_glp1r_comb_plots[[current_comb_name]] = umap_plot
}

#### need to plot them 1 by 1 to select zoom for each !
names(all_glp1r_comb_plots)

# "Glp1r_Pomc"
current_comb_name = "Glp1r_Pomc"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(0.5,4.5)+ylim(1,5)#+ggtitle(current_comb_name)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 4.3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 200,units="mm")

#source
source_ext_figure8_current =current_plot$data
colnames(source_ext_figure8_current)[4] = current_comb_name
data.table::fwrite(source_ext_figure8_current ,paste0(results_path_extended_figure8,"source_figure8_",current_comb_name,".txt"),sep="\t")

# "Glp1r_Pomc_Anxa2"
current_comb_name = "Glp1r_Pomc_Anxa2"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(0.5,4.5)+ylim(1,5)#+ggtitle(current_comb_name)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 4.3)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 200,units="mm")
#source
source_ext_figure8_current =current_plot$data
colnames(source_ext_figure8_current)[4] = current_comb_name
data.table::fwrite(source_ext_figure8_current ,paste0(results_path_extended_figure8,"source_figure8_",current_comb_name,".txt"),sep="\t")

# "Glp1r_Ghrh"
current_comb_name = "Glp1r_Ghrh"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3.5)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")
#source
source_ext_figure8_current =current_plot$data
colnames(source_ext_figure8_current)[4] = current_comb_name
data.table::fwrite(source_ext_figure8_current ,paste0(results_path_extended_figure8,"source_figure8_",current_comb_name,".txt"),sep="\t")

# "Glp1r_Sst"
current_comb_name = "Glp1r_Sst"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3.5)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")
#source
source_ext_figure8_current =current_plot$data
colnames(source_ext_figure8_current)[4] = current_comb_name
data.table::fwrite(source_ext_figure8_current ,paste0(results_path_extended_figure8,"source_figure8_",current_comb_name,".txt"),sep="\t")

# "Glp1r_Sst_Unc13c"
current_comb_name = "Glp1r_Sst_Unc13c"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3.5)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")
#source
source_ext_figure8_current =current_plot$data
colnames(source_ext_figure8_current)[4] = current_comb_name
data.table::fwrite(source_ext_figure8_current ,paste0(results_path_extended_figure8,"source_figure8_",current_comb_name,".txt"),sep="\t")

# "Glp1r_Tbx19_Anxa2"
current_comb_name = "Glp1r_Tbx19_Anxa2"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3.5)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")
#source
source_ext_figure8_current =current_plot$data
colnames(source_ext_figure8_current)[4] = current_comb_name
data.table::fwrite(source_ext_figure8_current ,paste0(results_path_extended_figure8,"source_figure8_",current_comb_name,".txt"),sep="\t")

# "Glp1r_Nkx2"
current_comb_name = "Glp1r_Nkx2-4_Trh"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(-6,0)+ylim(-11,-2)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3.5)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 260,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 260,units="mm")
#source
source_ext_figure8_current =current_plot$data
colnames(source_ext_figure8_current)[4] = current_comb_name
data.table::fwrite(source_ext_figure8_current ,paste0(results_path_extended_figure8,"source_figure8_",current_comb_name,".txt"),sep="\t")

# "Glp1r_Oxt"
current_comb_name = "Glp1r_Oxt"
current_plot = all_glp1r_comb_plots[[current_comb_name]]+xlim(-3,0)+ylim(6.5,9.5)#+ggtitle(current_comb_name)
rasterized_plot = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = 3.5)
rasterized_plot

ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.png"),
       plot = rasterized_plot, "png",dpi=450,width=220,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure8,current_comb_name,"_zoom_umap.pdf"),
       plot = rasterized_plot, "pdf",dpi=450,width=220,height = 200,units="mm")

##########
### source data
##########

#source
source_ext_figure8_bg_plot =bg_plot$data
# get all scores:
glp1r_gene_combinations = list(c("Glp1r","Ghrh"),c("Glp1r","Pomc"),c("Glp1r","Pomc","Anxa2"),c("Glp1r","Sst"),c("Glp1r","Sst","Unc13c"),
                               c("Glp1r","Tbx19","Anxa2"),c("Glp1r","Nkx2-4","Trh"),c("Glp1r","Oxt"))
for( i in 1:length(glp1r_gene_combinations)){
  current_comb = glp1r_gene_combinations[[i]]
  current_comb_name = paste0(current_comb,collapse = "_")
  source_ext_figure8_bg_plot$tmp_score = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = current_comb)
  colnames(source_ext_figure8_bg_plot)[colnames(source_ext_figure8_bg_plot) == "tmp_score"] = current_comb_name
}

data.table::fwrite(source_ext_figure8_bg_plot ,paste0(results_path_extended_figure8,"source_ext_figure8.txt"),sep="\t")


