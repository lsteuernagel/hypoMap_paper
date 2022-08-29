##########
### Load & Prepare
##########

#set path
results_path_extended_figure3 = "figure_outputs/figure_extended_3/"
system(paste0("mkdir -p ",results_path_extended_figure3))

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
### Figure 1b
##########

# DiscretePalette(5, palette = "alphabet")

# prepare column for campbell cluster names (with others NA)
campbell_names= names(table(hypoMap_v2_seurat@meta.data$Author_CellType[hypoMap_v2_seurat@meta.data$Dataset=="CampbellDropseq"]))[table(hypoMap_v2_seurat@meta.data$Author_CellType[hypoMap_v2_seurat@meta.data$Dataset=="CampbellDropseq"]) > 5]
campbell_names = campbell_names[!grepl("NA_",campbell_names)]

hypoMap_v2_seurat@meta.data$campbell_anno_col=NA
hypoMap_v2_seurat@meta.data$campbell_anno_col[hypoMap_v2_seurat@meta.data$Dataset=="CampbellDropseq"] = hypoMap_v2_seurat@meta.data$Author_CellType[hypoMap_v2_seurat@meta.data$Dataset=="CampbellDropseq"]
hypoMap_v2_seurat@meta.data$campbell_anno_col[! hypoMap_v2_seurat@meta.data$campbell_anno_col %in% campbell_names] =NA
hypoMap_v2_seurat@meta.data$campbell_anno_col = gsub("_Neurons[0-9]","",hypoMap_v2_seurat@meta.data$campbell_anno_col)

# plot
# hypoMap_v2_seurat_original_order_cells = hypoMap_v2_seurat@meta.data$Cell_ID
# hypoMap_v2_seurat@meta.data = hypoMap_v2_seurat@meta.data[order(hypoMap_v2_seurat@meta.data$campbell_anno_col,na.last = FALSE),]
# hypoMap_v2_seurat@meta.data$campbell_anno_col = factor(hypoMap_v2_seurat@meta.data$campbell_anno_col,levels = levels(factor(hypoMap_v2_seurat@meta.data$campbell_anno_col)))
# Idents(hypoMap_v2_seurat) = "campbell_anno_col"
campbell_anno_plot=DimPlot(hypoMap_v2_seurat,group.by = "campbell_anno_col",reduction = paste0("umap_","scvi"),pt.size = seurat_pt_size,raster = F,cols = getOkabeItoPalette(62),
                           label = TRUE,label.size = 5.5,repel = TRUE,order = TRUE,na.value = bg_col,raster.dpi = c(rasterize_px,rasterize_px))+
  NoLegend()+NoAxes()+ggtitle("Campbell ARH celltypes")

# change order -- old: stopped working properly at some point when scattermore stopped respecting the order in campbell_anno_plot$data
plot_data = campbell_anno_plot$data
campbell_anno_plot$data = plot_data[order(plot_data$campbell_anno_col,na.last = FALSE),]
#!!!!! make sure cattermore is < 0.8 else it casues problems with sorting!
campbell_anno_plot_r = scUtils::rasterize_ggplot(campbell_anno_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
campbell_anno_plot_r

ggsave(filename = paste0(results_path_extended_figure3,"campbell_annotations.png"),
       plot = campbell_anno_plot_r, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure3,"campbell_annotations.pdf"),
       plot = campbell_anno_plot_r, "pdf",dpi=600,width=300,height = 300,units="mm")

# source files
ext_figure3_sourcedata =  campbell_anno_plot_r$data %>% dplyr::mutate(Cell_ID = rownames(campbell_anno_plot_r$data)) #%>% dplyr::select(Cell_ID, campbell_anno_col),by="Cell_ID")

##########
### VIP inlay plot
##########

#plot
Idents(hypoMap_v2_seurat) <- "C66"
seurat_vip_subset = subset(hypoMap_v2_seurat,subset = C66 %in% c("C66-31","C66-29") & umapscvi_1 > 4 & umapscvi_1 < 9 & umapscvi_2 < -2 & umapscvi_2 > -9)
#seurat_vip_subset = subset(seurat_vip_subset,subset = umapscvi_1 > 1 & umapscvi_2 < -0.5)
seurat_vip_subset@meta.data$C185_named_wo = stringr::str_remove(seurat_vip_subset@meta.data$C185_named,pattern = "C185-[0-9]+: ")
vip_small_plot=Seurat::DimPlot(seurat_vip_subset,group.by = "C185_named_wo",label = TRUE,label.size = 7,repel = TRUE,
                               order = TRUE,na.value = bg_col,raster = TRUE,raster.dpi = c(1536,1536),pt.size = 3, cols = getOkabeItoPalette(7))+
  NoLegend()+NoAxes()+ggtitle("SCN neurons reference map")
#vip_small_plot = rasterize_ggplot(vip_small_plot,pixel_raster = 1536,pointsize = 1.8)
vip_small_plot

ggsave(filename = paste0(results_path_extended_figure3,"vip_small_plot.png"),
       plot = vip_small_plot, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure3,"vip_small_plot.pdf"),
       plot = vip_small_plot, "pdf",dpi=450,width=200,height = 200,units="mm")

# source files
ext_figure3_sourcedata =  dplyr::left_join(ext_figure3_sourcedata,vip_small_plot$data %>% dplyr::mutate(Cell_ID = rownames(vip_small_plot$data)) %>% dplyr::select(Cell_ID, C185_named_scn = C185_named_wo),by="Cell_ID")

data.table::fwrite(ext_figure3_sourcedata,paste0(results_path_extended_figure3,"source_ext_figure3_umap_data.txt"),sep="\t")


# 
# ##########
# ### Calculate Neuron subset UMAP
# ##########
# 
# hypoMap_v2_neurons = subset(hypoMap_v2_seurat, C2_named == "C2-1: Neurons")
# 
# hypoMap_v2_neurons = RunUMAP(hypoMap_v2_neurons,
#                              reduction = "scvi",
#                              seed.use= 1234567,
#                              dims=1:ncol(hypoMap_v2_seurat@reductions[["scvi"]]@cell.embeddings),
#                              reduction.name=paste0("umap_scvi_neurons"),
#                              reduction.key = paste0("umap_scvi_neurons"),
#                              verbose=TRUE,
#                              n.neighbors =25,
#                              return.model = TRUE)
# 
# DimPlot(hypoMap_v2_neurons,group.by = "C185_named",reduction = "umap_scvi_neurons",raster.dpi = c(2048,2048),pt.size = 1.5)+NoLegend()+NoAxes()
# DimPlot(hypoMap_v2_neurons,group.by = "C185",reduction = "umap_scvi_neurons",label=TRUE,label.size = 3,raster.dpi = c(2048,2048),pt.size = 1.5)+NoLegend()+NoAxes()
# 
# FeaturePlot(hypoMap_v2_neurons,features = "Slc32a1",reduction = "umap_scvi_neurons",raster.dpi = c(2048,2048),pt.size = 1.5)+NoLegend()+NoAxes()
# 
# # save
# #saveRDS(hypoMap_v2_neurons,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/hypoMap_v2_neurons.rds")
# #hypoMap_v2_neurons = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/hypoMap_v2_neurons.rds")
# 
# ##########
# ### Figure Xa:neuron subset UMAP
# ##########
# 
# # factor:
# hypoMap_v2_neurons@meta.data$C66_named = factor(hypoMap_v2_neurons@meta.data$C66_named,levels = unique(hypoMap_v2_neurons@meta.data$C66_named))
# 
# neuron_umap_plot_c66 = DimPlot(hypoMap_v2_neurons,group.by = "C66_named",reduction = paste0("umap_scvi_neurons"),cols=getOkabeItoPalette(50),
#                            label = TRUE,label.size = 4,raster = TRUE,pt.size = seurat_pt_size,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
#   theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# neuron_umap_plot_c66
# 
# 
# ggsave(filename = paste0(results_path_extended_figure3,"neuron_umap_C66.png"),
#        plot = neuron_umap_plot_c66, "png",dpi=600,width=300,height = 300,units="mm")
# ggsave(filename = paste0(results_path_extended_figure3,"neuron_umap_C66.pdf"),
#        plot = neuron_umap_plot_c66, "pdf",dpi=600,width=300,height = 300,units="mm")
# 
# 
# ##########
# ### Figure Xb:neuron subset UMAP
# ##########
# 
# hypoMap_v2_neurons@meta.data$C185_named = factor(hypoMap_v2_neurons@meta.data$C185_named,levels = unique(hypoMap_v2_neurons@meta.data$C185_named))
# 
# neuron_umap_plot_C185 = DimPlot(hypoMap_v2_neurons,group.by = "C185_named",reduction = paste0("umap_scvi_neurons"),cols=getOkabeItoPalette(130),
#                                label = TRUE,label.size = 3,raster = TRUE,pt.size = seurat_pt_size,repel = TRUE,raster.dpi = c(rasterize_px,rasterize_px))+NoAxes()+NoLegend()+
#   theme(text = element_text(size=text_size))+ggtitle("Cell types on hypothalamus reference map")
# neuron_umap_plot_C185
# 
# 
# ggsave(filename = paste0(results_path_extended_figure3,"neuron_umap_C185.png"),
#        plot = neuron_umap_plot_C185, "png",dpi=600,width=300,height = 300,units="mm")
# ggsave(filename = paste0(results_path_extended_figure3,"neuron_umap_C185.pdf"),
#        plot = neuron_umap_plot_C185, "pdf",dpi=600,width=300,height = 300,units="mm")
# 
# ##########
# ### source data
# ##########
# 
# # source: 
# source_ext_figure3 =neuron_umap_plot_c66$data
# source_ext_figure3$C185_named = neuron_umap_plot_C185$data$C185_named
# data.table::fwrite(source_ext_figure3,paste0(results_path_extended_figure3,"source_ext_figure3_umap_data.txt"),sep="\t")



