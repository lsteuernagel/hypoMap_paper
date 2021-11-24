##########
### Load & Prepare
##########

#set path
results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_supplementary_7/"
system(paste0("mkdir -p ",results_path))
# load everything required
source("scripts/paper_figures_new/load_data.R")
source("utils.R")

# path with output files
data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_input/"

# subsample ids for plotting
subsample_ids = data.table::fread(paste0(data_path,"_subsampled_Cell_IDs_neuronMap.txt"),data.table = F,header = F)[,1]

## load mapped object
query_snseq_neurons = readRDS(paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/nucseq_neurons_map.rds"))

## plotting
rasterize_point_size = 2.2
rasterize_pixels = 2048
colorvec = RColorBrewer::brewer.pal(9, "Blues")
colorvec[1] =  "#dedede"

##########
### Figure 7: Glp1r
##########


### Glp1r expression in neuron map
glp1r_neuron_map = FeaturePlot(neuron_map_seurat,features = "Glp1r",reduction =paste0("umap_scvi"), combine = TRUE,order = TRUE,pt.size = 0.2)+
  NoAxes() +scale_color_gradientn(colours = colorvec)  #+ scale_color_gradient(low="lightgrey",high = colorvec[9])
  #scale_color_gradientn(colours = colorvec) #+ NoLegend() #+scale_color_gradient(low="lightgrey",high=max_color)
glp1r_neuron_map = rasterize_ggplot(glp1r_neuron_map,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
glp1r_neuron_map

#save
ggsave(filename = paste0(results_path,"glp1r_neuron_map.png"),
       plot = glp1r_neuron_map, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path,"glp1r_neuron_map.pdf"),
       plot = glp1r_neuron_map, "pdf",dpi=600,width=350,height = 300,units="mm")


### Glp1r expression in nucseq
glp1r_snseq =FeaturePlot(query_snseq_neurons,features = "Glp1r",reduction =paste0("umap_scvi"), combine = TRUE,order = TRUE,pt.size = 0.2)+
  NoAxes() +scale_color_gradientn(colours = colorvec)   #+ scale_color_gradientn(colours = colorvec) #+ NoLegend() #+scale_color_gradient(low="lightgrey",high=max_color)
glp1r_snseq = rasterize_ggplot(glp1r_snseq,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
glp1r_snseq

#save
ggsave(filename = paste0(results_path,"glp1r_snseq.png"),
       plot = glp1r_snseq, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path,"glp1r_snseq.pdf"),
       plot = glp1r_snseq, "pdf",dpi=600,width=350,height = 300,units="mm")

### glp1r in bacTRAp
# USE unlabelled plot from figure 3 script !
system(paste0("cp ",gsub("figure_supplementary_7","figure_3",results_path),"bacTRAP_glp1r_rbo_plot.pdf"," ",results_path))


# selected clusters !
anno_df = neuron_map_seurat@misc$annotations
selected_cluster_K169_named = c("Ttr.Pomc.Tcf7l2.Tbx3.HY1","Sstr1.Sst.Il1rapl2.Otp.HY1","Ghrh.Gsx1.Hmx2.HY1","Crabp1.Nfib.Six3.Hmx2.HY1","Trh.Nkx2-4.Gsx1.Hmx2.HY1","Ebf3.Oxt")
# make a subset 
subset_seurat_glp1r = subset(neuron_map_seurat, K169_named %in% selected_cluster_K169_named)
# only keep labels for celltypes with mapped neuron
subset_seurat_glp1r@meta.data$label_col = subset_seurat_glp1r@meta.data$K169_named
neuron_map_seurat@meta.data$label_col = as.factor(neuron_map_seurat@meta.data$K169_named)
neuron_map_seurat@meta.data$label_col[!neuron_map_seurat@meta.data$K169_named %in% unique(subset_seurat_glp1r@meta.data$label_col)] = NA
# make plot using mapscvi:
selected_clusters_plot = mapscvi::plot_query_labels(query_seura_object=subset_seurat_glp1r,reference_seurat=neuron_map_seurat,label_col="label_col",label_col_query="label_col", 
                                          overlay = TRUE,query_pt_size = 0.1,labelonplot = TRUE,label.size=6,repel=TRUE)+ggtitle("Selected clusters")
selected_clusters_plot = rasterize_ggplot(selected_clusters_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
selected_clusters_plot

# #save
# ggsave(filename = paste0(results_path,"selected_glp1r_clusters_plot.png"),
#        plot = selected_clusters_plot, "png",dpi=600,width=350,height = 300,units="mm")
# ggsave(filename = paste0(results_path,"selected_glp1r_clusters_plot"),
#        plot = selected_clusters_plot, "pdf",dpi=600,width=350,height = 300,units="mm")


# use percentages from RNAscope
ish_quantification = data.table::fread(paste0(gsub("figure_supplementary_7","figure_6",results_path),"pct_expressed_cells_clusters_ISH.txt"),data.table = FALSE)
ish_quantification_summary = ish_quantification %>% group_by(Experiment) %>% dplyr::summarise(mean_expressed_cells = mean(total_2_pct))
ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Ghrh"] = "Ghrh.Gsx1.Hmx2.HY1"
ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Oxt"] = "Ebf3.Oxt"
ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Pomc"] = "Ttr.Pomc.Tcf7l2.Tbx3.HY1"
ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Pomc Anxa2"] = "Anxa2.Pomc.Tcf7l2.Tbx3.HY1"
ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Sst Unc13c"] = "Sstr1.Sst.Il1rapl2.Otp.HY1"
ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Tbx19 Anxa2"] = "Crabp1.Nfib.Six3.Hmx2.HY1"
ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Trh Nkx2-4"] = "Trh.Nkx2-4.Gsx1.Hmx2.HY1"
#ish_quantification_summary$celltype[ish_quantification_smmary$Experiment == "Npw Nkx2-4"] = "Ebf1.Nkx2-4.Gsx1.Hmx2.HY1"

neuron_map_seurat@meta.data = neuron_map_seurat@meta.data[,!colnames(neuron_map_seurat@meta.data) %in% c("Experiment","mean_expressed_cells","new_label_col")]
neuron_map_seurat@meta.data = dplyr::left_join(neuron_map_seurat@meta.data,ish_quantification_summary, by=c("K169_named"="celltype"))
rownames(neuron_map_seurat@meta.data) = neuron_map_seurat@meta.data$Cell_ID

colorvec = RColorBrewer::brewer.pal(9, "Blues")
#colorvec[1] =  "lightgrey"
neuron_map_seurat@meta.data$new_label_col = NA
neuron_map_seurat@meta.data$new_label_col[!is.na(neuron_map_seurat@meta.data$mean_expressed_cells)] = neuron_map_seurat@meta.data$K169_named[!is.na(neuron_map_seurat@meta.data$mean_expressed_cells)]
Idents(neuron_map_seurat) = "new_label_col"
ish_pct_umap_plot = FeaturePlot(neuron_map_seurat,features = "mean_expressed_cells",label = TRUE,label.size = 5,repel = TRUE)+NoAxes()+
  scale_color_gradientn(colours = colorvec,na.value = "lightgrey",limits=c(0,100))+ggtitle("Pct of expressing cells in ISH")
ish_pct_umap_plot

#save
ggsave(filename = paste0(results_path,"ish_pct_umap_plot.png"),
       plot = ish_pct_umap_plot, "png",dpi=600,width=350,height = 300,units="mm")
ggsave(filename = paste0(results_path,"ish_pct_umap_plot.pdf"),
       plot = ish_pct_umap_plot, "pdf",dpi=600,width=350,height = 300,units="mm")
