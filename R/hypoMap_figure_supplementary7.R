##########
### Load & Prepare
##########

#set path
results_path = "figure_outputs/figure_supplementary_7/"
system(paste0("mkdir -p ",results_path))
# load everything required
source("R/load_data.R")
source("R/plot_functions.R")
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

## load mapped object
query_snseq_neurons = readRDS(paste0(large_data_path,"nucseq_neurons_map.rds"))

## plotting
rasterize_point_size = 2.2
rasterize_pixels = 2048
cols_for_feature_plot = c("#dedede","#0b3ebd") # "#0b3ebd"


##########
###  DEGs in all clusters
##########

Idents(query_snseq_neurons) = "predicted_K169_named"
all_clusters = unique(query_snseq_neurons@meta.data$predicted_K169_named)
all_conditionGenes_list = list()
for(i in 1:length(all_clusters)){
  current_cluster = all_clusters[i]
  message(current_cluster)
  cells_adlib = length(query_snseq_neurons@meta.data$predicted_K169_named[query_snseq_neurons@meta.data$predicted_K169_named == current_cluster & query_snseq_neurons@meta.data$Diet=="adlib"])
  cells_fasting = length(query_snseq_neurons@meta.data$predicted_K169_named[query_snseq_neurons@meta.data$predicted_K169_named == current_cluster & query_snseq_neurons@meta.data$Diet=="fast"])
  if(cells_adlib >= min_cells & cells_fasting >= min_cells ){
    conditionGenes_current = Seurat::FindMarkers(query_snseq_neurons, ident.1 = "adlib",ident.2 = "fast" , group.by = "Diet", 
                                                 subset.ident = current_cluster,min.pct = 0.1,logfc.threshold = 0.25,max.cells.per.ident = 5000)
    conditionGenes_current$gene = rownames(conditionGenes_current) # add gene name 
    conditionGenes_current$pct_diff = conditionGenes_current$pct.1 - conditionGenes_current$pct.2
    conditionGenes_current$current_cluster = current_cluster
    all_conditionGenes_list[[current_cluster]] = conditionGenes_current
  }
}

# rbind
all_conditionGenes = do.call(rbind,all_conditionGenes_list)

# save !
conditionGenes_all_file = paste0(results_path,"all_clusters_fasting_DEG.txt")
data.table::fwrite(all_conditionGenes,conditionGenes_all_file,sep="\t")


##########
### plot on UMAP
##########

#all_conditionGenes = data.table::fread(conditionGenes_all_file,data.table = FALSE)

all_conditionGenes_filtered = all_conditionGenes[all_conditionGenes$p_val_adj < 0.01,] # filter pval
# how many express Fos:
fos_degs = all_conditionGenes[all_conditionGenes$gene == "Fos",]

n_cells_per_cluster = query_snseq_neurons@meta.data %>% dplyr::group_by(predicted_K169_named) %>% dplyr::count(name="n_cells_per_cluster")
n_DEG_per_cluster = all_conditionGenes_filtered %>% dplyr::group_by(current_cluster) %>% dplyr::count(name="n_DEG_per_cluster") %>% 
  dplyr::full_join(n_cells_per_cluster,by=c("current_cluster"="predicted_K169_named"))
n_DEG_per_cluster$n_DEG_per_cluster[is.na(n_DEG_per_cluster$n_DEG_per_cluster)] = 0

ggplot(n_DEG_per_cluster,aes(n_DEG_per_cluster,n_cells_per_cluster))+geom_point()

n_DEG_per_cluster$deg_score = n_DEG_per_cluster$n_DEG_per_cluster / (n_DEG_per_cluster$n_cells_per_cluster - 0)

temp_meta = dplyr::left_join(query_snseq_neurons@meta.data,n_DEG_per_cluster,by=c("predicted_K169_named"="current_cluster"))
rownames(temp_meta) = temp_meta$Cell_ID
query_snseq_neurons@meta.data = temp_meta
query_snseq_neurons@meta.data$n_DEG_per_cluster_zscore = (query_snseq_neurons@meta.data$n_DEG_per_cluster - mean(query_snseq_neurons@meta.data$n_DEG_per_cluster)) / sd(query_snseq_neurons@meta.data$n_DEG_per_cluster)

### Glp1r expression in neuron map
ndeg_sn_seq = FeaturePlot(query_snseq_neurons,features = "deg_score",reduction =paste0("umap_scvi"),cols = cols_for_feature_plot,order = TRUE,pt.size = 0.2)+
  NoAxes() 
ndeg_sn_seq = rasterize_ggplot(ndeg_sn_seq,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
ndeg_sn_seq

# save
ggsave(filename = paste0(results_path,"snseq_n_deg_fasting_umap.png"),
       plot = ndeg_sn_seq, "png",dpi=400,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path,"snseq_n_deg_fasting_umap.pdf"),
       plot = ndeg_sn_seq, "pdf",dpi=400,width=300,height =300,units="mm")

## plot z_score version (but I think above plot is better for figure)
ndeg_sn_seq = FeaturePlot(query_snseq_neurons,features = "n_DEG_per_cluster_zscore",reduction =paste0("umap_scvi"),cols = cols_for_feature_plot,order = TRUE,pt.size = 0.2)+
  NoAxes() 
ndeg_sn_seq = rasterize_ggplot(ndeg_sn_seq,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
ndeg_sn_seq


##########
### Plot fos and ieg
##########

# subset tto campbell
campbell_diet = subset(neuron_map_seurat,subset = Dataset=="Campbell" & Diet %in% c("Normal","Fasted" ))

#  feature plot
Idents(campbell_diet) <- "K31_named"
fos_plot = Seurat::FeaturePlot(campbell_diet,"Fos",split.by = "Diet",label = TRUE,label.size = 3.5,repel = TRUE,pt.size = 0.5,
                               keep.scale="feature",order = TRUE,combine = FALSE,cols = cols_for_feature_plot)
fos_plot_fasting = fos_plot[[1]]+NoAxes()+theme(panel.border = element_blank())+ylab("") +ggplot2::ggtitle("fasted") #+ scale_color_gradientn(colours = colorvec) #+scale_color_gradient(low = "lightgrey",high = "#8c390a")
fos_plot_adlib = fos_plot[[2]]+NoAxes()+theme(panel.border = element_blank())+ylab("") +ggplot2::ggtitle("adlib")#+ scale_color_gradientn(colours = colorvec) #+scale_color_gradient(low = "lightgrey",high = "#8c390a")

fos_plot_adlib = rasterize_ggplot(fos_plot_adlib,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
fos_plot_adlib
fos_plot_fasting = rasterize_ggplot(fos_plot_fasting,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
fos_plot_fasting

# save
ggsave(filename = paste0(results_path,"fos_campbell_adlib.png"),
       plot = fos_plot_adlib, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"fos_campbell_adlib.pdf"),
       plot = fos_plot_adlib, "pdf",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"fos_campbell_fasting.png"),
       plot = fos_plot_fasting, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"fos_campbell_fasting.pdf"),
       plot = fos_plot_fasting, "pdf",dpi=450,width=200,height = 200,units="mm")



