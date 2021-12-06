
##########
### Load & Prepare
##########

#set path
results_path = "figure_outputs/figure_supplementary_6/"
system(paste0("mkdir -p ",results_path))

# load required functions
require(mapscvi)
source("R/utility_functions.R")
source("R/plot_functions.R")

# where to find large data objects (need to change to local dir !)
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

# load seurat objects via large_data_path
load_required_files(large_data_path = large_data_path)

## load mapped object full
query_snseq_all = readRDS(paste0(large_data_path,"nucseq_all_map.rds"))

## plotting
rasterize_point_size = 2.2
rasterize_pixels = 2048

##########
### Original UMAP
##########

original_umap = data.table::fread(paste0("data_inputs/nucSeq_originalUMAP.txt"),data.table = FALSE)
rownames(original_umap) = query_snseq_all@meta.data$Cell_ID
original_umap_dimred = CreateDimReducObject(embeddings = as.matrix(original_umap),key = "original_umap")

# add back in
query_snseq_all@reductions[["original_umap"]] = original_umap_dimred
# make plot
original_full_umap  = DimPlot(query_snseq_all,group.by = "Cluster_IDs",reduction = "original_umap",label = TRUE,label.size = 3.5,repel = TRUE,raster = F,pt.size = 0.2)+NoAxes()+NoLegend()+
  theme(text = element_text(size=20))+ggtitle("NucSeq data - independent analysis")+guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
original_full_umap = rasterize_ggplot(original_full_umap,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
original_full_umap

# save
ggsave(filename = paste0(results_path,"snseq_original_full_umap.png"),
       plot = original_full_umap, "png",dpi=400,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path,"snseq_original_full_umap.pdf"),
       plot = original_full_umap, "pdf",dpi=400,width=300,height =300,units="mm")


##########
### Original UMAP
##########

projected_full_umap  = DimPlot(query_snseq_all,group.by = "predicted_Curated_Class",reduction = "umap_scvi",label = TRUE,label.size = 7,repel = TRUE,raster = F,pt.size = 0.2)+NoAxes()+NoLegend()+
  theme(text = element_text(size=20))+ggtitle("NucSeq data - all celltypes")+guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
projected_full_umap = rasterize_ggplot(projected_full_umap,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
projected_full_umap

# save
ggsave(filename = paste0(results_path,"snseq_projected_full_umap.png"),
       plot = projected_full_umap, "png",dpi=400,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path,"snseq_projected_full_umap.pdf"),
       plot = projected_full_umap, "pdf",dpi=400,width=300,height =300,units="mm")


##########
### Mapping probability & cluster representation -- not included at the moment !!
##########

# Make one with mapping probability vs representation of cluster in sc and sn seq
# probably could be added to this supplementary figure (5)
# build on earlier stuff included in mapscvi !?
# do for query neurons!

#mapscvi::compare_clustering(query_seura_object = query_snseq_neurons,)
mapscvi::plot_propagation_stats(query_seura_object = query_snseq_neurons,reference_seurat = neuron_map_seurat,
                                label_col_query = "predicted_K169_named",label_col = "K169_named") # show plot
propagation_stats_sn= mapscvi::plot_propagation_stats(query_seura_object = query_snseq_neurons,reference_seurat = neuron_map_seurat,
                                                      label_col_query = "predicted_K169_named",label_col = "K169_named",return_data = TRUE) # get data

# compare relative log2 ratio and mapping probability (avergaged per cluster !?)
probability_per_cluster = query_snseq_neurons@meta.data[,c("predicted_K169_named","prediction_probability")] %>% dplyr::group_by(predicted_K169_named) %>%
  dplyr::summarise(avg_prob = mean (prediction_probability))
# join
probability_per_cluster = dplyr::left_join(probability_per_cluster,propagation_stats_sn[,c("group","reference_count","query_count","reference_pct","query_pct","relative_log2_ratio")], by=c("predicted_K169_named"="group"))
ggplot(probability_per_cluster[probability_per_cluster$query_count>30,],aes(relative_log2_ratio,avg_prob,color=query_pct))+geom_point(size=2)+
  scale_color_gradient2(low="white",mid=colorvec[4],high = colorvec[9],midpoint = 0.1)+
  geom_vline(xintercept = 0,color="grey60")+theme(text=element_text(size=20))+
  geom_smooth(method = "lm",color = "grey60")

cor(probability_per_cluster$avg_prob[probability_per_cluster$query_count>30],probability_per_cluster$relative_log2_ratio[probability_per_cluster$query_count>30])

ggplot(probability_per_cluster,aes(relative_log2_ratio,query_pct,color=avg_prob))+geom_point()+
  scale_color_gradient(low=colorvec[2],high = colorvec[9])+
  geom_vline(xintercept = 0,color="grey60")

