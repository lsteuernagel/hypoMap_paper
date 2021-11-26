##########
### Load & Prepare
##########

#set path
results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_supplementary_5/"
system(paste0("mkdir -p ",results_path))
# load everything required
source("load_data.R")

# path with output files
data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_input/"

## load mapped object
query_snseq_neurons = readRDS(paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/nucseq_neurons_map.rds"))
## load mapped object
query_snseq_all = readRDS(paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/nucseq_all_map.rds"))

## plotting
rasterize_point_size = 2.2
rasterize_pixels = 2048

##########
### Original UMAP
##########

original_umap = data.table::fread(paste0(data_path,"nucSeq_originalUMAP.txt"),data.table = FALSE)
rownames(original_umap) = query_snseq_all@meta.data$Cell_ID
original_umap_dimred = CreateDimReducObject(embeddings = as.matrix(original_umap),key = "original_umap")

# add back in
query_snseq_all@reductions[["original_umap"]] = original_umap_dimred
# make plot
original_full_umap  = DimPlot(query_snseq_all,group.by = "Cluster_IDs",reduction = "original_umap",label = TRUE,label.size = 5,repel = TRUE,raster = F,pt.size = 0.2)+NoAxes()+NoLegend()+
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
### Mapping probability & cluster representation
##########

# per request from brian:
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


##########
### TODO: optionally correlations catter plot --> at the moment not included!
##########

# TODO: import valid genes from Figure 3 script !!!
# shared_genes = ....

# # load marker genes sc seq
# sc_seq_markers_K169 = neuron_map_seurat@misc$markers_comparisons_all[neuron_map_seurat@misc$markers_comparisons_all$p_val_adj<0.01,]
# 
# # load marker genes sn seq
# output_file_name = paste0(data_path,"sn_seq_mapped_neurons_K169_markers_2_sampleID.txt")
# sn_seq_markers_K169 = data.table::fread(output_file_name,data.table = F)
# sn_seq_markers_K169$specificity = (sn_seq_markers_K169$pct.1 / sn_seq_markers_K169$pct.2) * sn_seq_markers_K169$avg_logFC
# sn_seq_markers_K169$avg_log2FC = sn_seq_markers_K169$avg_logFC
# 
# 
# # defien cluster
# current_cluster= "K169-99"# "K169-99"#"K169-1"  #"K169-53"  #"K169-32"  #"K169-1"  #"K169-141"
# current_cluster= "K169-1"# "K169-99"#"K169-1"  #"K169-53"  #"K169-32"  #"K169-1"  #"K169-141"
# 
# neuron_map_seurat@misc$annotations$cluster_name[neuron_map_seurat@misc$annotations$cluster_id == current_cluster]
# 
# # define marker genes
# sn_Seq_based_markers = sn_seq_markers_K169$gene[sn_seq_markers_K169$p_val_adj<padj_cut & sn_seq_markers_K169$avg_log2FC>fc_min &
#                                                   sn_seq_markers_K169$cluster==current_cluster & sn_seq_markers_K169$pct.2<0.5]
# sc_Seq_based_markers = sc_seq_markers_K169$gene[sc_seq_markers_K169$p_val_adj<padj_cut & sc_seq_markers_K169$avg_logFC>fc_min &
#                                                   sc_seq_markers_K169$cluster==current_cluster & sc_seq_markers_K169$pct.2<0.5]
# union_markers = union(sn_Seq_based_markers,sc_Seq_based_markers)
# #union_markers = sn_Seq_based_markers
# union_markers = union_markers[union_markers %in% shared_genes]
# 
# # get mean
# mean_expr_sc = data.frame(mean_expr = neuron_map_mean_expression[union_markers,current_cluster],gene = union_markers)
# mean_expr_sn = data.frame(mean_expr = snseq_mean_expression[union_markers,current_cluster],gene = union_markers)
# 
# #### compare mean between the two:
# compare_mean_between = dplyr::left_join(mean_expr_sc,mean_expr_sn,by="gene",suffix=c("_sc","_sn"))
# cor(compare_mean_between$mean_expr_sc,compare_mean_between$mean_expr_sn,use = "pairwise.complete.obs",method="pearson")
# 
# # add marker info:
# compare_mean_between$classification[compare_mean_between$gene %in% sn_Seq_based_markers & !compare_mean_between$gene %in% sc_Seq_based_markers] = "Marker in sc-seq"
# compare_mean_between$classification[!compare_mean_between$gene %in% sn_Seq_based_markers & compare_mean_between$gene %in% sc_Seq_based_markers] = "Marker in sn-seq"
# compare_mean_between$classification[compare_mean_between$gene %in% sn_Seq_based_markers & compare_mean_between$gene %in% sc_Seq_based_markers] = "Marker in both"
# 
# #plot
# ggplot(compare_mean_between,aes(mean_expr_sc,mean_expr_sn))+geom_point(aes(color=classification),size=1)+theme(text = element_text(size=20))+
#   geom_smooth(method = "lm",color="grey40")
# 





