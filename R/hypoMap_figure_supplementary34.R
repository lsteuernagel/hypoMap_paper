##########
### Load & Prepare
##########

#set path
results_path = "figure_outputs/figure_supplementary_3_4/"
system(paste0("mkdir -p ",results_path))
# load everything required
source("R/load_data.R")
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

# neuron metric results
neurons_metrics = data.table::fread(paste0("data_inputs/hypothalamusMapNeurons_v4_comparison_457fc60c3c4f1911bcbc6c5d46127037.txt"),data.table = F)

text_size =20

##########
### Helper functions
##########

##### Helper function for per celltype comparison:

# function to get input for line plots!
# expects columns as embeddings and names_to_check to be in colnames! and cell_ids as rownames!
per_celltype_summary_selected = function(input_values, celltype_ids, names_to_check){
  # get df of celltype ids
  celltype_df = plyr::ldply(celltype_ids, data.frame)
  colnames(celltype_df) = c("celltype","Cell_ID")
  # get df of values of selected embeddings
  values_df = as.data.frame(input_values[,names_to_check])
  colnames(values_df) = names(names_to_check)
  # set rownames as new column for join
  values_df$Cell_ID =rownames(values_df)
  # join with celltypes
  values_df = dplyr::left_join(values_df,celltype_df,by="Cell_ID")
  # convert to long and summarise per celltype
  values_df_long = values_df %>% tidyr::gather(-celltype,-Cell_ID,key="reduction",value="score")
  values_df_long$score = as.numeric(values_df_long$score)
  values_df_summary = values_df_long %>% dplyr::filter(!is.na(celltype))  %>%  dplyr::group_by(reduction,celltype) %>% dplyr::summarise(median_score = median(score))
  values_df_summary
  return(values_df_summary)
}

##########
### Supplemental Figure 3: 
##########

##### load evaluation results
evaluation_mixing_probNorm_all = data.table::fread(paste0(large_data_path,"evaluation_mixing_probNorm.Batch_ID.20000.100.123467_all.txt"),data.table = F)
rownames(evaluation_mixing_probNorm_all) = evaluation_mixing_probNorm_all$Cell_ID
evaluation_purity_knn_all = data.table::fread(paste0("data_inputs/evaluation_purity_knn.24.20.123467_all.txt"),data.table = F)
rownames(evaluation_purity_knn_all) = evaluation_purity_knn_all$mapped_celltype
evaluation_entropy_knn_all = data.table::fread(paste0(large_data_path,"evaluation_entropy_knn_all.Batch_ID.20.123467_all.txt"),data.table = F)
rownames(evaluation_entropy_knn_all) = evaluation_entropy_knn_all$Cell_ID

##### Define which integration results should be compared
# add information
neurons_metrics$assay=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[a-zA-Z]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
neurons_metrics$ndim=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")%>% as.numeric()
neurons_metrics$features=neurons_metrics$reduction %>% stringr::str_extract(pattern="\\.\\..*\\.features") %>% stringr::str_replace(pattern = "\\.\\..*\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.features",replacement = "")
neurons_metrics$method = stringr::str_extract(neurons_metrics$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
neurons_metrics$features_ngenes = stringr::str_extract(neurons_metrics$reduction,pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
# filter
top_from_each_method=neurons_metrics %>% dplyr::filter(ndim %in% c("30","50","60","80","90",NA)) %>% dplyr::mutate(rank_sum = purity_score+mixing_score,rank_diff = abs(purity_score-mixing_score)) %>%
  dplyr::mutate(rank_col =rank_sum-rank_diff) %>% dplyr::arrange(desc(rank_col)) %>% dplyr::group_by(method) %>% dplyr::filter(silhouette_score_norm>50 | silhouette_score_norm ==max(silhouette_score_norm)) %>% dplyr::top_n(n = 1)
names_to_check = top_from_each_method$reduction
names(names_to_check) = paste0("best_",top_from_each_method$method)

## add 2 more examples:
names_to_check = c(names_to_check,"scvi_cov"="scVI_1_300_0.025_3_256_gene_zinb_cov2..scVI..90..features_RNA.log.vst.split_Batch_ID.features.750_73e0c7a2fe2b9ab7fe2f1b1afa2980d6")
names_to_check = c(names_to_check,"harmony_lowdim"="harmony_8_Batch_ID_2_0.1_130_PCA..RNA..30..RNA.log.vst.split_Batch_ID.features.750_9625657cf22d373e46bdde89031e6f00")

##### Comparison line plots

# rf mxing
mixing_prob_perCelltype = per_celltype_summary_selected(evaluation_mixing_probNorm_all, mapped_celltypes, names_to_check)
mixing_prob_perCelltype_lineplot = ggplot(mixing_prob_perCelltype,aes(x=celltype,y=median_score,group=reduction,color=reduction))+geom_line(size=0.7)+geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 0.5))+
  xlab("Celltype")+ylab("Median score")+ggtitle("Mixing score (RF) on selected celltypes")
mixing_prob_perCelltype_lineplot
#theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 0.5))+ggtitle("Mean mixing vs mean purity")+xlim(0,100)+ylim(0,100)#+ theme_bw()

# knn entropy mixing
mixing_entropy_perCelltype = per_celltype_summary_selected(evaluation_entropy_knn_all, mapped_celltypes, names_to_check)
mixing_entropy_perCelltype_lineplot = ggplot(mixing_entropy_perCelltype,aes(x=celltype,y=median_score,group=reduction,color=reduction))+geom_line(size=0.7)+geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 0.5))+
  xlab("Celltype")+ylab("Median score")+ggtitle("Mixing score (knn) on selected celltypes")
mixing_entropy_perCelltype_lineplot

# knn purity:
purity_perCelltype = evaluation_purity_knn_all[,names_to_check]
colnames(purity_perCelltype) =names(names_to_check)
purity_perCelltype$celltype = rownames(purity_perCelltype)
purity_perCelltype = purity_perCelltype %>% tidyr::gather(-celltype,key="reduction",value="median_score")
purity_perCelltype_lineplot = ggplot(purity_perCelltype,aes(x=celltype,y=median_score,group=reduction,color=reduction))+geom_line(size=0.7)+geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 0.5))+
  xlab("Celltype")+ylab("Median score")+ggtitle("Purity score (knn) on selected celltypes")
purity_perCelltype_lineplot

### save results
ggsave(filename = paste0(results_path,"mixing_rf_perCelltype_methods.png"),
       plot = mixing_prob_perCelltype_lineplot, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path,"mixing_rf_perCelltype_methods.pdf"),
       plot = mixing_prob_perCelltype_lineplot, "pdf",dpi=400,width=300,height =200,units="mm")

ggsave(filename = paste0(results_path,"mixing_entropy_perCelltype_methods.png"),
       plot = mixing_entropy_perCelltype_lineplot, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path,"mixing_entropy_perCelltype_methods.pdf"),
       plot = mixing_entropy_perCelltype_lineplot, "pdf",dpi=400,width=300,height =200,units="mm")

ggsave(filename = paste0(results_path,"purity_perCelltype_methods.png"),
       plot = purity_perCelltype_lineplot, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path,"purity_perCelltype_methods.pdf"),
       plot = purity_perCelltype_lineplot, "pdf",dpi=400,width=300,height =200,units="mm")


##########
### Supplemental Figure 4: UMAPs of best methods
##########

## make plots
rasterize_point_size = 1.5
rasterize_pixels = 1536

## CAREFUL THIS RUNS SOME TIME TO CALC ALL UMAPS!

# add reductions:
integration_res_path = paste0(large_data_path,"best_reductions_per_method/")
for( i in 1:length(names_to_check)){
  neuron_map_seurat = add_reduction_seurat(neuron_map_seurat,integration_name=names_to_check[i],new_name=names(names_to_check)[i],
                                           integration_res_path,max_dim=200,calc_umap=TRUE,k_param_umap=30,overwrite =F,overwrite2=F)
}

# umap by Dataset or Batch_ID
p <-list()
for(i in 1:length(names_to_check)) {
  p[[i]] <- DimPlot(neuron_map_seurat,group.by = "Dataset",reduction = paste0("umap_",names(names_to_check)[i]),label=F)+ NoAxes() + NoLegend()+ggtitle(names(names_to_check)[i])
  p[[i]] <- rasterize_ggplot( p[[i]],pixel_raster = rasterize_pixels, pointsize =  rasterize_point_size)
}
compareMthods_cp_batch = cowplot::plot_grid(plotlist = p,ncol=4)
compareMthods_cp_batch

#save
ggsave(filename = paste0(results_path,"umap_dataset_methods.png"),
       plot = compareMthods_cp_batch, "png",dpi=400,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path,"umap_dataset_methods.pdf"),
       plot = compareMthods_cp_batch, "pdf",dpi=400,width=300,height =200,units="mm")





