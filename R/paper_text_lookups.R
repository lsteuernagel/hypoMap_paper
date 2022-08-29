## get cluster annotation overview:
anno_df = hypoMap_v2_seurat@misc$annotation_result %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})
anno_df$pasted = paste0(anno_df$cluster_id,": ",anno_df$cluster_name)
edgelist = hypoMap_v2_seurat@misc$clustering_edgelist

# helper
stringr::str_match(scUtils::find_children("C66-4",edgelist),"C185-[0-9]+") %>% na.omit() %>% as.character()

## get percentage per major cluster
table(hypoMap_v2_seurat@meta.data$C2_named) / sum(table(hypoMap_v2_seurat@meta.data$C2_named))
table(hypoMap_v2_seurat@meta.data$C7_named) / sum(table(hypoMap_v2_seurat@meta.data$C7_named))

table(hypoMap_v2_seurat@meta.data$Author_Class_Curated) / sum(table(hypoMap_v2_seurat@meta.data$Author_Class_Curated))

## how many neronal clusters
length(stringr::str_match(scUtils::find_children("C2-1",edgelist),"C66-[0-9]+") %>% na.omit() %>% as.character())
length(stringr::str_match(scUtils::find_children("C2-1",edgelist),"C25-[0-9]+") %>% na.omit() %>% as.character())
length(stringr::str_match(scUtils::find_children("C2-1",edgelist),"C185-[0-9]+") %>% na.omit() %>% as.character())

# gnrh1
sort(table(hypoMap_v2_seurat@meta.data$C66_named))

# regiona anno
region_stat = hypoMap_v2_seurat@meta.data %>% dplyr::group_by(C286,C286_named,Region_summarized) %>% dplyr::count()

# non - neurons
table(hypoMap_v2_seurat@meta.data$C7_named[hypoMap_v2_seurat@meta.data$C2 == "C2-2"]) / sum(table(hypoMap_v2_seurat@meta.data$C7_named[hypoMap_v2_seurat@meta.data$C2 == "C2-2"]))
length(stringr::str_match(scUtils::find_children("C7-5",edgelist),"C185-[0-9]+") %>% na.omit() %>% as.character())
length(stringr::str_match(scUtils::find_children("C66-52",edgelist),"C185-[0-9]+") %>% na.omit() %>% as.character())
stringr::str_match(scUtils::find_children("C25-18",edgelist),"C66-[0-9]+") %>% na.omit() %>% as.character()

### nucseq distribution
length(table(dowsett_subset@meta.data$C185)[table(dowsett_subset@meta.data$C185) >= 30])
under_rep_clusters = table(dowsett_subset@meta.data$C286_named)[table(dowsett_subset@meta.data$C286_named)< 30]
table(region_stat$Region_summarized[region_stat$C286_named %in% names(under_rep_clusters)])

##
cor_summary_per_class = data.table::fread(paste0(results_path_figure3,"summary_per_class.txt"),data.table = F)
cor_summary_per_class

# cor cluster
per_dataset_cor_df_print = data.table::fread(paste0(results_path_figure3,"per_cluster_correlations.txt"),data.table = F)
per_dataset_cor_df_print[per_dataset_cor_df_print$C185 %in% c("C185-11","C185-12"),]

# fasting:
# ---> see figure script itself

# romanov:
romanovpredicted_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/mapped_data/romanov_predicted_seurat.rds")

length(table(romanovpredicted_seurat@meta.data$C185_named_predicted)[table(romanovpredicted_seurat@meta.data$C185_named_predicted) >= 3])
sort(table(romanovpredicted_seurat@meta.data$C25_named_predicted) / sum(table(romanovpredicted_seurat@meta.data$C25_named_predicted)))
sort(table(romanovpredicted_seurat@meta.data$C185_named_predicted))# / sum(table(romanovpredicted_seurat@meta.data$C185_named_predicted)))
table(romanovpredicted_seurat@meta.data$C185_named_predicted[romanovpredicted_seurat@meta.data$Author_CellType == "Neurons_GABA 13 (Galanin)"])
sort(table(romanovpredicted_seurat@meta.data$C185_named_predicted[romanovpredicted_seurat@meta.data$Author_CellType == "Neurons_GABA 7 (Pomc+/-)"]))
cbind(FetchData(romanovpredicted_seurat,vars = "Pomc",cells = romanovpredicted_seurat@meta.data$Cell_ID[romanovpredicted_seurat@meta.data$Author_CellType == "Neurons_GABA 7 (Pomc+/-)"] ), romanovpredicted_seurat@meta.data$C185_named_predicted[romanovpredicted_seurat@meta.data$Author_CellType == "Neurons_GABA 7 (Pomc+/-)"])


## bacTRAP
rbo_result_bacTRAP=data.table::fread(paste0(results_path_figure5,"bacTRAP_rbo_enrichment.txt"),data.table = F)
a1=rbo_result_bacTRAP[rbo_result_bacTRAP$signature_name == "pnoc",] %>% dplyr::arrange(desc(rbo)) %>% dplyr::left_join(anno_df,by=c("cluster"="cluster_id"))
rbo_result_bacTRAP[rbo_result_bacTRAP$signature_name == "pnoc",] %>% dplyr::arrange(desc(rbo)) %>% dplyr::left_join(anno_df,by=c("cluster"="cluster_id")) %>% head(n=20)

pnoc_pcts = scUtils::gene_pct_cluster(seurat_object = hypoMap_v2_seurat,col_name = "C185_named",genes = "Pnoc")
pnoc_pcts$cluster_name = rownames(pnoc_pcts)
pnoc_pcts[pnoc_pcts$cluster_name %in% c("C185-117: Npy.Sst.GABA-4","C185-118: Otp.Sst.GABA-4"),]
# Otp.Sst.GABA-4
pnoc_pcts = scUtils::gene_pct_cluster(seurat_object = hypoMap_v2_seurat,col_name = "C286_named",genes = "Pnoc")
pnoc_pcts$cluster_name = rownames(pnoc_pcts)
pnoc_pcts[pnoc_pcts$cluster_name %in% c("C286-158: Vgll3.Tbx3.GABA-1","C286-159: Crabp1.Sytl4.Tbx3.GABA-1"),]


## mousebrain discussion
leaf_level_column = "C185"
leaf_level = 6
# make data for first heatmap with percentages per dataset
heatmap_data = hypoMap_v2_seurat@meta.data %>% dplyr::select(Cell_ID,Dataset,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Dataset) %>% #dplyr::filter(predicted_Campbell!="NA") 
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label")) 
  tidyr::spread(key = Dataset,value=presence) %>% as.data.frame()

length(heatmap_data[,"Mousebrainorg10x"][heatmap_data[,"Mousebrainorg10x"] >= 1])

(185 - length(heatmap_data[,"Mousebrainorg10x"][is.na(heatmap_data[,"Mousebrainorg10x"])])) / 185
length(heatmap_data[,"Mousebrainorg10x"][heatmap_data[,"Mousebrainorg10x"] >= 2]) / 179


## astrocytes: C66-54: Lgals3.Astrocytes
lgals_astros = subset(hypoMap_v2_seurat,subset = C66_named == "C66-54: Lgals3.Astrocytes")

table(lgals_astros@meta.data$Diet)



