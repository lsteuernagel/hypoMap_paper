##########
### Load & Prepare
##########

# I don't want to have excel screw up numbers in scientifc notation:
#options(scipen = 999)
fwrite_scipen = 999 # I'll leave it one for some sc-seq tables because they are difficult to read else

results_path_tables = "table_outputs/"
system(paste0("mkdir -p ",results_path_tables))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
#require(openxlsx)
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
### Overview of required Tables
##########

# - [ ]  dataset overview
# - [x]  neuron signatures
# - [x]  integration results across methods (neurons)
# - [x]  integration results scVI
# - [x]  integration results scVI full vs neurons
# - [ ]  hypomap cluster overview
# - include dataset contribution pcts
# - [ ]  hypomap cluster region prediction
# - [ ]  cluster tree edgelist
# - [ ]  top marker genes (add deseq pval)
# - top marker genes (siblings) ? (add deseq pval)
# - [ ]  correlations per gene
# - [ ]  correlations per cluster
# - [ ]  degs IEGs (add deseq pval)
# - [ ]  degs fasting (add deseq pval)
# - [ ]  nucseq gobp Agrp fasting
# - [ ]  bacTRAP signature AgRP
# - [ ]  bacTRAP signature Pomc
# - [ ]  bacTRAP signature Pomc-lepr
# - [ ]  bacTRAP signature Pomc-Glp1r
# - [ ]  bacTRAP signature Glp1r
# - [ ]  bacTRAP signature Pnoc
# - [ ]  bacTRAp enrichment results (all in one)
# - [ ]  rnascope glp1r
# - [ ]  rnascope Pnoc

##########
### Existing Supplementary tables 1 + 2
##########

## Table with dataset overview
# load table 1
# dataset_overview = data.table::fread(paste0(results_path_tables,"supplementary_table_1_dataset_overview.csv"))
dataset_overview = data.table::fread(paste0("data_inputs/dataset_overview_inputV2.csv"))
dataset_overview$Name[dataset_overview$Name == "Lam"] = "Dowsett"
dataset_overview = dataset_overview %>% dplyr::arrange((Name))
#  and add more metadata ?
# ncells_after_qc , n_batches
dataset_stats = hypoMap_v2_seurat@meta.data %>% dplyr::group_by(Dataset) %>% dplyr::add_count(name="ncells")   %>%
   dplyr::distinct(Dataset, ncells,Batch_ID) %>% dplyr::group_by(Dataset) %>% dplyr::add_count(name="nbatches") %>%
  dplyr::distinct(Dataset, ncells,nbatches)

# join
dataset_overview = dplyr::bind_cols(dataset_overview, dataset_stats)#, by=c("Name"="Dataset"))

# finalize =
dataset_overview = dataset_overview %>% dplyr::select(Dataset,SRA_Acc = `SRA Acc`,GEO_Acc = `GEO Acc`,publication = paper,year,Process = Count,Technology ,Region=Subregion,ncells,nbatches)
dataset_overview$used_for_integration_method_eval = "n"
dataset_overview$used_for_integration_method_eval[dataset_overview$Dataset %in% c("CampbellDropseq","ChenDropseq" ,"Flynn10x" ,"Kim10x" ,"KimDev10x", "LeeDropseq" ,"Mickelsen10x" ,"Moffit10x" ,"Mousebrainorg10x","RomanovDev10x","RossiDropseq","Wen10x","wenDropseq")] = "y"

# save
data.table::fwrite(dataset_overview,paste0(results_path_tables,"supplementary_table_1_dataset_overview.csv"))

##########
### Existing Supplementary tables 2
##########

## Table with purity signature overview
#quick fix, don't run again
# neuron_signatures = data.table::fread(paste0("data_inputs/","supplementary_table_2_neuron_signatures.csv"))
# curated_neuron_celltype_signatures = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/documentation/curated_neuron_celltype_signatures.rds")
# curated_neuron_celltype_signatures_df = data.frame(name= names(curated_neuron_celltype_signatures),genes = sapply(curated_neuron_celltype_signatures,paste0,collapse = ";"))
# neuron_signatures = neuron_signatures %>% dplyr::select(-marker_genes,-n_cells) %>% dplyr::left_join(curated_neuron_celltype_signatures_df,by="name")
# data.table::fwrite(neuron_signatures,paste0("data_inputs/","supplementary_table_2_neuron_signatures.csv"))
# call_celltypes = sapply(jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scIntegration/data/hypothalamus_celltype_signatures.json"),unlist)
# missing_celltypes = call_celltypes[!names(call_celltypes) %in% neuron_signatures$name]
# missing_celltypes_df = data.frame(name= names(missing_celltypes),genes = sapply(missing_celltypes,paste0,collapse = ";"),signature_length= sapply(missing_celltypes,length))
# all_signatures = dplyr::bind_rows(neuron_signatures,missing_celltypes_df)
# data.table::fwrite(all_signatures,paste0("data_inputs/","supplementary_table_2_all_signatures.csv"))

all_signatures = data.table::fread(paste0("data_inputs/","supplementary_table_2_all_signatures.csv"))
data.table::fwrite(all_signatures,paste0(results_path_tables,"supplementary_table_2_all_signatures.csv"))

##########
### Metrics methods comparison
##########

# REMOVE --> SOURCE DATA

# # for Methods eval:
# results_path_supplementary_figure1 = "figure_outputs/figure_supplementary_1/"
# neurons_metrics_curated = data.table::fread(paste0(results_path_supplementary_figure1,"neurons_metrics_curated.csv"),data.table = FALSE)
# neurons_metrics_curated = neurons_metrics_curated %>% dplyr::select(assay,method,mixing_score,purity_score,reduction,ndim,features_ngenes)
# # save:
# data.table::fwrite(neurons_metrics_curated,paste0(results_path_tables,"supplementary_table_3_integration_methods_metrics_overview.csv"))#,scipen = fwrite_scipen)

##########
### Metrics scvi specifically
##########

# REMOVE --> SOURCE DATA

# results_path_supplementary_figure2 = "figure_outputs/figure_supplementary_2/"
# complete_evaluation_result = data.table::fread("data_inputs/complete_evaluation_results_updated.txt",data.table = FALSE)
# complete_evaluation_result = complete_evaluation_result[,c(2,1,3:ncol(complete_evaluation_result))] %>% dplyr::select(-scvi_params) %>% as.data.frame()
# # save:
# data.table::fwrite(complete_evaluation_result,paste0(results_path_tables,"supplementary_table_4_scvi_metrics_overview.csv"))#,scipen = fwrite_scipen)
# 
# full_vs_neurons_evaluation_result = data.table::fread("data_inputs/full_vs_neurons_evaluation_results.txt",data.table = F)
# full_vs_neurons_evaluation_result = full_vs_neurons_evaluation_result %>% dplyr::select(-scvi_params) %>% as.data.frame()
# # save:
# data.table::fwrite(full_vs_neurons_evaluation_result,paste0(results_path_tables,"supplementary_table_5_scvi_neurons_metrics_overview.csv"))#,scipen = fwrite_scipen)

##########
### HypoMap Cluster overview tables
##########

## Table with all clusters in HypoMap

# regiona anno
region_stat = hypoMap_v2_seurat@meta.data %>% dplyr::group_by(C286,Region_summarized) %>% dplyr::count() %>% dplyr::select(-n)
hypo_map_cluster_overview = dplyr::left_join(hypoMap_v2_seurat@misc$annotation_result,region_stat,by=c("cluster_id"="C286")) %>% 
  dplyr::group_by(clusterlevel) %>% dplyr::mutate(pct_of_all = round(ncells / sum(ncells), 4) * 100) %>% dplyr::ungroup() #%>% dplyr::mutate(name_pasted = paste0())
hypo_map_cluster_overview = hypo_map_cluster_overview %>% dplyr::select(cluster_id,cluster_name = clean_names_withID,clusterlevel,Region_summarized,ncells,pct_of_all)

# save:
data.table::fwrite(hypo_map_cluster_overview,paste0(results_path_tables,"supplementary_table_3_hypomap_cluster_overview.csv"))

## Table with mrtree edgelist for TREE
edgelist_clusters_hypoMap = hypoMap_v2_seurat@misc$clustering_edgelist %>% dplyr::select(-isLeaf,-level,-totalCount,nchildren = count)
# save:
data.table::fwrite(edgelist_clusters_hypoMap,paste0(results_path_tables,"supplementary_table_4_hypomap_cluster_edgelist.csv"))#,scipen = fwrite_scipen)


##########
### top marker genes (add deseq pval)
##########

## table with hypomap cluster marker genes
hypomap_marker_genes = hypoMap_v2_seurat@misc$marker_genes_all

# maybe restrict to high specificity to cut down size ?
hypomap_marker_genes = hypomap_marker_genes %>% dplyr::filter(specificity > 1 & p_val_adj < 0.00001) 
# restrict to top 50 per clusters
old_order = unique(hypomap_marker_genes$cluster_id)
hypomap_marker_genes2 = hypomap_marker_genes %>% dplyr::group_by(cluster_id) %>% dplyr::top_n(n = 50,wt = specificity*(-log10(p_val_adj+1e-200))) %>% dplyr::arrange(desc(specificity*(-log10(p_val_adj+1e-200)))) #%>% dplyr::arr
# reorder
hypomap_marker_genes2$cluster_id = factor(hypomap_marker_genes2$cluster_id,levels = hypoMap_v2_seurat@misc$annotation_result$cluster_id)
hypomap_marker_genes2 = hypomap_marker_genes2 %>% dplyr::arrange(cluster_id)
hypomap_marker_genes2$cluster_id = as.character(hypomap_marker_genes2$cluster_id)
# save:
data.table::fwrite(hypomap_marker_genes2,paste0(results_path_tables,"supplementary_table_5_hypomap_marker_genes.csv"))


##########
### top marker genes (siblings) ? (add deseq pval)
##########


## table with hypomap cluster marker genes
hypomap_sibling_marker_genes = hypoMap_v2_seurat@misc$marker_genes_siblings

# maybe restrict to high specificity to cut down size ?
hypomap_sibling_marker_genes = hypomap_sibling_marker_genes %>% dplyr::filter(specificity > 1 & p_val_adj < 0.00001) 
# restrict to top 50 per clusters
old_order = unique(hypomap_sibling_marker_genes$cluster_id)
hypomap_sibling_marker_genes2 = hypomap_sibling_marker_genes %>% dplyr::group_by(cluster_id) %>% dplyr::top_n(n = 50,wt = specificity*(-log10(p_val_adj+1e-200))) %>% dplyr::arrange(desc(specificity*(-log10(p_val_adj+1e-200)))) #%>% dplyr::arr
# reorder
hypomap_sibling_marker_genes2$cluster_id = factor(hypomap_sibling_marker_genes2$cluster_id,levels = hypoMap_v2_seurat@misc$annotation_result$cluster_id)
hypomap_sibling_marker_genes2 = hypomap_sibling_marker_genes2 %>% dplyr::arrange(cluster_id)
hypomap_sibling_marker_genes2$cluster_id = as.character(hypomap_sibling_marker_genes2$cluster_id)
# save:
data.table::fwrite(hypomap_sibling_marker_genes,paste0(results_path_tables,"supplementary_table_6_hypomap_sibling_marker_genes.csv"))

##########
### hypomap per_dataset_contribution
##########

# add per cluster pct
leaf_level_column = "C185_named"
per_dataset_contribution = hypoMap_v2_seurat@meta.data %>% dplyr::select(Cell_ID,Dataset,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Dataset) %>% #dplyr::filter(predicted_Campbell!="NA") 
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label")) 
  tidyr::spread(key = Dataset,value=presence) %>% as.data.frame()
# save:
data.table::fwrite(per_dataset_contribution,paste0(results_path_tables,"supplementary_table_7_hypomap_dataset_contribution_C185.csv"))#,scipen = fwrite_scipen)

##########
### hypomap cluster region prediction
##########

scores_per_target_level_region_all_short = data.table::fread(paste0("data_inputs/","aba_enrichment_scores_top10_C286.txt"),data.table = F)
# save:
data.table::fwrite(scores_per_target_level_region_all_short,paste0(results_path_tables,"supplementary_table_8_hypomap_region_prediction_C286.csv"))#,scipen = fwrite_scipen)

##########
### correlations per gene
##########

results_path_figure3 = "figure_outputs/figure_4/"
per_gene_cor = data.table::fread(paste0(results_path_figure4,"per_gene_correlations_with_class.txt"),data.table = FALSE)
data.table::fwrite(per_gene_cor,paste0(results_path_tables,"supplementary_table_9_hypomap_nucSeq_per_gene_correlation.csv"))


##########
### correlations per cluster
##########

# REMOVE --> SOURCE DATA

# per_dataset_cor_df = data.table::fread(paste0(results_path_figure3,"per_cluster_correlations.txt"),data.table = FALSE) %>% as.data.frame()
# data.table::fwrite(per_dataset_cor_df,paste0(results_path_tables,"supplementary_table_13_hypomap_nucSeq_per_cluster_correlation.csv"))

##########
### fasting IEGS  
##########

results_path_figure5 = "figure_outputs/figure_5/"

## table witH IEGs
activation_per_cluster_p = data.table::fread(paste0(results_path_figure5,"activation_genes_per_cluster.txt"),data.table = FALSE)
activation_per_cluster_p = activation_per_cluster_p %>% dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_name,cluster_id),by=c("current_cluster"="cluster_id")) %>% 
  dplyr::filter(p_val_bonferroni < 0.05 & abs(avg_log2FC) > 0.1) %>%
  dplyr::select(cluster_id = current_cluster,cluster_name = cluster_name, gene, avg_log2FC, p_val_bonferroni,p_val_bonferroni_negbinom, avg_log2FC, pct.1, pct.2,pct_diff) 
# save:
data.table::fwrite(activation_per_cluster_p,paste0(results_path_tables,"supplementary_table_10_nucseq_IEGs_fasting.csv"))

##########
### fasting DEGs 
##########

# REMOVE --> SOURCE DATA --> I KEEP THIS because it is useful Supplemental information

## table witH DEGs
all_clusters_fasting_DEG = data.table::fread("figure_outputs/figure_extended_7/all_clusters_fasting_DEG.txt",data.table = F)
all_clusters_fasting_DEG = all_clusters_fasting_DEG %>% dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_name,cluster_id),by=c("current_cluster"="cluster_id"))
all_clusters_fasting_DEG = all_clusters_fasting_DEG %>% dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::select(cluster_id = current_cluster,cluster_name = cluster_name, gene, avg_log2FC, p_val_adj, p_val_adj_negbinom, avg_log2FC, pct.1, pct.2,pct_diff) %>% 
  dplyr::group_by(gene) %>% dplyr::add_count(name = "N_DEG")
# save:
data.table::fwrite(all_clusters_fasting_DEG,paste0(results_path_tables,"supplementary_table_11_nucseq_DEG_fasting.csv"))

##########
### GO enrichment Agrp
##########

# REMOVE --> SOURCE DATA
# 
# ## table with go enrichment results for Agrp ?
# agrp_fasting_go_enrichment_simplified = data.table::fread(paste0(results_path_figure4,"/agrp_fasting_go_enrichment_simplified.txt"),data.table = F)
# #agrp_fasting_go_enrichment_simplified = agrp_fasting_go_enrichment_simplified %>% dplyr::select(ID)
# data.table::fwrite(agrp_fasting_go_enrichment_simplified,paste0(results_path_tables,"supplementary_table_16_nucseq_gobp_Agrp_fasting.csv"))


##########
### ampbell vs nucseq
##########


# REMOVE --> SOURCE DATA

# agrp_sn_vs_campbell_DEG = data.table::fread( paste0(results_path_figure4,"agrp_sn_vs_campbell_DEG.txt"),data.table = F)
# agrp_sn_vs_campbell_DEG = agrp_sn_vs_campbell_DEG %>% dplyr::select(gene,avg_log2FC_sn,avg_log2FC_sc,p_val_adj_sn,p_val_adj_negbinom_sn,p_val_adj_sc,p_val_adj_negbinom_sc,regulated)
# # save:
# data.table::fwrite(agrp_sn_vs_campbell_DEG,paste0(results_path_tables,"supplementary_table_17_nucseq_vs_campbell_DEG_agrp_fasting.csv"))       


##########
### GO enrichment overlap genes
##########


# REMOVE --> SOURCE DATA

# ## table with go enrichment 
# results_path_supplementary_figure4 = "figure_outputs/figure_supplementary_4/"
# overlap_genes_go_enrichment_simplified = data.table::fread(paste0(results_path_supplementary_figure4,"overlap_genes_go_enrichment_simplified.txt"),data.table = F)
# #agrp_fasting_go_enrichment_simplified = agrp_fasting_go_enrichment_simplified %>% dplyr::select(ID)
# data.table::fwrite(overlap_genes_go_enrichment_simplified,paste0(results_path_tables,"supplementary_table_18_nucseq_gobp_overlap_genes.csv"))

##########
### Romanov overview projection
##########

results_path_figure6 = "figure_outputs/figure_6/"

romanov_cells_projected = data.table::fread(paste0(results_path_figure6,"romanov_cells_projected.txt"),data.table = F)
data.table::fwrite(romanov_cells_projected,paste0(results_path_tables,"supplementary_table_12_romanov_projected.csv"))

##########
### bacTRAp results
##########

# This is basically a prettified version of the input DEseq2 tables for Figure 3. 

# run on all files:
start_idx_for_supplemental_file = 13
bacTRAP_files = list.files("data_inputs/",pattern = "bacTRAP_deseq2")
for(i in 1:length(bacTRAP_files)){
  current_res = data.table::fread(paste0("data_inputs/",bacTRAP_files[i]),data.table = FALSE) 
  print(colnames(current_res))
  colnames(current_res)[colnames(current_res)=="row"] = "ensembl_gene_id"
  current_res =  current_res %>%dplyr::select(gene = external_gene_name,ensembl_gene_id,log2FoldChange,padj)
  current_res = current_res[!is.na(current_res$log2FoldChange),]
  current_res$padj[is.na(current_res$padj)] = 1
  current_res$log2FoldChange = round(current_res$log2FoldChange,5)
  #current_res$padj = round(current_res$padj,7)
  current_name = stringr::str_to_title(gsub("bacTRAP_deseq2_|\\.csv","",bacTRAP_files[i]))
  assign(x = paste0("result_bacTRAP_",current_name),value = current_res)
  message(paste0(results_path_tables,"supplementary_table_",(start_idx_for_supplemental_file+i-1),"_bacTRAP_",current_name,".csv"))
  data.table::fwrite(current_res,paste0(results_path_tables,"supplementary_table_",(start_idx_for_supplemental_file+i-1),"_bacTRAP_",current_name,".csv"))#,scipen = fwrite_scipen)
}
print(start_idx_for_supplemental_file+i)

##########
### nucseq Cluster overview tables
##########

# ## Table with all projected HypoMap clusters in nucSeq
# # name, id, n_cells projected, region, N-deg fasting, correlation with HypoMap
# cell_cluster_map =query_snseq_neurons@meta.data %>% dplyr::select("Cell_ID",tidyselect::matches("K[0-9]+_pruned")) %>% tidyr::gather(-Cell_ID,key="clusterlevel",value="cluster_id") # group all levels ing long format
# cell_cluster_map$clusterlevel = gsub("predicted_|_pruned","",cell_cluster_map$clusterlevel)
# nucseq_cluster_overview = cell_cluster_map %>% dplyr::group_by(cluster_id) %>% dplyr::add_count(name="ncells_nuc") %>% dplyr::distinct(cluster_id,clusterlevel,ncells_nuc) %>% # count cells
#   dplyr::group_by(clusterlevel) %>% dplyr::mutate(pct_of_all_nuc = round(ncells_nuc / sum(ncells_nuc), 4) * 100) %>% dplyr::ungroup() %>% # calc pct
#   dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_id,cluster_name,likely_region,pct_of_all_ref = pct_of_all),by="cluster_id") %>% # add hypomap infos
#   dplyr::filter(!is.na(cluster_id)) %>% dplyr::mutate(cluster_id = factor(cluster_id,levels = hypo_map_cluster_overview$cluster_id[hypo_map_cluster_overview$cluster_id %in% cell_cluster_map$cluster_id])) %>%
#   dplyr::arrange(cluster_id) %>% # sort similar to hypoMap overview
#   dplyr::select(projected_cluster_id = cluster_id, projected_cluster_name = cluster_name, clusterlevel, ncells_nuc, pct_of_all_nuc, pct_of_all_ref, likely_region) # reorder and name columns 
# # add correlations
# per_cluster_correlations = data.table::fread("figure_outputs/figure_4/per_cluster_correlations.txt",data.table = FALSE)
# per_cluster_correlations = per_cluster_correlations[,2:ncol(per_cluster_correlations)]
# colnames(per_cluster_correlations)[3:length(colnames(per_cluster_correlations))] = paste0(colnames(per_cluster_correlations)[3:length(colnames(per_cluster_correlations))],"_cor")
# nucseq_cluster_overview = nucseq_cluster_overview %>% dplyr::left_join(per_cluster_correlations,by=c("projected_cluster_id"="K169_pruned"))
# # add N degs in fasting
# all_clusters_fasting_DEG = data.table::fread("figure_outputs/figure_supplementary_7/all_clusters_fasting_DEG.txt")
# all_clusters_fasting_DEG_stat = all_clusters_fasting_DEG %>% dplyr::filter(p_val_adj < 0.01) %>% dplyr::group_by(current_cluster) %>% dplyr::count(name = "n_degs_fasting")
# nucseq_cluster_overview = nucseq_cluster_overview %>% dplyr::left_join(all_clusters_fasting_DEG_stat,by=c("projected_cluster_name"="current_cluster"))
# # save:
# data.table::fwrite(nucseq_cluster_overview,paste0(results_path_tables,"supplementary_table_13_nucseq_cluster_overview.csv"))

##########
### bacTRAp result
##########

# REMOVE --> SOURCE DATA --> I still keep it!

rbo_result_bacTRAP = data.table::fread(paste0(results_path_figure6,"bacTRAP_rbo_enrichment.txt"),data.table = F)
rbo_result_bacTRAP = rbo_result_bacTRAP %>% dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_name,cluster_id),by=c("cluster"="cluster_id")) %>%
  dplyr::select(cluster_id =cluster,cluster_name,rbo_score = rbo,bacTRAP_signature = signature_name) %>% dplyr::arrange(desc(rbo_score)) %>% dplyr::arrange(bacTRAP_signature)

data.table::fwrite(rbo_result_bacTRAP,paste0(results_path_tables,"supplementary_table_19_bacTRAP_enrichment.csv"))


##########
### ISH results Glp1r
##########


# REMOVE --> SOURCE DATA

# This is the input table for Figure 6 which is based on the manual counting of images.
# pct_expressed_cells_clusters_ISH = data.table::fread("figure_outputs/figure_6/pct_expressed_cells_glp1r_ISH.txt",data.table = F)
# data.table::fwrite(pct_expressed_cells_clusters_ISH,paste0(results_path_tables,"supplementary_table_27_rnascope_glp1r_result.csv"))#,scipen = fwrite_scipen)

##########
### ISH results Pnoc
##########

# REMOVE --> SOURCE DATA

# pct_expressed_cells_pnoc_sst_ISH = data.table::fread("figure_outputs/figure_6/pct_expressed_cells_pnoc_sst_ISH.txt",data.table = F)
# data.table::fwrite(pct_expressed_cells_pnoc_sst_ISH,paste0(results_path_tables,"supplementary_table_28_rnascope_pnoc_sst_result.csv"))#,scipen = fwrite_scipen)
# 
# 
# # This is the input table for Figure 6 which is based on the manual counting of images.
# # pct_expressed_cells_clusters_ISH = data.table::fread("figure_outputs/figure_6/pct_expressed_cells_clusters_ISH.txt")
# pct_expressed_cells_pnoc_crabp1_ISH = data.table::fread("figure_outputs/figure_6/pct_expressed_cells_pnoc_crabp1_ISH.txt",data.table = F)
# data.table::fwrite(pct_expressed_cells_pnoc_crabp1_ISH,paste0(results_path_tables,"supplementary_table_29_rnascope_pnoc_crabp1_result.csv"))#,scipen = fwrite_scipen)

##########
### metadata
##########

cell_labels = hypoMap_v2_seurat@meta.data %>% dplyr::select(Cell_ID, Dataset, Author_CellType,Class_Harmonized = Author_Class_Curated,C2,C7,C25,C66,C185,C286,C465)
#cell_labels = hypoMap_v2_seurat@meta.data %>% dplyr::select(Cell_ID, Dataset, Author_CellType,Class_Harmonized = Author_Class_Curated,C2_named,C7_named,C25_named,C66_named,C185_named,C286_named,C465_named)
data.table::fwrite(cell_labels,paste0(results_path_tables,"supplementary_table_20_all_cell_labels.csv"))#,scipen = fwrite_scipen)


##########
### methods table with reagents etc
##########

methods_table_reagents = data.table::fread("data_inputs/methods_table_reagents.txt",data.table = F)
data.table::fwrite(methods_table_reagents,paste0(results_path_tables,"supplementary_table_21_reagents.csv"))#,scipen = fwrite_scipen)


##########
### join in xlsx
##########

supplementary_files = list.files("table_outputs/",pattern = ".csv")
supp_table_list = list()
dictionary_list = list()
supp_table_list[["0_dictionary"]] = data.frame()# dummy to fill in later
for(i in 1:length(supplementary_files)){
  current_table = data.table::fread(paste0("table_outputs/",supplementary_files[i]),data.table = FALSE) 
  current_name = gsub("supplementary_table_|\\.csv","",supplementary_files[i])
  current_colnames = colnames(current_table)
  df_for_dictionary = data.frame(table = c(current_name,current_colnames), explanation = NA)
  df_for_dictionary = rbind(df_for_dictionary,data.frame(table =NA,explanation =NA))
  dictionary_list[[current_name]] = df_for_dictionary
  supp_table_list[[current_name]] = current_table
}
# make dictionary 
dictionary_list = dictionary_list[base::order(as.numeric(stringr::str_extract(names(dictionary_list),pattern = "[0-9]+")))]
dictionary_df = do.call("rbind",dictionary_list)
supp_table_list[["0_dictionary"]] = dictionary_df
# order numerically
supp_table_list = supp_table_list[base::order(as.numeric(stringr::str_extract(names(supp_table_list),pattern = "[0-9]+")))]
names(supp_table_list)
# write to one file
#openxlsx::write.xlsx(x = supp_table_list,file = "table_outputs/supplementary_tables_preliminary.xlsx",colNames=TRUE,rowNames=FALSE,keepNA=FALSE)

# shorten colnames
names(supp_table_list)[sapply(names(supp_table_list),nchar) > 31] = sapply(names(supp_table_list)[sapply(names(supp_table_list),nchar) > 31],substr,start=1,stop=31)
#write
WriteXLS::WriteXLS(x = supp_table_list,ExcelFileName = "table_outputs/supplementary_tables_preliminary.xlsx",col.names=TRUE)

###### I manually created the file "supplementary_tables.xlsx" that includes the full column descriptions (not only skeleton) !!







