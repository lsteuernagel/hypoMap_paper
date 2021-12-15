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
require(openxlsx)
source("R/utility_functions.R")
source("R/plot_functions.R")

# where to find large data objects (need to change to local dir !)
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

# load seurat objects via large_data_path
load_required_files(large_data_path = large_data_path)

##########
### Overview of required Tables
##########

## Table with dataset overview
# load table 1 and add more metadata

## Table with purity signature overview

## Table with all clusters in HypoMap
# name, id, n_cells, region, ????

## Table with all projected HypoMap clusters in nucSeq
# name, id, n_cells projected, region, N-deg fasting, correlation with HypoMap

## Table with mrtree edgelist for TREE

## table with cluster marker genes
# standard columns ?
# maybe restrict to high specificty to cut down size ?

## table witH DEGs
# details ?

## table with go enrichment results for Agrp ?

## tables with DEseq2 bacTRAP results

## tables with ISH numbers


##########
### Existing Supplementary tables 1 + 2
##########

## Table with dataset overview
# load table 1
dataset_overview = data.table::fread(paste0(results_path_tables,"supplementary_table_1_dataset_overview.csv"))
#  and add more metadata ?

## Table with purity signature overview
neuron_signatures = data.table::fread(paste0(results_path_tables,"supplementary_table_2_neuron_signatures.csv"))

##########
### Metrics
##########

neurons_metrics_curated = data.table::fread(paste0(results_path_figure1,"neurons_metrics_curated.csv"),data.table = FALSE)
neurons_metrics_curated = neurons_metrics_curated %>% dplyr::select(assay,method,mixing_score,purity_score,reduction,ndim,features_ngenes)
# save:
data.table::fwrite(neurons_metrics_curated,paste0(results_path_tables,"supplementary_table_3_integration_metrics_neurons_overview.csv"),scipen = fwrite_scipen)

full_metrics_curated = data.table::fread(paste0(results_path_supplementary_figure2,"full_metrics_curated.csv"),data.table = FALSE)
full_metrics_curated = full_metrics_curated %>% dplyr::select(assay,method,mixing_score,purity_score,reduction,ndim,features_ngenes)
# save:
data.table::fwrite(full_metrics_curated,paste0(results_path_tables,"supplementary_table_4_integration_metrics_full_overview.csv"),scipen = fwrite_scipen)


##########
### HypoMap Cluster overview tables
##########

## Table with all clusters in HypoMap
# name, id, n_cells, region, ????
region_anno =  neuron_map_seurat@meta.data %>% dplyr::distinct(suggested_region_curated,K169_pruned) %>% dplyr::select(cluster_id = K169_pruned,likely_region = suggested_region_curated)
hypo_map_cluster_overview = dplyr::left_join(neuron_map_seurat@misc$annotations,region_anno,by="cluster_id") %>% 
  dplyr::group_by(clusterlevel) %>% dplyr::mutate(pct_of_all = round(ncells / sum(ncells), 4) * 100) %>% dplyr::ungroup()
# save:
data.table::fwrite(hypo_map_cluster_overview,paste0(results_path_tables,"supplementary_table_5_hypomap_cluster_overview.csv"))

## Table with mrtree edgelist for TREE
edgelist_clusters_hypoMap = neuron_map_seurat@misc$mrtree_edgelist %>% dplyr::select(-isLeaf,-level,-totalCount,nchildren = count)
# save:
data.table::fwrite(edgelist_clusters_hypoMap,paste0(results_path_tables,"supplementary_table_6_hypomap_cluster_edgelist.csv"),scipen = fwrite_scipen)

##########
### bacTRAp results
##########

# This is basically a prettified version of the input DEseq2 tables for Figure 3. 
# run on all files:
start_idx_for_supplemental_file = 7
bacTRAP_files = list.files("data_inputs/",pattern = "bacTRAP_deseq2")
for(i in 1:length(bacTRAP_files)){
  current_res = data.table::fread(paste0("data_inputs/",bacTRAP_files[i]),data.table = FALSE) 
  colnames(current_res)[colnames(current_res)=="row"] = "ensembl_gene_id"
  current_res =  current_res %>%dplyr::select(gene = external_gene_name,ensembl_gene_id,log2FoldChange,padj)
  current_res = current_res[!is.na(current_res$log2FoldChange),]
  current_res$padj[is.na(current_res$padj)] = 1
  current_res$log2FoldChange = round(current_res$log2FoldChange,5)
  current_res$padj = round(current_res$padj,7)
  current_name = stringr::str_to_title(gsub("bacTRAP_deseq2_|\\.csv","",bacTRAP_files[i]))
  assign(x = paste0("result_bacTRAP_",current_name),value = current_res)
  message(paste0(results_path_tables,"supplementary_table_",(start_idx_for_supplemental_file+i-1),"_bacTRAP_",current_name,".csv"))
  data.table::fwrite(current_res,paste0(results_path_tables,"supplementary_table_",(start_idx_for_supplemental_file+i-1),"_bacTRAP_",current_name,".csv"),scipen = fwrite_scipen)
}
print(start_idx_for_supplemental_file+i-1)

##########
### nucseq Cluster overview tables
##########

## Table with all projected HypoMap clusters in nucSeq
# name, id, n_cells projected, region, N-deg fasting, correlation with HypoMap
cell_cluster_map =query_snseq_neurons@meta.data %>% dplyr::select("Cell_ID",tidyselect::matches("K[0-9]+_pruned")) %>% tidyr::gather(-Cell_ID,key="clusterlevel",value="cluster_id") # group all levels ing long format
cell_cluster_map$clusterlevel = gsub("predicted_|_pruned","",cell_cluster_map$clusterlevel)
nucseq_cluster_overview = cell_cluster_map %>% dplyr::group_by(cluster_id) %>% dplyr::add_count(name="ncells_nuc") %>% dplyr::distinct(cluster_id,clusterlevel,ncells_nuc) %>% # count cells
  dplyr::group_by(clusterlevel) %>% dplyr::mutate(pct_of_all_nuc = round(ncells_nuc / sum(ncells_nuc), 4) * 100) %>% dplyr::ungroup() %>% # calc pct
  dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_id,cluster_name,likely_region,pct_of_all_ref = pct_of_all),by="cluster_id") %>% # add hypomap infos
  dplyr::filter(!is.na(cluster_id)) %>% dplyr::mutate(cluster_id = factor(cluster_id,levels = hypo_map_cluster_overview$cluster_id[hypo_map_cluster_overview$cluster_id %in% cell_cluster_map$cluster_id])) %>%
  dplyr::arrange(cluster_id) %>% # sort similar to hypoMap overview
  dplyr::select(projected_cluster_id = cluster_id, projected_cluster_name = cluster_name, clusterlevel, ncells_nuc, pct_of_all_nuc, pct_of_all_ref, likely_region) # reorder and name columns 
# add correlations
per_cluster_correlations = data.table::fread("figure_outputs/figure_4/per_cluster_correlations.txt")[,2:ncol(per_cluster_correlations)]
colnames(per_cluster_correlations)[3:length(colnames(per_cluster_correlations))] = paste0(colnames(per_cluster_correlations)[3:length(colnames(per_cluster_correlations))],"_cor")
nucseq_cluster_overview = nucseq_cluster_overview %>% dplyr::left_join(per_cluster_correlations,by=c("projected_cluster_id"="K169_pruned"))
# add N degs in fasting
all_clusters_fasting_DEG = data.table::fread("figure_outputs/figure_supplementary_7/all_clusters_fasting_DEG.txt")
all_clusters_fasting_DEG_stat = all_clusters_fasting_DEG %>% dplyr::filter(p_val_adj < 0.01) %>% dplyr::group_by(current_cluster) %>% dplyr::count(name = "n_degs_fasting")
nucseq_cluster_overview = nucseq_cluster_overview %>% dplyr::left_join(all_clusters_fasting_DEG_stat,by=c("projected_cluster_name"="current_cluster"))
# save:
data.table::fwrite(nucseq_cluster_overview,paste0(results_path_tables,"supplementary_table_13_nucseq_cluster_overview.csv"))

##########
### Markers & DEGs
##########

## table with hypomap cluster marker genes
hypomap_marker_genes = neuron_map_seurat@misc$markers_comparisons_all
hypomap_marker_genes = hypomap_marker_genes %>% dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_name,cluster_id),by=c("cluster_1"="cluster_id"))
# maybe restrict to high specificity to cut down size ?
hypomap_marker_genes = hypomap_marker_genes %>% dplyr::filter(specificity > 1 & p_val_adj < 0.001) %>% 
  dplyr::select(cluster_id = cluster_1,cluster_name, gene, specificity, p_val_adj, avg_log2FC =  avg_logFC, pct.1, pct.2)
# save:
data.table::fwrite(hypomap_marker_genes,paste0(results_path_tables,"supplementary_table_14_hypomap_marker_genes.csv"))


## table with nucseq cluster marker genes
nucseq_marker_genes = data.table::fread("data_inputs/sn_seq_mapped_neurons_K169_markers_2_sampleID.txt")
nucseq_marker_genes = nucseq_marker_genes %>% dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_name,cluster_id),by=c("cluster"="cluster_id"))
# maybe restrict to high specificity to cut down size ?
nucseq_marker_genes = nucseq_marker_genes %>% dplyr::mutate(specificity = avg_logFC* (pct.1/pct.2)) %>% dplyr::filter(specificity > 1 & p_val_adj < 0.001) %>% 
  dplyr::select(cluster_id = cluster,cluster_name, gene, specificity, p_val_adj, avg_log2FC = avg_logFC, pct.1, pct.2) %>%
  dplyr::filter(!is.na(cluster_id)) 
nucseq_marker_genes$cluster_id = factor(nucseq_marker_genes$cluster_id,levels = unique(hypomap_marker_genes$cluster_id[hypomap_marker_genes$cluster_id %in% nucseq_marker_genes$cluster_id]))
nucseq_marker_genes = nucseq_marker_genes %>% dplyr::arrange(cluster_id) # sort similar to hypoMap overview
# save:
data.table::fwrite(nucseq_marker_genes,paste0(results_path_tables,"supplementary_table_15_nucseq_marker_genes.csv"))


## table witH DEGs
all_clusters_fasting_DEG = data.table::fread("figure_outputs/figure_supplementary_7/all_clusters_fasting_DEG.txt")
all_clusters_fasting_DEG = all_clusters_fasting_DEG %>% dplyr::left_join(hypo_map_cluster_overview %>% dplyr::select(cluster_name,cluster_id),by=c("current_cluster"="cluster_name"))
all_clusters_fasting_DEG = all_clusters_fasting_DEG %>% dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::select(cluster_id,cluster_name = current_cluster, gene, avg_log2FC, p_val_adj, avg_log2FC, pct.1, pct.2,pct_diff)
# save:
data.table::fwrite(all_clusters_fasting_DEG,paste0(results_path_tables,"supplementary_table_16_nucseq_DEG_fasting.csv"))

## table with go enrichment results for Agrp ?
agrp_fasting_go_enrichment_simplified = data.table::fread("figure_outputs/figure_5/agrp_fasting_go_enrichment_simplified.txt")
#agrp_fasting_go_enrichment_simplified = agrp_fasting_go_enrichment_simplified %>% dplyr::select(ID)
data.table::fwrite(agrp_fasting_go_enrichment_simplified,paste0(results_path_tables,"supplementary_table_17_nucseq_gobp_Agrp_fasting.csv"))


##########
### ISH results
##########

# This is the input table for Figure 6 which is based on the manual counting of images.
pct_expressed_cells_clusters_ISH = data.table::fread("figure_outputs/figure_6/pct_expressed_cells_clusters_ISH.txt")
data.table::fwrite(pct_expressed_cells_clusters_ISH,paste0(results_path_tables,"supplementary_table_18_rnascope_glp1r_result.csv"),scipen = fwrite_scipen)

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
openxlsx::write.xlsx(x = supp_table_list,file = "table_outputs/supplementary_tables_preliminary.xlsx",colNames=TRUE,rowNames=FALSE,keepNA=FALSE)

###### I manually created the file "supplementary_tables.xlsx" that includes the full column descriptions (not only skeleton) !!







