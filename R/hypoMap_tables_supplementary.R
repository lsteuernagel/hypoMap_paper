##########
### Load & Prepare
##########

results_path_tables = "table_outputs/"
system(paste0("mkdir -p ",results_path_tables))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
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
### Cluster overview tables
##########

## Table with all clusters in HypoMap
# name, id, n_cells, region, ????
region_anno =  neuron_map_seurat@meta.data %>% dplyr::distinct(suggested_region_curated,K169_pruned) %>% dplyr::select(cluster_id = K169_pruned,likely_region = suggested_region_curated)
hypo_map_cluster_overview = dplyr::left_join(neuron_map_seurat@misc$annotations,region_anno,by="cluster_id") %>% 
  dplyr::group_by(clusterlevel) %>% dplyr::mutate(pct_of_all = round(ncells / sum(ncells), 4) * 100) %>% dplyr::ungroup()
# save:
data.table::fwrite(hypo_map_cluster_overview,paste0(results_path_tables,"supplementary_table_3_hypomap_cluster_overview.csv"))

## Table with mrtree edgelist for TREE
edgelist_clusters_hypoMap = neuron_map_seurat@misc$mrtree_edgelist %>% dplyr::select(-isLeaf,-level,-totalCount,nchildren = count)
# save:
data.table::fwrite(edgelist_clusters_hypoMap,paste0(results_path_tables,"supplementary_table_4_hypomap_cluster_edgelist.csv"))


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
data.table::fwrite(nucseq_cluster_overview,paste0(results_path_tables,"supplementary_table_5_nucseq_cluster_overview.csv"))


