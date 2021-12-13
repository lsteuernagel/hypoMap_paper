
##########
### Define paths
##########

small_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/hypoMap_paper/data_inputs/"
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"
# define some relevant files for loading of data
key="hypothalamusMapNeurons_v4"
scHarmonize_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/"

##########
### move seurat objects into one folder
##########

## see the romanov_scvi RMD in /scHarmonize/scripts/mapscvi or the mapscvi vignette
#query_romanov_neurons = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypothalamus/romanov/mapped_data/query_romanov_neurons.rds")
system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypothalamus/romanov/mapped_data/query_romanov_neurons.rds ",large_data_path))

# nuc seq all
system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/nucseq_all_map.rds ",large_data_path))

# nuc seq neurons
system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/nucseq_neurons_map.rds ",large_data_path))

# map full
system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/harmonization_results/hypothalamus_full_map.rds ",large_data_path))

# map neurons
system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_objects/hypothalamus_neurons_reference.rds ",large_data_path))

# temp
# system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/agrp_fasting_all.txt ","figure_outputs/figure_5/"))
# system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/agrp_fasting_all_campbell.txt ","figure_outputs/figure_5/"))
# system(command = paste0("cp ","/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/Sox14_fasting_all.txt ","figure_outputs/figure_5/"))
# system(command = paste0("cp ", "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/hk2_fasting_all.txt ","figure_outputs/figure_5/"))
# system(command = paste0("cp ", "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/global_fasting_all.txt ","figure_outputs/figure_5/"))


##########
### Existing Supplementary tables
##########

results_path_tables = "table_outputs/"
system(paste0("mkdir -p ",results_path_tables))

##### Move to collect data:
# load the manually created supplementary table
results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/paper_results/"
dataset_overview = data.table::fread(paste0(results_path,"Supplementary_Tables/supplementalTable_1_dataset_overview.csv"))
dataset_overview = dataset_overview %>% dplyr::select(-Age)
# save in supplement
data.table::fwrite(dataset_overview,paste0(results_path_tables,"supplementary_table_1_dataset_overview.csv"))

# curated neuron signatures
curated_neuron_celltypes = data.table::fread(paste0(results_path,"Supplementary_Tables/supplementalTable_2_purity_neuron_signatures.csv"))
# restrict to selected ones and remove col
curated_neuron_celltypes_restrict = curated_neuron_celltypes %>% dplyr::filter(selected_for_purity) %>% dplyr::select(-selected_for_purity)
# save in supplement
data.table::fwrite(curated_neuron_celltypes_restrict,paste0(results_path_tables,"supplementary_table_2_neuron_signatures.csv"))

## regions:
# all_target_regions_result =  data.table::fread(paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_objects/region_prediction_results_suggested_hypothalamus_neurons_reference.txt"))
# data.table::fwrite(sn_seq_markers_K169,paste0(small_data_path,"sn_seq_mapped_neurons_K169_markers_2_sampleID.txt"),sep="\t")

##########
### move reductions
##########

system(paste0("mkdir -p ",paste0(large_data_path,"best_reductions_per_method/")))

# add reductions:
key="hypothalamusMapNeurons_v4"
scHarmonize_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/"
for( i in 1:length(names_to_check)){
  integration_res_path = paste0(scHarmonize_path,key,"/","integration_results","/")
  origin_path = paste0(integration_res_path,list.files(path=integration_res_path,pattern = names_to_check[i],recursive = TRUE)[1])
  system(command = paste0("cp ",origin_path," ",large_data_path,"best_reductions_per_method/"))
 }

##########
### sn seq markers results
##########

# load marker genes sn seq
output_file_name = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/","sn_seq_mapped_neurons_K169_markers_2_sampleID.txt")
sn_seq_markers_K169 = data.table::fread(output_file_name,data.table = F)
data.table::fwrite(sn_seq_markers_K169,paste0(small_data_path,"sn_seq_mapped_neurons_K169_markers_2_sampleID.txt"),sep="\t")

##########
### ISH results
##########

#### load ish results
ish_quantification = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/harmonization_results/hypothalamus_neurons_reference/paper_results/ish_quantification_glp1r_updated.csv",data.table=FALSE)
data.table::fwrite(ish_quantification,paste0(small_data_path,"ish_quantification_glp1r_updated.csv"),sep="\t")

##########
### metrics from integration
##########

# neuron metric results
neurons_metrics = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapNeurons_v4/documentation/hypothalamusMapNeurons_v4_comparison_457fc60c3c4f1911bcbc6c5d46127037.txt",data.table = F)
data.table::fwrite(neurons_metrics,paste0(small_data_path,"hypothalamusMapNeurons_v4_comparison_457fc60c3c4f1911bcbc6c5d46127037.txt"),sep="\t")
### load comparison data full
full_metrics = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/documentation/hypothalamusMapFull_v4_comparison_8af8a1cd950067bb6859cbc1225c818d.txt",data.table = F)
data.table::fwrite(full_metrics,paste0(small_data_path,"hypothalamusMapFull_v4_comparison_8af8a1cd950067bb6859cbc1225c818d.txt"),sep="\t")


##########
### celltype mapping etc.
##########

## relevant files
mapped_celltypes_file=paste0(scHarmonize_path,key,"/","seurat_objects","/",key,"_mapped_celltypes.rds")
id_file_name = paste0(scHarmonize_path,key,"/","seurat_objects","/",key,"_subsampled_Cell_IDs.txt")

# get mapped celltypes
print(mapped_celltypes_file)
mapped_celltypes =readRDS(mapped_celltypes_file)
saveRDS(mapped_celltypes,paste0(small_data_path,"mapped_celltypes_neuronMap.rds"))

## subsample ids as start for plotting
subsample_ids_df = data.table::fread(id_file_name,data.table = F,header = F)
data.table::fwrite(subsample_ids_df,paste0(small_data_path,"_subsampled_Cell_IDs_neuronMap.txt"),sep="\t")

##########
### full evaluation results
##########

## original metadata
seurat_metasmall_data_path = paste0(scHarmonize_path,key,"/","seurat_objects","/",key,"_processed_metadata.txt")
meta_data = data.table::fread(seurat_metasmall_data_path,data.table = FALSE)
# just use from finals seurat object ?

## subsample ids as start for plotting
id_file_name = paste0(scHarmonize_path,key,"/","seurat_objects","/",key,"_subsampled_Cell_IDs.txt")
subsample_ids = data.table::fread(id_file_name,data.table = F,header = F)[,1]

# evaluation res file names
evaluation_file_path = paste0(scHarmonize_path,key,"/evaluation_results/")
evaluation_mixing_probNorm_file = paste0(evaluation_file_path,"evaluation_mixing_probNorm.Batch_ID.20000.100.123467_all.txt")
evaluation_purity_knn_file = paste0(evaluation_file_path,"evaluation_purity_knn.24.20.123467_all.txt") 
evaluation_entropy_knn_all_file = paste0(evaluation_file_path,"evaluation_entropy_knn_all.Batch_ID.20.123467_all.txt") 

##### load evaluation results
evaluation_mixing_probNorm_all = data.table::fread(evaluation_mixing_probNorm_file,data.table = F)
rownames(evaluation_mixing_probNorm_all) = subsample_ids
evaluation_purity_knn_all = data.table::fread(evaluation_purity_knn_file,data.table = F)
rownames(evaluation_purity_knn_all) = names(mapped_celltypes)
evaluation_entropy_knn_all = data.table::fread(evaluation_entropy_knn_all_file,data.table = F)
rownames(evaluation_entropy_knn_all) = meta_data$Cell_ID

# save results with IDs
#### I Save these in a different folde rbceause they are large!!
evaluation_mixing_probNorm_all = cbind(Cell_ID=rownames(evaluation_mixing_probNorm_all),evaluation_mixing_probNorm_all)
data.table::fwrite(evaluation_mixing_probNorm_all,paste0(large_data_path,"evaluation_mixing_probNorm.Batch_ID.20000.100.123467_all.txt"),sep="\t")
evaluation_entropy_knn_all = cbind(Cell_ID=rownames(evaluation_entropy_knn_all),evaluation_entropy_knn_all)
data.table::fwrite(evaluation_entropy_knn_all,paste0(large_data_path,"evaluation_entropy_knn_all.Batch_ID.20.123467_all.txt"),sep="\t")
# save in regular folder
evaluation_purity_knn_all = cbind(mapped_celltype=rownames(evaluation_purity_knn_all),evaluation_purity_knn_all)
data.table::fwrite(evaluation_purity_knn_all,paste0(small_data_path,"evaluation_purity_knn.24.20.123467_all.txt"),sep="\t")

##########
### bacTRAP signature file
##########

# signature files is generated via all_bacTRAP_signatures.R from bacTRAP_map project 
bacTRAP_signature_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/signatures/bacTRAP_signatures_rbo_2207.rds"
bacTRAP_signatures = readRDS(bacTRAP_signature_file)
saveRDS(bacTRAP_signatures,paste0(small_data_path,"bacTRAP_signatures_rbo.rds"))

##########
### original UMAP
##########

snuc_hypo_master = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/snuc_hypo_master_211101.RDS")
original_umap = snuc_hypo_master@reductions$umap@cell.embeddings
data.table::fwrite(original_umap,paste0(small_data_path,"nucSeq_originalUMAP.txt"),sep="\t")


