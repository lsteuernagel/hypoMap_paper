# map romanov using old package version

library(mapscvi)
require(dplyr)
require(ggplot2)
require(Seurat)
source("R/utility_functions.R")
source("R/plot_functions.R")

# specify query data:
romanov_original_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypothalamus/romanov/hypothalamus_romanov_seurat_051020.rds") # seurat object to load

# prepare:
romanov_prepared_seurat = prepare_query(romanov_original_seurat,suffix="romanov")

# load var features
var_features = unlist(jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/feature_set.json"))

# predict
romanovpredicted_seurat = predict_query(query_seurat_object = romanov_prepared_seurat,
                                        model_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/hypoMap_harmonized_scVI_model/",
                                        query_reduction="scvi",
                                        var_names=var_features,
                                        max_epochs = 30,
                                        assay="RNA",
                                        use_reticulate = FALSE,
                                        global_seed=12345)


# load  reference map
#hypoMap_harmonized_curated = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/hypoMap_harmonized_curated.rds")
# load seurat objects via large_data_path
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/"
load_required_files(large_data_path = large_data_path)

# make label_vec
reference_reduction = "scvi"
reference_map_metadata = hypoMap_v2_seurat@meta.data
label_vec=as.character(reference_map_metadata$C465)
ref_umap = hypoMap_v2_seurat@reductions[[paste0("umap_",reference_reduction)]]

# project
romanovpredicted_seurat = project_query(query_seurat_object = romanovpredicted_seurat,
                                        reference_map_reduc = hypoMap_v2_seurat@reductions[[reference_reduction]],
                                        reference_map_umap = hypoMap_v2_seurat@reductions[[paste0("umap_",reference_reduction)]],
                                        label_vec =label_vec,
                                        use_projectUMAP = FALSE,
                                        n_neighbors = 25,
                                        global_seed=12345)
# rename:
romanovpredicted_seurat@meta.data$C465_predicted = romanovpredicted_seurat@meta.data$predicted
romanovpredicted_seurat@meta.data$C465_prediction_probability = romanovpredicted_seurat@meta.data$prediction_probability

# add other clusters levels of tree:
other_levels = mapscvi::reconstruct_levels(hypoMap_v2_seurat@misc$clustering_edgelist,input_annotation = romanovpredicted_seurat@meta.data$C465_predicted,level_prefix = "C[0-9]+",result_prefix="predicted_")
# throws an error ....

# predict others
prediction_C25_named= mapscvi::propagate_labels_prob(neighbors_object = romanovpredicted_seurat@neighbors$query_ref_nn,
                                                                                   label_vec = hypoMap_v2_seurat@meta.data$C25_named,
                                                                                   query_seurat_object = romanovpredicted_seurat,
                                                                                   k.param = 25)
romanovpredicted_seurat@meta.data$C25_named_predicted = prediction_C25_named$predicted
romanovpredicted_seurat@meta.data$C25_named_predicted_prob = prediction_C25_named$prediction_probability

##
prediction_C185_named = mapscvi::propagate_labels_prob(neighbors_object = romanovpredicted_seurat@neighbors$query_ref_nn,
                                               label_vec = hypoMap_v2_seurat@meta.data$C185_named,
                                               query_seurat_object = romanovpredicted_seurat,
                                               k.param = 25)
romanovpredicted_seurat@meta.data$C185_named_predicted = prediction_C185_named$predicted
romanovpredicted_seurat@meta.data$C185_named_predicted_prob = prediction_C185_named$prediction_probability


# save
saveRDS(romanovpredicted_seurat,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/mapped_data/romanov_predicted_seurat.rds")

romanovpredicted_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/mapped_data/romanov_predicted_seurat.rds")

length(table(romanovpredicted_seurat@meta.data$C185_named_predicted))

# # make plot
# romanov_mapped_plot = mapscvi::plot_query_labels(query_seura_object=romanovpredicted_seurat,
#                                                  reference_seurat=hypoMap_v2_seurat,
#                                                  label_col="C25_named",
#                                                  label_col_query = "C25_predicted",
#                                                  overlay = TRUE,
#                                                  bg_col = "grey90",
#                                                  query_pt_size = 0.6,
#                                                  labelonplot = TRUE,
#                                                  label.size=5,
#                                                  repel=TRUE)+ggtitle("Romanov et al. mapped on HypoMap")
# 
# romanov_mapped_plot_R = rasterize_ggplot(romanov_mapped_plot,pixel_raster = rasterize_pixels,pointsize = 3.4)
# romanov_mapped_plot_R



