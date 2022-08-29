##########
### Load & Prepare
##########

#set path
results_path_extended_figure5 = "figure_outputs/figure_extended_5/"
system(paste0("mkdir -p ",results_path_extended_figure5))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
source("R/utility_functions.R")
source("R/plot_functions.R")

# load seurat objects via large_data_path
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/"
load_required_files(large_data_path = large_data_path)

# load colors 
load_colors()
getLongPalette = colorRampPalette(long_palette_strong)
getOkabeItoPalette = colorRampPalette(short_palette)
color_value_vector =unlist(jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scHarmonization/data/region_prediction_mapping_colors.json"))

## plotting
load_plot_params()

##########
### Prepare heatmaps
##########


leaf_level_column = "C185"
leaf_level = 6

# http://yulab-smu.top/treedata-book/chapter2.html

# make data for first heatmap with percentages per dataset
heatmap_data = hypoMap_v2_seurat@meta.data %>% dplyr::select(Cell_ID,Dataset,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Dataset) %>% #dplyr::filter(predicted_Campbell!="NA") 
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label")) 
  tidyr::spread(key = Dataset,value=presence) %>% as.data.frame()
heatmap_matrix = as.matrix(heatmap_data[,2:ncol(heatmap_data)])
rownames(heatmap_matrix) = heatmap_data[,leaf_level_column]
heatmap_matrix[is.na(heatmap_matrix)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix)),Dataset = colnames(heatmap_matrix))
data.table::fwrite(colnames_overview,paste0(results_path_extended_figure5,"circular_tree_colnames.txt"),sep="\t")
colnames(heatmap_matrix) = colnames_overview$number

# make data for second heatmap with n cells
heatmap_data2 = hypoMap_v2_seurat@meta.data %>% dplyr::select(Cell_ID,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column)) %>% #dplyr::filter(predicted_Campbell!="NA") 
  dplyr::count(name = "ncells")  %>% dplyr::ungroup()  %>% dplyr::mutate(ncells_pct = ncells / sum(ncells)*100)  %>% as.data.frame()
heatmap_matrix2 = as.matrix(heatmap_data2[,"ncells_pct",drop=FALSE])
colnames(heatmap_matrix2) = "%"
rownames(heatmap_matrix2) = heatmap_data2[,leaf_level_column]
heatmap_matrix2[is.na(heatmap_matrix2)] = 0

# make data for third heatmap with regions
heatmap_data3 = hypoMap_v2_seurat@meta.data %>% dplyr::select(Cell_ID,Region_summarized,!!sym(leaf_level_column)) %>%
  dplyr::group_by(!!sym(leaf_level_column),Region_summarized) %>% dplyr::count() %>% dplyr::group_by(!!sym(leaf_level_column)) %>%
  dplyr::top_n(n = 1,wt = n) %>% ungroup() %>%
  dplyr::distinct(!!sym(leaf_level_column),Region_summarized,.keep_all=TRUE) %>% as.data.frame()
heatmap_matrix3 = as.matrix(heatmap_data3[,"Region_summarized",drop=F])
rownames(heatmap_matrix3) = heatmap_data3[,leaf_level_column]
colnames(heatmap_matrix3) = "R"

## make annotation data frame
anno_df = hypoMap_v2_seurat@misc$annotation_result %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})


##########
### tree new - neurons
##########

root = "C2-1"

edgelist = hypoMap_v2_seurat@misc$clustering_edgelist
edgelist_neurons = edgelist[edgelist$to %in% scUtils::find_children(nodes = root,edges = edgelist),]

# plot cluster tree:
tree_color = "grey70"
circular_tree = plot_cluster_tree(edgelist = edgelist_neurons,
                                  leaf_level=leaf_level,
                                  anno_df = anno_df ,
                                  metadata=hypoMap_v2_seurat@meta.data,
                                  label_size = 3.5, 
                                  show_genes = TRUE,
                                  vjust_label = -0.25,
                                  edge_color = tree_color, 
                                  node_color = tree_color)
circular_tree = rotate_tree(circular_tree, -90)
circular_tree

# plot tree with heatmap 1
circular_tree_heat = add_heatmap(circular_tree=circular_tree,
                                 heatmap_matrix = heatmap_matrix,
                                 heatmap_colors=c(bg_col,"darkred"),
                                 scale_limits = c(0,100),
                                 heatmap_colnames =TRUE, 
                                 legend_title = "Pct Dataset",
                                 matrix_offset = 0.2,
                                 matrix_width =0.4,
                                 colnames_angle=0,
                                 legend_text_size = 5,
                                 hjust_colnames=1,
                                 na_color = "white",
                                 heatmap_text_size=2)
# plot tree with heatmaps
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix2,
                                 heatmap_colors=c(bg_col,"#17056e"),
                                 scale_limits = c(0,2),
                                 heatmap_colnames =TRUE, 
                                 legend_title = "Cells",
                                 matrix_offset = 1.8,
                                 matrix_width =0.05,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=3, 
                                 na_color = "white",
                                 heatmap_text_size=3)
#circular_tree_heat

# plot tree with heatmap2
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix3,
                                 heatmap_colors=c(bg_col,"darkred"),
                                 heatmap_colnames =TRUE, 
                                 legend_title = "Pct Dataset",
                                 matrix_offset = 2.05,
                                 matrix_width =0.05,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=4, 
                                 na_color = "grey80",
                                 heatmap_text_size=3)
# change colors for regions:
circular_tree_heat = circular_tree_heat + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")

# show
circular_tree_heat

#store:
ggsave(filename = paste0(results_path_extended_figure5,"circular_tree_heat_neuron.png"),
       plot = circular_tree_heat, "png",dpi=600,width=400,height = 400,units="mm")

ggsave(filename = paste0(results_path_extended_figure5,"circular_tree_heat_neuron.pdf"),
       plot = circular_tree_heat, "pdf",dpi=600,width=400,height = 400,units="mm")

# source:
neuron_clusters = scUtils::find_children("C2-1",hypoMap_v2_seurat@misc$clustering_edgelist)
figure5_neurons_tree_heatmap_sourcedata = heatmap_data %>% bind_cols(heatmap_data2 %>% dplyr::select(-C185))  %>% bind_cols(heatmap_data3 %>% dplyr::select(-C185,-n)) %>% dplyr::filter(C185 %in% neuron_clusters)
data.table::fwrite(figure5_neurons_tree_heatmap_sourcedata,paste0(results_path_extended_figure5,"source_figure4e_neuron_tree_heatmap.txt"),sep="\t")

data.table::fwrite(edgelist_neurons,paste0(results_path_extended_figure5,"source_figure4e_neuron_tree_edgelist.txt"),sep="\t")

data.table::fwrite(anno_df[anno_df$cluster_id %in% c("C2-1",neuron_clusters),],paste0(results_path_extended_figure5,"source_figure4e_neuron_tree_annotations.txt"),sep="\t")

##########
### tree new - non-neurons
##########

root = "C2-2"

edgelist = hypoMap_v2_seurat@misc$clustering_edgelist
edgelist_nonneurons = edgelist[edgelist$to %in% scUtils::find_children(nodes = root,edges = edgelist),]


# plot cluster tree:
tree_color = "grey70"
circular_tree = plot_cluster_tree(edgelist = edgelist_nonneurons,
                                  leaf_level=leaf_level,
                                  anno_df = anno_df ,
                                  metadata=hypoMap_v2_seurat@meta.data,
                                  label_size = 3.5, 
                                  show_genes = TRUE,
                                  vjust_label = -0.25,
                                  edge_color = tree_color, 
                                  node_color = tree_color)
circular_tree = rotate_tree(circular_tree, -90)
circular_tree


# plot tree with heatmap 1
circular_tree_nonneuron = add_heatmap(circular_tree=circular_tree,
                                      heatmap_matrix = heatmap_matrix,
                                      heatmap_colors=c(bg_col,"darkred"),
                                      scale_limits = c(0,100),
                                      heatmap_colnames =TRUE, 
                                      legend_title = "Pct Dataset",
                                      matrix_offset = 0.2,
                                      matrix_width =0.4,
                                      colnames_angle=0,
                                      legend_text_size = 5,
                                      hjust_colnames=1,
                                      na_color = "white",
                                      heatmap_text_size=2)
# plot tree with heatmaps
circular_tree_nonneuron = add_heatmap(circular_tree=circular_tree_nonneuron,
                                      heatmap_matrix = heatmap_matrix2,
                                      heatmap_colors=c(bg_col,"#17056e"),
                                      scale_limits = c(0,2),
                                      heatmap_colnames =TRUE, 
                                      legend_title = "Cells",
                                      matrix_offset = 1.8,
                                      matrix_width =0.05,
                                      colnames_angle=0,
                                      legend_text_size = 8,
                                      hjust_colnames=3, 
                                      na_color = "white",
                                      heatmap_text_size=3)

# show
circular_tree_nonneuron

#store:
ggsave(filename = paste0(results_path_extended_figure5,"circular_tree_heat_nonneuron.png"),
       plot = circular_tree_nonneuron, "png",dpi=600,width=400,height = 400,units="mm")

ggsave(filename = paste0(results_path_extended_figure5,"circular_tree_heat_nonneuron.pdf"),
       plot = circular_tree_nonneuron, "pdf",dpi=600,width=400,height = 400,units="mm")

# source:
nonneuron_clusters = scUtils::find_children("C2-2",hypoMap_v2_seurat@misc$clustering_edgelist)
figure5_nonneurons_tree_heatmap_sourcedata = heatmap_data %>% bind_cols(heatmap_data2 %>% dplyr::select(-C185))  %>% bind_cols(heatmap_data3 %>% dplyr::select(-C185,-n)) %>% dplyr::filter(C185 %in% nonneuron_clusters)
data.table::fwrite(figure5_nonneurons_tree_heatmap_sourcedata,paste0(results_path_extended_figure5,"source_figure4e_nonneuron_tree_heatmap.txt"),sep="\t")

data.table::fwrite(edgelist_nonneurons,paste0(results_path_extended_figure5,"source_figure4e_nonneuron_tree_edgelist.txt"),sep="\t")

data.table::fwrite(anno_df[anno_df$cluster_id %in% c("C2-2",nonneuron_clusters),],paste0(results_path_extended_figure5,"source_figure4e_nonneuron_tree_annotations.txt"),sep="\t")

