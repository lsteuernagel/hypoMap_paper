##########
### Load & Prepare
##########

#set path
results_path_extended_figure4 = "figure_outputs/figure_extended_4/"
system(paste0("mkdir -p ",results_path_extended_figure4))

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

## plotting
load_plot_params()


##########
### circular_tree BIG
##########

leaf_level_column = "C465"
leaf_level =8

# http://yulab-smu.top/treedata-book/chapter2.html

## make annotation data frame
anno_df = hypoMap_v2_seurat@misc$annotation_result %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

# plot cluster tree:
tree_color = "grey70"
circular_tree_big = plot_cluster_tree(edgelist = hypoMap_v2_seurat@misc$clustering_edgelist,
                                  leaf_level=leaf_level,
                                  anno_df = anno_df ,
                                  metadata=hypoMap_v2_seurat@meta.data,
                                  label_size = 2.5,
                                  label_size_tip = 1.75,
                                  show_genes = TRUE,
                                  vjust_label = -0.25,
                                  edge_color = tree_color, 
                                  node_color = tree_color)
circular_tree_big = rotate_tree(circular_tree_big, -90)
circular_tree_big

#store:
ggsave(filename = paste0(results_path_extended_figure4,"circular_tree_big.png"),
       plot = circular_tree_big, "png",dpi=600,width=400,height = 400,units="mm")

ggsave(filename = paste0(results_path_extended_figure4,"circular_tree_big.pdf"),
       plot = circular_tree_big, "pdf",dpi=600,width=400,height = 400,units="mm")

# source
data.table::fwrite(hypoMap_v2_seurat@misc$clustering_edgelist,paste0(results_path_extended_figure4,"source_ext_figure4_tree_edgelist.txt"),sep="\t")
data.table::fwrite(anno_df,paste0(results_path_extended_figure4,"source_ext_figure4_tree_annotations.txt"),sep="\t")





