
##########
### Load & Prepare
##########

results_path_figure2 = "figure_outputs/figure_2/"
system(paste0("mkdir -p ",results_path_figure2))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
source("R/utility_functions.R")
source("R/plot_functions.R")

# where to find large data objects (need to change to local dir !)
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/"

# load seurat objects via large_data_path
load_required_files(large_data_path = large_data_path)
#hypoMap_v2_seurat@meta.data$Author_Class_Curated[hypoMap_v2_seurat@meta.data$Author_Class_Curated=="Differentiating"] = "Dividing"

# load colors 
load_colors()
getLongPalette = colorRampPalette(long_palette_strong)
getOkabeItoPalette = colorRampPalette(short_palette)
color_value_vector =unlist(jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scHarmonization/data/region_prediction_mapping_colors.json"))


## plotting
load_plot_params()

##########
### tree new
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
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames.txt"),sep="\t")
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

# plot cluster tree:
tree_color = "grey70"
circular_tree = plot_cluster_tree(edgelist = hypoMap_v2_seurat@misc$clustering_edgelist,
                                  leaf_level=leaf_level,
                                  anno_df = anno_df ,
                                  metadata=hypoMap_v2_seurat@meta.data,
                                  label_size = 2.5, 
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
                                 matrix_offset = 2.2,
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
                                 matrix_offset = 2.5,
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
ggsave(filename = paste0(results_path_figure2,"circular_tree_heat_label.png"),
       plot = circular_tree_heat, "png",dpi=600,width=400,height = 400,units="mm")

ggsave(filename = paste0(results_path_figure2,"circular_tree_heat_label.pdf"),
       plot = circular_tree_heat, "pdf",dpi=600,width=400,height = 400,units="mm")

#### source data tree

figure2_tree_heatmap_sourcedata = heatmap_data %>% bind_cols(heatmap_data2 %>% dplyr::select(-C185))  %>% bind_cols(heatmap_data3 %>% dplyr::select(-C185,-n))
data.table::fwrite(figure2_tree_heatmap_sourcedata,paste0(results_path_figure2,"source_figure2_a_tree_heatmap.txt"),sep="\t")

data.table::fwrite(hypoMap_v2_seurat@misc$clustering_edgelist,paste0(results_path_figure2,"source_figure2_a_tree_edgelist.txt"),sep="\t")

data.table::fwrite(anno_df,paste0(results_path_figure2,"source_figure2_a_tree_annotations.txt"),sep="\t")


##########
### dotplot
##########

## make annotation data frame
anno_df = hypoMap_v2_seurat@misc$annotation_result %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

anno_df_dotplot = anno_df

target_col="C66_named"
cluster_to_include = as.character(unique(hypoMap_v2_seurat@meta.data[hypoMap_v2_seurat@meta.data$C2_named=="C2-1: Neurons",target_col]))
#cluster_to_include = cluster_to_include[1:10]
cluster_to_include_wo = stringr::str_extract(cluster_to_include,"C66-[0-9]+")

# use just the top feature ?
# features_dotplot = hypoMap_v2_seurat@misc$marker_genes_all %>% dplyr::filter(specificity > 5 & cluster_id %in% cluster_to_include_wo) %>%
#   dplyr::filter(!grepl("Rik|Gm",gene)) %>%
#   dplyr::arrange(desc(specificity)) %>% dplyr::group_by(cluster_id) %>% dplyr::top_n(n = 1,wt = specificity)

anno_df_dotplot = anno_df[anno_df$cluster_id %in% cluster_to_include_wo,]
anno_df_dotplot$gene_dotplot = anno_df_dotplot$first_cluster_name
# use names:
anno_df_dotplot$gene_dotplot[! anno_df_dotplot$gene_dotplot %in% rownames(hypoMap_v2_seurat@assays$RNA@counts)] = NA

# need to reorder factor level in seurat
# order by number --> then it also is ordered by tree
unordered = unique(hypoMap_v2_seurat@meta.data[,target_col])
unordered = as.numeric(stringr::str_remove(stringr::str_extract(unordered,"-[0-9]+"),"-"))
names(unordered) = unique(hypoMap_v2_seurat@meta.data[,target_col])
hypoMap_v2_seurat@meta.data[,target_col] = factor(as.character(hypoMap_v2_seurat@meta.data[,target_col]),levels = names(sort(unordered)))
# also reorder anno_df_dotplot
rownames(anno_df_dotplot) = anno_df_dotplot$cluster_id
anno_df_dotplot = anno_df_dotplot[match(stringr::str_extract(names(sort(unordered)),"C66-[0-9]+"),rownames(anno_df_dotplot)),] %>% na.omit()

# plot
Idents(hypoMap_v2_seurat) = target_col
dotplot_neurons = Seurat::DotPlot(object = hypoMap_v2_seurat,
                                  features = unique(anno_df_dotplot$gene_dotplot[!is.na(anno_df_dotplot$gene_dotplot)]),
                                  idents= as.character(cluster_to_include),
                                  scale = FALSE,
                                  cluster.idents = F)
dotplot_neurons2 = dotplot_neurons + guides(color=guide_colourbar('Avg. Exp.'),size = guide_legend("Pct. Exp.")) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 22,angle = 90, vjust = 0.35, hjust=0.75),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_gradient(low = cols_for_feature_plot[1],high = cols_for_feature_plot[2],limits = c(0,4), oob = scales::squish)
dotplot_neurons2

# save
ggsave(filename = paste0(results_path_figure2,"dotplot_neurons.png"),
       plot = dotplot_neurons2, "png",dpi=400,width=400,height = 300,units="mm")
ggsave(filename = paste0(results_path_figure2,"dotplot_neurons.pdf"),
       plot = dotplot_neurons2, "pdf",dpi=400,width=400,height =300,units="mm")
ggsave(filename = paste0(results_path_figure2,"dotplot_neurons_500.pdf"),
       plot = dotplot_neurons2, "pdf",dpi=400,width=500,height =300,units="mm")
##source data
source_data_dotplot_neurons = dotplot_neurons$data
data.table::fwrite(source_data_dotplot_neurons,paste0(results_path_figure2,"source_figure2_b_dotplot_neurons.txt"),sep="\t")


##########
### dotplot
##########

## make annotation data frame
anno_df = hypoMap_v2_seurat@misc$annotation_result %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

anno_df_dotplot = anno_df

target_col="C66_named"
cluster_to_include = unique(hypoMap_v2_seurat@meta.data[hypoMap_v2_seurat@meta.data$C2_named=="C2-2: Non-Neurons",target_col])
#cluster_to_include = cluster_to_include[1:10]
cluster_to_include_wo = stringr::str_extract(cluster_to_include,"C66-[0-9]+")

#use just the top feature ?
features_dotplot = hypoMap_v2_seurat@misc$marker_genes_all %>% dplyr::filter(specificity > 5 & cluster_id %in% cluster_to_include_wo) %>%
  dplyr::filter(!grepl("Rik|Gm",gene)) %>%
  dplyr::arrange(desc(specificity)) %>% dplyr::group_by(cluster_id) %>% dplyr::top_n(n = 1,wt = specificity)

anno_df_dotplot = anno_df[anno_df$cluster_id %in% cluster_to_include_wo,]
anno_df_dotplot$gene_dotplot = anno_df_dotplot$first_cluster_name

anno_df_dotplot = dplyr::left_join(anno_df_dotplot,features_dotplot[,1:2],by=c("cluster_id"="cluster_id"))
anno_df_dotplot$gene_dotplot[anno_df_dotplot$cluster_id %in% c("C66-51","C66-52","C66-57","C66-61","C66-62","C66-63","C66-64")] = anno_df_dotplot$gene[anno_df_dotplot$cluster_id %in% c("C66-51","C66-52","C66-57","C66-61","C66-62","C66-63","C66-64")] 
anno_df_dotplot$gene_dotplot[anno_df_dotplot$cluster_id %in% c("C66-51")] = "Col23a1"
# need to reorder factor level in seurat
# order by number --> then it also is ordered by tree
unordered = unique(hypoMap_v2_seurat@meta.data[,target_col])
unordered = as.numeric(stringr::str_remove(stringr::str_extract(unordered,"-[0-9]+"),"-"))
names(unordered) = unique(hypoMap_v2_seurat@meta.data[,target_col])
hypoMap_v2_seurat@meta.data[,target_col] = factor(as.character(hypoMap_v2_seurat@meta.data[,target_col]),levels = names(sort(unordered)))
# also reorder anno_df_dotplot
rownames(anno_df_dotplot) = anno_df_dotplot$cluster_id
anno_df_dotplot = anno_df_dotplot[match(stringr::str_extract(names(sort(unordered)),"C66-[0-9]+"),rownames(anno_df_dotplot)),] %>% na.omit()

# plot
Idents(hypoMap_v2_seurat) = target_col
dotplot_non_neuron = Seurat::DotPlot(object = hypoMap_v2_seurat,
                                     features = anno_df_dotplot$gene_dotplot[!is.na(anno_df_dotplot$gene_dotplot)],
                                     idents= cluster_to_include,
                                     scale = FALSE,
                                     cluster.idents = F)#+facet_grid( ~ C23_named)
dotplot_non_neuron2 = dotplot_non_neuron + guides(color=guide_colourbar('Avg. Exp.'),size = guide_legend("Pct. Exp.")) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 25,angle = 90, vjust = 0.35, hjust=0.75),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_gradient(low = cols_for_feature_plot[1],high = cols_for_feature_plot[2],limits = c(0,4), oob = scales::squish)

dotplot_non_neuron2

# save
ggsave(filename = paste0(results_path_figure2,"dotplot_non_neurons.png"),
       plot = dotplot_non_neuron2, "png",dpi=400,width=400,height = 300,units="mm")
ggsave(filename = paste0(results_path_figure2,"dotplot_non_neurons.pdf"),
       plot = dotplot_non_neuron2, "pdf",dpi=400,width=400,height =300,units="mm")

##source data
source_data_dotplot_nonneurons = dotplot_non_neuron$data
data.table::fwrite(source_data_dotplot_nonneurons,paste0(results_path_figure2,"source_figure2_c_dotplot_nonneurons.txt"),sep="\t")
