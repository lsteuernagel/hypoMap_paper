##########
### Load & Prepare
##########

results_path_extended_figure6 = "figure_outputs/figure_extended_6/"
system(paste0("mkdir -p ",results_path_extended_figure6))


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
### plot selected cell type
##########

anno_df = hypoMap_v2_seurat@misc$annotation_result

# male stacked vilin plot of marker genes for selected cluster 

## selected celltypes:

neuron_celltypes = hypoMap_v2_seurat@meta.data$C185[hypoMap_v2_seurat@meta.data$C2=="C2-1"]
celltypes= c("C185-115","C185-48")#c("C185-65","C185-47","C185-22","C185-103")

all_vln_plots = list()
for(i in 1:length(celltypes)){
  print(i)
  celltype_to_use = celltypes[i]
  marker_genes = hypoMap_v2_seurat@misc$marker_genes_all[hypoMap_v2_seurat@misc$marker_genes_all$cluster_id==celltype_to_use,] 
  marker_genes$log_pval_adj = (-log10(marker_genes$p_val_adj))
  marker_genes$log_pval_adj[marker_genes$log_pval_adj > 300] = 300
  marker_genes$score = marker_genes$specificity * marker_genes$log_pval_adj
  marker_genes = marker_genes %>% top_n(n = 5,wt = score)
  hypoMap_v2_seurat_celltype_subset = subset(hypoMap_v2_seurat,subset = C185 == celltype_to_use)
  Idents(hypoMap_v2_seurat_celltype_subset) = "Dataset"
  fullname= anno_df$clean_names_withID[anno_df$cluster_id == celltype_to_use ]
  #color by Dataset:
  # stacked_vln = VlnPlot(hypoMap_v2_seurat_celltype_subset,features = marker_genes$gene,stack=TRUE,split.by = "Dataset",cols=getOkabeItoPalette(18))+
  #   theme(text = element_text(size=text_size)) %>% suppressWarnings()
  # color by gene
  stacked_vln = Seurat::VlnPlot(hypoMap_v2_seurat_celltype_subset,features = marker_genes$gene,stack=TRUE,cols=getOkabeItoPalette(7))+
    theme(text = element_text(size=text_size),axis.text.x =  element_text(size=14),axis.text.y =  element_text(size=20)) + ylab("Dataset") +
    ggtitle(fullname)+NoLegend() %>% suppressWarnings()
  
  #manually add colors:
  #stacked_vln$layers[[1]]$aes_params$colour <- cols_for_feature_plot[2]
  stacked_vln$layers[[1]]$aes_params$fill <- cols_for_feature_plot[2]

  all_vln_plots[[celltype_to_use]] = stacked_vln
  
}

## save
for(i in 1:length(all_vln_plots)){
  print(i)
  current_name = names(all_vln_plots)[i]
  current_plot = all_vln_plots[[current_name]]
 
  #current_plot_r = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)
  # and save:
  current_name = gsub("-","_",current_name)
  ggsave(filename = paste0(results_path_extended_figure6,current_name,"_top5_features.png"),
         plot = current_plot, "png",dpi=600,width=400,height = 300,units="mm")
  ggsave(filename = paste0(results_path_extended_figure6,current_name,"_top5_features.pdf"),
         plot = current_plot, "pdf",dpi=600,width=400,height = 300,units="mm")
  
}

# source
extended_figure6_pomc = all_vln_plots$`C185-48`$data
data.table::fwrite(extended_figure6_pomc,paste0(results_path_extended_figure6,"source_ext_figure6_a_pomc.txt"),sep="\t")
extended_figure6_agrp = all_vln_plots$`C185-115`$data
data.table::fwrite(extended_figure6_agrp,paste0(results_path_extended_figure6,"source_ext_figure6_b_agrp.txt"),sep="\t")

