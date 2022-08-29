
##########
### Load & Prepare
##########

results_path_figure6 = "figure_outputs/figure_6/"
system(paste0("mkdir -p ",results_path_figure6))

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
hypoMap_v2_seurat@meta.data$Author_Class_Curated[hypoMap_v2_seurat@meta.data$Author_Class_Curated=="Differentiating"] = "Dividing"


# load colors 
load_colors()
getLongPalette = colorRampPalette(long_palette_strong)
getOkabeItoPalette = colorRampPalette(short_palette)

## plotting
load_plot_params()

##########
### Figure 1c
##########

# load
romanovpredicted_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/mapped_data/romanov_predicted_seurat.rds")

compare_clustering_romanov =mapscvi::compare_clustering(romanovpredicted_seurat,clustering_1 = "Author_CellType",clustering_2 = "C185_named_predicted" ,
                                                        min_cells = 0,min_pct = 0,return_data = TRUE)
data.table::fwrite(compare_clustering_romanov,paste0(results_path_figure6,"compare_clustering_romanov.txt"),sep="\t")
data.table::fwrite(romanovpredicted_seurat@meta.data[,c("Cell_ID","C25_named_predicted","C185_named_predicted","Author_CellType")],paste0(results_path_figure6,"compare_clustering_romanov_full.txt"),sep="\t")


## romanov_full_overview
romanov_full_overview = romanovpredicted_seurat@meta.data %>% dplyr::select(Cell_ID,Author_CellType,C25_named_predicted,C185_named_predicted,C185_named_predicted_prob)
romanov_full_overview = cbind(romanov_full_overview,romanovpredicted_seurat@reductions$umap_scvi@cell.embeddings[,1:2])
data.table::fwrite(romanov_full_overview,paste0(results_path_figure6,"romanov_cells_projected.txt"),sep="\t")


# make plot
romanov_mapped_plot = mapscvi::plot_query_labels(query_seura_object=romanovpredicted_seurat,
                                                 reference_seurat=hypoMap_v2_seurat,
                                                 label_col="C25_named",
                                                 label_col_query = "C25_predicted",
                                                 overlay = TRUE,
                                                 bg_col = "grey90",
                                                 query_pt_size = 0.6,
                                                 labelonplot = TRUE,
                                                 cols_plot = getOkabeItoPalette(25),
                                                 label.size=5,
                                                 repel=TRUE)+ggtitle("Romanov et al. mapped on HypoMap")

romanov_mapped_plot_R = rasterize_ggplot(romanov_mapped_plot,pixel_raster = rasterize_pixels,pointsize = 3.4)
romanov_mapped_plot_R


ggsave(filename = paste0(results_path_figure6,"romanov_mapped_C25.png"),
       plot = romanov_mapped_plot_R, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_figure6,"romanov_mapped_C25.pdf"),
       plot = romanov_mapped_plot_R, "pdf",dpi=600,width=300,height = 300,units="mm")

## clustering prob
# romanovpredicted_seurat@meta.data$dummy=NA
# rownames(romanovpredicted_seurat@reductions$umap_scvi@cell.embeddings) = romanovpredicted_seurat@meta.data$Cell_ID
# romanov_prob_plot = FeaturePlot(romanovpredicted_seurat,features = c("C25_named_predicted_prob"), combine = TRUE,raster = TRUE,order=F,reduction = "umap_scvi",
#                                 cols = cols_for_feature_plot,raster.dpi = c(rasterize_px,rasterize_px),pt.size = 3.4)+
#   scale_color_gradient(low=cols_for_feature_plot[1],high = cols_for_feature_plot[2],limits =c(0,1))+ggtitle("")+NoAxes()
# romanov_prob_plot


plot_data = cbind(hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings,color_col = NA)
plot_data2 = cbind(romanovpredicted_seurat@reductions$umap_scvi@cell.embeddings,color_col = romanovpredicted_seurat@meta.data$C25_named_predicted_prob)
plot_data = as.data.frame(rbind(plot_data,plot_data2))
romanov_prob_plot = ggplot2::ggplot(plot_data,aes(x = umapscvi_1, y = umapscvi_2,color=color_col))+
  #ggplot2::scale_fill_manual(values = cols_for_feature_plot  ,na.value = bg_col)+
  ggplot2::scale_color_gradient(low = "white",high = cols_for_feature_plot[2] ,na.value = bg_col,name = "Probability")+
  scattermore::geom_scattermore(pixels=c(2048,2048),pointsize=3)+
  cowplot::theme_cowplot()+NoAxes() #+ guides(color=guide_legend(title="Probability"),)

#p1r=rasterize_ggplot(p1,pixel_raster = 2048,pixel_raster_y = 2048,pointsize = 2.4)
#p1r
romanov_prob_plot

# save
ggsave(filename = paste0(results_path_figure6,"romanov_mapped_prob.png"),
       plot = romanov_prob_plot, "png",dpi=600,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_figure6,"romanov_mapped_prob.pdf"),
       plot = romanov_prob_plot, "pdf",dpi=600,width=300,height = 300,units="mm")

# source

# need to build it for this one
umap_ref = hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings %>% as.data.frame()
umap_ref$Cell_ID = rownames(umap_ref)
romanov_data = bind_cols(romanovpredicted_seurat@reductions$umap_scvi@cell.embeddings,romanovpredicted_seurat@meta.data %>% dplyr::select(C25_named_predicted,C25_named_predicted_prob))
romanov_data$Cell_ID = rownames(romanov_data)
# source
source_figure6_a_b = bind_rows(romanov_data,umap_ref)
data.table::fwrite(source_figure6_a_b,paste0(results_path_figure6,"source_figure6_ab_romanov.txt"),sep="\t")


##########
### Make bacTRAP Signatures
##########

fc_cut = 0.5
padj_cut = 0.01
all_signatures_bacTRAP= list()

# convert list to ordered named list of foldchanges
convert_df_toVec = function(df){
  vec = as.numeric(df[,1])
  names(vec) = as.character(df[,2])
  vec = sort(vec,decreasing = TRUE)
  return(vec)
}

# TODO: restrict to protein coding ?
# library(biomaRt)
# mart <- useMart(dataset="mmusculus_gene_ensembl",biomart='ensembl',host="nov2020.archive.ensembl.org")
# mouse_biotype = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','gene_biotype'),mart = mart)
# mouse_protein_coding_genes = unique(mouse_biotype$external_gene_name[mouse_biotype$gene_biotype=="protein_coding"])

# run on all files:
bacTRAP_files = list.files("data_inputs/",pattern = "bacTRAP_deseq2")
for(i in 1:length(bacTRAP_files)){
  current_res = data.table::fread(paste0("data_inputs/",bacTRAP_files[i]),data.table = FALSE)
  current_signature = current_res[current_res$log2FoldChange >fc_cut & current_res$padj<padj_cut & !is.na(current_res$external_gene_name),]
  #current_signature = current_signature[current_signature$external_gene_name %in% mouse_protein_coding_genes, ]
  current_name = gsub("bacTRAP_deseq2_|\\.csv","",bacTRAP_files[i])
  all_signatures_bacTRAP[[current_name]] = convert_df_toVec(df = current_signature[,c("log2FoldChange","external_gene_name")]) 
}

##########
### Prepare markers
##########

# first filter
marker_table = hypoMap_v2_seurat@misc$marker_genes_all
marker_table = marker_table[marker_table$p_val_adj<0.001 & marker_table$specificity> 0.5,]
marker_table$specificity[marker_table$specificity>1000] = 1000

# subset to C286
marker_table = marker_table[grepl("C286",marker_table$cluster_id),]

# split into list of dfs
marker_list = base::split(marker_table[,c("avg_log2FC","gene")],f=marker_table[,"cluster_id"]) # or: avg_log2fc or fc_mast
marker_list = lapply(marker_list, convert_df_toVec)

all_genes_for_subset = as.character(unique(marker_table[,"gene"]))

##########
### RBO enrichment
##########

# this code requires the rbo function from the gesper package:
# https://rdrr.io/bioc/gespeR/src/R/gespeR-concordance.R 
# https://www.bioconductor.org/packages/release/bioc/html/gespeR.html 

# Here I am using an adapted version (to avoid installing all gesper dependencies)

# parameters
all_p = c(0.98)# c(0.94,0.96,0.98,0.99)
k="auto"
aggregated_w=0.999
assay="RNA"

signature_list = all_signatures_bacTRAP

bacTRAP_signatures_rbo = list()
for(k in 1:length(signature_list)){
  signature_name = names(signature_list)[k]
  message("Running signature: ",signature_name)
  signature = signature_list[[k]]
  # delete genes that are not in marker table
  message(" Subsetting signature from ",length(signature)," to ",length(which(names(signature) %in% all_genes_for_subset)))
  signature = signature[names(signature) %in% all_genes_for_subset]
  
  for(i in 1:length(all_p)){
    current_p = all_p[i]
    # auto k
    #aggregated_w = 0.999
    aggregated_prob = cumsum(sapply(1:max(length(signature),length(all_genes_for_subset)),function(x,current_p){return((1-current_p)*current_p^(x-1))},current_p=current_p))
    auto_k = which(aggregated_prob==min(aggregated_prob[aggregated_prob>aggregated_w]))
    message(" Cumulative probability of ",aggregated_w," with current_p = ",current_p," reached at k = ",auto_k,".")
    message(" Setting k to ",auto_k)
    k=auto_k
    # run
    signature_rbo_l=lapply(marker_list,FUN=rbo2,t=signature,p=current_p,k=k,side="top", uneven.lengths = TRUE) # rbo2 is based on gesper package and loaded from utility_functions.R
    signature_rbo <- data.frame(rbo=matrix(unlist(signature_rbo_l), nrow=length(signature_rbo_l), byrow=TRUE),stringsAsFactors=FALSE)
    signature_rbo$cluster = names(signature_rbo_l)
    signature_rbo$signature_name = signature_name
    signature_rbo$current_p = current_p
    
    # store
    bacTRAP_signatures_rbo[[paste0(signature_name,"_",current_p)]] =  signature_rbo
  }
}

rbo_result_bacTRAP = do.call(rbind,bacTRAP_signatures_rbo)

#save
output_file_name=paste0(results_path_figure6,"bacTRAP_rbo_enrichment.txt")
data.table::fwrite(rbo_result_bacTRAP,file=output_file_name,sep="\t")

rbo_result_bacTRAP = data.table::fread(paste0(results_path_figure6,"bacTRAP_rbo_enrichment.txt"),data.table = F)

##########
### plot with barplot
##########

min_rbo = 0.047

signature_names = unique(rbo_result_bacTRAP$signature_name)
signature_genes = lapply(stringr::str_split(signature_names,"_"),stringr::str_to_title)
names(signature_genes) = signature_names
signature_genes[["pomc_vglut2"]] = c("Pomc","Slc17a6")

umap_with_rbo_list=list()
highlight_list = list(
  "agrp"=c("C286-176","C286-177","C286-178"),
  "glp1r"=c("C286-81","C286-82","C286-75","C286-76","C286-77","C286-175","C286-130","C286-174","C286-181"),
  "pnoc"=c("C286-183","C286-181","C286-179","C286-159", "C286-158","C286-171","C286-170","C286-122"),
  "pomc_glp1r"=c("C286-75","C286-76","C286-77"),
  "pomc_lepr"=c("C286-75","C286-76","C286-77"),
  "pomc" =c("C286-75","C286-76","C286-77")
)

for(i in 1:length(signature_names)){
  print(i)
  current_name = signature_names[i]
  current_genes = signature_genes[[current_name]]
  current_clusters_highlight = highlight_list[[current_name]]
  print(current_clusters_highlight)
  current_result = rbo_result_bacTRAP[rbo_result_bacTRAP$signature_name == current_name & rbo_result_bacTRAP$rbo > min_rbo,]
  coordinate_centers = get_coordinates(current_result$cluster,label_column = "C286",seurat_object = hypoMap_v2_seurat)
  current_result = dplyr::left_join(current_result,coordinate_centers,by=c("cluster"="label"))
  hypoMap_v2_seurat@meta.data$temp_label = NA
  hypoMap_v2_seurat@meta.data$temp_label[hypoMap_v2_seurat@meta.data$C286 %in% current_clusters_highlight] = hypoMap_v2_seurat@meta.data$C286[hypoMap_v2_seurat@meta.data$C286 %in% current_clusters_highlight]
  Idents(hypoMap_v2_seurat) = "temp_label"
  if(length(current_genes)==1){
    umap_plot =FeaturePlot(hypoMap_v2_seurat,features = current_genes,raster = F,order=TRUE,cols = cols_for_feature_plot,label = TRUE,label.size = 4,repel = TRUE)+
      ggtitle(paste0("Positive cells ",current_name," + bacTRAP enrichment"))+NoAxes()
  }else{
    hypoMap_v2_seurat@meta.data$tmp_score = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = current_genes)
    umap_plot =FeaturePlot(hypoMap_v2_seurat,features = "tmp_score",raster = F,order=TRUE,cols = c(bg_col,"#c96410"),label = TRUE,label.size = 4,repel = TRUE)+
      ggtitle(paste0("Double positive cells ",current_name," + bacTRAP enrichment"))+NoAxes()
  }
  
  umap_with_rbo = barplots_on_umap(scatter_plot = umap_plot,
                                   data_barplot = current_result,
                                   max_height = 1.75,
                                   below_center_y = 0.2,
                                   max_width = 0.875,
                                   max_value_display = NULL,
                                   color = "#009E73", # bluishgreen from okabe ito
                                   color_background = "grey50",
                                   border_color = "black",
                                   alpha_background = 0,
                                   value_col = "rbo",
                                   cluster_col = "cluster",
                                   scatter_1_col = "umapscvi_1",
                                   scatter_2_col= "umapscvi_2")
  save_label = umap_with_rbo$layers[[2]]
  umap_with_rbo$layers[[length(umap_with_rbo$layers)+1]] = save_label
  umap_with_rbo$layers[[2]] = NULL
  
  umap_with_rbo_list[[current_name]] = umap_with_rbo
  
}

scUtils::rasterize_ggplot(umap_with_rbo_list[["pomc_glp1r"]],pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)

# example
## make annotation data frame
anno_df = hypoMap_v2_seurat@misc$annotation_result %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

scUtils::rasterize_ggplot(umap_with_rbo_list[["glp1r"]],pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)
rbo_result_bacTRAP %>% dplyr::filter(signature_name %in% c("glp1r")) %>% dplyr::arrange(desc(rbo)) %>% head(n=10) %>% dplyr::left_join(anno_df,by=c("cluster"="cluster_id")) %>% dplyr::select(signature_name,cluster,cluster_name,rbo)

scUtils::rasterize_ggplot(umap_with_rbo_list[["pnoc"]],pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)
rbo_result_bacTRAP %>% dplyr::filter(signature_name %in% c("pnoc")) %>% dplyr::arrange(desc(rbo)) %>% head(n=10) %>% dplyr::left_join(anno_df,by=c("cluster"="cluster_id")) %>% dplyr::select(signature_name,cluster,cluster_name,rbo)

# scUtils::rasterize_ggplot(umap_with_rbo_list[["pomc_gck"]],pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)
# scUtils::rasterize_ggplot(umap_with_rbo_list[["pomc_vglut2"]],pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)
# p1= FeaturePlot(hypoMap_v2_seurat,features = "Slc17a6",raster = F,order=TRUE,cols = cols_for_feature_plot)+NoLegend()
# scUtils::rasterize_ggplot(p1,pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)
# #  %>% dplyr::top_n(n = 10,wt = rbo)  failed
# rbo_result_bacTRAP %>% dplyr::filter(signature_name %in% c("pomc_vglut2")) %>% dplyr::arrange(desc(rbo)) %>% head(n=10) %>% dplyr::left_join(anno_df,by=c("cluster"="cluster_id")) %>% dplyr::select(signature_name,cluster,cluster_name,rbo)

##########
### overview for Pnoc section of text
##########

# need:
# - rbo scores
# - pct of Pnoc in sc
# - cluster name
# - region anno

pnoc_overview = hypoMap_v2_seurat@meta.data %>% dplyr::select(C286,C286_named,Region_summarized) %>% dplyr::distinct(C286,C286_named,Region_summarized)
pnoc_pcts = scUtils::gene_pct_cluster(hypoMap_v2_seurat,col_name = "C286",genes = "Pnoc")
pnoc_pcts$C286 = rownames(pnoc_pcts)
pnoc_overview = dplyr::left_join(pnoc_overview,pnoc_pcts,by=c("C286"="C286"))
pnoc_overview = dplyr::left_join(pnoc_overview,rbo_result_bacTRAP[rbo_result_bacTRAP$signature_name == "pnoc",c("cluster","rbo")],by=c("C286"="cluster"))

pnoc_overview_arc = pnoc_overview %>% dplyr::filter(Region_summarized == "Arcuate hypothalamic nucleus") %>% dplyr::arrange(desc(rbo))

##########
### save
##########

for(i in 1:length(signature_names)){
  print(i)
  current_name = signature_names[i]
  current_plot = umap_with_rbo_list[[current_name]]
  current_plot_r = scUtils::rasterize_ggplot(current_plot,pixel_raster =rasterize_pixels,pointsize = rasterize_point_size)
  # and save:
  ggsave(filename = paste0(results_path_figure6,current_name,"gene_and_rbo_barplot.png"),
         plot = current_plot_r, "png",dpi=600,width=330,height = 300,units="mm")
  ggsave(filename = paste0(results_path_figure6,current_name,"gene_and_rbo_barplot.pdf"),
         plot = current_plot_r, "pdf",dpi=600,width=330,height = 300,units="mm")
  
}

# source data 1:
# bactRAP results
source_figure6_rbores = rbo_result_bacTRAP %>% dplyr::select(-current_p) %>% dplyr::filter(signature_name %in% c("agrp","pomc","pomc_glp1r","pomc_lepr")) %>% tidyr::spread(key=signature_name,value=rbo)
data.table::fwrite(source_figure6_rbores,paste0(results_path_figure6,"source_figure6_cf_rbo_result.txt"),sep="\t")

# source data 2:
# umap coords and gene expression (just re run)
source_figure6_umap_data = hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings  %>% as.data.frame()  %>% dplyr::mutate(Cell_ID = rownames(hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings))
source_figure6_umap_data$Agrp = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = "Agrp")
source_figure6_umap_data$Pomc = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = "Pomc")
source_figure6_umap_data$Pomc_Glp1r = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = c("Glp1r","Pomc"))
source_figure6_umap_data$Pomc_Lepr = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = c("Lepr","Pomc"))

#source_figure6_umap_data = bind_cols(source_figure6_umap_data,FetchData(hypoMap_v2_seurat,c("Agrp","Pomc"))) 
data.table::fwrite(source_figure6_umap_data,paste0(results_path_figure6,"source_figure6_cf_umap_data.txt"),sep="\t")


##########
### source data UMAPs for Figure 7 (see also Figure 7 script)
##########

results_path_figure7 = "figure_outputs/figure_7/"

# source data 1:
# bactRAP results
source_figure7_rbores = rbo_result_bacTRAP %>% dplyr::select(-current_p)  %>% dplyr::filter(signature_name %in% c("glp1r","pnoc")) %>% tidyr::spread(key=signature_name,value=rbo)
data.table::fwrite(source_figure7_rbores,paste0(results_path_figure7,"source_figure7_ad_rbo_result.txt"),sep="\t")

# source data 2:
# umap coords and gene expression (just re run)
source_figure7_umap_data = hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings  %>% as.data.frame()  %>% dplyr::mutate(Cell_ID = rownames(hypoMap_v2_seurat@reductions$umap_scvi@cell.embeddings))
source_figure7_umap_data$Glp1r = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = "Glp1r")
source_figure7_umap_data$Pnoc = scUtils::CalculateMultScore(seurat_object = hypoMap_v2_seurat,features = "Pnoc")

data.table::fwrite(source_figure7_umap_data,paste0(results_path_figure7,"source_figure7_ad_umap_data.txt"),sep="\t")
