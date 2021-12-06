
##########
### Load & Prepare
##########

#set path
results_path = "figure_outputs/figure_3/"
system(paste0("mkdir -p ",results_path))

# load required functions
require(mapscvi)
source("R/utility_functions.R")
source("R/plot_functions.R")

# where to find large data objects (need to change to local dir !)
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

# load seurat objects via large_data_path
load_required_files(large_data_path = large_data_path)

######### TODO: 
# TODO: discuss this with Paul regarding geo etc. and whether I am using all the right results

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
library(biomaRt)
mart <- useMart(dataset="mmusculus_gene_ensembl",biomart='ensembl',host="nov2020.archive.ensembl.org")
mouse_biotype = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','gene_biotype'),mart = mart)
mouse_protein_coding_genes = unique(mouse_biotype$external_gene_name[mouse_biotype$gene_biotype=="protein_coding"])

# run on all files:
bacTRAP_files = list.files("data_inputs/",pattern = "bacTRAP_deseq2")
for(i in 1:length(bacTRAP_files)){
  current_res = data.table::fread(paste0("data_inputs/",bacTRAP_files[i]),data.table = FALSE)
  current_signature = current_res[current_res$log2FoldChange >fc_cut & current_res$padj<padj_cut & !is.na(current_res$external_gene_name),]
  current_signature = current_signature[current_signature$external_gene_name %in% mouse_protein_coding_genes, ]
  current_name = gsub("bacTRAP_deseq2_|\\.csv","",bacTRAP_files[i])
  all_signatures_bacTRAP[[current_name]] = convert_df_toVec(df = current_signature[,c("log2FoldChange","external_gene_name")]) 
}

##########
### Prepare markers
##########

# first filter
marker_table = neuron_map_seurat@misc$markers_comparisons_all
marker_table = marker_table[marker_table$p_val_adj<0.001 & marker_table$specificity> 0.5,]
marker_table$specificity[marker_table$specificity>1000] = 1000

# subset to K169_pruned
marker_table = marker_table[grepl("K169",marker_table$cluster_1),]

# split into list of dfs
marker_list = base::split(marker_table[,c("avg_logFC","gene")],f=marker_table[,"cluster_1"]) # or: avg_log2fc or fc_mast
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
output_file_name=paste0(results_path,"bacTRAP_rbo_enrichment.txt")
data.table::fwrite(rbo_result_bacTRAP,file=output_file_name)

##########
### plotting helper function
##########

# This function creates ggplots similar to Suerat's Feature plot to plot the rbo score per cell

plot_rbo = function(rbo_result,map_seurat,clusterlevel,colorvec,na_color,center_cutoff =NULL,label_col_name = "cluster_name",label_size=5,remove_grep_text = "",
                    relevant_clusters=NULL,text_color="black",text_size=5,plot_max=NULL,point_size=0.2,nudge_x=0){
  require(ggplot2)
  require(scales)
  # gett coordinates and named column from seurat metadata 
  plot_cluster = cbind(map_seurat@meta.data[,c("Cell_ID",paste0(clusterlevel,"_named"))],as.data.frame(map_seurat@reductions$umap_scvi@cell.embeddings))
  colnames(plot_cluster)[2] = "cluster_name"
  # add the rbo result with the pruned id and rbo values
  plot_cluster=left_join(plot_cluster,rbo_result,by=c("cluster_name"="cluster_name"))
  # get idx of column that will be used to name cluster on the plot
  col_labels = colnames(plot_cluster)
  label_idx = which(col_labels == label_col_name)
  if(is.null(plot_max)){plot_max = max(0.1,max(plot_cluster[,5],na.rm = TRUE))}
  # arrange plot_cluster by rbo
  plot_cluster = plot_cluster %>% dplyr::arrange(rbo)
  # make initial plot without labels
  p=ggplot(plot_cluster,aes_string(x=col_labels[3],y=col_labels[4],color=col_labels[5]))+geom_point(size=point_size)+ 
    scale_color_gradient(low = colorvec[1],high = colorvec[2],limits = c(0, plot_max),na.value = na_color,oob=squish)#+ # +ggtitle("RBO per cluster")
  #scale_fill_gradientn(colours = colorvec,limits = c(0, plot_max),na.value = colorvec[1],oob=squish)
  # optionally add labels 
  if(!is.null(center_cutoff)){
    # calculate centers
    label_centers = plot_cluster %>% dplyr::group_by(!!sym(label_col_name)) %>% dplyr::summarise_at(vars(matches("umap")), mean)
    label_centers=as.data.frame(label_centers)
    # if not prespecified clusters are given: use all
    if(is.null(relevant_clusters)){relevant_clusters=label_centers[,label_col_name]}
    # mark cluster with a low score
    score_to_low_clusters = unique(plot_cluster[plot_cluster$rbo<center_cutoff ,label_col_name])
    label_centers[label_centers[,label_col_name] %in% score_to_low_clusters,label_col_name] = NA
    label_centers[!label_centers[,label_col_name] %in% unique(relevant_clusters),label_col_name] = NA
    label_centers[grepl("problematic",label_centers[,label_col_name]),label_col_name] = NA
    label_centers[,label_col_name] = gsub(remove_grep_text,"",label_centers[,label_col_name] )
    # update plot
    p = p + ggrepel::geom_text_repel(data=label_centers,mapping=aes_string(x=col_labels[3],y=col_labels[4],label=label_col_name),
                                     inherit.aes = F , size=label_size,color=text_color,min.segment.length = 0, seed = 42,nudge_x = nudge_x)#,
  }
  p = p + theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"))+NoAxes()+ggtitle("RBO per cluster")#,axis.line=element_line(color="grey40",size = 1)
  #output
  p
}

##########
### RBO results plots
##########

# need annotations
anno_df = neuron_map_seurat@misc$annotations

require(RColorBrewer)
# graph params
clusterlevel = "K169"
colorvec = RColorBrewer::brewer.pal(9, "Blues")
colorvec[1] =  "#dedede"
cols_for_feature_plot = c("#dedede","#0b3ebd") # "#0b3ebd"
na_color = "#dedede"
text_color="black"
text_size = 15
label_size = 7
point_size = 0.3
# others
max_scale= 0.3
nudge_x=4
nudge_x_klabel=2.5
current_p=0.98
remove_grep_text = ""
# raster
rasterize_pixels = 1536
rasterize_point_size = 1.5

all_signature_names =  unique(rbo_result_bacTRAP$signature_name)#c("bacTRAP_agrp_ctrl","bacTRAP_pomc" ,"bacTRAP_pomc_lepr","bacTRAP_pomc_glp1r","bacTRAP_pnoc_cd","bacTRAP_glp1r")# unique(rbo_result_bacTRAP$signature_name)
all_signature_names
center_cutoff_vector = c(0.15,0.06,0.04,0.15,0.15,0.15)
all_signature_plots = list()
all_signature_plots_unlabelled = list()
for(i in 1:length(all_signature_names)){
  current_signature_name = all_signature_names[i]
  center_cutoff = center_cutoff_vector[i] # has to be the same order as all_signature_names, then this is used to select a cutoff for labels
  rbo_result_current_sig = rbo_result_bacTRAP %>% dplyr::filter(signature_name==current_signature_name & current_p==current_p & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
    dplyr::arrange(desc(rbo))
  rbo_result_current_sig=left_join(rbo_result_current_sig,anno_df[,c("cluster_id","cluster_name","ncells")],by=c("cluster"="cluster_id"))
  # make plot w/ labels
  relevant_clusters = anno_df$cluster_name[anno_df$clusterlevel=="K169"]
  # customize labels for heterogeneous bacTRAPs
  if(current_signature_name=="glp1r"){## glp1r
    arh_clusters = c("K169-32","K169-82","K169-116","K169-17","K169-39","K169-86","K169-58","K169-89","K169-149","K169-18","K169-1","K169-68","K169-75","K169-12","K169-44","K169-105","K169-70","K169-60","K169-154","K169-21","K169-35","K169-125","K169-130","K169-28","K169-24","K169-72","K169-48","K169-50","K169-130","K169-172","K169-6","K169-45")
    relevant_clusters = anno_df$cluster_name[anno_df$cluster_id %in%  c("K169-40","K169-144","K169-24","K169-128",arh_clusters)]
  }
  if(current_signature_name=="pnoc"){
    arh_clusters = c("K169-32","K169-82","K169-116","K169-17","K169-39","K169-86","K169-58","K169-89","K169-149","K169-18","K169-1","K169-68","K169-75","K169-12","K169-44","K169-105","K169-70","K169-60","K169-154","K169-21","K169-35","K169-125","K169-130","K169-28","K169-24","K169-72","K169-48","K169-50","K169-130","K169-172","K169-6","K169-45")
    relevant_clusters = anno_df$cluster_name[anno_df$cluster_id  %in%  c("K169-139","K169-107",arh_clusters)]
  }
  all_signature_plots[[current_signature_name]] =plot_rbo(rbo_result_current_sig,neuron_map_seurat,clusterlevel,cols_for_feature_plot,na_color=na_color,center_cutoff=center_cutoff,
                                                          label_col_name="cluster_name",label_size=label_size,remove_grep_text=remove_grep_text,relevant_clusters=relevant_clusters,
                                                          text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=5,plot_max = max_scale)
  all_signature_plots[[current_signature_name]] = rasterize_ggplot(all_signature_plots[[current_signature_name]],pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
  # make plot w/o labels
  all_signature_plots_unlabelled[[current_signature_name]] =plot_rbo(rbo_result_current_sig,neuron_map_seurat,clusterlevel,cols_for_feature_plot,na_color=na_color,center_cutoff=NULL,
                                                                     label_size=label_size,remove_grep_text=remove_grep_text,text_color=text_color,
                                                                     text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
  all_signature_plots_unlabelled[[current_signature_name]] = rasterize_ggplot(all_signature_plots_unlabelled[[current_signature_name]],pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
  
}

names(all_signature_plots)
all_signature_plots$pnoc

all_signature_plots_unlabelled$glp1r

##########
### save all files
##########

for(i in 1:length(all_signature_names)){
  current_signature_name = all_signature_names[i]
  message(current_signature_name)
  current_plot = all_signature_plots[[current_signature_name]]
  # and save:
  ggsave(filename = paste0(results_path,current_signature_name,"_rbo_plot.png"),
         plot = current_plot, "png",dpi=600,width=330,height = 300,units="mm")
  ggsave(filename = paste0(results_path,current_signature_name,"_rbo_plot.pdf"),
         plot = current_plot, "pdf",dpi=600,width=330,height = 300,units="mm")
  
  current_plot_unlabelled = all_signature_plots_unlabelled[[current_signature_name]]
  # save non_labelled as well
  ggsave(filename = paste0(results_path,current_signature_name,"_rbo_nolabel_plot.png"),
         plot = current_plot_unlabelled, "png",dpi=600,width=330,height = 300,units="mm")
  ggsave(filename = paste0(results_path,current_signature_name,"_rbo_nolabel_plot.pdf"),
         plot = current_plot_unlabelled, "pdf",dpi=600,width=330,height = 300,units="mm")
  
}

##########
### EXPORT
##########
# 
# list_rbo_plots = list( p_agrp = p_agrp  , p_pomc = p_pomc, p_pomc_lepr = p_pomc_lepr, p_pomc_glp1r= p_pomc_glp1r,p_pnoc = p_pnoc, p_glp1r  = p_glp1r)
# 
# list_rbo_plots$p_agrp
# 
# saveRDS(list_rbo_plots,paste0(results_path,"Figure_3_plots.rds"))
# 
