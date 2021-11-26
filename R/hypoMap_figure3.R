
##########
### Load & Prepare
##########

#set path
results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_3/"
system(paste0("mkdir -p ",results_path))

# load everything required
source("load_data.R")
require(mapscvi)

# subsample ids for plotting
subsample_ids = data.table::fread(paste0(data_path,"_subsampled_Cell_IDs_neuronMap.txt"),data.table = F,header = F)[,1]

######### TODO: 
# TODO: Currently loading the deseq2 results directly! Move into paper results folder as well !
# TODO: discuss this with Paul regarding geo etc.
# TODO loads rbo2 finction from scHarmonize
source("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scHarmonize/harmonization/gesper_rbo_functions.R")

##########
### Make bacTRAP Signatures
##########

fc_cut = 0.5
padj_cut = 0.01
all_signatures_bacTRAP= list()

#corinna glp1r:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2021-04-corinna-traps/release/dereportr/hypo/deseq_diff/deseq2_diff.csv"
corinna_deseq2_diff_bac_hypo = as.data.frame(fread(paste0(datafile),data.table = F))
bacTRAP_glp1r = corinna_deseq2_diff_bac_hypo[corinna_deseq2_diff_bac_hypo$log2FoldChange>fc_cut & corinna_deseq2_diff_bac_hypo$padj<padj_cut & !is.na(corinna_deseq2_diff_bac_hypo$external_gene_name),]
all_signatures_bacTRAP[["bacTRAP_glp1r"]] = convert_df_toVec(df = bacTRAP_glp1r[,c("log2FoldChange","external_gene_name")])

#alex pnoc:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2019-02-alex-trap-nfcore/release/trapdiff/bactrap/de.rds"
alex_deseq2_diff_bac_pnoc = as.data.frame(readRDS(paste0(datafile)))
bacTRAP_pnoc = alex_deseq2_diff_bac_pnoc[alex_deseq2_diff_bac_pnoc$log2FoldChange>fc_cut & alex_deseq2_diff_bac_pnoc$padj<padj_cut & alex_deseq2_diff_bac_pnoc$comparison %in% c("ip_input_cd")  & !is.na(alex_deseq2_diff_bac_pnoc$external_gene_name) & !is.na(alex_deseq2_diff_bac_pnoc$padj) ,]
all_signatures_bacTRAP[["bacTRAP_pnoc_cd"]] = convert_df_toVec(df = bacTRAP_pnoc[,c("log2FoldChange","external_gene_name")])

#nasim pomc lepr:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-02-nasim-bactrap/release/deseq/de/lepr_deseq2_diff.xlsx"
nasim_deseq2_diff_bac_pomc_lepr = as.data.frame(readxl::read_xlsx(paste0(datafile),sheet = 1))
bacTRAP_pomc_lepr = nasim_deseq2_diff_bac_pomc_lepr[nasim_deseq2_diff_bac_pomc_lepr$log2FoldChange>fc_cut & nasim_deseq2_diff_bac_pomc_lepr$padj<padj_cut & !is.na(nasim_deseq2_diff_bac_pomc_lepr$external_gene_name) & !is.na(nasim_deseq2_diff_bac_pomc_lepr$padj) & nasim_deseq2_diff_bac_pomc_lepr$gene_biotype == "protein_coding",]
all_signatures_bacTRAP[["bacTRAP_pomc_lepr"]] = convert_df_toVec(df = bacTRAP_pomc_lepr[,c("log2FoldChange","external_gene_name")])

#nasim pomc glp1r:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-02-nasim-bactrap/release/deseq/de/glp1r_deseq2_diff.xlsx"
nasim_deseq2_diff_bac_pomc_glp1r =  as.data.frame(readxl::read_xlsx(paste0(datafile),sheet = 1))
bacTRAP_pomc_glp1r = nasim_deseq2_diff_bac_pomc_glp1r[nasim_deseq2_diff_bac_pomc_glp1r$log2FoldChange>fc_cut & nasim_deseq2_diff_bac_pomc_glp1r$padj<padj_cut & !is.na(nasim_deseq2_diff_bac_pomc_glp1r$external_gene_name) & !is.na(nasim_deseq2_diff_bac_pomc_glp1r$padj) & nasim_deseq2_diff_bac_pomc_glp1r$gene_biotype == "protein_coding",]
all_signatures_bacTRAP[["bacTRAP_pomc_glp1r"]] = convert_df_toVec(df = bacTRAP_pomc_glp1r[,c("log2FoldChange","external_gene_name")])

#nasim pomc vglut:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-02-nasim-bactrap/release/deseq/de/vglut_deseq2_diff.xlsx"
nasim_deseq2_diff_bac_pomc_vglut = as.data.frame(readxl::read_xlsx(paste0(datafile),sheet = 1))
bacTRAP_pomc_vglut= nasim_deseq2_diff_bac_pomc_vglut[nasim_deseq2_diff_bac_pomc_vglut$log2FoldChange>fc_cut & nasim_deseq2_diff_bac_pomc_vglut$padj<padj_cut  & !is.na(nasim_deseq2_diff_bac_pomc_vglut$external_gene_name) & !is.na(nasim_deseq2_diff_bac_pomc_vglut$padj) & nasim_deseq2_diff_bac_pomc_vglut$gene_biotype == "protein_coding",]
all_signatures_bacTRAP[["bacTRAP_pomc_vglut"]] = convert_df_toVec(df = bacTRAP_pomc_vglut[,c("log2FoldChange","external_gene_name")])

# #nasim pomc vglut:
# datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-02-nasim-bactrap/release/deseq/de/vglut_deseq2_diff.xlsx"
# nasim_deseq2_diff_bac_pomc_vglut = as.data.frame(readxl::read_xlsx(paste0(datafile),sheet = 1))
# bacTRAP_pomc_vglut= nasim_deseq2_diff_bac_pomc_vglut[nasim_deseq2_diff_bac_pomc_vglut$log2FoldChange>fc_cut & nasim_deseq2_diff_bac_pomc_vglut$padj<padj_cut  & !is.na(nasim_deseq2_diff_bac_pomc_vglut$external_gene_name) & !is.na(nasim_deseq2_diff_bac_pomc_vglut$padj) & nasim_deseq2_diff_bac_pomc_vglut$gene_biotype == "protein_coding",]
# all_signatures_bacTRAP[["bacTRAP_pomc_vglut"]] = convert_df_toVec(df = bacTRAP_pomc_vglut[,c("log2FoldChange","external_gene_name")])

#nasim pomc gck:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-02-nasim-bactrap/release/deseq/de/gck_deseq2_diff.xlsx"
nasim_deseq2_diff_bac_pomc_gck = as.data.frame(readxl::read_xlsx(paste0(datafile),sheet = 1))
bacTRAP_pomc_gck= nasim_deseq2_diff_bac_pomc_gck[nasim_deseq2_diff_bac_pomc_gck$log2FoldChange>fc_cut & nasim_deseq2_diff_bac_pomc_gck$padj<padj_cut  & !is.na(nasim_deseq2_diff_bac_pomc_gck$external_gene_name) & !is.na(nasim_deseq2_diff_bac_pomc_gck$padj) & nasim_deseq2_diff_bac_pomc_gck$gene_biotype == "protein_coding",]
all_signatures_bacTRAP[["bacTRAP_pomc_gck"]] = convert_df_toVec(df = bacTRAP_pomc_gck[,c("log2FoldChange","external_gene_name")])

#kasia agrp:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-05-kasia-bactrap/release/trapdiff/fasted_vs_control/de.rds"
kasia_bac_agrp=  as.data.frame(readRDS(paste0(datafile)))
bacTRAP_agrp_ctrl= kasia_bac_agrp[kasia_bac_agrp$log2FoldChange>fc_cut & kasia_bac_agrp$padj<padj_cut & kasia_bac_agrp$comparison =="ip_input_control"  & !is.na(kasia_bac_agrp$external_gene_name) & !is.na(kasia_bac_agrp$padj),]
all_signatures_bacTRAP[["bacTRAP_agrp_ctrl"]] = convert_df_toVec(df = bacTRAP_agrp_ctrl[,c("log2FoldChange","external_gene_name")])

#alai Pomc agrp:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2021-05-alai-bactrap/release/trapdiff/dre_cre/de.rds"
alai_bac_pomc =  as.data.frame(readRDS(paste0(datafile)))
bacTRAP_Pomc = alai_bac_pomc[alai_bac_pomc$log2FoldChange>fc_cut & alai_bac_pomc$padj<padj_cut & alai_bac_pomc$comparison =="ip_input_cre"  & !is.na(alai_bac_pomc$external_gene_name) & !is.na(alai_bac_pomc$padj),]
all_signatures_bacTRAP[["bacTRAP_pomc"]] = convert_df_toVec(df = bacTRAP_Pomc[,c("log2FoldChange","external_gene_name")])

sapply(all_signatures_bacTRAP,length)
# save as backup
saveRDS(all_signatures_bacTRAP,paste0(data_path,"bacTRAP_signatures_rbo.rds"))


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
marker_list = base::split(marker_table[,c("specificity","gene")],f=marker_table[,"cluster_1"]) # or: avg_log2fc or fc_mast
# convert list to ordered named list of foldchanges
convert_df_toVec = function(df){
  vec = as.numeric(df[,1])
  names(vec) = as.character(df[,2])
  vec = sort(vec,decreasing = TRUE)
  return(vec)
}
marker_list = lapply(marker_list, convert_df_toVec)

all_genes_for_subset = as.character(unique(marker_table[,"gene"]))

##########
### RBO enrichment
##########

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
    signature_rbo_l=lapply(marker_list,FUN=rbo2,t=signature,p=current_p,k=k,side="top", uneven.lengths = TRUE)
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
# data.table::fwrite(rbo_result_bacTRAP,file=output_file_name)

##########
### plotting helper function
##########

# This function creates ggplots similar to Suerat's Feature plot to plot the rbo score per cell

plot_rbo = function(rbo_result,map_seurat,clusterlevel,colorvec,center_cutoff =NULL,label_col_name = "cluster_name",label_size=5,remove_grep_text = "",
                    relevant_clusters=NULL,text_color="black",text_size=5,plot_max=NULL,point_size=0.2,nudge_x=0,rasterize_plot=FALSE,rasterize_pixels = 1536,rasterize_point_size = 1.5){
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
    scale_color_gradientn(colours = colorvec,limits = c(0, plot_max),na.value = colorvec[1],oob=squish)#+ # +ggtitle("RBO per cluster")
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
  # rasterize
  if(rasterize_plot){
    p = rasterize_ggplot(p,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
  }
  #output
  p
}

##########
### RBO results plots
##########

# need annotations
anno_df = neuron_map_seurat@misc$annotations

### downsample for plotting
neuron_map_seurat_downsampled =  subset(neuron_map_seurat,cells = subsample_ids)

require(RColorBrewer)

# graph params
clusterlevel = "K169"
colorvec = RColorBrewer::brewer.pal(9, "Blues")
colorvec[1] = "#dedede"
text_color="black"
text_size = 20
label_size = 7
point_size = 0.3
# others
max_scale= 0.3
nudge_x=4.5
nudge_x_klabel=2.5
current_p=0.98
remove_grep_text = ""

all_signature_names = c("bacTRAP_agrp_ctrl","bacTRAP_pomc" ,"bacTRAP_pomc_lepr","bacTRAP_pomc_glp1r","bacTRAP_pnoc_cd","bacTRAP_glp1r")# unique(rbo_result_bacTRAP$signature_name)
center_cutoff_vector = c(0.2,0.1,0.2,0.2,0.06,0.06)
all_signature_plots = list()
all_signature_plots_unlabelled = list()
for(i in 1:length(all_signature_names)){
  current_signature_name = all_signature_names[i]
  center_cutoff = center_cutoff_vector[i] # has to be the same order as all_signature_names, then this is used to select a cutoff for labels
  rbo_result_current_sig = rbo_result_bacTRAP %>% dplyr::filter(signature_name==current_signature_name & current_p==current_p & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
    dplyr::arrange(desc(rbo))
  rbo_result_current_sig=left_join(rbo_result_current_sig,anno_df[,c("cluster_id","cluster_name","ncells")],by=c("cluster"="cluster_id"))
  # make plot w/ labels
  all_signature_plots[[current_signature_name]] =plot_rbo(rbo_result_current_sig,neuron_map_seurat_downsampled,clusterlevel,colorvec,center_cutoff=center_cutoff,
                                                          label_col_name="cluster_name",label_size=label_size,remove_grep_text=remove_grep_text,
                                                          text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=5,plot_max = max_scale)
  # make plot w/o labels
  all_signature_plots_unlabelled[[current_signature_name]] =plot_rbo(rbo_result_current_sig,neuron_map_seurat_downsampled,clusterlevel,colorvec,center_cutoff=NULL,
                                                                     label_size=label_size,remove_grep_text=remove_grep_text,text_color=text_color,
                                                                     text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
}

names(all_signature_plots)
all_signature_plots$bacTRAP_pnoc_cd

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
### Old:
##########


## pomc lepr
rbo_result_pomc_lepr = rbo_result_bacTRAP %>% dplyr::filter(signature_name=="bacTRAP_pomc_lepr" & current_p==0.98 & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
  dplyr::arrange(desc(rbo))
rbo_result_pomc_lepr=left_join(rbo_result_pomc_lepr,anno_df[,c("cluster_id","Map_CellType","ncells")],by=c("cluster"="cluster_id"))
p_pomc_lepr=plot_rbo(rbo_result_pomc_lepr,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.2,label_col_name="cluster_name",label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                     text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=5,plot_max = max_scale)
p_pomc_lepr

# and save:
ggsave(filename = paste0(results_path,"pomc_lepr_rbo_plot.png"),
       plot = p_pomc_lepr, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_lepr_rbo_plot.pdf"),
       plot = p_pomc_lepr, "pdf",dpi=600,width=330,height = 300,units="mm")
#without label
p_pomc_lepr_unlabelled=plot_rbo(rbo_result_pomc_lepr,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=NULL,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                                text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
p_pomc_lepr_unlabelled
ggsave(filename = paste0(results_path,"pomc_lepr_rbo_nolabel_plot.png"),
       plot = p_pomc_lepr_unlabelled, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_lepr_rbo_nolabel_plot.pdf"),
       plot = p_pomc_lepr_unlabelled, "pdf",dpi=600,width=330,height = 300,units="mm")
#with K label
p_pomc_lepr_klabel =plot_rbo(rbo_result_pomc_lepr,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.2,label_col_name="cluster",label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                             text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x_klabel,plot_max = max_scale)
p_pomc_lepr_klabel
ggsave(filename = paste0(results_path,"pomc_lepr_rbo_klabel_plot.png"),
       plot = p_pomc_lepr_klabel, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_lepr_rbo_klabel_plot.pdf"),
       plot = p_pomc_lepr_klabel, "pdf",dpi=600,width=330,height = 300,units="mm")


## pomc glp1r
rbo_result_pomc_glp1r = rbo_result_bacTRAP %>% dplyr::filter(signature_name=="bacTRAP_pomc_glp1r" & current_p==0.98 & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
  dplyr::arrange(desc(rbo))
rbo_result_pomc_glp1r=left_join(rbo_result_pomc_glp1r,anno_df[,c("cluster_id","Map_CellType","ncells")],by=c("cluster"="cluster_id"))
p_pomc_glp1r=plot_rbo(rbo_result_pomc_glp1r,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.2,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                      text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=4.5,plot_max = max_scale)
p_pomc_glp1r
# and save:
ggsave(filename = paste0(results_path,"pomc_glp1r_rbo_plot.png"),
       plot = p_pomc_glp1r, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_glp1r_rbo_plot.pdf"),
       plot = p_pomc_glp1r, "pdf",dpi=600,width=330,height = 300,units="mm")
#without label
p_pomc_glp1r_unlabelled=plot_rbo(rbo_result_pomc_glp1r,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=NULL,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                                 text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
p_pomc_glp1r_unlabelled
ggsave(filename = paste0(results_path,"pomc_glp1r_rbo_nolabel_plot.png"),
       plot = p_pomc_glp1r_unlabelled, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_glp1r_rbo_nolabel_plot.pdf"),
       plot = p_pomc_glp1r_unlabelled, "pdf",dpi=600,width=330,height = 300,units="mm")
#with K label
p_pomc_glp1r_klabel =plot_rbo(rbo_result_pomc_glp1r,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.2,label_col_name="cluster",label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                              text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x_klabel,plot_max = max_scale)
p_pomc_glp1r_klabel
ggsave(filename = paste0(results_path,"pomc_glp1r_rbo_klabel_plot.png"),
       plot = p_pomc_glp1r_klabel, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_glp1r_rbo_klabel_plot.pdf"),
       plot = p_pomc_glp1r_klabel, "pdf",dpi=600,width=330,height = 300,units="mm")


## glp1r
arh_clusters = c("K169-32","K169-82","K169-116","K169-17","K169-39","K169-86","K169-58","K169-89","K169-149","K169-18","K169-1","K169-68","K169-75","K169-12","K169-44","K169-105","K169-70","K169-60","K169-154","K169-21","K169-35","K169-125","K169-130","K169-28","K169-24","K169-72","K169-48","K169-50","K169-130","K169-172","K169-6","K169-45")
relevant_clusters_glp1r = anno_df$Map_CellType[anno_df$cluster_id %in%  c("K169-40","K169-144","K169-24","K169-128",arh_clusters)]
rbo_result_glp1r = rbo_result_bacTRAP %>% dplyr::filter(signature_name=="bacTRAP_glp1r" & current_p==0.98 & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
  dplyr::arrange(desc(rbo))
rbo_result_glp1r=left_join(rbo_result_glp1r,anno_df[,c("cluster_id","Map_CellType","ncells")],by=c("cluster"="cluster_id"))
p_glp1r=plot_rbo(rbo_result_glp1r,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.06,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                 text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=3.25,plot_max = max_scale,relevant_clusters = relevant_clusters_glp1r)
p_glp1r
# and save:
ggsave(filename = paste0(results_path,"glp1r_rbo_plot.png"),
       plot = p_glp1r, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"glp1r_rbo_plot.pdf"),
       plot = p_glp1r, "pdf",dpi=600,width=330,height = 300,units="mm")
#without label
p_glp1r_unlabelled=plot_rbo(rbo_result_glp1r,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=NULL,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                            text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
p_glp1r_unlabelled
ggsave(filename = paste0(results_path,"glp1r_rbo_nolabel_plot.png"),
       plot = p_glp1r_unlabelled, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"glp1r_rbo_nolabel_plot.pdf"),
       plot = p_glp1r_unlabelled, "pdf",dpi=600,width=330,height = 300,units="mm")
#with K label
p_glp1r_klabel =plot_rbo(rbo_result_glp1r,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.05,label_col_name="cluster",label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                         text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x_klabel,plot_max = max_scale,relevant_clusters = c("K169-40","K169-144","K169-24","K169-128",arh_clusters))
p_glp1r_klabel
ggsave(filename = paste0(results_path,"glp1r_rbo_klabel_plot.png"),
       plot = p_glp1r_klabel, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"glp1r_rbo_klabel_plot.pdf"),
       plot = p_glp1r_klabel, "pdf",dpi=600,width=330,height = 300,units="mm")


## pnoc
arh_clusters = c("K169-32","K169-82","K169-116","K169-17","K169-39","K169-86","K169-58","K169-89","K169-149","K169-18","K169-1","K169-68","K169-75","K169-12","K169-44","K169-105","K169-70","K169-60","K169-154","K169-21","K169-35","K169-125","K169-130","K169-28","K169-24","K169-72","K169-48","K169-50","K169-130","K169-172","K169-6","K169-45")
relevant_clusters_pnoc = anno_df$Map_CellType[anno_df$cluster_id %in%  c("K169-139","K169-107",arh_clusters)]
rbo_result_pnoc= rbo_result_bacTRAP %>% dplyr::filter(signature_name=="bacTRAP_pnoc_cd" & current_p==0.98 & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
  dplyr::arrange(desc(rbo))
rbo_result_pnoc=left_join(rbo_result_pnoc,anno_df[,c("cluster_id","Map_CellType","ncells")],by=c("cluster"="cluster_id"))
p_pnoc=plot_rbo(rbo_result_pnoc,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.04,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=3.25,plot_max = max_scale,relevant_clusters = relevant_clusters_pnoc)
p_pnoc
# and save:
ggsave(filename = paste0(results_path,"pnoc_rbo_plot.png"),
       plot = p_pnoc, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pnoc_rbo_plot.pdf"),
       plot = p_pnoc, "pdf",dpi=600,width=330,height = 300,units="mm")
#without label
p_pnoc_unlabelled=plot_rbo(rbo_result_pnoc,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=NULL,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                           text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
p_pnoc_unlabelled
ggsave(filename = paste0(results_path,"pnoc_rbo_nolabel_plot.png"),
       plot = p_pnoc_unlabelled, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pnoc_rbo_nolabel_plot.pdf"),
       plot = p_pnoc_unlabelled, "pdf",dpi=600,width=330,height = 300,units="mm")
#with K label
p_pnoc_klabel =plot_rbo(rbo_result_pnoc,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.035,label_col_name="cluster",label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                        text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=2,plot_max = max_scale, relevant_clusters = c("K169-139","K169-107",arh_clusters))
p_pnoc_klabel
ggsave(filename = paste0(results_path,"pnoc_rbo_klabel_plot.png"),
       plot = p_pnoc_klabel, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pnoc_rbo_klabel_plot.pdf"),
       plot = p_pnoc_klabel, "pdf",dpi=600,width=330,height = 300,units="mm")

## agrp
rbo_result_agrp= rbo_result_bacTRAP %>% dplyr::filter(signature_name=="bacTRAP_agrp_ctrl" & current_p==0.98 & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
  dplyr::arrange(desc(rbo))
rbo_result_agrp=left_join(rbo_result_agrp,anno_df[,c("cluster_id","Map_CellType","ncells")],by=c("cluster"="cluster_id"))
p_agrp=plot_rbo(rbo_result_agrp,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.2,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=4,plot_max = max_scale)
p_agrp
# and save:
ggsave(filename = paste0(results_path,"agrp_rbo_plot.png"),
       plot = p_agrp, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"agrp_rbo_plot.pdf"),
       plot = p_agrp, "pdf",dpi=600,width=330,height = 300,units="mm")
#without label
p_agrp_unlabelled=plot_rbo(rbo_result_agrp,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=NULL,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                           text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
p_agrp_unlabelled
ggsave(filename = paste0(results_path,"agrp_rbo_nolabel_plot.png"),
       plot = p_agrp_unlabelled, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"agrp_rbo_nolabel_plot.pdf"),
       plot = p_agrp_unlabelled, "pdf",dpi=600,width=330,height = 300,units="mm")
#with K label
p_agrp_klabel =plot_rbo(rbo_result_agrp,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.2,label_col_name="cluster",label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                        text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x_klabel,plot_max = max_scale)
p_agrp_klabel
ggsave(filename = paste0(results_path,"agrp_rbo_klabel_plot.png"),
       plot = p_agrp_klabel, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"agrp_rbo_klabel_plot.pdf"),
       plot = p_agrp_klabel, "pdf",dpi=600,width=330,height = 300,units="mm")

## pomc
rbo_result_pomc= rbo_result_bacTRAP %>% dplyr::filter(signature_name=="bacTRAP_pomc" & current_p==0.98 & grepl(clusterlevel,cluster)) %>% # & rbo > 0.03 
  dplyr::arrange(desc(rbo))
rbo_result_pomc=left_join(rbo_result_pomc,anno_df[,c("cluster_id","Map_CellType","ncells")],by=c("cluster"="cluster_id"))
p_pomc=plot_rbo(rbo_result_pomc,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.1,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=3.75,plot_max = max_scale)
p_pomc
# and save:
ggsave(filename = paste0(results_path,"pomc_rbo_plot.png"),
       plot = p_pomc, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_rbo_plot.pdf"),
       plot = p_pomc, "pdf",dpi=600,width=330,height = 300,units="mm")
#without label
p_pomc_unlabelled=plot_rbo(rbo_result_pomc,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=NULL,label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                           text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x,plot_max = max_scale)
p_pomc_unlabelled
ggsave(filename = paste0(results_path,"pomc_rbo_nolabel_plot.png"),
       plot = p_pomc_unlabelled, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_rbo_nolabel_plot.pdf"),
       plot = p_pomc_unlabelled, "pdf",dpi=600,width=330,height = 300,units="mm")
#with K label
p_pomc_klabel =plot_rbo(rbo_result_pomc,neuron_map_seurat,clusterlevel,colorvec,center_cutoff=0.15,label_col_name="cluster",label_size=label_size,remove_grep_text="Slc32a1\\.|Slc17a6\\.",
                        text_color=text_color,text_size=text_size,point_size = point_size,nudge_x=nudge_x_klabel,plot_max = max_scale)
p_pomc_klabel
ggsave(filename = paste0(results_path,"pomc_rbo_klabel_plot.png"),
       plot = p_pomc_klabel, "png",dpi=600,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path,"pomc_rbo_klabel_plot.pdf"),
       plot = p_pomc_klabel, "pdf",dpi=600,width=330,height = 300,units="mm")
### and save with ggsave pdf
# ggsave(filename = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/analysis_results/cers/kimVMH_cers1to6_genes.pdf",
#        plot = p_all_kimVMH, "pdf",dpi=300,width=300,height = 200,units="mm")

##########
### Pnoc
##########

pnoc_stats = cbind(neuron_map_seurat@meta.data[,c("Cell_ID","Dataset","K169_named")] ,FetchData(neuron_map_seurat,vars = "Pnoc")) %>% group_by(K169_named) %>%
  dplyr::summarise(mean_pnoc= mean(Pnoc) )

rbo_result_pnoc_with_expression = dplyr::left_join(rbo_result_pnoc,pnoc_stats,by=c("clean_names"="K169_named"))

cellshigh=neuron_map_seurat@meta.data$Cell_ID[neuron_map_seurat@meta.data$K169_named=="Slc32a1.Six6.Arx.Penk"]
DimPlot(neuron_map_seurat,group.by = "K329_named",reduction = paste0("umap_","scvi"),label = F,label.size = 3,repel = TRUE,cells.highlight = cellshigh,sizes.highlight = 0.1)+NoAxes()+NoLegend()

FeaturePlot(neuron_map_seurat,features = c("Trh"),reduction = paste0("umap_","scvi"),order = TRUE)+NoAxes()#+NoLegend()

table(neuron_map_seurat@meta.data$Dataset[neuron_map_seurat@meta.data$K98_named == "Slc17a6.Nrn1.Sim1.Trh"])

#bacTRAp_signatures$bacTRAP_pnoc_cd[names(bacTRAp_signatures$bacTRAP_pnoc_cd) %in% markers_comparisons_all$gene[markers_comparisons_all$cluster_1=="K169-161" & markers_comparisons_all$specificity>5]]

##########
### EXPORT
##########

list_rbo_plots = list( p_agrp = p_agrp  , p_pomc = p_pomc, p_pomc_lepr = p_pomc_lepr, p_pomc_glp1r= p_pomc_glp1r,p_pnoc = p_pnoc, p_glp1r  = p_glp1r)

list_rbo_plots$p_agrp

saveRDS(list_rbo_plots,paste0(results_path,"Figure_3_july/","Figure_3_plots.rds"))

