
##########
### Load & Prepare
##########

results_path_figure4 = "figure_outputs/figure_4/"
system(paste0("mkdir -p ",results_path_figure4))

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
### Plot nucseq on rest
##########

hypoMap_v2_seurat@meta.data$Dowsett10xnuc_anno_col=NA
hypoMap_v2_seurat@meta.data$Dowsett10xnuc_anno_col[hypoMap_v2_seurat@meta.data$Dataset=="Dowsett10xnuc"] = as.character(hypoMap_v2_seurat@meta.data$C25_named[hypoMap_v2_seurat@meta.data$Dataset=="Dowsett10xnuc"])
hypoMap_v2_seurat@meta.data$Dowsett10xnuc_anno_col = factor(hypoMap_v2_seurat@meta.data$Dowsett10xnuc_anno_col,levels = unique(hypoMap_v2_seurat@meta.data$Dowsett10xnuc_anno_col))
#hypoMap_v2_seurat@meta.data$Dowsett10xnuc_anno_col = factor(hypoMap_v2_seurat@meta.data$Dowsett10xnuc_anno_col,levels = levels(as.factor(hypoMap_v2_seurat@meta.data$C185)))
# plot
# getLongPalette or getOkabeItoPalette
nucseq_in_hypoMap_plot=DimPlot(hypoMap_v2_seurat,group.by = "Dowsett10xnuc_anno_col",reduction = paste0("umap_","scvi"),pt.size = 0.2,raster = F,cols = getOkabeItoPalette(25),
                               label = TRUE,label.size = 5,repel = TRUE,order = TRUE,na.value = bg_col,raster.dpi = c(rasterize_px,rasterize_px))+
  NoLegend()+NoAxes()+ggtitle("Nucseq clusters")
# change order
plot_data = nucseq_in_hypoMap_plot$data
nucseq_in_hypoMap_plot$data = plot_data[order(plot_data$Dowsett10xnuc_anno_col,na.last = FALSE),]
# i rasterize afterwards so that I can use standard ggplot functions to set color palettes etc
nucseq_in_hypoMap_plot = rasterize_ggplot(nucseq_in_hypoMap_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
nucseq_in_hypoMap_plot

ggsave(filename = paste0(results_path_figure4,"nucseq_in_hypoMap.png"),
       plot = nucseq_in_hypoMap_plot, "png",dpi=450,width=300,height = 300,units="mm")
ggsave(filename = paste0(results_path_figure4,"nucseq_in_hypoMap.pdf"),
       plot = nucseq_in_hypoMap_plot, "pdf",dpi=450,width=300,height = 300,units="mm")

##source data
nucseq_in_hypoMap_data = nucseq_in_hypoMap_plot$data
nucseq_in_hypoMap_data$Cell_ID = rownames(nucseq_in_hypoMap_data)
data.table::fwrite(nucseq_in_hypoMap_data,paste0(results_path_figure4,"source_figure4_a_umap_nucseq.txt"),sep="\t")



###########################
########## Dynamic range of expression 


##########
### Cluster marker genes
##########

# padj_cut = 0.001
# fc_min = 0.25
markers_C185 = hypoMap_v2_seurat@misc$marker_genes_all[hypoMap_v2_seurat@misc$marker_genes_all$p_val_adj < 0.001 & 
                                                         hypoMap_v2_seurat@misc$marker_genes_all$specificity > 0.25 & 
                                                         stringr::str_extract(hypoMap_v2_seurat@misc$marker_genes_all$cluster_id,"C[0-9]+") == "C185",]

# which cluster seem relevant (have a sufficient number of cells mapped in sn-seq data to calculate reliable averages)
min_cells = 30
clusters_to_check = names(table(hypoMap_v2_seurat@meta.data$C185[hypoMap_v2_seurat@meta.data$Dataset=="Dowsett10xnuc"]))[table(hypoMap_v2_seurat@meta.data$C185[hypoMap_v2_seurat@meta.data$Dataset=="Dowsett10xnuc"])>=min_cells]
length(clusters_to_check)

##########
### gene-centric comparison:

##########
### global gene mean calculation
##########

# which clusters to exclude fromsc-seq comparison
not_scseq_datasets = c("Dowsett10xnuc","Affinati10x","Rupp10x")

# build metacell expression matrix
build_cell = function(cells,expr_matrix){
  if(length(cells)>1){
    return(rowMeans(expr_matrix[,expr_matrix@Dimnames[[2]] %in% cells]))
  }else{
    return(expr_matrix[,expr_matrix@Dimnames[[2]] %in% cells])
  }
}
st=Sys.time()
cluster_list_sc_seq = split(hypoMap_v2_seurat@meta.data$Cell_ID[hypoMap_v2_seurat@meta.data$C185 %in% clusters_to_check &  ! hypoMap_v2_seurat@meta.data$Dataset %in% not_scseq_datasets],
                                                                f =hypoMap_v2_seurat@meta.data$C185[hypoMap_v2_seurat@meta.data$C185 %in% clusters_to_check & !hypoMap_v2_seurat@meta.data$Dataset %in% not_scseq_datasets] )
mean_expr_sc = lapply(cluster_list_sc_seq,build_cell,expr_matrix = hypoMap_v2_seurat@assays$RNA@data)#seurat_object_harmonized
message("Time used to summarise expression matrix : ",Sys.time()-st)
mean_expression_sc = as.matrix(do.call(cbind,mean_expr_sc))

# for snseq
st=Sys.time()
cluster_list_sn_seq = split(hypoMap_v2_seurat@meta.data$Cell_ID[hypoMap_v2_seurat@meta.data$C185 %in% clusters_to_check & hypoMap_v2_seurat@meta.data$Dataset == "Dowsett10xnuc"],
                            f =hypoMap_v2_seurat@meta.data$C185[hypoMap_v2_seurat@meta.data$C185 %in% clusters_to_check & hypoMap_v2_seurat@meta.data$Dataset == "Dowsett10xnuc"] )
mean_expr_sn = lapply(cluster_list_sn_seq,build_cell,expr_matrix = hypoMap_v2_seurat@assays$RNA@data)#seurat_object_harmonized
message("Time used to summarise expression matrix : ",Sys.time()-st)
mean_expression_sn = as.matrix(do.call(cbind,mean_expr_sn))

##########
### Calculate highest percentage per cluster to subset genes
##########

# minimum values for genes to keep
min_pct_genes = 0.2
min_mean_genes = 0.2

## add a max pct value
Idents(hypoMap_v2_seurat) = "C185"
all_clusterstats = list()
for(cluster in unique(hypoMap_v2_seurat@meta.data$C185)){
  message(cluster)
  cells1 = hypoMap_v2_seurat@meta.data$Cell_ID[hypoMap_v2_seurat@meta.data$C185 == cluster]
  all_clusterstats[[cluster]] = get_expression_stats(object=hypoMap_v2_seurat, cells.1=cells1,features = NULL,  thresh.min = 0)
  all_clusterstats[[cluster]]$cluster = cluster
}
all_clusterstats_df = do.call(rbind,all_clusterstats)
all_clusterstats_df_summarized = all_clusterstats_df %>% dplyr::filter(pct.1 > 0) %>% dplyr::group_by(gene) %>% dplyr::filter(pct.1 == max(pct.1)) %>% distinct(gene,.keep_all = TRUE)
genes_to_keep = all_clusterstats_df_summarized$gene[all_clusterstats_df_summarized$mean.1>min_mean_genes | all_clusterstats_df_summarized$pct.1 > min_pct_genes]

##########
### Calculate per gene correlations
##########

# # need to subset to same genes 
# shared_genes = sort(intersect(rownames(neuron_map_mean_expression),rownames(snseq_mean_expression)))
# shared_genes = shared_genes[shared_genes %in% all_genes_to_keep] 
# 
# # and same order
# neuron_map_mean_expression_subset = neuron_map_mean_expression[genes_to_keep,]
# snseq_mean_expression_subset = snseq_mean_expression[genes_to_keep,]

all_cors=sapply(genes_to_keep,function(gene,mat1,mat2){
  suppressWarnings(cor(mat1[gene,],mat2[gene,],use="pairwise.complete.obs",method ="pearson"))
},mat1=mean_expression_sc,mat2=mean_expression_sn)
all_cors_spearman=sapply(genes_to_keep,function(gene,mat1,mat2){
  suppressWarnings(cor(mat1[gene,],mat2[gene,],use="pairwise.complete.obs",method ="spearman"))
},mat1=mean_expression_sc,mat2=mean_expression_sn)
per_gene_cor = data.frame(gene = names(all_cors), pearson = all_cors,spearman = all_cors_spearman)

## add more info to gene_cors

padj_cut = 0.001
fc_min = 0.25

## simple classificationfor markers:
per_gene_cor$marker_gene = FALSE
per_gene_cor$marker_gene[per_gene_cor$gene %in% markers_C185$gene ] = TRUE

# save gene cors!

data.table::fwrite(per_gene_cor,paste0(results_path_figure4,"per_gene_correlations.txt"),sep = "\t")

#per_gene_cor = data.table::fread(paste0(results_path_figure4,"per_gene_correlations.txt"),data.table = FALSE)

##########
### Heatmaps per gene 'class'
##########

## we defined gene classes based on KEGG, GO and Ingenuity that can be loaded from a json file provided in this repo
# to make this scipt shroter I have removed the full creation procedure.
class_list = jsonlite::read_json("data_inputs/gene_classes.json")
names(class_list) = paste0(names(class_list),"s")

# add markers from above to list
# class_list[["Markers"]] = per_gene_cor$gene[per_gene_cor$both_marker]
# class_list[["Unique sc-seq"]] = per_gene_cor$gene[per_gene_cor$sc_marker]
# class_list[["Unique Nuc-seq"]] = per_gene_cor$gene[per_gene_cor$sn_marker]

# make a combined receptor class
class_list[["enzyme"]] = unique(c(class_list[["peptidase"]],class_list[["Gphosphatase"]],class_list[["kinase"]],class_list[["enzyme"]]))

# reduce class lists to the ones we want to show in the final map:
class_list_subset = class_list[names(class_list) %in% c("ligand-dependent nuclear receptors", "transmembrane receptors", "G-protein coupled receptors",
                                                        "neuropeptide/hormones","transcription regulators","translation regulators",
                                                        "ion channels","growth factors","Markers" ,"Unique sc-seq","Unique Nuc-seq")]

# add to per_gene_cor and save
class_list_subset_df <- purrr::map_df(class_list_subset, ~as.data.frame(.x), .id="class")
class_list_subset_df = class_list_subset_df %>% tidyr::gather(key="key",value="gene",-class) %>% dplyr::filter(!is.na(gene)) %>% dplyr::select(-key)  %>%
  dplyr::distinct(class,gene)
rownames(class_list_subset_df) = 1:nrow(class_list_subset_df)
class_list_subset_df$val = TRUE
class_list_subset_df_wide = class_list_subset_df %>% tidyr::spread(key=class,value = val)
class_list_subset_df_wide[is.na(class_list_subset_df_wide)] = FALSE
per_gene_cor_with_class = dplyr::left_join(per_gene_cor,class_list_subset_df_wide,by=c("gene"="gene"))
per_gene_cor_with_class[is.na(per_gene_cor_with_class)] = FALSE
data.table::fwrite(per_gene_cor_with_class,paste0(results_path_figure4,"per_gene_correlations_with_class.txt"),sep = "\t")

#per_gene_cor_with_class = data.table::fread(paste0(results_path_figure4,"per_gene_correlations_with_class.txt"),data.table = FALSE)
                                   
# get density per class
n_density = 300
list_df = lapply(class_list_subset,function(genes,per_gene_cor,subset_genes,min_genes = 10,n_density = n_density){
  genes = unlist(genes)
  genes = genes[genes %in% subset_genes & genes %in% unique(per_gene_cor$gene)]
  if(length(genes)>min_genes){
    y = density(per_gene_cor$pearson[per_gene_cor$gene %in% genes] %>% na.omit(),n=n_density,from = -1,to = 1)$y
  }else{
    y=NULL 
  }
  return(y)
}, per_gene_cor = per_gene_cor,subset_genes = per_gene_cor$gene,n_density=n_density)

# make a dataframe for plotting with ggplot2
steps = seq(-1,1,(1/(n_density*0.5)))
density_per_class_df = as.data.frame(cbind(steps=steps[2:length(steps)],do.call(cbind,list_df)))
density_per_class_df_long = density_per_class_df %>% tidyr::gather(key="class","density",-steps)

# add classes to split heatmap by
density_per_class_df_long$group = "gene classes"
density_per_class_df_long$group[density_per_class_df_long$class %in% c("Markers" ,"Unique sc-seq","Unique Nuc-seq")] = "celltype markers"
density_per_class_df_long$group = factor(density_per_class_df_long$group,levels=(c("celltype markers","gene classes")))

# neuropeptide/hormone
# add numbers per class
gene_per_class = sapply(class_list_subset,function(genes,subset_genes){length(genes[genes %in% subset_genes])},subset_genes = per_gene_cor$gene)
gene_per_class = data.frame(class = names(gene_per_class),n_genes = gene_per_class)
density_per_class_df_long = density_per_class_df_long %>% dplyr::left_join(gene_per_class,by="class")

# reorder based on maximal pearson value, see below solution for one that worked in this case
density_per_class_df_long = density_per_class_df_long %>% dplyr::group_by(class) %>% dplyr::mutate(max_pearson = steps[density == max(density)])
# density_per_class_df_long = density_per_class_df_long %>% dplyr::group_by(class) %>% dplyr::mutate(q95_closest = (density - quantile(density,probs = 0.95)) == min(density - quantile(density,probs = 0.95))) #%>% 
#   dplyr::mutate(max_pearson = steps[q95_closest == density])
density_per_class_df_long = density_per_class_df_long %>% dplyr::arrange((max_pearson))
density_per_class_df_long$class_ordered = factor(density_per_class_df_long$class,levels = unique(density_per_class_df_long$class))
#or:
#density_per_class_df_long$class = forcats::fct_reorder(.f = density_per_class_df_long$class, .x = density_per_class_df_long$max_pearson, .fun = mean)

# ormanually
density_per_class_df_long$class_ordered = factor(density_per_class_df_long$class,levels = rev(c("neuropeptide/hormones" ,"G-protein coupled receptors", "transmembrane receptors","growth factors" ,"ion channels" , "ligand-dependent nuclear receptors","transcription regulators","translation regulators" )))


## make heatmap 
require(scales)
class_heatmap = ggplot2::ggplot(density_per_class_df_long[density_per_class_df_long$steps> (-0.3),], aes(x = class_ordered, y = steps , fill = density)) +
  geom_tile() +
  geom_text(aes(x = class, y = 1.11,label=n_genes), hjust =1,size=5)+ # add cell numbers
  ggplot2::scale_fill_gradient(low="white",high="#ff1414",na.value = bg_col,
                               limits=c(0,max(density_per_class_df_long$density)), oob=squish) +
  xlab("Gene class")+ylab("Pearson correlation")+ guides(fill=guide_legend(title="Density"))+
  scale_fill_gradient(low = "white",high = reference_sc_color)+
  guides(fill = guide_colorbar()) + # ensure continous colorbar
  geom_hline(yintercept = 0,color="grey60",linetype = "dashed")+
  coord_flip()+facet_grid(group ~ ., scales = "free",space = "free_y",drop=TRUE)  + 
  theme(text = element_text(size=25), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing = unit(2, "lines"),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white"))+ 
  scale_y_continuous(breaks=seq(0,1,0.25))
# and show
class_heatmap

#store:
ggsave(filename = paste0(results_path_figure4,"geneclass_cor_heatmap.png"),
       plot = class_heatmap, "png",dpi=600,width=330,height = 200,units="mm")

ggsave(filename = paste0(results_path_figure4,"geneclass_cor_heatmap.pdf"),
       plot = class_heatmap, "pdf",dpi=600,width=330,height = 200,units="mm")

##source data
class_heatmap_data = class_heatmap$data
data.table::fwrite(class_heatmap_data,paste0(results_path_figure4,"source_figure4_b_heatmap_data.txt"),sep="\t")
# also add the gene correlations ? --> no keep in supplemenatry tables

### calculate table with stats per class for paper Text
summary_per_class_list = lapply(class_list_subset,function(genes,per_gene_cor,subset_genes,min_genes = 10,n_density = n_density){
  genes = genes[genes %in% subset_genes & genes %in% unique(per_gene_cor$gene)]
  vals = per_gene_cor$pearson[per_gene_cor$gene %in% genes] %>% na.omit()
  tmp = data.frame(mean = mean(vals), sd = sd(vals), sem = sd(vals)/sqrt(length(vals)))
  return(tmp)
}, per_gene_cor = per_gene_cor,subset_genes = per_gene_cor$gene,n_density=n_density)
summary_per_class = do.call(rbind,summary_per_class_list)
summary_per_class$class = rownames(summary_per_class)
data.table::fwrite(summary_per_class,paste0(results_path_figure4,"summary_per_class.txt"),sep="\t")
# specifically for bimodal:
# vals = per_gene_cor$pearson[per_gene_cor$gene %in% unlist(class_list_subset$`transmembrane receptor`)]
# length(vals[vals<0.6]) / length(vals)
# length(vals[vals>0.6])  / length(vals)

##########
### cluster-centric comparison:

##########
### Per cluster and Dataset: RSQ and correlation
##########

# how many cells should be there per cluster AND dataset to include averaged value
min_cells_with_dataset = 10
# miniumum average expression value:
min_expr = 0.1
# remove uncorrelated genes
include_genes = unique(per_gene_cor$gene[per_gene_cor$pearson>0.3])
# ensure that only genes occuring in both are compared:
# use shared_genes which are the same as all genes in the per_gene_cor df as defined above!!

#### Need to get mean per dataset and cluster
datasets = unique(hypoMap_v2_seurat@meta.data$Dataset)
datasets = datasets[! datasets %in% c("Affinati10x","Rupp10x")]
datasets_mean_expression_list =list()
for(dset in datasets){
  message(dset)
  st=Sys.time()
  cluster_list_dset = split(hypoMap_v2_seurat@meta.data$Cell_ID[hypoMap_v2_seurat@meta.data$C185 %in% clusters_to_check & hypoMap_v2_seurat@meta.data$Dataset %in% dset],
                            f =hypoMap_v2_seurat@meta.data$C185[hypoMap_v2_seurat@meta.data$C185 %in% clusters_to_check & hypoMap_v2_seurat@meta.data$Dataset %in% dset] 
  )
  # filter n cells
  cluster_list_dset = cluster_list_dset[sapply(cluster_list_dset,length) >= min_cells_with_dataset]
  mean_expr = suppressMessages(lapply(cluster_list_dset,build_cell,expr_matrix = hypoMap_v2_seurat@assays$RNA@data))#seurat_object_harmonized
  message("Time used to summarise expression matrix : ",Sys.time()-st)
  datasets_mean_expression_list[[dset]] = as.matrix(do.call(cbind,mean_expr))
}

# init
cor_method="pearson"
mean_lm_res = list()
mean_adj_rsq = numeric()
per_dataset_cor_res=list()
global_cor_res =numeric()
n_markers =numeric()
# loop over all clusters
for(i in 1:length(clusters_to_check)){
  current_cluster= clusters_to_check[i]
  #message("  ",current_cluster)
  
  # get cluster markers
  #sn_Seq_based_markers = sn_seq_markers_C185$gene[sn_seq_markers_C185$p_val_adj<padj_cut & sn_seq_markers_C185$avg_log2FC>fc_min & sn_seq_markers_C185$cluster==current_cluster & sn_seq_markers_C185$gene %in% per_gene_cor$gene]
  #sc_Seq_based_markers = sc_seq_markers_C185$gene[sc_seq_markers_C185$p_val_adj<padj_cut & sc_seq_markers_C185$avg_logFC>fc_min & sc_seq_markers_C185$cluster==current_cluster & sc_seq_markers_C185$gene %in% per_gene_cor$gene]
  #union_markers = union(sn_Seq_based_markers,sc_Seq_based_markers)
  #union_markers = sn_Seq_based_markers
  union_markers = markers_C185$gene[markers_C185$cluster_id == current_cluster]
  n_markers = c(n_markers,length(union_markers))
  names(n_markers)[i] = current_cluster
  ##### RSQ between sn-seq and sc-seq (split up by Dataset)
  # calculate mean expression 
  mean_expressions= do.call(cbind,lapply(datasets_mean_expression_list, function(x,cluster,genes){if(cluster %in% colnames(x)){x[genes,cluster]}},cluster=current_cluster,genes = union_markers))
  #mean_expressions= cbind(sn_seq = snseq_mean_expression[union_markers,current_cluster],mean_expressions)
  
  # filter (TODO: do I want to filter here?)
  # mean_expressions_keep = apply(mean_expressions,1,function(x,min_expr){if(length(x[x>min_expr])>1){return(TRUE)}else{return(FALSE)}},min_expr=min_expr)
  # mean_expressions = mean_expressions[mean_expressions_keep,]
  # construct formula:
  predictors = colnames(mean_expressions)[!colnames(mean_expressions) %in% c("gene","Dowsett10xnuc")]
  #as.formula(paste0("snseq ~ ",paste0(predictors,collapse=" + ")))
  # fit lm
  if(length(union_markers)>3){
    if(length(predictors)>=1){
      mean_lm_res[[current_cluster]] = lm(formula = as.formula(paste0("Dowsett10xnuc ~ ",paste0(predictors,collapse=" + "))),data = as.data.frame(mean_expressions))
      mean_adj_rsq[current_cluster] = summary(mean_lm_res[[current_cluster]])$adj.r.squared
      # cor per cluster
      dataset_cor_vec = cor(mean_expressions[,"Dowsett10xnuc"],mean_expressions[,predictors],use = "pairwise.complete.obs",method=cor_method)
      if(length(dataset_cor_vec)>1){
        names(dataset_cor_vec)= predictors#colnames(mean_expressions[,2:ncol(mean_expressions)])
        dataset_cor_vec = dataset_cor_vec[match(names(datasets_mean_expression_list),names(dataset_cor_vec))] # include missing datasets with NA
        names(dataset_cor_vec) = names(datasets_mean_expression_list)
      }else{
        save = dataset_cor_vec
        dataset_cor_vec = rep(NA,length(names(datasets_mean_expression_list)))
        names(dataset_cor_vec) = names(datasets_mean_expression_list)
        dataset_cor_vec[predictors] = save
      }
      per_dataset_cor_res[[current_cluster]] =dataset_cor_vec
    }
    ## correlation across all
    mean_expr_sc = data.frame(mean_expr = mean_expression_sc[,current_cluster],gene = rownames(mean_expression_sc))
    mean_expr_sn = data.frame(mean_expr = mean_expression_sn[,current_cluster],gene = rownames(mean_expression_sn))
    compare_mean_between = dplyr::left_join(mean_expr_sc,mean_expr_sn,by="gene",suffix=c("_sc","_sn")) %>% dplyr::filter(gene %in% union_markers)
    global_cor_res[current_cluster] = cor(compare_mean_between$mean_expr_sc,compare_mean_between$mean_expr_sn,use = "pairwise.complete.obs",method=cor_method)
  }
}

## summarize:
ndigits = 5
# per_dataset_cor_res_sq = lapply(per_dataset_cor_res,function(x,ndigits){round(x^2,digits = ndigits)},ndigits=ndigits)
per_dataset_cor_df =as.data.frame(do.call(rbind, per_dataset_cor_res))  #as.data.frame(do.call(rbind, per_dataset_cor_res_sq))
#per_dataset_cor_df$mean_adj_rsq = round(mean_adj_rsq,digits = ndigits)
per_dataset_cor_df$global_cor = round(global_cor_res,digits = ndigits) #round(global_cor_res^2,digits = ndigits)
per_dataset_cor_df$C185 = rownames(per_dataset_cor_df)
per_dataset_cor_df$n_markers = n_markers
# add paired anno
#per_dataset_cor_df$C185_named = mapscvi::add_paired_annotation(input_annotation = per_dataset_cor_df$C185,reference_annotations = hypoMap_v2_seurat@meta.data[,c("C185","C185_named")])

# add Rupp and Affinati as all NA columns
per_dataset_cor_df$Rupp10x = NA
per_dataset_cor_df$Affinati10x = NA

# save gene cors!
per_dataset_cor_df_print = per_dataset_cor_df[,c(18,17,19,21,1:14,20,15,16)]
data.table::fwrite(per_dataset_cor_df_print,paste0(results_path_figure4,"per_cluster_correlations.txt"),sep = "\t")

#per_dataset_cor_df = data.table::fread(paste0(results_path_figure4,"per_cluster_correlations.txt"),data.table = FALSE)

##########
### Make heatmap - currently use treeplot below!
##########

## heatmap
# per_dataset_cor_df = per_dataset_cor_df %>% dplyr::select(-mean_adj_rsq )
# per_dataset_cor_long = per_dataset_cor_df %>% dplyr::arrange(desc(global_cor_res)) %>% tidyr::gather("Dataset","RSQ",-C185,-C185_named) 
# per_dataset_cor_long$Dataset = factor(per_dataset_cor_long$Dataset,levels = c("global_cor","Kim","Campbell","Chen","wenDropSeq","Lee_Idol","Mickelsen","Rossi","Moffit","Flynn","RomanovDev","Mousebrainorg","wen10x","kimDev"))
# 
# per_dataset_cor_long$C185_named = factor(per_dataset_cor_long$C185_named , levels = per_dataset_cor_df$C185_named[order(per_dataset_cor_df$global_cor)])
# 
# # make heatmap 
# require(scales)
# dff_heatmap = ggplot(per_dataset_cor_long, aes(x = C185_named, y = Dataset , fill = RSQ)) +
#   geom_tile() +
#   # facet_grid(sample~., scales="free_y") +
#   #theme(axis.text.y = element_blank())+ 
#   ggplot2::scale_fill_gradient(low="#1a1a82",high="#ff1414",na.value = "white",
#                                limits=c(0,1), oob=squish) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# # and show
# dff_heatmap

##### cell numbers

# sc_seq_cell_numbers = data.frame(cluster_id = names(cluster_list_sc_seq), ncells_scseq = sapply(cluster_list_sc_seq,length))
# sn_seq_cell_numbers = data.frame(cluster_id = names(cluster_list_sn_seq), ncells_snseq = sapply(cluster_list_sn_seq,length))
# 
# cell_numbers = left_join(sc_seq_cell_numbers,sn_seq_cell_numbers,by="cluster_id")
# cell_numbers$nucseq_ratio = cell_numbers$ncells_snseq /  cell_numbers$ncells_scseq
# cell_numbers = left_join(cell_numbers,heatmap_data[,c("global_cor","C185","n_markers")],by=c("cluster_id"="C185"))
# 
# cor(cell_numbers[,2:ncol(cell_numbers)])

##########
### Make Tree heatmap 
##########

leaf_level_column = "C185"
leaf_level = 6

# make data for first heatmap with correlations
datasets = unique(hypoMap_v2_seurat@meta.data$Dataset)#[ unique(hypoMap_v2_seurat@meta.data$Dataset) != "Dowsett10xnuc"]
heatmap_data = per_dataset_cor_df

heatmap_matrix1 = as.matrix(heatmap_data[,"n_markers"]) # might want to change order
colnames(heatmap_matrix1) = "M"
rownames(heatmap_matrix1) = heatmap_data$C185

heatmap_matrix2 = as.matrix(heatmap_data[,"global_cor"]) # might want to change order
colnames(heatmap_matrix2) = "Cor"
rownames(heatmap_matrix2) = heatmap_data$C185

heatmap_matrix3 = as.matrix(heatmap_data[,datasets]) # might want to change order
rownames(heatmap_matrix3) = heatmap_data$C185

# numbers
colnames_overview = data.frame(number = 1:ncol(heatmap_matrix3),Dataset = colnames(heatmap_matrix3))
data.table::fwrite(colnames_overview,paste0(results_path_figure4,"circular_tree_colnames.txt"),sep="\t")
colnames(heatmap_matrix3) = colnames_overview$number

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
circular_tree_heat_cor = add_heatmap(circular_tree=circular_tree,
                                 heatmap_matrix = heatmap_matrix1,
                                 heatmap_colors=c(bg_col,"darkred"),
                                 scale_limits = c(0,400),
                                 heatmap_colnames =TRUE, 
                                 legend_title = "N Markers",
                                 matrix_offset = 0.2,
                                 matrix_width =0.05,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=0.75, 
                                 na_color = "white",
                                 heatmap_text_size=4)
# global cor
circular_tree_heat_cor = add_heatmap(circular_tree=circular_tree_heat_cor,
                                 heatmap_matrix = heatmap_matrix2,
                                 heatmap_colors= c("#1a1a82","#14fa14"),
                                 scale_limits = c(0,1),
                                 heatmap_colnames =TRUE, 
                                 legend_title = "Correlation",
                                 matrix_offset = 0.5,
                                 matrix_width =0.05,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=1, 
                                 na_color = "white",
                                 heatmap_text_size=4)
# per dataset cor
circular_tree_heat_cor = add_heatmap(circular_tree=circular_tree_heat_cor,
                                 heatmap_matrix = heatmap_matrix3,
                                 heatmap_colors= c("#1a1a82","#ff1414"),
                                 scale_limits = c(0,1),
                                 heatmap_colnames =TRUE, 
                                 legend_title = "Cor Dataset",
                                 matrix_offset = 0.9,
                                 matrix_width =0.4,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=2.5, 
                                 na_color = "white",
                                 heatmap_text_size=2.5)

circular_tree_heat_cor

#save:
ggsave(filename = paste0(results_path_figure4,"circular_tree_correlation.png"),
       plot = circular_tree_heat_cor, "png",dpi=600,width=400,height = 400,units="mm")
ggsave(filename = paste0(results_path_figure4,"circular_tree_correlation.pdf"),
       plot = circular_tree_heat_cor, "pdf",dpi=600,width=400,height = 400,units="mm")


figure4_tree_heatmap_sourcedata = heatmap_data
data.table::fwrite(figure4_tree_heatmap_sourcedata,paste0(results_path_figure4,"source_figure4_c_tree_heatmap.txt"),sep="\t")

data.table::fwrite(hypoMap_v2_seurat@misc$clustering_edgelist,paste0(results_path_figure4,"source_figure4_c_tree_edgelist.txt"),sep="\t")

data.table::fwrite(anno_df,paste0(results_path_figure4,"source_figure4_c_tree_annotations.txt"),sep="\t")




