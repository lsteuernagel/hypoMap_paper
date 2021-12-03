
##########
### Load & Prepare
##########

results_path = "figure_outputs/figure_4/"
system(paste0("mkdir -p ",results_path))

# set path for additional data
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"

# load functions
source("R/utility_functions.R")
source("R/plot_functions.R")
library(mapscvi) # please install the mapscvi package that was used to project the nucseq data and provides additional visualization functions for the projected data

# load everything required
source("R/load_data.R")

## load nucseq mapped object
query_snseq_neurons = readRDS(paste0(large_data_path,"nucseq_neurons_map.rds"))

# colors
reference_color = "#cc2118"
query_sn_color = "#302ac9"
bg_col = "#dedede"

rasterize_point_size = 2.2
rasterize_pixels = 2048

##########
### Plot mapping results
##########

## with labels on hypomap background:
mapped_query_snseq_neurons_plot = mapscvi::plot_query_labels(query_seura_object=query_snseq_neurons,reference_seurat=neuron_map_seurat,label_col="K31_named",
                                                             label_col_query = "predicted_K31_named",overlay = TRUE,bg_col = bg_col,
                                                             query_pt_size = 0.05,labelonplot = TRUE,label.size=5)
mapped_query_snseq_neurons_plot = rasterize_ggplot(mapped_query_snseq_neurons_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
mapped_query_snseq_neurons_plot

ggsave(filename = paste0(results_path,"mapping_snseq.png"),
       plot = mapped_query_snseq_neurons_plot, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"mapping_snseq.pdf"),
       plot = mapped_query_snseq_neurons_plot, "pdf",dpi=450,width=200,height = 200,units="mm")


## mapping quality:
#colnames(query_snseq_neurons@meta.data)
quality_query_snseq_neurons_plot=Seurat::FeaturePlot(query_snseq_neurons,features = "prediction_probability",reduction = paste0("umap_","scvi"),label = FALSE,raster=FALSE)+
  NoAxes()+ggtitle("Mapping probability")
quality_query_snseq_neurons_plot = rasterize_ggplot(quality_query_snseq_neurons_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
quality_query_snseq_neurons_plot

ggsave(filename = paste0(results_path,"mapping_prob_snseq.png"),
       plot = quality_query_snseq_neurons_plot, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"mapping_prob_snseq.pdf"),
       plot = quality_query_snseq_neurons_plot, "pdf",dpi=450,width=200,height = 200,units="mm")


##########
### Sankey plots
##########

# we selected the VMH and the SCN neurons as examples

# create cluster overview using mapscvi:
overview_clustering = mapscvi::compare_clustering(query_snseq_neurons,"predicted_K169_named","Cluster_IDs",min_cells = 10,min_pct = 0.1,return_data=TRUE)
data.table::fwrite(overview_clustering,paste0(results_path,"overview_clustering_nucSeq.txt"),sep="\t")

#### SCN plots
# sankey
clustering_2_filter = c("C32-C1ql3/Rgs16","C40-Vip")
clustering_1_filter = overview_clustering$clustering_1[overview_clustering$clustering_2 %in% clustering_2_filter]
sankey_scn = mapscvi::plot_sankey_comparison(overview_clustering,clustering_1_filter = clustering_1_filter,clustering_2_filter = clustering_2_filter,
                                             text_size=20, col1 = reference_color, col2 = query_sn_color)
sankey_scn
# get the numbers:
sankey_scn_numbers = mapscvi::plot_sankey_comparison(overview_clustering,clustering_1_filter = clustering_1_filter,clustering_2_filter = clustering_2_filter,text_size=20, col1 = reference_color, col2 = query_sn_color,return_data = TRUE)
sum(sankey_scn_numbers$edges$n) / ncol(query_snseq_neurons) # total number of SCN neurons in this analysis (with the sankey)
length(query_snseq_neurons@meta.data$Cell_ID[query_snseq_neurons@meta.data$predicted_K169_named %in% sankey_scn_numbers$edges$clustering_1]) # this is all cells mapped to the likely SCN clusters
# are the clustering_1_filter really scn:
table(neuron_map_seurat@meta.data$suggested_region_curated[neuron_map_seurat@meta.data$K169_named %in% clustering_1_filter],neuron_map_seurat@meta.data$K169_named[neuron_map_seurat@meta.data$K169_named %in% clustering_1_filter])
unique(neuron_map_seurat@meta.data$other_likely_regions[neuron_map_seurat@meta.data$K169_named %in% clustering_1_filter])


# orientation umap:
clustering_2_filter = c("C32-C1ql3/Rgs16","C40-Vip")
cellsh = query_snseq_neurons@meta.data$Cell_ID[query_snseq_neurons@meta.data$Cluster_IDs %in% sankey_scn_numbers$edges$clustering_2 & 
                                                 query_snseq_neurons@meta.data$predicted_K169_named %in% sankey_scn_numbers$edges$clustering_1]
scn_dimplot = DimPlot(query_snseq_neurons,cells.highlight = cellsh,sizes.highlight = 0.15)+NoLegend()+NoAxes()
scn_dimplot = rasterize_ggplot(scn_dimplot,pixel_raster = 1024,pointsize = 1.1)
scn_dimplot

# save sankeys
library(webshot)
# https://stackoverflow.com/questions/65158327/how-can-i-save-a-networkd3sankeynetwork-into-a-static-image-automatically-vi
networkD3::saveNetwork(sankey_scn, paste0(results_path,"sankey_scn_neurons.html"))
# convert it
# need: webshot::install_phantomjs()
webshot::webshot(paste0(results_path,"sankey_scn_neurons.html"),file=paste0(results_path,"sankey_scn_neurons.png"), vwidth = 1000, vheight = 900)
webshot::webshot(paste0(results_path,"sankey_scn_neurons.html"),file=paste0(results_path,"sankey_scn_neurons.pdf"), vwidth = 1000, vheight = 900)

# save dimplot
ggsave(filename = paste0(results_path,"scn_dimplot.png"),
       plot = scn_dimplot, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"scn_dimplot.pdf"),
       plot = scn_dimplot, "pdf",dpi=450,width=200,height = 200,units="mm")


#### VMH plots
# sankey
clustering_1_filter = c("Fgf10.Gpr149.Cd40.Fezf1.HY2")
clustering_2_filter = overview_clustering$clustering_2[overview_clustering$clustering_1=="Fgf10.Gpr149.Cd40.Fezf1.HY2"]
sankey_vmh = mapscvi::plot_sankey_comparison(overview_clustering,clustering_1_filter = clustering_1_filter,clustering_2_filter = clustering_2_filter,text_size=20, col1 = reference_color, col2 = query_sn_color)
sankey_vmh

## if we would like to check on a lower level:
# overview_clustering_2 = mapscvi::compare_clustering(query_snseq_neurons,"predicted_K329_named","Cluster_IDs",min_cells = 10,min_pct = 0.1,return_data=TRUE)
# clustering_1_filter = c("Fam19a1.Gpc5.Fgf10.Gpr149.Cd40.Fezf1.HY2","Rbms2.Fgf10.Gpr149.Cd40.Fezf1.HY2","Tac1.Gpc5.Fgf10.Gpr149.Cd40.Fezf1.HY2")
# clustering_2_filter = overview_clustering_2$clustering_2[overview_clustering_2$clustering_1 %in% clustering_1_filter]
# sankey_vmh = mapscvi::plot_sankey_comparison(overview_clustering_2,clustering_1_filter = clustering_1_filter,clustering_2_filter = clustering_2_filter,text_size=20, col1 = reference_color, col2 = query_sn_color)
# sankey_vmh

length(query_snseq_neurons@meta.data$Cell_ID[query_snseq_neurons@meta.data$predicted_K169_named %in% sankey_scn_numbers$edges$clustering_1]) # this is all cells mapped to the likely SCN clusters

# orientation umap:
clustering_2_filter = overview_clustering$clustering_2[overview_clustering$clustering_1=="Fgf10.Gpr149.Cd40.Fezf1.HY2"]
cellsh = query_snseq_neurons@meta.data$Cell_ID[query_snseq_neurons@meta.data$Cluster_IDs %in% clustering_2_filter]
vmh_dimplot = DimPlot(query_snseq_neurons,cells.highlight = cellsh,sizes.highlight = 0.15)+NoLegend()+NoAxes()
vmh_dimplot = rasterize_ggplot(vmh_dimplot,pixel_raster = 1024,pointsize = 1.1)
vmh_dimplot

# save sankeys
networkD3::saveNetwork(sankey_vmh, paste0(results_path,"sankey_vmh_neurons.html"))
# convert it
webshot::webshot(paste0(results_path,"sankey_vmh_neurons.html"),file=paste0(results_path,"sankey_vmh_neurons.png"), vwidth = 1000, vheight = 900)
webshot::webshot(paste0(results_path,"sankey_vmh_neurons.html"),file=paste0(results_path,"sankey_vmh_neurons.pdf"), vwidth = 1000, vheight = 900)

# save dimplot
ggsave(filename = paste0(results_path,"vmh_dimplot.png"),
       plot = vmh_dimplot, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"vmh_dimplot.pdf"),
       plot = vmh_dimplot, "pdf",dpi=450,width=200,height = 200,units="mm")


###########################
########## Dynamic range of expression 


##########
### Cluster marker genes
##########

# load marker genes sc seq
sc_seq_markers_K169 = neuron_map_seurat@misc$markers_comparisons_all[neuron_map_seurat@misc$markers_comparisons_all$p_val_adj<0.01,]

# load marker genes sn seq
sn_seq_markers_K169 = data.table::fread(paste0("data_inputs/sn_seq_mapped_neurons_K169_markers_2_sampleID.txt"),data.table = F)
sn_seq_markers_K169$specificity = (sn_seq_markers_K169$pct.1 / sn_seq_markers_K169$pct.2) * sn_seq_markers_K169$avg_logFC
sn_seq_markers_K169$avg_log2FC = sn_seq_markers_K169$avg_logFC

# which cluster seem relevant (have a sufficient number of cells mapped in sn-seq data to calculate reliable averages)
min_cells = 30
clusters_to_check = names(table(query_snseq_neurons@meta.data$predicted_K169_pruned))[table(query_snseq_neurons@meta.data$predicted_K169_pruned)>=min_cells]
length(clusters_to_check)

##########
### gene-centric comparison:

##########
### global gene mean calculation
##########

# build metacell expression matrix
build_cell = function(cells,expr_matrix){
  if(length(cells)>1){
    return(rowMeans(expr_matrix[,expr_matrix@Dimnames[[2]] %in% cells]))
  }else{
    return(expr_matrix[,expr_matrix@Dimnames[[2]] %in% cells])
  }
}
st=Sys.time()
cluster_list_neuron_map = split(neuron_map_seurat@meta.data$Cell_ID[neuron_map_seurat@meta.data$K169_pruned %in% clusters_to_check],f =neuron_map_seurat@meta.data$K169_pruned[neuron_map_seurat@meta.data$K169_pruned %in% clusters_to_check] )
mean_expr = lapply(cluster_list_neuron_map,build_cell,expr_matrix = neuron_map_seurat@assays$RNA@data)#seurat_object_harmonized
message("Time used to summarise expression matrix : ",Sys.time()-st)
neuron_map_mean_expression = as.matrix(do.call(cbind,mean_expr))

# for snseq
st=Sys.time()
cluster_list_snseq = split(query_snseq_neurons@meta.data$Cell_ID[query_snseq_neurons@meta.data$predicted_K169_pruned %in% clusters_to_check],f =query_snseq_neurons@meta.data$predicted_K169_pruned[query_snseq_neurons@meta.data$predicted_K169_pruned %in% clusters_to_check] )
mean_expr_sn = lapply(cluster_list_snseq,build_cell,expr_matrix = query_snseq_neurons@assays$RNA@data)#seurat_object_harmonized
message("Time used to summarise expression matrix : ",Sys.time()-st)
snseq_mean_expression = as.matrix(do.call(cbind,mean_expr_sn))

##########
### Calculate highest percentage per cluster to subset genes
##########

# minimum values for genes to keep
min_pct_genes = 0.2
min_mean_genes = 0.2

## add a max pct value
Idents(neuron_map_seurat) = "K169_pruned"
all_clusterstats_sc = list()
for(cluster in unique(neuron_map_seurat@meta.data$K169_pruned)){
  message(cluster)
  cells1 = neuron_map_seurat@meta.data$Cell_ID[neuron_map_seurat@meta.data$K169_pruned == cluster]
  all_clusterstats_sc[[cluster]] = get_expression_stats(object=neuron_map_seurat, cells.1=cells1,features = NULL,  thresh.min = 0)
  all_clusterstats_sc[[cluster]]$cluster = cluster
}
all_clusterstats_sc_df = do.call(rbind,all_clusterstats_sc)
all_clusterstats_sc_df_summarized = all_clusterstats_sc_df %>% dplyr::filter(pct.1 > 0) %>% dplyr::group_by(gene) %>% dplyr::filter(pct.1 == max(pct.1)) %>% distinct(gene,.keep_all = TRUE)
sc_genes_to_keep = all_clusterstats_sc_df_summarized$gene[all_clusterstats_sc_df_summarized$mean.1>min_mean_genes | all_clusterstats_sc_df_summarized$pct.1 > min_pct_genes]

# run the same for sn
Idents(query_snseq_neurons) = "predicted_K169_pruned"
all_clusterstats_sn = list()
for(cluster in unique(query_snseq_neurons@meta.data$predicted_K169_pruned)){
  message(cluster)
  cells1 = query_snseq_neurons@meta.data$Cell_ID[query_snseq_neurons@meta.data$predicted_K169_pruned == cluster]
  if(length(cells1)>=5){
    all_clusterstats_sn[[cluster]] = get_expression_stats(object=query_snseq_neurons, cells.1=cells1,features = NULL,  thresh.min = 0)
    all_clusterstats_sn[[cluster]]$cluster = cluster
  }
}
all_clusterstats_sn_df = do.call(rbind,all_clusterstats_sn)
all_clusterstats_sn_df_summarized = all_clusterstats_sn_df %>% dplyr::filter(pct.1 > 0) %>% dplyr::group_by(gene) %>% dplyr::filter(pct.1 == max(pct.1)) %>% distinct(gene,.keep_all = TRUE)
sn_genes_to_keep = all_clusterstats_sn_df_summarized$gene[all_clusterstats_sn_df_summarized$mean.1>min_mean_genes | all_clusterstats_sn_df_summarized$pct.1 > min_pct_genes]

all_genes_to_keep = base::union(sc_genes_to_keep,sn_genes_to_keep)

# make a table that can be saved 
all_clusterstats_both_summarized = dplyr::full_join(all_clusterstats_sc_df_summarized,all_clusterstats_sn_df_summarized, by =c("gene"="gene"),suffix=c("_sc","_sn"))
# data.table::fwrite(all_clusterstats_both_summarized,paste0(results_path,"max_expression_in_any_cluster_per_gene.txt"),sep = "\t")

##########
### Calculate per gene correlations
##########

# need to subset to same genes 
shared_genes = sort(intersect(rownames(neuron_map_mean_expression),rownames(snseq_mean_expression)))
shared_genes = shared_genes[shared_genes %in% all_genes_to_keep] 

# and same order
neuron_map_mean_expression_subset = neuron_map_mean_expression[shared_genes,]
snseq_mean_expression_subset = snseq_mean_expression[shared_genes,]

all_cors=sapply(shared_genes,function(gene,mat1,mat2){
  suppressWarnings(cor(mat1[gene,],mat2[gene,],use="pairwise.complete.obs",method ="pearson"))
},mat1=neuron_map_mean_expression_subset,mat2=snseq_mean_expression_subset)
all_cors_spearman=sapply(shared_genes,function(gene,mat1,mat2){
  suppressWarnings(cor(mat1[gene,],mat2[gene,],use="pairwise.complete.obs",method ="spearman"))
},mat1=neuron_map_mean_expression_subset,mat2=snseq_mean_expression_subset)
per_gene_cor = data.frame(gene = names(all_cors), pearson = all_cors,spearman = all_cors_spearman)

## add more info to gene_cors

padj_cut = 0.001
fc_min = 0.25

## simple classificationfor markers:
per_gene_cor$sn_marker = FALSE
per_gene_cor$sn_marker[per_gene_cor$gene %in% sn_seq_markers_K169$gene[sn_seq_markers_K169$p_val_adj<padj_cut & sn_seq_markers_K169$avg_log2FC > fc_min & sn_seq_markers_K169$cluster %in% clusters_to_check]] = TRUE
per_gene_cor$sc_marker = FALSE
per_gene_cor$sc_marker[per_gene_cor$gene %in% sc_seq_markers_K169$gene[sc_seq_markers_K169$p_val_adj<padj_cut & sc_seq_markers_K169$avg_logFC > fc_min & sc_seq_markers_K169$cluster %in% clusters_to_check]] = TRUE
per_gene_cor$both_marker = FALSE
per_gene_cor$both_marker[per_gene_cor$sn_marker & per_gene_cor$sc_marker] = TRUE
per_gene_cor$sn_marker[per_gene_cor$both_marker] = FALSE
per_gene_cor$sc_marker[per_gene_cor$both_marker] = FALSE

## add max value
all_clusterstats_both_summarized_for_join = all_clusterstats_both_summarized %>% 
  dplyr::mutate(max_pct_any = base::pmax(pct.1_sc,pct.1_sn,na.rm = TRUE), max_mean_any =  base::pmax(pct.1_sc,mean.1_sn,na.rm = TRUE)) %>%
  dplyr::select(gene,max_pct_any,max_mean_any)
per_gene_cor = dplyr::left_join(per_gene_cor,all_clusterstats_both_summarized_for_join,by=c("gene"="gene"))
# save gene cors!
data.table::fwrite(per_gene_cor,paste0(results_path,"per_gene_correlations.txt"),sep = "\t")

#per_gene_cor = data.table::fread(paste0(results_path,"per_gene_correlations.txt"),data.table = FALSE)

##########
### Heatmaps per gene 'class'
##########

## we defined gene classes based on KEGG, GO and Ingenuity that can be loaded from a json file provided in this repo
# to make this scipt shroter I have removed the full creation procedure.
class_list = jsonlite::read_json("data_inputs/gene_classes.json")
names(class_list) = paste0(names(class_list),"s")

# add markers from above to list
class_list[["Markers"]] = per_gene_cor$gene[per_gene_cor$both_marker]
class_list[["Unique sc-seq"]] = per_gene_cor$gene[per_gene_cor$sc_marker]
class_list[["Unique Nuc-seq"]] = per_gene_cor$gene[per_gene_cor$sn_marker]

# make a combined receptor class
class_list[["enzyme"]] = unique(c(class_list[["peptidase"]],class_list[["Gphosphatase"]],class_list[["kinase"]],class_list[["enzyme"]]))

# reduce class lists to the ones we want to show in the final map:
class_list_subset = class_list[names(class_list) %in% c("ligand-dependent nuclear receptors", "transmembrane receptors", "G-protein coupled receptors",
                                                        "neuropeptide/hormones","transcription regulators","translation regulators",
                                                        "ion channels","growth factors","Markers" ,"Unique sc-seq","Unique Nuc-seq")]
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
density_per_class_df_long = density_per_class_df_long %>% dplyr::arrange((max_pearson))
density_per_class_df_long$class_ordered = factor(density_per_class_df_long$class,levels = unique(density_per_class_df_long$class))
#or:
#density_per_class_df_long$class = forcats::fct_reorder(.f = density_per_class_df_long$class, .x = density_per_class_df_long$max_pearson, .fun = mean)

# rename classes
# density_per_class_df_long$class[density_per_class_df_long$class=="markers_both"]="Markers"
# density_per_class_df_long$class[density_per_class_df_long$class=="markers_sn"]="Unique Nuc-seq"
# density_per_class_df_long$class[density_per_class_df_long$class=="markers_sc"]="Unique sc-seq"
# density_per_class_df_long$class[density_per_class_df_long$class=="neuropeptide/hormones"]="neuropeptides/hormones"


## make heatmap 
require(scales)
class_heatmap = ggplot2::ggplot(density_per_class_df_long[density_per_class_df_long$steps> (-0.3),], aes(x = class_ordered, y = steps , fill = density)) +
  geom_tile() +
  geom_text(aes(x = class, y = 1.11,label=n_genes), hjust =1,size=5)+ # add cell numbers
  ggplot2::scale_fill_gradient(low="white",high="#ff1414",na.value = "grey80",
                               limits=c(0,max(density_per_class_df_long$density)), oob=squish) +
  xlab("Gene class")+ylab("Pearson correlation")+ guides(fill=guide_legend(title="Density"))+
  scale_fill_gradient(low = "white",high = reference_color)+
  geom_hline(yintercept = 0,color="grey60",linetype = "dashed")+
  coord_flip()+facet_grid(group ~ ., scales = "free",space = "free_y",drop=TRUE)  + 
  theme(text = element_text(size=25), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing = unit(2, "lines"),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white"))
# and show
class_heatmap

#store:
ggsave(filename = paste0(results_path,"geneclass_cor_heatmap.png"),
       plot = class_heatmap, "png",dpi=600,width=330,height = 200,units="mm")

ggsave(filename = paste0(results_path,"geneclass_cor_heatmap.pdf"),
       plot = class_heatmap, "pdf",dpi=600,width=330,height = 200,units="mm")


### calculate table with stats per class for paper Text
summary_per_class_list = lapply(class_list_subset,function(genes,per_gene_cor,subset_genes,min_genes = 10,n_density = n_density){
  genes = genes[genes %in% subset_genes & genes %in% unique(per_gene_cor$gene)]
  vals = per_gene_cor$pearson[per_gene_cor$gene %in% genes] %>% na.omit()
  tmp = data.frame(mean = mean(vals), sd = sd(vals), sem = sd(vals)/sqrt(length(vals)))
  return(tmp)
}, per_gene_cor = per_gene_cor,subset_genes = per_gene_cor$gene,n_density=n_density)
summary_per_class = do.call(rbind,summary_per_class_list)

# specifically for bimodal:
vals = per_gene_cor$pearson[per_gene_cor$gene %in% unlist(class_list_subset$`transmembrane receptor`)]
length(vals[vals<0.6]) / length(vals)
length(vals[vals>0.6])  / length(vals)

##########
### gene characteristics -- we excluded this for now
##########
# 
# ## which genes do we want to check ?
# genes_of_interest = per_gene_cor$gene
# 
# ## query biomart to obtain relevant characteristics
# # I use nov2020.archive.ensembl.org because it is still on 38 mm genome.
# library(biomaRt)
# mart <- useMart(dataset="mmusculus_gene_ensembl",biomart='ensembl',host="nov2020.archive.ensembl.org") # ,host="feb2021.archive.ensembl.org"
# #all_attr = listAttributes(mart)
# all_attributes = c('ensembl_gene_id', 'external_gene_name',"gene_biotype","start_position","end_position","ensembl_transcript_id","external_transcript_name","transcript_start", "transcript_end","strand","chromosome_name",
#                    "transcript_biotype","transcript_count","transcript_length","transcription_start_site","transcript_count","percentage_gene_gc_content") # ,"ncoils","tmhmm","transcript_exon_intron"
# mouse_gene_info = getBM(attributes = all_attributes,filters='external_gene_name',values=genes_of_interest,mart = mart)
# 
# # ideas
# # gene length
# # shortest transcript
# # longest transcript
# # transcript count
# # protein-coding transcript count
# # percentage_gene_gc_content
# # "ncoils","tmhmm"
# 
# # define characteristics
# mouse_gene_info_charac = mouse_gene_info %>% dplyr::mutate(
#   gene_length = as.numeric(end_position) - as.numeric(start_position),
# ) %>% dplyr::group_by(ensembl_gene_id) %>% dplyr::mutate(
#   shortest_transcript = transcript_length[transcript_length == min(transcript_length)][1],
#   longest_transcript = transcript_length[transcript_length == max(transcript_length)][1],
# ) %>% dplyr::select(ensembl_gene_id,external_gene_name,gene_biotype,gene_length,transcript_count,shortest_transcript,longest_transcript,percentage_gene_gc_content) %>%
#   dplyr::distinct(ensembl_gene_id,external_gene_name,.keep_all = TRUE)
# 
# ## correlate numeric charac
# mouse_gene_info_charac_with_cor = dplyr::left_join(mouse_gene_info_charac, per_gene_cor %>% dplyr::select(gene,pearson_cor = pearson),by=c("external_gene_name"="gene"))
# 
# # remove some outliers
# mouse_gene_info_charac_with_cor$gene_length[mouse_gene_info_charac_with_cor$gene_length>1000000] = 1000000
# mouse_gene_info_charac_with_cor$shortest_transcript[mouse_gene_info_charac_with_cor$shortest_transcript>15000] = 15000
# mouse_gene_info_charac_with_cor$longest_transcript[mouse_gene_info_charac_with_cor$longest_transcript>15000] = 1500
# # example plot
# ggplot(mouse_gene_info_charac_with_cor,aes(pearson_cor,gene_length))+geom_point(size=0.3)
# # correlate
# cor(mouse_gene_info_charac_with_cor$pearson_cor,mouse_gene_info_charac_with_cor[,c("gene_length","transcript_count","shortest_transcript","longest_transcript", "percentage_gene_gc_content")],method="pearson",use = "pairwise.complete.obs")
# ## by biotype classes:
# table(mouse_gene_info_charac_with_cor$gene_biotype)
# sort(tapply(mouse_gene_info_charac_with_cor$pearson_cor,INDEX = mouse_gene_info_charac_with_cor$gene_biotype,FUN = mean),decreasing = TRUE)
# # generally to low



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
datasets = unique(neuron_map_seurat@meta.data$Dataset)
datasets_mean_expression_list =list()
for(dset in datasets){
  message(dset)
  st=Sys.time()
  cluster_list_dset = split(neuron_map_seurat@meta.data$Cell_ID[neuron_map_seurat@meta.data$K169_pruned %in% clusters_to_check & neuron_map_seurat@meta.data$Dataset %in% dset],
                            f =neuron_map_seurat@meta.data$K169_pruned[neuron_map_seurat@meta.data$K169_pruned %in% clusters_to_check & neuron_map_seurat@meta.data$Dataset %in% dset] 
  )
  # filter n cells
  cluster_list_dset = cluster_list_dset[sapply(cluster_list_dset,length) >= min_cells_with_dataset]
  mean_expr = suppressMessages(lapply(cluster_list_dset,build_cell,expr_matrix = neuron_map_seurat@assays$RNA@data))#seurat_object_harmonized
  message("Time used to summarise expression matrix : ",Sys.time()-st)
  datasets_mean_expression_list[[dset]] = as.matrix(do.call(cbind,mean_expr))
}

# init
cor_method="pearson"
mean_lm_res = list()
mean_adj_rsq = numeric()
per_dataset_cor_res=list()
global_cor_res =numeric()
# loop over all clusters
for(i in 1:length(clusters_to_check)){
  current_cluster= clusters_to_check[i]
  #message("  ",current_cluster)
  # get cluster markers
  sn_Seq_based_markers = sn_seq_markers_K169$gene[sn_seq_markers_K169$p_val_adj<padj_cut & sn_seq_markers_K169$avg_log2FC>fc_min & sn_seq_markers_K169$cluster==current_cluster & sn_seq_markers_K169$gene %in% per_gene_cor$gene]
  sc_Seq_based_markers = sc_seq_markers_K169$gene[sc_seq_markers_K169$p_val_adj<padj_cut & sc_seq_markers_K169$avg_logFC>fc_min & sc_seq_markers_K169$cluster==current_cluster & sc_seq_markers_K169$gene %in% per_gene_cor$gene]
  union_markers = union(sn_Seq_based_markers,sc_Seq_based_markers)
  union_markers = sn_Seq_based_markers
  
  ##### RSQ between sn-seq and sc-seq (split up by Dataset)
  # calculate mean expression 
  mean_expressions= do.call(cbind,lapply(datasets_mean_expression_list, function(x,cluster,genes){if(cluster %in% colnames(x)){x[genes,cluster]}},cluster=current_cluster,genes = union_markers))
  mean_expressions= cbind(sn_seq = snseq_mean_expression[union_markers,current_cluster],mean_expressions)
  # filter (TODO: do I want to filter here?)
  # mean_expressions_keep = apply(mean_expressions,1,function(x,min_expr){if(length(x[x>min_expr])>1){return(TRUE)}else{return(FALSE)}},min_expr=min_expr)
  # mean_expressions = mean_expressions[mean_expressions_keep,]
  # construct formula:
  predictors = colnames(mean_expressions)[!colnames(mean_expressions) %in% c("gene","sn_seq")]
  #as.formula(paste0("snseq ~ ",paste0(predictors,collapse=" + ")))
  # fit lm
  if(length(union_markers)>3){
    if(length(predictors)>=1){
      mean_lm_res[[current_cluster]] = lm(formula = as.formula(paste0("sn_seq ~ ",paste0(predictors,collapse=" + "))),data = as.data.frame(mean_expressions))
      mean_adj_rsq[current_cluster] = summary(mean_lm_res[[current_cluster]])$adj.r.squared
      # cor per cluster
      dataset_cor_vec = cor(mean_expressions[,1],mean_expressions[,2:ncol(mean_expressions)],use = "pairwise.complete.obs",method=cor_method)
      if(length(dataset_cor_vec)>1){
        names(dataset_cor_vec)=colnames(mean_expressions[,2:ncol(mean_expressions)])
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
    mean_expr_sc = data.frame(mean_expr = neuron_map_mean_expression[,current_cluster],gene = rownames(neuron_map_mean_expression))
    mean_expr_sn = data.frame(mean_expr = snseq_mean_expression[,current_cluster],gene = rownames(snseq_mean_expression))
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
per_dataset_cor_df$K169_pruned = rownames(per_dataset_cor_df)
# add paired anno
per_dataset_cor_df$K169_named = add_paired_annotation(input_annotation = per_dataset_cor_df$K169_pruned,reference_annotations = neuron_map_seurat@meta.data[,c("K169_pruned","K169_named")])

# save gene cors!
per_dataset_cor_df_print = per_dataset_cor_df[,c(16,15,14,1:13)]
data.table::fwrite(per_dataset_cor_df_print,paste0(results_path,"per_cluster_correlations.txt"),sep = "\t")

##########
### Make heatmap - currently use treeplot below!
##########

## heatmap
# per_dataset_cor_df = per_dataset_cor_df %>% dplyr::select(-mean_adj_rsq )
# per_dataset_cor_long = per_dataset_cor_df %>% dplyr::arrange(desc(global_cor_res)) %>% tidyr::gather("Dataset","RSQ",-K169_pruned,-K169_named) 
# per_dataset_cor_long$Dataset = factor(per_dataset_cor_long$Dataset,levels = c("global_cor","Kim","Campbell","Chen","wenDropSeq","Lee_Idol","Mickelsen","Rossi","Moffit","Flynn","RomanovDev","Mousebrainorg","wen10x","kimDev"))
# 
# per_dataset_cor_long$K169_named = factor(per_dataset_cor_long$K169_named , levels = per_dataset_cor_df$K169_named[order(per_dataset_cor_df$global_cor)])
# 
# # make heatmap 
# require(scales)
# dff_heatmap = ggplot(per_dataset_cor_long, aes(x = K169_named, y = Dataset , fill = RSQ)) +
#   geom_tile() +
#   # facet_grid(sample~., scales="free_y") +
#   #theme(axis.text.y = element_blank())+ 
#   ggplot2::scale_fill_gradient(low="#1a1a82",high="#ff1414",na.value = "white",
#                                limits=c(0,1), oob=squish) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# # and show
# dff_heatmap

##########
### Make Tree heatmap 
##########

# make data for first heatmap with correlations
heatmap_data = per_dataset_cor_df
heatmap_matrix = as.matrix(heatmap_data[,"global_cor"]) # might want to change order
colnames(heatmap_matrix) = "Cor"
rownames(heatmap_matrix) = heatmap_data$K169_pruned
heatmap_matrix2 = as.matrix(heatmap_data[,datasets]) # might want to change order
rownames(heatmap_matrix2) = heatmap_data$K169_pruned

circular_tree_cor = plot_cluster_tree(edgelist = neuron_map_seurat@misc$mrtree_edgelist,
                                      heatmap_matrix=heatmap_matrix,heatmap_matrix2 = heatmap_matrix2,
                                      leaf_level=6,metadata=neuron_map_seurat@meta.data,
                                      label_size = 2, show_genes = TRUE, legend_title_1 = "Cor",legend_title_2 = "Cor",
                                      matrix_offset = 0.2, matrix_width =0.15,matrix_width_2 = 0.5,heatmap_colnames = TRUE,
                                      manual_off_second = 1.5,legend_text_size = 8,heatmap_text_size = 2,colnames_angle=0,hjust_colnames=0.5,
                                      heatmap_colors =c("#1a1a82","#ff1414")) + ggplot2::scale_fill_gradient(low="#1a1a82",high="#ff1414",na.value = "grey90",
                                                                                                             limits=c(0,1), oob=squish)
#circular_tree_cor
require(ggtree)
circular_tree_cor_rotated = rotate_tree(circular_tree_cor, -90)
circular_tree_cor_rotated

#save:
ggsave(filename = paste0(results_path,"circular_tree_correlation.png"),
       plot = circular_tree_cor, "png",dpi=600,width=400,height = 400,units="mm")

ggsave(filename = paste0(results_path,"circular_tree_correlation.pdf"),
       plot = circular_tree_cor, "pdf",dpi=600,width=400,height = 400,units="mm")






