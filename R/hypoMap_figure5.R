
##########
### Load & Prepare
##########


results_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/figure_5/"
system(paste0("mkdir -p ",results_path))

# load everything required
source("load_data.R")
require(mapscvi)
## load mapped object
query_snseq_neurons = readRDS(paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/nucseq_neurons_map.rds"))

query_sn_color = "#302ac9"
fasting_color = "#28a7c9"
adlib_color = "#772ac9"
bg_col = "lightgrey"
colorvec = RColorBrewer::brewer.pal(9, "Blues")
colorvec[1] =  "#dedede"

rasterize_point_size = 2.2
rasterize_pixels = 2048

##########
### Plot fos and ieg
##########

#  feature plot
Idents(query_snseq_neurons) <- "predicted_K31_named"
fos_plot = Seurat::FeaturePlot(query_snseq_neurons,"Fos",split.by = "Diet",label = TRUE,label.size = 3.5,repel = TRUE,pt.size = 0.4,
                               keep.scale="feature",order = TRUE,combine = FALSE)
fos_plot_adlib = fos_plot[[1]]+NoAxes()+theme(panel.border = element_blank())+ scale_color_gradientn(colours = colorvec) +ggplot2::ggtitle("adlib")#+scale_color_gradient(low = "lightgrey",high = "#8c390a")
fos_plot_fasting = fos_plot[[2]]+NoAxes()+theme(panel.border = element_blank())+ scale_color_gradientn(colours = colorvec) +ggplot2::ggtitle("fasted")#+scale_color_gradient(low = "lightgrey",high = "#8c390a")

fos_plot_adlib = rasterize_ggplot(fos_plot_adlib,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
fos_plot_adlib
fos_plot_fasting = rasterize_ggplot(fos_plot_fasting,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
fos_plot_fasting

# save
ggsave(filename = paste0(results_path,"fos_plot_adlib.png"),
       plot = fos_plot_adlib, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"fos_plot_adlib.pdf"),
       plot = fos_plot_adlib, "pdf",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"fos_plot_fasting.png"),
       plot = fos_plot_fasting, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path,"fos_plot_fasting.pdf"),
       plot = fos_plot_fasting, "pdf",dpi=450,width=200,height = 200,units="mm")

##########
### Immediate early genes
##########

# load geneshot based list
immediate_early_genes = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/paper_results/immediate_early_genes.txt",data.table = F)

# get the gene occurence in the dataset ( at least 100 cells to be somewhat relevant )
ieg_expression = FetchData(query_snseq_neurons,vars = immediate_early_genes$mouse_symbol)
ieg_expression[ieg_expression>0] = 1
ieg_occurence = data.frame(occ = colSums(ieg_expression), mouse_symbol = names(colSums(ieg_expression)))
immediate_early_genes = dplyr::left_join(immediate_early_genes,ieg_occurence,by = "mouse_symbol")
# filter to a miniumu occurence (I use 200 for now)
ieg_set = immediate_early_genes$mouse_symbol[immediate_early_genes$occ > 300]
# filter to relevant IEG and add 1700016P03Rik
#ieg_set = c(immediate_early_genes_filtered$mouse_symbol,"1700016P03Rik")# NO!!!!
# remove some gene manually e.g. JUND (reverese effect?), Crebzf (extremely high count of cells) , SH3GL3 (arc interaction !? but high counts)
ieg_set = ieg_set[! ieg_set %in% c("Jund","Crebzf","Sh3gl3")]

## how to define this set best ? Only use Agrp changing genes ? only use gloabl genes ?
# --> less genes means less strongly adjusted pvalues when testing !

## calculate agrp fcs:
Idents(query_snseq_neurons) <- "predicted_K98_pruned"
ieg_in_agrp_fcs = Seurat::FindMarkers(query_snseq_neurons, ident.1 = "adlib",ident.2 = "fast" , group.by = "Diet", 
                                      subset.ident = "K98-4",min.pct = 0,logfc.threshold = 0,features = ieg_set)
ieg_in_agrp_fcs$gene = rownames(ieg_in_agrp_fcs)
#ieg_in_agrp_fcs$padj =  p.adjust(ieg_in_agrp_fcs$p_val, method = "fdr")
# which genes are significant:
activation_genes = ieg_in_agrp_fcs$gene[ieg_in_agrp_fcs$p_val_adj<0.05]
# manually add two other well known genes:
activation_genes = c(activation_genes,"Nr4a3","Egr1")
activation_genes

# example plot
# Seurat::FeaturePlot(query_snseq_neurons,"Btg1",split.by = "Diet",label = TRUE,label.size = 3.5,repel = TRUE,pt.size = 0.4,
#                    keep.scale="feature",order = TRUE)

##########
### Immediate early genes in Agrp neurons - Violin plot
##########

Idents(query_snseq_neurons) <- "predicted_K98_pruned"
p <- Seurat::VlnPlot(query_snseq_neurons,features = activation_genes,split.by = "Diet",idents = c("K98-4"),cols=c(adlib_color,fasting_color),
                     ncol =4,same.y.lims = TRUE,combine = FALSE) 
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(axis.title.x = element_text(size = 10),axis.ticks.x =  element_blank(),axis.text.x=element_blank())+
    NoLegend()+xlab(NULL)+ylab(NULL) #xlab("Agrp neurons")+ylab("")
}
p[[8]] = cowplot::plot_grid(legend)
activation_genes_agrp_vlnPlot = cowplot::plot_grid(plotlist = p,ncol = 4)

# save
ggsave(filename = paste0(results_path,"activation_genes_agrp_vlnPlot.png"),
       plot = activation_genes_agrp_vlnPlot, "png",dpi=450,width=250,height = 200,units="mm")
ggsave(filename = paste0(results_path,"activation_genes_agrp_vlnPlot.pdf"),
       plot = activation_genes_agrp_vlnPlot, "pdf",dpi=450,width=250,height = 200,units="mm")

##########
### Immediate early genes quantification across celtypes
##########

# how to best quantify the response
# fisher summed pvalues ? 
# take mean across cell and make only one test (more elegant regarding pvalues, but what about baseline differences ?)


# calculate marker statistics but require 10% occurence in either group
min_pct = 0.1
source("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/2020-11-jule-single-cell-rnaseq/functions/findDEG.R")
activation_per_cluster = FindAll_DEGs(query_snseq_neurons,group_var = "Diet",idents_name = "predicted_K169_named",max.cells.per.ident = 5000,
                                      features = activation_genes,logfc.threshold=0,pval_filter=Inf,min.pct = min_pct)
# make a manual fdr adjustement
activation_per_cluster$p_val_adjusted = p.adjust(activation_per_cluster$p_val,method = "hochberg") # manual fdr
head(activation_per_cluster)

## summarise:
# define a column for pvalue aggregation
activation_per_cluster$pvalue_for_calc = activation_per_cluster$p_val
pval_ieg_max = 0.05
fc_ieg_min = 0.1
# summarise per cluster
activation_per_cluster_stat = activation_per_cluster %>% group_by(Identity) %>% 
  dplyr::summarise( fisher_sig_pval_down = -2 * sum(log(pvalue_for_calc[pvalue_for_calc<pval_ieg_max & avg_log2FC < (-fc_ieg_min)])),  # summarise tests
                    fisher_sig_pval_up = -2 * sum(log(pvalue_for_calc[pvalue_for_calc<pval_ieg_max & avg_log2FC > fc_ieg_min])), 
                    mean_fc_up = mean(avg_log2FC[pvalue_for_calc<pval_ieg_max & avg_log2FC > fc_ieg_min]),
                    mean_fc_down = mean(avg_log2FC[pvalue_for_calc<pval_ieg_max  & avg_log2FC < (-fc_ieg_min)]),
                    n_down_sig = sum(pvalue_for_calc<pval_ieg_max  & avg_log2FC < (-fc_ieg_min)), # count genes
                    n_up_sig = sum(pvalue_for_calc<pval_ieg_max  & avg_log2FC > fc_ieg_min)) %>%
  dplyr::filter(!is.nan(mean_fc_up) | !is.nan(mean_fc_down)) %>% # general filter for anything significant
  dplyr::filter((n_down_sig >= 1 & fisher_sig_pval_down>20) | (n_up_sig >= 1 & fisher_sig_pval_up>20)) %>% # reduce to a bit more relevant
  dplyr::arrange(desc(fisher_sig_pval_down))
n_cells_per_cluster = query_snseq_neurons@meta.data %>% dplyr::group_by(predicted_K169_named) %>% dplyr::count(name="n_cells_cluster")
activation_per_cluster_stat = dplyr::left_join(activation_per_cluster_stat,n_cells_per_cluster,by=c("Identity"="predicted_K169_named"))
activation_per_cluster_stat

activation_per_cluster_stat$Identity = factor(activation_per_cluster_stat$Identity,levels = rev(activation_per_cluster_stat$Identity))
activation_celltype_barplot = ggplot(activation_per_cluster_stat,aes(x=Identity,y=fisher_sig_pval_down,fill= -1*mean_fc_down))+geom_bar(stat="identity")+
  geom_text(aes(x = Identity, y = max(fisher_sig_pval_down)*1.1,label=n_cells_cluster), hjust =1,size=6)+ # add cell numbers
  scale_fill_gradient(low="grey80",high=fasting_color,limits=c(0,max(abs(activation_per_cluster_stat$mean_fc_down))))+
  ylim(c(0,max(activation_per_cluster_stat$fisher_sig_pval_down)*1.1)) + # extend the y axis a bit so that the n cell labels doN#t overlap with the highest plot!
  coord_flip()+ theme(text = element_text(size=20),axis.title.x = element_text(size = 20))+ylab("Summed adjusted p-value")+xlab(NULL)+ 
  labs(fill='Mean log2-Foldchange')
activation_celltype_barplot

# save
ggsave(filename = paste0(results_path,"activation_celltype_barplot.png"),
       plot = activation_celltype_barplot, "png",dpi=450,width=250,height = 200,units="mm")
ggsave(filename = paste0(results_path,"activation_celltype_barplot.pdf"),
       plot = activation_celltype_barplot, "pdf",dpi=450,width=250,height = 200,units="mm")

# TODO: OR just use table !?
#write.table(activation_per_cluster_stat,file="",quote = FALSE,row.names = FALSE)


##########
### Transcriptional changes in Agrp
##########

# get all DEG of Agrp
conditionGenes_Agrp_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/agrp_fasting_all.txt"
if(!file.exists(conditionGenes_Agrp_file)){
  Idents(query_snseq_neurons) <- "predicted_K98_pruned"
  conditionGenes_Agrp = Seurat::FindMarkers(query_snseq_neurons, ident.1 = "adlib",ident.2 = "fast" , group.by = "Diet", subset.ident = "K98-4",min.pct = 0,logfc.threshold = 0)
  conditionGenes_Agrp$gene = rownames(conditionGenes_Agrp) # add gene name 
  conditionGenes_Agrp$pct_diff = conditionGenes_Agrp$pct.1 - conditionGenes_Agrp$pct.2
  data.table::fwrite(conditionGenes_Agrp,conditionGenes_Agrp_file,sep="\t")
}else{
  # load and save above code which runs for 8 minutes !
  conditionGenes_Agrp = data.table::fread(conditionGenes_Agrp_file,data.table = F)
}

#filter
conditionGenes_Agrp_filtered = conditionGenes_Agrp[conditionGenes_Agrp$p_val_adj < 0.05,] # filter pval
conditionGenes_Agrp_filtered %>% dplyr::arrange((pct_diff),(avg_log2FC))

### make volcano plot
volcano_df = conditionGenes_Agrp
volcano_df$color = "not regulated"
volcano_df$color[volcano_df$avg_log2FC< (-0.35) & volcano_df$p_val_adj < 0.001 ] = "up-regulated in fasting"
volcano_df$color[volcano_df$avg_log2FC > 0.35 & volcano_df$p_val_adj < 0.001 ] = "down-regulated in fasting"

# define gene sets that are labelled:
up_regulated_genes = conditionGenes_Agrp_filtered$gene[conditionGenes_Agrp_filtered$avg_log2FC< (-1) & conditionGenes_Agrp_filtered$p_val_adj < 0.0001 ]
up_regulated_genes = up_regulated_genes[!grepl("Rps|Rpl|Nduf|ENSMUS",up_regulated_genes)]
up_regulated_genes = c(up_regulated_genes,"Fam107b","Lepr","Zbtb16","Tnik","Asic2","Acvr1c","Ncoa2","Vgf")
down_regulated_genes = conditionGenes_Agrp_filtered$gene[conditionGenes_Agrp_filtered$avg_log2FC > 1.25 & conditionGenes_Agrp_filtered$p_val_adj < 0.0001 ]
down_regulated_genes = down_regulated_genes[!grepl("Rps|Rpl|Nduf|ENSMUS",down_regulated_genes)]
down_regulated_genes = c(down_regulated_genes,"Agrp","Sst","Cartpt","Npy")
volcano_df$label = NA
volcano_df$label[volcano_df$gene %in% c(down_regulated_genes,up_regulated_genes)] = volcano_df$gene[volcano_df$gene %in% c(down_regulated_genes,up_regulated_genes)]

agrp_fasting_volcano = ggplot(volcano_df,aes(-1*avg_log2FC,-log10(p_val_adj),label=label,color=color))+geom_point(size=0.3)+
  ggrepel::geom_text_repel(max.overlaps=100,size=6)+
  scale_color_manual(values = c("up-regulated in fasting"=fasting_color,"down-regulated in fasting"=adlib_color,"not regulated"="grey"))+
  theme(text = element_text(size=30))+labs(color='Condition')
agrp_fasting_volcano

#save
ggsave(filename = paste0(results_path,"agrp_fasting_volcano.png"),
       plot = agrp_fasting_volcano, "png",dpi=450,width=250,height = 200,units="mm")
ggsave(filename = paste0(results_path,"agrp_fasting_volcano.pdf"),
       plot = agrp_fasting_volcano, "pdf",dpi=450,width=250,height = 200,units="mm")

##########
### Go enrichment using clusterProfiler
##########

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)
source("utils")

# which gene to enrich for
agrp_upregulated_genes_names = conditionGenes_Agrp_filtered$gene[conditionGenes_Agrp_filtered$avg_log2FC<(-0.35) & conditionGenes_Agrp_filtered$p_val_adj<0.001]

# need to map to entrez
library(biomaRt)
mart <- useMart(dataset="mmusculus_gene_ensembl",biomart='ensembl',host="feb2021.archive.ensembl.org")
up_regulated_genes_ids = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','entrezgene_id'),
                               filters = "external_gene_name",values =agrp_upregulated_genes_names,mart = mart)
# entrez_ids = up_regulated_genes

# need to map background to entrex
cellsh = query_snseq_neurons@meta.data$Cell_ID[query_snseq_neurons@meta.data$predicted_K98_pruned == "K98-4"]
agrp_baseline_expression = get_expression_stats(query_snseq_neurons,cells.1 = cellsh)
agrp_expressed_genes = agrp_baseline_expression$gene[agrp_baseline_expression$mean.1 > 0.02 & agrp_baseline_expression$pct.1 > 0.1]
background_genes_ids = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','entrezgene_id'),
                             filters = "external_gene_name",values = agrp_expressed_genes,mart = mart)

#data.table::fwrite(background_genes_ids,file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/test_objects/agrp_expressed_genes.txt",sep="\t")

# run enrichment on GO BP
agrp_fasting_go_enrichment = clusterProfiler::enrichGO(gene  =  as.character(up_regulated_genes_ids$entrezgene_id%>%na.omit()) ,
                                                       universe      = as.character(background_genes_ids$entrezgene_id%>%na.omit()),
                                                       OrgDb         = org.Mm.eg.db,
                                                       ont           = "BP",
                                                       pAdjustMethod = "BH",
                                                       pvalueCutoff  = 0.01,
                                                       qvalueCutoff  = 0.05,
                                                       readable      = TRUE)
## simplify terms with strict cutoff
agrp_fasting_go_enrichment_simplified = clusterProfiler::simplify(agrp_fasting_go_enrichment,cutoff=0.5)

#get result
agrp_fasting_go_enrichment_simplified_res=agrp_fasting_go_enrichment_simplified@result

# make dotplot
go_bp_agrp_enrich_dotplot = enrichplot::dotplot(agrp_fasting_go_enrichment_simplified, showCategory=20) + 
  ggtitle("Dotplot for ORA")+theme(axis.text.y = element_text(size=20))+
  scale_color_gradient(low=fasting_color,high ="#a1bdc4")
go_bp_agrp_enrich_dotplot

# save
ggsave(filename = paste0(results_path,"go_bp_agrp_enrich_dotplot.png"),
       plot = go_bp_agrp_enrich_dotplot, "png",dpi=450,width=250,height = 200,units="mm")
ggsave(filename = paste0(results_path,"go_bp_agrp_enrich_dotplot.pdf"),
       plot = go_bp_agrp_enrich_dotplot, "pdf",dpi=450,width=250,height = 200,units="mm")

##########
### Transcriptional changes in Agrp of Campbell vs NucSeq
##########

### Comparison with Campbell Agrp changes:

campbell_agrp_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/agrp_fasting_all_campbell.txt"
if(!file.exists(campbell_agrp_file)){
  campbell_diet = subset(neuron_map_seurat,subset = Dataset=="Campbell")
  Idents(campbell_diet) <- "K98_pruned"
  conditionGenes_campbell = Seurat::FindMarkers(campbell_diet, ident.1 = "Normal",ident.2 = "Fasted" , group.by = "Diet", subset.ident = "K98-4",min.pct = 0,logfc.threshold = 0)
  conditionGenes_campbell$gene = rownames(conditionGenes_campbell) # add gene name 
  conditionGenes_campbell$pct_diff = conditionGenes_campbell$pct.1 - conditionGenes_campbell$pct.2
  data.table::fwrite(conditionGenes_campbell,campbell_agrp_file,sep="\t")
}else{
  conditionGenes_campbell = data.table::fread(campbell_agrp_file,data.table = F) 
}

# campbell_diet2 = subset(campbell_diet, subset = Diet %in% c("Normal","Fasted"))
# FeaturePlot(campbell_diet2,features = "Fos",split.by = "Diet")

#### Compare with sc-seq data:
agrp_sn_vs_sc = dplyr::full_join(conditionGenes_Agrp,conditionGenes_campbell,suffix=c("_sn","_sc"),by=c("gene"="gene"))
colnames(agrp_sn_vs_sc)

#### Define categories
sc_p_cut = 0.001
min_fc = 0.2

##### invert FCs and make simple classification:
agrp_sn_vs_sc_plot = agrp_sn_vs_sc

agrp_sn_vs_sc_plot$avg_log2FC_sn = agrp_sn_vs_sc_plot$avg_log2FC_sn*(-1)
agrp_sn_vs_sc_plot$avg_log2FC_sc = agrp_sn_vs_sc_plot$avg_log2FC_sc*(-1)
agrp_sn_vs_sc_plot$regulated = "not in both"
agrp_sn_vs_sc_plot$regulated[ ( agrp_sn_vs_sc_plot$avg_log2FC_sn > min_fc & agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut ) | 
                                (agrp_sn_vs_sc_plot$avg_log2FC_sc  > min_fc & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut)] = "up in fasting - one"

agrp_sn_vs_sc_plot$regulated[agrp_sn_vs_sc_plot$avg_log2FC_sn>min_fc & agrp_sn_vs_sc_plot$avg_log2FC_sc>min_fc & 
                               agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut] = "up in fasting - both"

agrp_sn_vs_sc_plot$regulated[ ( agrp_sn_vs_sc_plot$avg_log2FC_sn< (-1*min_fc) & agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut ) | 
                                (agrp_sn_vs_sc_plot$avg_log2FC_sc  < (-1*min_fc) & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut)] = "down in fasting - one"

agrp_sn_vs_sc_plot$regulated[agrp_sn_vs_sc_plot$avg_log2FC_sn< (-1*min_fc) & agrp_sn_vs_sc_plot$avg_log2FC_sc  < (-1*min_fc) & 
                               agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut] = "down in fasting - both"

agrp_sn_vs_sc_plot$regulated = factor(agrp_sn_vs_sc_plot$regulated,levels = c("not in both","down in fasting - one", "up in fasting - one","down in fasting - both","up in fasting - both"))

agrp_Campvell_vs_sn_plot = ggplot(agrp_sn_vs_sc_plot,aes(x=avg_log2FC_sc,y=avg_log2FC_sn))+
  geom_point(alpha=0.6,aes(color=regulated))+geom_smooth(method="lm",color="grey60")+
  scale_color_manual(values = c("up in fasting - both" = fasting_color,"not in both" = "grey80","down in fasting - both" =adlib_color,
                                "up in fasting - one" = "#a3bbc2", "down in fasting - one" = "#ae9bc2"))+
  xlab("log2FC sc-seq")+ylab("log2FC sn-seq")+theme(text=element_text(size=20))
agrp_Campvell_vs_sn_plot

# save
ggsave(filename = paste0(results_path,"agrp_Campvell_vs_sn_plot.png"),
       plot = agrp_Campvell_vs_sn_plot, "png",dpi=450,width=250,height = 200,units="mm")
ggsave(filename = paste0(results_path,"agrp_Campvell_vs_sn_plot.pdf"),
       plot = agrp_Campvell_vs_sn_plot, "pdf",dpi=450,width=250,height = 200,units="mm")


##########
### Example genes for last section of Figure
##########

# Zbtb16, Fam107b, Vgf, Sv2c, 1700016P03Rik
genes_to_use =c("Zbtb16", "Fam107b", "Vgf", "Sv2c", "1700016P03Rik")

## Option 1: just a nice feature plot of a gene of choice
Idents(query_snseq_neurons) <- "predicted_K31_named"
example_gene_plot = Seurat::FeaturePlot(query_snseq_neurons,"1700016P03Rik",split.by = "Diet",label = TRUE,label.size = 3.5,repel = TRUE,pt.size = 0.4,
                                        keep.scale="feature",order = TRUE,combine = FALSE)
example_gene_plot[[1]] = example_gene_plot[[1]]+NoAxes()+theme(panel.border = element_blank())+ggplot2::ggtitle("adlib")#+scale_color_gradient(low = "lightgrey",high = "#8c390a")
example_gene_plot[[2]] = example_gene_plot[[2]]+NoAxes()+theme(panel.border = element_blank())+ggplot2::ggtitle("fasted")#+scale_color_gradient(low = "lightgrey",high = "#8c390a")
example_gene_plot = cowplot::plot_grid(plotlist = example_gene_plot)
example_gene_plot

## Option 2: some violin plot
# use all other cells and top 5 from enrichement plus pomc plus cck
classes_to_use = c("Pomc.Tcf7l2.Tbx3.HY1","Agrp.Npy.Otp.HY1","Sox14.Lef1.Hmx2.HY1","Hk2.HY14.HY10.HY2","Tmem114.Lef1.Hmx2.HY1")
query_snseq_neurons@meta.data$custom_annotation = "other celltypes"
query_snseq_neurons@meta.data$custom_annotation[query_snseq_neurons@meta.data$predicted_K98_named %in% classes_to_use] = query_snseq_neurons@meta.data$predicted_K98_named[query_snseq_neurons@meta.data$predicted_K98_named %in% classes_to_use]
query_snseq_neurons@meta.data$custom_annotation = factor(query_snseq_neurons@meta.data$custom_annotation,
                                                         levels=rev(c("Agrp.Npy.Otp.HY1","Hk2.HY14.HY10.HY2","Tmem114.Lef1.Hmx2.HY1","Sox14.Lef1.Hmx2.HY1","Pomc.Tcf7l2.Tbx3.HY1","other celltypes")))

# make plot
Idents(query_snseq_neurons) <- "custom_annotation"
interesting_genes_violin_plot = Seurat::VlnPlot(query_snseq_neurons,features = genes_to_use,split.by = "Diet",cols=c(adlib_color,fasting_color),
                                                same.y.lims = TRUE,combine = TRUE,stack = TRUE,adjust = 0.75,flip=TRUE,pt.size = 1) 

#Seurat::DotPlot(query_snseq_neurons,features = genes_to_use,split.by = "Diet",cols=c(adlib_color,fasting_color))

# save
ggsave(filename = paste0(results_path,"interesting_genes_violin_plot.png"),
       plot = interesting_genes_violin_plot, "png",dpi=450,width=250,height = 200,units="mm")
ggsave(filename = paste0(results_path,"interesting_genes_violin_plot.pdf"),
       plot = interesting_genes_violin_plot, "pdf",dpi=450,width=250,height = 200,units="mm")


##########
### Currently not part of the Figure but relevant analysis:

##########
### Transcriptional changes in other activated celltypes
##########

# get all DEG of Hk2
# 	K169-123
conditionGenes_Hk2_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/hk2_fasting_all.txt"
if(!file.exists(conditionGenes_Hk2_file)){
  Idents(query_snseq_neurons) <- "predicted_K169_named"
  conditionGenes_Hk2 = Seurat::FindMarkers(query_snseq_neurons, ident.1 = "adlib",ident.2 = "fast" , group.by = "Diet", subset.ident = "Hk2.HY14.HY10.HY2",min.pct = 0,logfc.threshold = 0)
  conditionGenes_Hk2$gene = rownames(conditionGenes_Hk2) # add gene name 
  conditionGenes_Hk2$pct_diff = conditionGenes_Hk2$pct.1 - conditionGenes_Hk2$pct.2
  data.table::fwrite(conditionGenes_Hk2,conditionGenes_Hk2_file,sep="\t")
}else{
  # load and save above code which runs for 8 minutes !
  conditionGenes_Hk2 = data.table::fread(conditionGenes_Hk2_file,data.table = F)
}
conditionGenes_Hk2_filtered = conditionGenes_Hk2[conditionGenes_Hk2$p_val_adj < 0.05,] # filter pval

# get all DEG of Lef1/Sox14
conditionGenes_Sox14_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/Sox14_fasting_all.txt"
if(!file.exists(conditionGenes_Sox14_file)){
  Idents(query_snseq_neurons) <- "predicted_K169_named"
  conditionGenes_Sox14 = Seurat::FindMarkers(query_snseq_neurons, ident.1 = "adlib",ident.2 = "fast" , group.by = "Diet", subset.ident = "Sox14.Lef1.Hmx2.HY1",min.pct = 0,logfc.threshold = 0)
  conditionGenes_Sox14$gene = rownames(conditionGenes_Sox14) # add gene name 
  conditionGenes_Sox14$pct_diff = conditionGenes_Sox14$pct.1 - conditionGenes_Sox14$pct.2
  data.table::fwrite(conditionGenes_Sox14,conditionGenes_Sox14_file,sep="\t")
}else{
  # load and save above code which runs for 8 minutes !
  conditionGenes_Sox14 = data.table::fread(conditionGenes_Sox14_file,data.table = F)
}
conditionGenes_Sox14_filtered = conditionGenes_Sox14[conditionGenes_Sox14$p_val_adj < 0.05,] # filter pval

##########
### Comparison with global
##########

## global changes
# load global changes
# get all DEG of Hk2
conditionGenes_global_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/mapdata/global_fasting_all.txt"
if(!file.exists(conditionGenes_global_file)){
  conditionGenes_global = Seurat::FindMarkers(query_snseq_neurons, ident.1 = "adlib",ident.2 = "fast" , group.by = "Diet",min.pct = 0,logfc.threshold = 0.1,max.cells.per.ident = 10000)
  conditionGenes_global$gene = rownames(conditionGenes_global) # add gene name 
  conditionGenes_global$pct_diff = conditionGenes_global$pct.1 - conditionGenes_global$pct.2
  data.table::fwrite(conditionGenes_global,conditionGenes_global_file,sep="\t")
}else{
  # load and save above code which runs for 8 minutes !
  conditionGenes_global = data.table::fread(conditionGenes_global_file,data.table = FALSE)
}
conditionGenes_global_filtered = conditionGenes_global[conditionGenes_global$p_val_adj < 0.05,] # filter pval


# flag genes that not changing globally
agrp_vs_global = dplyr::full_join(conditionGenes_Agrp,conditionGenes_global,by=c("gene"="gene"),suffix=c("_agrp","_global"))
agrp_vs_global$avg_log2FC_global[is.na(agrp_vs_global$avg_log2FC_global)] = 0 # set na to 0 (was below 0.1)
agrp_vs_global$p_val_adj_global[is.na(agrp_vs_global$p_val_adj_global)] = 1 # set na to 1
agrp_vs_global$fc_diff = agrp_vs_global$avg_log2FC_agrp - agrp_vs_global$avg_log2FC_global #
agrp_vs_global$significant_agrp = FALSE
agrp_vs_global$significant_agrp[agrp_vs_global$p_val_adj_agrp < 0.001] = TRUE
agrp_vs_global$label = NA
agrp_vs_global$label[abs(agrp_vs_global$fc_diff) > 0.75 & (agrp_vs_global$p_val_adj_agrp < 0.001 | agrp_vs_global$p_val_adj_global < 0.001)] =
  agrp_vs_global$gene[abs(agrp_vs_global$fc_diff) > 0.75 & (agrp_vs_global$p_val_adj_agrp < 0.001 | agrp_vs_global$p_val_adj_global < 0.001)]

ggplot2::ggplot(agrp_vs_global,aes(-1*avg_log2FC_global,-1*avg_log2FC_agrp,color=significant_agrp,label=label))+geom_point(size=0.3)+
  ggrepel::geom_text_repel(max.overlaps  = 100,size=5)+scale_color_manual(values = c("TRUE"=fasting_color,"FALSE"=adlib_color))+theme(text=element_text(size=20))


# flag genes that not changing globally
sox14_vs_global = dplyr::full_join(conditionGenes_Sox14,conditionGenes_global,by=c("gene"="gene"),suffix=c("_sox14","_global"))
sox14_vs_global$avg_log2FC_global[is.na(sox14_vs_global$avg_log2FC_global)] = 0 # set na to 0 (was below 0.1)
sox14_vs_global$p_val_adj_global[is.na(sox14_vs_global$p_val_adj_global)] = 1 # set na to 1
sox14_vs_global$fc_diff = sox14_vs_global$avg_log2FC_sox14 - sox14_vs_global$avg_log2FC_global #
sox14_vs_global$significant_sox14 = FALSE
sox14_vs_global$significant_sox14[sox14_vs_global$p_val_adj_sox14 < 0.001] = TRUE
sox14_vs_global$label = NA
sox14_vs_global$label[abs(sox14_vs_global$fc_diff) > 0.75 & (sox14_vs_global$p_val_adj_sox14 < 0.001 | sox14_vs_global$p_val_adj_global < 0.001)] = sox14_vs_global$gene[abs(sox14_vs_global$fc_diff) > 0.75 & (sox14_vs_global$p_val_adj_sox14 < 0.001 | sox14_vs_global$p_val_adj_global < 0.001)]
sox14_vs_global$label[sox14_vs_global$gene %in% c("Sox14","Lef1")] = sox14_vs_global$gene[sox14_vs_global$gene %in% c("Sox14","Lef1")]

ggplot2::ggplot(sox14_vs_global,aes(-1*avg_log2FC_global,-1*avg_log2FC_sox14,color=significant_sox14,label=label))+geom_point(size=0.3)+
  ggrepel::geom_text_repel(max.overlaps  = 100,size=5)+scale_color_manual(values = c("TRUE"=fasting_color,"FALSE"=adlib_color))+theme(text=element_text(size=20))

# flag genes that not changing globally
Hk2_vs_global = dplyr::full_join(conditionGenes_Hk2,conditionGenes_global,by=c("gene"="gene"),suffix=c("_Hk2","_global"))
Hk2_vs_global$avg_log2FC_global[is.na(Hk2_vs_global$avg_log2FC_global)] = 0 # set na to 0 (was below 0.1)
Hk2_vs_global$p_val_adj_global[is.na(Hk2_vs_global$p_val_adj_global)] = 1 # set na to 0 (was below 0.1)
Hk2_vs_global$fc_diff = Hk2_vs_global$avg_log2FC_Hk2 - Hk2_vs_global$avg_log2FC_global #
Hk2_vs_global$significant_Hk2 = FALSE
Hk2_vs_global$significant_Hk2[Hk2_vs_global$p_val_adj_Hk2 < 0.001] = TRUE
Hk2_vs_global$label = NA
Hk2_vs_global$label[abs(Hk2_vs_global$fc_diff) > 0.75 & (Hk2_vs_global$p_val_adj_Hk2 < 0.001 | Hk2_vs_global$p_val_adj_global < 0.001)] = Hk2_vs_global$gene[abs(Hk2_vs_global$fc_diff) > 0.75 & (Hk2_vs_global$p_val_adj_Hk2 < 0.001 | Hk2_vs_global$p_val_adj_global < 0.001)]

ggplot2::ggplot(Hk2_vs_global,aes(-1*avg_log2FC_global,-1*avg_log2FC_Hk2,color=significant_Hk2,label=label))+geom_point(size=0.3)+
  ggrepel::geom_text_repel(max.overlaps  = 100)+scale_color_manual(values = c("TRUE"=fasting_color,"FALSE"=adlib_color))

##########
### Which genes are shared between celltypes
##########

intersect(conditionGenes_Hk2_filtered$gene,conditionGenes_Agrp_filtered$gene)
intersect(conditionGenes_Sox14_filtered$gene,conditionGenes_Agrp_filtered$gene)

conditionGenes_Hk2


