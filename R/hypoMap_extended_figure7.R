##########
### Load & Prepare
##########

results_path_extended_figure7 = "figure_outputs/figure_extended_7/"
system(paste0("mkdir -p ",results_path_extended_figure7))


# load required functions
require(dplyr)
require(ggplot2)a
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
### Make subset for Campbell
##########

rasterize_point_size_inc = 5.4

# subset to dowsett data:
campbell_subset = hypoMap_v2_seurat@meta.data[hypoMap_v2_seurat@meta.data$Dataset=="CampbellDropseq","Cell_ID"]
campbell_subset = subset(hypoMap_v2_seurat,cells=campbell_subset)

#DimPlot(campbell_subset,group.by = "C286",label=TRUE,label.size = 2)+NoLegend()+NoAxes()

# make agrp subset
campbell_subset_agrp = subset(campbell_subset,subset = C66 == "C66-46" & umapscvi_1 < -2 & umapscvi_2 < 9)

##########
### Plot fos and ieg
##########

rasterize_point_size_inc = 4.4
cols_for_feature_plot_edit = c("grey80",cols_for_feature_plot[2])

# subplot A: highlight the agrp cluster on UMAP of dowsett data
cellsh = hypoMap_v2_seurat@meta.data$Cell_ID[hypoMap_v2_seurat@meta.data$C66=="C66-46"]
subplot_A = DimPlot(campbell_subset, cells.highlight = cellsh,sizes.highlight = 0.1,cols.highlight = "#D55E00",na.value = bg_col,label = F)+NoLegend()+NoAxes()
subplot_A = rasterize_ggplot(subplot_A,pixel_raster = rasterize_pixels,pointsize = 2.7)
subplot_A

# subplot B: show agrp expression in Agrp
subplot_B = FeaturePlot(campbell_subset_agrp,"Agrp",order = TRUE,cols = cols_for_feature_plot_edit)+NoAxes()
subplot_B = rasterize_ggplot(subplot_B,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size_inc)
subplot_B

# subplot c/d: fos in NCD / Fasting
subplot_CD = FeaturePlot(campbell_subset_agrp,"Fos",split.by = "Diet",keep.scale = "all",pt.size = 1,order = TRUE,combine = FALSE,cols = cols_for_feature_plot_edit)
fos_plot_adlib = subplot_CD[[3]]+NoAxes()+theme(panel.border = element_blank()) #+ggplot2::ggtitle("adlib") #+ scale_color_gradientn(colours = colorvec) #+scale_color_gradient(low = "lightgrey",high = "#8c390a")
fos_plot_fasting = subplot_CD[[1]]+NoAxes()+theme(panel.border = element_blank()) #+ggplot2::ggtitle("fasted")#+ scale_color_gradientn(colours = colorvec) #+scale_color_gradient(low = "lightgrey",high = "#8c390a")
fos_plot_adlib = rasterize_ggplot(fos_plot_adlib,pixel_raster = rasterize_pixels,pointsize = 6.5)
fos_plot_adlib
fos_plot_fasting = rasterize_ggplot(fos_plot_fasting,pixel_raster = rasterize_pixels,pointsize = 6.5)
fos_plot_fasting

# save
#a
ggsave(filename = paste0(results_path_extended_figure7,"agrp_in_campbell_umap.png"),
       plot = subplot_A, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure7,"agrp_in_campbell_umap.pdf"),
       plot = subplot_A, "pdf",dpi=450,width=200,height = 200,units="mm")
#b
ggsave(filename = paste0(results_path_extended_figure7,"agrp_in_campbell_subset.pdf"),
       plot = subplot_B, "pdf",dpi=450,width=200,height = 200,units="mm")
# cd
ggsave(filename = paste0(results_path_extended_figure7,"fos_plot_adlib_campbell.png"),
       plot = fos_plot_adlib, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure7,"fos_plot_adlib_campbell.pdf"),
       plot = fos_plot_adlib, "pdf",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure7,"fos_plot_fasting_campbell.png"),
       plot = fos_plot_fasting, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure7,"fos_plot_fasting_campbell.pdf"),
       plot = fos_plot_fasting, "pdf",dpi=450,width=200,height = 200,units="mm")

## source data
source_ext_figure7_a_umap_all = subplot_A$data
source_ext_figure7_a_umap_all$Cell_ID = rownames(source_ext_figure7_a_umap_all)
data.table::fwrite(source_ext_figure7_a_umap_all,paste0(results_path_extended_figure7,"source_ext_figure7_a_umap_all.txt"),sep="\t")

source_ext_figure7_a_umaps_small = subplot_B$data %>% dplyr::mutate(Cell_ID = rownames(subplot_B$data )) %>% 
  dplyr::left_join(fos_plot_adlib$data  %>% dplyr::mutate(Cell_ID = rownames(fos_plot_adlib$data )) %>% dplyr::select(Cell_ID, Fos_adlib =Fos) %>% dplyr::mutate(Fos_adlib = as.numeric(Fos_adlib)) ,by="Cell_ID" ) %>% 
  dplyr::left_join(fos_plot_fasting$data  %>% dplyr::mutate(Cell_ID = rownames(fos_plot_fasting$data )) %>% dplyr::select(Cell_ID, Fos_fasting =Fos) %>% dplyr::mutate(Fos_fasting = as.numeric(Fos_fasting)) ,by="Cell_ID" ) 

data.table::fwrite(source_ext_figure7_a_umaps_small,paste0(results_path_extended_figure7,"source_ext_figure7_a_umaps_small.txt"),sep="\t")


##########
###  Make dowsett subset
##########

# subset to dowsett data:
dowsett_subset = hypoMap_v2_seurat@meta.data[hypoMap_v2_seurat@meta.data$Dataset=="Dowsett10xnuc","Cell_ID"]
dowsett_subset = subset(hypoMap_v2_seurat,cells=dowsett_subset)

# modify agrp
agrp_cluster_to_use = "C286-178"
dowsett_subset@meta.data$C286_modified = dowsett_subset@meta.data$C286
dowsett_subset@meta.data$C286_modified[dowsett_subset@meta.data$C286_modified %in% c("C286-176","C286-177","C286-178")] = agrp_cluster_to_use

#DimPlot(dowsett_subset,group.by = "C286",label=TRUE,label.size = 2)+NoLegend()+NoAxes()

##########
###  DEGs in all clusters
##########

# !!! run some time

min_cells = 10 # ?
target_cluster = "C286_modified"
Idents(dowsett_subset) = target_cluster
all_clusters = unique(dowsett_subset@meta.data[,target_cluster])
all_conditionGenes_list = list()
for(i in 1:length(all_clusters)){
  current_cluster = all_clusters[i]
  message(current_cluster)
  cells_fasting= length(dowsett_subset@meta.data[dowsett_subset@meta.data[,target_cluster] == current_cluster & dowsett_subset@meta.data$Diet=="Fasted",target_cluster])
  cells_adlib = length(dowsett_subset@meta.data[dowsett_subset@meta.data[,target_cluster] == current_cluster & dowsett_subset@meta.data$Diet=="Normal chow",target_cluster])
  if(cells_adlib >= min_cells & cells_fasting >= min_cells ){
    conditionGenes_current = Seurat::FindMarkers(dowsett_subset, ident.1 = "Fasted",ident.2 = "Normal chow" , group.by = "Diet", 
                                                 subset.ident = current_cluster,min.pct = 0.1,logfc.threshold = 0.25,max.cells.per.ident = 10000)
    conditionGenes_current$gene = rownames(conditionGenes_current) # add gene name 
    conditionGenes_current$pct_diff = conditionGenes_current$pct.1 - conditionGenes_current$pct.2
    conditionGenes_current$current_cluster = current_cluster
    all_conditionGenes_list[[current_cluster]] = conditionGenes_current
  }
}

# rbind
all_conditionGenes = do.call(rbind,all_conditionGenes_list)

all_conditionGenes_filtered = all_conditionGenes[all_conditionGenes$p_val_adj < 0.01,] # filter pval

##########
###  load negbinom version (ran via slurm)
##########


conditionGenes_all_file = paste0(results_path_extended_figure7,"dowsett_all_clusters_fasting_DEG_","negbinom",".txt")
all_conditionGenes_negbinom = data.table::fread(conditionGenes_all_file,data.table = F)


##########
###  save
##########

# save !
conditionGenes_all_file = paste0(results_path_extended_figure7,"all_clusters_fasting_DEG.txt")

# data.table::fwrite(all_conditionGenes_filtered,conditionGenes_all_file,sep="\t") # only save filtered version
#all_conditionGenes_filtered = data.table::fread(conditionGenes_all_file,data.table = F)

# join
all_conditionGenes_filtered_with_negbinom = dplyr::left_join(all_conditionGenes_filtered,
                                                    all_conditionGenes_negbinom %>% dplyr::select(p_val_adj_negbinom = p_val_adj,current_cluster,gene),by=c("gene"="gene","current_cluster"="current_cluster"))

data.table::fwrite(all_conditionGenes_filtered_with_negbinom,conditionGenes_all_file,sep="\t") # only save filtered version


##########
### plot on UMAP
##########

conditionGenes_all_file = paste0(results_path_extended_figure7,"all_clusters_fasting_DEG.txt")
all_conditionGenes = data.table::fread(conditionGenes_all_file,data.table = FALSE)

all_conditionGenes_filtered = all_conditionGenes[all_conditionGenes$p_val_adj < 0.01,] # filter pval

# how many express Fos:
fos_degs = all_conditionGenes[all_conditionGenes$gene == "Fos",]

target_col = "C286_modified"
n_cells_per_cluster = dowsett_subset@meta.data %>% dplyr::group_by(!!sym(target_col)) %>% dplyr::count(name="n_cells_per_cluster")
n_DEG_per_cluster = all_conditionGenes_filtered %>% dplyr::group_by(current_cluster) %>% dplyr::count(name="n_DEG_per_cluster") %>% 
  dplyr::full_join(n_cells_per_cluster,by=c("current_cluster"="C286_modified"))
n_DEG_per_cluster$n_DEG_per_cluster[is.na(n_DEG_per_cluster$n_DEG_per_cluster)] = 0

#ggplot(n_DEG_per_cluster,aes(n_DEG_per_cluster,n_cells_per_cluster))+geom_point()

n_DEG_per_cluster$deg_score = n_DEG_per_cluster$n_DEG_per_cluster / (n_DEG_per_cluster$n_cells_per_cluster - 0)

#dowsett_subset@meta.data = dowsett_subset@meta.data[,1:45]
temp_meta = dplyr::left_join(dowsett_subset@meta.data,n_DEG_per_cluster,by=c("C286_modified"="current_cluster"))
rownames(temp_meta) = temp_meta$Cell_ID
dowsett_subset@meta.data = temp_meta
dowsett_subset@meta.data$n_DEG_per_cluster_zscore = (dowsett_subset@meta.data$n_DEG_per_cluster - mean(dowsett_subset@meta.data$n_DEG_per_cluster)) / sd(dowsett_subset@meta.data$n_DEG_per_cluster)

# ### Glp1r expression in neuron map
# ndeg_sn_seq = FeaturePlot(dowsett_subset,features = "deg_score",reduction =paste0("umap_scvi"),cols = cols_for_feature_plot,order = TRUE,pt.size = 0.2)+
#   NoAxes() 
# ndeg_sn_seq = rasterize_ggplot(ndeg_sn_seq,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
# ndeg_sn_seq
# 
# ## plot z_score version (but I think above plot is better for figure)
# ndeg_sn_seq = FeaturePlot(dowsett_subset,features = "n_DEG_per_cluster_zscore",reduction =paste0("umap_scvi"),cols = cols_for_feature_plot,order = TRUE,pt.size = 0.2)+
#   NoAxes() 
# ndeg_sn_seq = rasterize_ggplot(ndeg_sn_seq,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
# ndeg_sn_seq

ndeg_sn_seq = FeaturePlot(dowsett_subset,features = "n_DEG_per_cluster",reduction =paste0("umap_scvi"),cols = cols_for_feature_plot,order = TRUE,pt.size = 0.2)+
  NoAxes() + theme(text = element_text(size=text_size))
ndeg_sn_seq = rasterize_ggplot(ndeg_sn_seq,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
ndeg_sn_seq

# save
ggsave(filename = paste0(results_path_extended_figure7,"snseq_n_deg_fasting_umap.png"),
       plot = ndeg_sn_seq, "png",dpi=400,width=330,height = 300,units="mm")
ggsave(filename = paste0(results_path_extended_figure7,"snseq_n_deg_fasting_umap.pdf"),
       plot = ndeg_sn_seq, "pdf",dpi=400,width=330,height =300,units="mm")

extended_figure7_c_data =ndeg_sn_seq$data
data.table::fwrite(extended_figure7_c_data,paste0(results_path_extended_figure7,"source_ext_figure7_c_DEGs.txt"),sep="\t")

##########
### Creb1 anaylsis
##########

# please see the additional script "tfbs_analysis.R"

extended_figure7_c_data =ndeg_sn_seq$data
data.table::fwrite(extended_figure7_c_data,paste0(results_path_extended_figure7,"source_ext_figure7_c_DEGs.txt"),sep="\t")


##########
### Overlapping genes
##########

# Calculate overlaps

n_DEG_per_cluster = all_conditionGenes_filtered %>% dplyr::group_by(current_cluster) %>% dplyr::count(name="n_DEG_per_cluster") %>% 
  dplyr::full_join(n_cells_per_cluster,by=c("current_cluster"="C286_modified"))

number_of_deg_clusters = length(n_DEG_per_cluster$current_cluster[n_DEG_per_cluster$n_DEG_per_cluster > 1] %>% na.omit())

overlap_DEG_per_cluster = all_conditionGenes_filtered %>% dplyr::group_by(gene) %>% dplyr::count(name="n_occ_DEG")  %>%
  dplyr::mutate(pct_occ = n_occ_DEG / number_of_deg_clusters)

#hist(overlap_DEG_per_cluster$pct_occ)

frequent_genes = overlap_DEG_per_cluster$gene[overlap_DEG_per_cluster$pct_occ > (1/5)]

# number of DEGs:
length(unique(overlap_DEG_per_cluster$gene))
# number of 50 % overlapping genes;
length(overlap_DEG_per_cluster$gene[overlap_DEG_per_cluster$pct_occ > (1/2)])
# number of 20 % overlapping genes;
length(overlap_DEG_per_cluster$gene[overlap_DEG_per_cluster$pct_occ > (1/5)])


##########
### Go enrichment using clusterProfiler
##########

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

# need to map to entrez
library(biomaRt)
mart <- useMart(dataset="mmusculus_gene_ensembl",biomart='ensembl',host="nov2020.archive.ensembl.org") # use 2020 release because better compatbile with most datasets  #or: feb2021.archive.ensembl.org
frequent_genes_ids = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','entrezgene_id'),
                               filters = "external_gene_name",values =frequent_genes,mart = mart)
# need to map background to entrex --> use all DEGs (1800)
background_genes_ids = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','entrezgene_id'),
                             filters = "external_gene_name",values = unique(overlap_DEG_per_cluster$gene),mart = mart)

# run enrichment on GO BP
overlap_genes_go_enrichment = clusterProfiler::enrichGO(gene  =  as.character(frequent_genes_ids$entrezgene_id%>%na.omit()) ,
                                                       universe      = as.character(background_genes_ids$entrezgene_id%>%na.omit()),
                                                       OrgDb         = org.Mm.eg.db,
                                                       ont           = "BP",
                                                       pAdjustMethod = "BH",
                                                       pvalueCutoff  = 0.01,
                                                       qvalueCutoff  = 0.05,
                                                       readable      = TRUE)
## simplify terms with strict cutoff
overlap_genes_go_enrichment_simplified = clusterProfiler::simplify(overlap_genes_go_enrichment,cutoff=0.5)

#get result
overlap_genes_go_enrichment_simplified_res=overlap_genes_go_enrichment_simplified@result

# save table
data.table::fwrite(overlap_genes_go_enrichment_simplified_res,file = paste0(results_path_extended_figure7,"overlap_genes_go_enrichment_simplified.txt"),sep="\t")

# make dotplot
go_bp_overlap_dotplot = enrichplot::dotplot(overlap_genes_go_enrichment_simplified, showCategory=20) + 
  ggtitle("Dotplot for ORA")+theme(axis.text.y = element_text(size=15))+
  scale_color_gradient(low=fasting_color,high ="#a1bdc4") # TODO: change color ?
go_bp_overlap_dotplot

# save
ggsave(filename = paste0(results_path_extended_figure7,"go_bp_overlap_enrich_dotplot.png"),
       plot = go_bp_overlap_dotplot, "png",dpi=450,width=250,height = 200,units="mm")
ggsave(filename = paste0(results_path_extended_figure7,"go_bp_overlap_enrich_dotplot.pdf"),
       plot = go_bp_overlap_dotplot, "pdf",dpi=450,width=250,height = 200,units="mm")

#source
source_ext_figure7_d_dotplot =go_bp_overlap_dotplot$data
data.table::fwrite(source_ext_figure7_d_dotplot ,paste0(results_path_extended_figure7,"source_ext_figure7_d_enrichment.txt"),sep="\t")

