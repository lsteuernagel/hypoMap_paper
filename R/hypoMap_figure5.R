
##########
### Load & Prepare
##########


results_path_figure5 = "figure_outputs/figure_5/"
system(paste0("mkdir -p ",results_path_figure5))

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
### Make subset objects for plots in this Figure !
##########

rasterize_point_size_inc = 5.4

# subset to dowsett data:
dowsett_subset = hypoMap_v2_seurat@meta.data[hypoMap_v2_seurat@meta.data$Dataset=="Dowsett10xnuc","Cell_ID"]
dowsett_subset = subset(hypoMap_v2_seurat,cells=dowsett_subset)

#DimPlot(dowsett_subset,group.by = "C286",label=TRUE,label.size = 2)+NoLegend()+NoAxes()

# make agrp subset
dowsett_subset_agrp = subset(dowsett_subset,subset = C66 == "C66-46" & umapscvi_1 < -2 & umapscvi_2 < 9)

##########
### Plot fos and ieg
##########

rasterize_point_size_inc = 5.4
cols_for_feature_plot_edit = c("grey80",cols_for_feature_plot[2])

# subplot A: highlight the agrp cluster on UMAP of dowsett data
cellsh = hypoMap_v2_seurat@meta.data$Cell_ID[hypoMap_v2_seurat@meta.data$C66=="C66-46"]
subplot_A = DimPlot(dowsett_subset, cells.highlight = cellsh,sizes.highlight = 0.1,cols.highlight = "#D55E00",na.value = bg_col,label = F)+NoLegend()+NoAxes()
subplot_A = rasterize_ggplot(subplot_A,pixel_raster = rasterize_pixels,pointsize = 3.2)
subplot_A

# subplot B: show agrp expression in Agrp
subplot_B = FeaturePlot(dowsett_subset_agrp,"Agrp",order = TRUE,cols = cols_for_feature_plot_edit)+NoAxes()
subplot_B = rasterize_ggplot(subplot_B,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size_inc)
subplot_B

# subplot c/d: fos in NCD / Fasting
subplot_CD = FeaturePlot(dowsett_subset_agrp,"Fos",split.by = "Diet",keep.scale = "all",pt.size = 1,order = TRUE,combine = FALSE,cols = cols_for_feature_plot_edit)
fos_plot_adlib = subplot_CD[[2]]+NoAxes()+theme(panel.border = element_blank()) #+ggplot2::ggtitle("adlib") #+ scale_color_gradientn(colours = colorvec) #+scale_color_gradient(low = "lightgrey",high = "#8c390a")
fos_plot_fasting = subplot_CD[[1]]+NoAxes()+theme(panel.border = element_blank())# +ggplot2::ggtitle("fasted")#+ scale_color_gradientn(colours = colorvec) #+scale_color_gradient(low = "lightgrey",high = "#8c390a")
fos_plot_adlib = rasterize_ggplot(fos_plot_adlib,pixel_raster = rasterize_pixels,pointsize = 10)
fos_plot_adlib
fos_plot_fasting = rasterize_ggplot(fos_plot_fasting,pixel_raster = rasterize_pixels,pointsize = 10)
fos_plot_fasting

# save
#a
ggsave(filename = paste0(results_path_figure5,"agrp_in_dowsett_umap.png"),
       plot = subplot_A, "png",dpi=450,width=200,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"agrp_in_dowsett_umap.pdf"),
       plot = subplot_A, "pdf",dpi=450,width=200,height = 200,units="mm")
#b
ggsave(filename = paste0(results_path_figure5,"agrp_in_subset.pdf"),
       plot = subplot_B, "pdf",dpi=450,width=220,height = 200,units="mm")
# cd
ggsave(filename = paste0(results_path_figure5,"fos_plot_adlib.png"),
       plot = fos_plot_adlib, "png",dpi=450,width=220,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"fos_plot_adlib.pdf"),
       plot = fos_plot_adlib, "pdf",dpi=450,width=220,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"fos_plot_fasting.png"),
       plot = fos_plot_fasting, "png",dpi=450,width=220,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"fos_plot_fasting.pdf"),
       plot = fos_plot_fasting, "pdf",dpi=450,width=220,height = 200,units="mm")

## source data
source_figure5_a_umap_all = subplot_A$data
source_figure5_a_umap_all$Cell_ID = rownames(source_figure5_a_umap_all)
data.table::fwrite(source_figure5_a_umap_all,paste0(results_path_figure5,"source_figure5_a_umap_all.txt"),sep="\t")

source_figure5_a_umaps_small = subplot_B$data %>% dplyr::mutate(Cell_ID = rownames(subplot_B$data )) %>% 
  dplyr::left_join(fos_plot_adlib$data  %>% dplyr::mutate(Cell_ID = rownames(fos_plot_adlib$data )) %>% dplyr::select(Cell_ID, Fos_adlib =Fos) %>% dplyr::mutate(Fos_adlib = as.numeric(Fos_adlib)) ,by="Cell_ID" ) %>% 
  dplyr::left_join(fos_plot_fasting$data  %>% dplyr::mutate(Cell_ID = rownames(fos_plot_fasting$data )) %>% dplyr::select(Cell_ID, Fos_fasting =Fos) %>% dplyr::mutate(Fos_fasting = as.numeric(Fos_fasting)) ,by="Cell_ID" ) 

data.table::fwrite(source_figure5_a_umaps_small,paste0(results_path_figure5,"source_figure5_a_umaps_small.txt"),sep="\t")


##########
### Immediate early genes
##########

# load geneshot based list
immediate_early_genes = data.table::fread(paste0("data_inputs/immediate_early_genes_wuetal.txt"),data.table = F)$gene

# get the gene occurence in the dataset ( at least 100 cells to be somewhat relevant )
ieg_expression = Seurat::FetchData(dowsett_subset,vars = immediate_early_genes)
ieg_expression[ieg_expression>0] = 1
ieg_occurence = data.frame(occ = colSums(ieg_expression), mouse_symbol = names(colSums(ieg_expression)))
# filter to a miniumu occurence (I use 200 for now)
ieg_set = ieg_occurence$mouse_symbol[ieg_occurence$occ > 300 & ieg_occurence$occ < 10000 ] # 
# filter to relevant IEG and add 1700016P03Rik
#ieg_set = c(immediate_early_genes_filtered$mouse_symbol,"1700016P03Rik")# NO!!!!

# remove some gene manually e.g. JUND (reverese effect?),
#ieg_set = ieg_set[! ieg_set %in% c("Jund","Crebzf","Sh3gl3")]

# remove some that are rather not true IEGs
#ieg_set = ieg_set[! ieg_set %in% c("Sgk1","Sh3gl3","Cxcl2")]
ieg_set = c(ieg_set,"1700016P03Rik")

# example plot
# Seurat::FeaturePlot(dowsett_subset,"Btg1",split.by = "Diet",label = TRUE,label.size = 3.5,repel = TRUE,pt.size = 0.4,
#                    keep.scale="feature",order = TRUE)


##########
### Immediate early genes quantification across celltypes
##########

# how to best quantify the response: use fisher summed pvalues ? or just count ?

# calculate marker statistics but require 10% occurence in either group
min_pct = 0.1
min_cells = 10

# we already know that Agrp nueornd are tricky!!!!
##.> see simple solution below
agrp_cluster_to_use = "C286-178"
dowsett_subset@meta.data$C286_modified = dowsett_subset@meta.data$C286
dowsett_subset@meta.data$C286_modified[dowsett_subset@meta.data$C286_modified %in% c("C286-176","C286-177","C286-178")] = agrp_cluster_to_use

### run
Idents(dowsett_subset) = "C286_modified"
all_clusters = unique(dowsett_subset@meta.data$C286_modified)
all_activationGenes_list = list()
all_activationGenes_negbinom_list = list()
for(i in 1:length(all_clusters)){
  current_cluster = all_clusters[i]
  message(current_cluster)
  cells_adlib = length(dowsett_subset@meta.data$C286_modified[dowsett_subset@meta.data$C286_modified == current_cluster & dowsett_subset@meta.data$Diet=="Normal chow"])
  cells_fasting = length(dowsett_subset@meta.data$C286_modified[dowsett_subset@meta.data$C286_modified == current_cluster & dowsett_subset@meta.data$Diet=="Fasted"])
  if(cells_adlib >= min_cells & cells_fasting >= min_cells ){
    activationGenes_current = Seurat::FindMarkers(dowsett_subset, ident.1 = "Fasted",ident.2 = "Normal chow" , group.by = "Diet",features = ieg_set,
                                                  subset.ident = current_cluster,min.pct = min_pct,logfc.threshold = 0.1,max.cells.per.ident = 5000)
    activationGenes_current_negbinom = Seurat::FindMarkers(dowsett_subset, ident.1 = "Normal chow",ident.2 = "Fasted" , group.by = "Diet",features = ieg_set,test.use="negbinom",
                                                           subset.ident = current_cluster,min.pct = min_pct,logfc.threshold = 0.1,max.cells.per.ident = 5000)
    #TODO: insert test_subset_with_negbinom function here!!!!!
    #activationGenes_current_negbinom = test_subset_with_negbinom(seurat_object=dowsett_subset, features_to_subset=ieg_set,cells.1=cells_fasting,cells.2=cells_adlib)
    if(nrow(activationGenes_current)>0){
      activationGenes_current$gene = rownames(activationGenes_current) # add gene name 
      activationGenes_current$pct_diff = activationGenes_current$pct.1 - activationGenes_current$pct.2
      activationGenes_current$p_val_fdr = p.adjust(activationGenes_current$p_val,method = "hochberg") # manual fdr
      activationGenes_current$p_val_bonferroni = p.adjust(activationGenes_current$p_val,method = "bonferroni") # manual fdr
      activationGenes_current$current_cluster = current_cluster
      all_activationGenes_list[[current_cluster]] = activationGenes_current
    }
    if(nrow(activationGenes_current_negbinom)>0){
      activationGenes_current_negbinom$gene = rownames(activationGenes_current_negbinom) # add gene name
      # activationGenes_current_negbinom$pct_diff = activationGenes_current_negbinom$pct.1 - activationGenes_current_negbinom$pct.2
      activationGenes_current_negbinom$current_cluster = current_cluster
      activationGenes_current_negbinom$p_val_bonferroni = p.adjust(activationGenes_current_negbinom$p_val,method = "bonferroni") # manual fdr
      all_activationGenes_negbinom_list[[current_cluster]] = activationGenes_current_negbinom
    }
    
    
  }
}

# rbind
activation_per_cluster = do.call(rbind,all_activationGenes_list)

## get negbinom pvals
activation_per_cluster_negbinom = do.call(rbind,all_activationGenes_negbinom_list)
activation_per_cluster_negbinom_sel = activation_per_cluster_negbinom %>% dplyr::select(gene,current_cluster,p_val_negbinom = p_val, p_val_bonferroni_negbinom = p_val_bonferroni, )
activation_per_cluster = dplyr::left_join(activation_per_cluster,activation_per_cluster_negbinom_sel,by=c("gene"="gene","current_cluster"="current_cluster"))

# get total number of expressed IEGs --> just use the 10% requirement from marker detection !
activation_per_cluster  = activation_per_cluster %>% dplyr::group_by(current_cluster) %>% dplyr::add_count(name = "expressed_genes")

# save table
data.table::fwrite(activation_per_cluster,file = paste0(results_path_figure5,"activation_genes_per_cluster.txt"),sep="\t")

#activation_per_cluster = data.table::fread(paste0(results_path_figure5,"activation_genes_per_cluster.txt"),data.table = FALSE)

## summarise:
# define a column for pvalue aggregation
activation_per_cluster$Identity = activation_per_cluster$current_cluster
activation_per_cluster$pvalue_for_calc = activation_per_cluster$p_val_bonferroni
pval_ieg_max = 0.05
fc_ieg_min = 0.1


# changed to simple number of genes:
activation_per_cluster_stat = activation_per_cluster %>% group_by(Identity) %>% 
  dplyr::summarise( mean_fc_up = mean(avg_log2FC[pvalue_for_calc<pval_ieg_max & avg_log2FC > fc_ieg_min]),
                    mean_fc_down = mean(avg_log2FC[pvalue_for_calc<pval_ieg_max  & avg_log2FC < (-fc_ieg_min)]),
                    n_down_sig = sum(pvalue_for_calc<pval_ieg_max  & avg_log2FC < (-fc_ieg_min)), # count genes
                    n_up_sig = sum(pvalue_for_calc<pval_ieg_max  & avg_log2FC > fc_ieg_min)) %>%
  dplyr::filter(!(is.nan(mean_fc_up) & is.nan(mean_fc_down))) %>% # general filter for anything significant
  #dplyr::filter( (n_down_sig - n_up_sig) >= 1 & n_down_sig >= 2) %>% # here I FILTER !!!!!
  dplyr::arrange(desc(n_up_sig),mean_fc_down)
n_cells_per_cluster = dowsett_subset@meta.data %>% dplyr::group_by(C286_modified) %>% dplyr::count(name="n_cells_cluster")
activation_per_cluster_stat = dplyr::left_join(activation_per_cluster_stat,n_cells_per_cluster,by=c("Identity"="C286_modified"))
activation_per_cluster_stat = dplyr::left_join(activation_per_cluster_stat,activation_per_cluster %>% dplyr::distinct(Identity,expressed_genes),by=c("Identity"="Identity"))
activation_per_cluster_stat

# adjust to percentage
activation_per_cluster_stat$pct_activated_genes = round(activation_per_cluster_stat$n_up_sig / activation_per_cluster_stat$expressed_genes,4) * 100 

# add names
activation_per_cluster_stat = dplyr::left_join(activation_per_cluster_stat,hypoMap_v2_seurat@misc$annotation_result %>% dplyr::select(Identity = cluster_id,Name = clean_names_withID),by="Identity")
activation_per_cluster_stat = activation_per_cluster_stat %>% dplyr::arrange(desc(pct_activated_genes))

# save:
#data.table::fwrite(activation_per_cluster_stat,file = paste0(results_path_figure5,"activation_stat_per_cluster.txt"),sep="\t")

# make barplot for Figure
neuron_clusters = scUtils::find_children("C2-1",edges = hypoMap_v2_seurat@misc$clustering_edgelist)
# I filter out clusters with only 1 activated gene !!!!
activation_per_cluster_stat$Name = factor(activation_per_cluster_stat$Name,levels = rev(activation_per_cluster_stat$Name))
activation_celltype_barplot = ggplot(activation_per_cluster_stat %>% dplyr::filter(n_up_sig > 2 & Identity %in% neuron_clusters),aes(x=Name,y=pct_activated_genes,fill= mean_fc_up))+geom_bar(stat="identity")+
  geom_text(aes(x = Name, y = max(pct_activated_genes)*1.22,label=n_cells_cluster), hjust =1,size=6)+ # add cell numbers
  geom_text(aes(x = Name, y = pct_activated_genes*1.05,label=n_up_sig), hjust =1,size=6)+ # add ieg numbers
  scale_fill_gradient(low=bg_col,high=fasting_color,limits=c(0,max(abs(activation_per_cluster_stat$mean_fc_up))))+
  ylim(c(0,max(activation_per_cluster_stat$pct_activated_genes)*1.22)) + # extend the y axis a bit so that the n cell labels doN#t overlap with the highest plot!
  # scale_y_continuous(breaks=c(2,4,6,8,10))+
  coord_flip()+ylab("Percent of activated IEGs")+xlab(NULL)+  #"Summed adjusted p-value"
  labs(fill='Mean log2FC')+
  theme_bw()+ theme(text = element_text(size=text_size),axis.title.x = element_text(size = text_size))
activation_celltype_barplot

# save
ggsave(filename = paste0(results_path_figure5,"activation_celltype_barplot.png"),
       plot = activation_celltype_barplot, "png",dpi=450,width=350,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"activation_celltype_barplot.pdf"),
       plot = activation_celltype_barplot, "pdf",dpi=450,width=350,height = 200,units="mm")

# source
source_figure5_c = activation_celltype_barplot$data
data.table::fwrite(source_figure5_c,paste0(results_path_figure5,"source_figure5_c_celltype_barplot.txt"),sep="\t")

##########
### Immediate early genes in Agrp neurons - Violin plot
##########

# top Agrp genes
activation_genes_plot_subset = unique(activation_per_cluster$gene[activation_per_cluster$pvalue_for_calc<0.05 & 
                                                                    activation_per_cluster$Identity %in% c(agrp_cluster_to_use) &
                                                                    activation_per_cluster$avg_log2FC > 0.35 ])
activation_genes_plot_subset = unique(c("Fos",activation_genes_plot_subset,"Sgk1"))
#activation_genes_plot_subset = c("Fos","1700016P03Rik","Gem","Btg2","Noct","Junb","Sgk1")
activation_genes_plot_subset

Idents(dowsett_subset) <- "C286_modified"
p <- Seurat::VlnPlot(dowsett_subset,features = activation_genes_plot_subset,split.by = "Diet",idents = c(agrp_cluster_to_use),cols=c(fasting_color,adlib_color),
                     ncol =4,same.y.lims = TRUE,combine = FALSE) 
legend <- cowplot::get_legend( p[[1]])
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(axis.title.x = element_text(size = 10),axis.ticks.x =  element_blank(),axis.text.x=element_blank())+
    NoLegend()+xlab(NULL)+ylab(NULL) #xlab("Agrp neurons")+ylab("")
}
p[[9]] = cowplot::plot_grid(legend)
activation_genes_agrp_vlnPlot = cowplot::plot_grid(plotlist = p,ncol = 4)
activation_genes_agrp_vlnPlot

# save
ggsave(filename = paste0(results_path_figure5,"activation_genes_agrp_vlnPlot.png"),
       plot = activation_genes_agrp_vlnPlot, "png",dpi=450,width=400,height = 250,units="mm")
ggsave(filename = paste0(results_path_figure5,"activation_genes_agrp_vlnPlot.pdf"),
       plot = activation_genes_agrp_vlnPlot, "pdf",dpi=450,width=400,height = 250,units="mm")

# source
source_figure5_b = cbind(p[[1]]$data[,c(3,1)],do.call(cbind,lapply(p[2:8],function(p){return(p$data[,c(1),drop=FALSE])})))
source_figure5_b$Cell_ID = rownames(source_figure5_b)
data.table::fwrite(source_figure5_b,paste0(results_path_figure5,"source_figure5_b_vln_iegs.txt"),sep="\t")

#### save as vertical version
# activation_genes_agrp_vlnPlot_vertical = cowplot::plot_grid(plotlist = p,ncol =2)
# activation_genes_agrp_vlnPlot_vertical
# 
# # save
# ggsave(filename = paste0(results_path_figure5,"activation_genes_agrp_vlnPlot_vertical.png"),
#        plot = activation_genes_agrp_vlnPlot_vertical, "png",dpi=450,width=150,height = 250,units="mm")
# ggsave(filename = paste0(results_path_figure5,"activation_genes_agrp_vlnPlot_vertical.pdf"),
#        plot = activation_genes_agrp_vlnPlot_vertical, "pdf",dpi=450,width=150,height = 250,units="mm")

##########
### Transcriptional changes in Agrp
##########

# get all DEG of Agrp
conditionGenes_Agrp_file = paste0(results_path_figure5,"agrp_fasting_all.txt")
conditionGenes_Agrp_negbinom_file = paste0(results_path_figure5,"agrp_fasting_negbinom_all.txt")
if(!file.exists(conditionGenes_Agrp_file)){
  Idents(dowsett_subset) <- "C286_modified"
  # I am not filtering the genes here to get values for all genes (for comparison with other datasets etc.) This does not have a strong effect on the Bonferroni correction done by Seurat!
  conditionGenes_Agrp = Seurat::FindMarkers(dowsett_subset,ident.1 = "Fasted",ident.2 = "Normal chow", group.by = "Diet", subset.ident = agrp_cluster_to_use,min.pct = 0,logfc.threshold = 0)
  conditionGenes_Agrp$gene = rownames(conditionGenes_Agrp) # add gene name 
  conditionGenes_Agrp$pct_diff = conditionGenes_Agrp$pct.1 - conditionGenes_Agrp$pct.2
  data.table::fwrite(conditionGenes_Agrp,conditionGenes_Agrp_file,sep="\t")
  
  # with deseq
  genes_negbinom_to_test = conditionGenes_Agrp$gene[conditionGenes_Agrp$pct.1 > 0.05 | conditionGenes_Agrp$pct.2 > 0.05]
  conditionGenes_Agrp_negbinom = Seurat::FindMarkers(dowsett_subset,ident.1 = "Fasted",ident.2 = "Normal chow", group.by = "Diet", subset.ident = agrp_cluster_to_use,min.pct = 0,logfc.threshold = 0,test.use = "negbinom")
  conditionGenes_Agrp_negbinom$gene = rownames(conditionGenes_Agrp_negbinom) # add gene name 
  conditionGenes_Agrp_negbinom$pct_diff = conditionGenes_Agrp_negbinom$pct.1 - conditionGenes_Agrp_negbinom$pct.2
  data.table::fwrite(conditionGenes_Agrp_negbinom,conditionGenes_Agrp_negbinom_file,sep="\t")
  
}else{
  # load and save above code which runs for 8 minutes !
  conditionGenes_Agrp = data.table::fread(conditionGenes_Agrp_file,data.table = F)
  conditionGenes_Agrp_negbinom = data.table::fread(conditionGenes_Agrp_negbinom_file,data.table = F)
}

conditionGenes_Agrp_with_negbinom = conditionGenes_Agrp  %>% dplyr::left_join(conditionGenes_Agrp_negbinom %>% dplyr::select(gene,p_val_adj_negbinom=p_val_adj),by="gene")

#filter
conditionGenes_Agrp_filtered = conditionGenes_Agrp[conditionGenes_Agrp$p_val_adj < 0.05,] # filter pval
#conditionGenes_Agrp_filtered %>% dplyr::arrange((pct_diff),(avg_log2FC))
conditionGenes_Agrp_filtered_with_negbinom = conditionGenes_Agrp_filtered %>% dplyr::left_join(conditionGenes_Agrp_negbinom %>% dplyr::select(gene,p_val_adj_negbinom=p_val_adj),by="gene")


### make volcano plot
volcano_df = conditionGenes_Agrp_filtered_with_negbinom
volcano_df$color = "not regulated"
volcano_df$color[volcano_df$avg_log2FC< (-0.35) & volcano_df$p_val_adj < 0.001 ] = "down-regulated in fasting"
volcano_df$color[volcano_df$avg_log2FC > 0.35 & volcano_df$p_val_adj < 0.001 ] = "up-regulated in fasting"

# save for sppl
# conditionGenes_Agrp_fasting_file = ""
# data.table::fwrite(conditionGenes_Agrp_filtered_with_negbinom,conditionGenes_Agrp_fasting_file,sep="\t")
# 
# conditionGenes_Agrp_filtered_with_negbinom = data.table::fread(conditionGenes_Agrp_fasting_file)

# define gene sets that are labelled:
up_regulated_genes = conditionGenes_Agrp_filtered$gene[conditionGenes_Agrp_filtered$avg_log2FC > 1 & conditionGenes_Agrp_filtered$p_val_adj < 0.0001 ]
up_regulated_genes = up_regulated_genes[!grepl("Rps|Rpl|Nduf|ENSMUS",up_regulated_genes)]
up_regulated_genes = c(up_regulated_genes,"Fam107b","Lepr","Zbtb16","Tnik","Asic2","Acvr1c","Ncoa2","Vgf" )
down_regulated_genes = conditionGenes_Agrp_filtered$gene[conditionGenes_Agrp_filtered$avg_log2FC < -1.25 & conditionGenes_Agrp_filtered$p_val_adj < 0.0001 ]
down_regulated_genes = down_regulated_genes[!grepl("Rps|Rpl|Nduf|ENSMUS|CT",down_regulated_genes)]
down_regulated_genes = c(down_regulated_genes,"Agrp","Sst","Cartpt","Npy")
volcano_df$label = NA
volcano_df$label[volcano_df$gene %in% c(down_regulated_genes,up_regulated_genes)] = volcano_df$gene[volcano_df$gene %in% c(down_regulated_genes,up_regulated_genes)]

agrp_fasting_volcano = ggplot(volcano_df[volcano_df$avg_log2FC !=0,],aes(avg_log2FC,-log10(p_val_adj),label=label,color=color))+
  # geom_point(size=0.3)+
  scattermore::geom_scattermore(pixels = c(2048,2048),pointsize = 4)+
  ggrepel::geom_text_repel(max.overlaps=100,size=6,show.legend = FALSE)+
  scale_color_manual(values = c("up-regulated in fasting"=fasting_color,"down-regulated in fasting"=adlib_color,"not regulated"=bg_col))+
  guides(colour = guide_legend(override.aes = list(size=7)))+
  xlab("log2FC")+ylab("-log10(adjusted pvalue)")+ggtitle("Changes in Agrp neurons")+
  xlim(c(-1*(max(abs(volcano_df$avg_log2FC))+0.1),max(abs(volcano_df$avg_log2FC))+0.1))+
  #ylim(c(0,-log10(0.00000000000001)))+
  theme_bw()+
  #cowplot::theme_cowplot()+
  theme(text = element_text(size=text_size))+labs(color='Condition')
agrp_fasting_volcano

#save
ggsave(filename = paste0(results_path_figure5,"agrp_fasting_volcano.png"),
       plot = agrp_fasting_volcano, "png",dpi=450,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"agrp_fasting_volcano.pdf"),
       plot = agrp_fasting_volcano, "pdf",dpi=450,width=300,height = 200,units="mm")

# source
source_figure5_d = agrp_fasting_volcano$data
data.table::fwrite(source_figure5_d,paste0(results_path_figure5,"source_figure5_d_volcano_agrp.txt"),sep="\t")

##########
### Go enrichment using clusterProfiler
##########

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

# which gene to enrich for
agrp_upregulated_genes_names = conditionGenes_Agrp_filtered$gene[conditionGenes_Agrp_filtered$avg_log2FC>(0.35) & conditionGenes_Agrp_filtered$p_val_adj<0.001]

# need to map to entrez
library(biomaRt)
mart <- useMart(dataset="mmusculus_gene_ensembl",biomart='ensembl',host="nov2020.archive.ensembl.org") # use 2020 release because better compatbile with most datasets  #or: feb2021.archive.ensembl.org
up_regulated_genes_ids = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','entrezgene_id'),
                               filters = "external_gene_name",values =agrp_upregulated_genes_names,mart = mart)
# entrez_ids = up_regulated_genes

# need to map background to entrex
cellsh = dowsett_subset@meta.data$Cell_ID[dowsett_subset@meta.data$C286_modified == agrp_cluster_to_use]
agrp_baseline_expression = get_expression_stats(dowsett_subset,cells.1 = cellsh)
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

# save table
data.table::fwrite(agrp_fasting_go_enrichment_simplified_res,file = paste0(results_path_figure5,"agrp_fasting_go_enrichment_simplified.txt"),sep="\t")

#agrp_fasting_go_enrichment_simplified_res = data.table::fread(file = paste0(results_path_figure5,"agrp_fasting_go_enrichment_simplified.txt"),data.table = F)

# make dotplot
go_bp_agrp_enrich_dotplot = enrichplot::dotplot(agrp_fasting_go_enrichment_simplified, showCategory=20) + 
  ggtitle("Dotplot for ORA")+theme(axis.text.y = element_text(size=15))+
  scale_color_gradient(low=fasting_color,high ="#a1bdc4")
go_bp_agrp_enrich_dotplot

# save
ggsave(filename = paste0(results_path_figure5,"go_bp_agrp_enrich_dotplot.png"),
       plot = go_bp_agrp_enrich_dotplot, "png",dpi=450,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"go_bp_agrp_enrich_dotplot.pdf"),
       plot = go_bp_agrp_enrich_dotplot, "pdf",dpi=450,width=300,height = 200,units="mm")

#source
source_figure5_d_dotplot =go_bp_agrp_enrich_dotplot$data
data.table::fwrite(source_figure5_d_dotplot,paste0(results_path_figure5,"source_figure5_d_enrichment_agrp.txt"),sep="\t")

##########
### Transcriptional changes in Agrp of Campbell vs NucSeq
##########

#load_required_files(large_data_path, overwrite_existing = TRUE, filenames = c("neuron_map_seurat" = "hypothalamus_neurons_map.rds" ))

campbell_diet = subset(hypoMap_v2_seurat,subset = Dataset=="CampbellDropseq")

# need to make modfied column in full object:
campbell_diet@meta.data$C286_modified = campbell_diet@meta.data$C286
campbell_diet@meta.data$C286_modified[campbell_diet@meta.data$C286_modified %in% c("C286-176","C286-177","C286-178")] = agrp_cluster_to_use

table(campbell_diet@meta.data$Diet)

### Comparison with Campbell Agrp changes:
campbell_agrp_file = paste0(results_path_figure5,"agrp_fasting_all_campbell.txt")
campbell_agrp_negbinom_file = paste0(results_path_figure5,"agrp_fasting_all_campbell_negbinom.txt")
if(!file.exists(campbell_agrp_file)){
  Idents(campbell_diet) <- "C286_modified"
  conditionGenes_campbell = Seurat::FindMarkers(campbell_diet, ident.1 = "Fasted",ident.2 = "Normal chow", group.by = "Diet", subset.ident = agrp_cluster_to_use,min.pct = 0,logfc.threshold = 0)
  conditionGenes_campbell$gene = rownames(conditionGenes_campbell) # add gene name 
  conditionGenes_campbell$pct_diff = conditionGenes_campbell$pct.1 - conditionGenes_campbell$pct.2
  data.table::fwrite(conditionGenes_campbell,campbell_agrp_file,sep="\t")
  # with negbinom
  conditionGenes_campbell_negbinom = Seurat::FindMarkers(campbell_diet,ident.1 = "Fasted",ident.2 = "Normal chow", group.by = "Diet", 
                                                         subset.ident = agrp_cluster_to_use,min.pct = 0,logfc.threshold = 0,test.use = "negbinom",min.cells.feature = 5)
  conditionGenes_campbell_negbinom$gene = rownames(conditionGenes_campbell_negbinom) # add gene name 
  conditionGenes_campbell_negbinom$pct_diff = conditionGenes_campbell_negbinom$pct.1 - conditionGenes_campbell_negbinom$pct.2
  data.table::fwrite(conditionGenes_campbell_negbinom,campbell_agrp_negbinom_file,sep="\t")
}else{
  # load and save above code which runs for 8 minutes !
  conditionGenes_campbell = data.table::fread(campbell_agrp_file,data.table = F)
  conditionGenes_campbell_negbinom = data.table::fread(campbell_agrp_negbinom_file,data.table = F)
}

#### Compare with sc-seq data:
conditionGenes_campbell_negbinom = conditionGenes_campbell  %>% dplyr::left_join(conditionGenes_campbell_negbinom %>% dplyr::select(gene,p_val_adj_negbinom=p_val_adj),by="gene")

agrp_sn_vs_sc = dplyr::full_join(conditionGenes_Agrp_with_negbinom,conditionGenes_campbell_negbinom,suffix=c("_sn","_sc"),by=c("gene"="gene"))
colnames(agrp_sn_vs_sc)

#### Define categories
sc_p_cut = 0.001
min_fc = 0.2

##### make simple classification:
agrp_sn_vs_sc_plot = agrp_sn_vs_sc

agrp_sn_vs_sc_plot$avg_log2FC_sn = agrp_sn_vs_sc_plot$avg_log2FC_sn#*(-1)
agrp_sn_vs_sc_plot$avg_log2FC_sc = agrp_sn_vs_sc_plot$avg_log2FC_sc#*(-1)
agrp_sn_vs_sc_plot$regulated = "not - both"
agrp_sn_vs_sc_plot$regulated[ ( agrp_sn_vs_sc_plot$avg_log2FC_sn > min_fc & agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut ) | 
                                (agrp_sn_vs_sc_plot$avg_log2FC_sc  > min_fc & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut)] = "up in fasting - one"

agrp_sn_vs_sc_plot$regulated[agrp_sn_vs_sc_plot$avg_log2FC_sn>min_fc & agrp_sn_vs_sc_plot$avg_log2FC_sc>min_fc & 
                               agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut] = "up in fasting - both"

agrp_sn_vs_sc_plot$regulated[ ( agrp_sn_vs_sc_plot$avg_log2FC_sn< (-1*min_fc) & agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut ) | 
                                (agrp_sn_vs_sc_plot$avg_log2FC_sc  < (-1*min_fc) & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut)] = "down in fasting - one"

agrp_sn_vs_sc_plot$regulated[agrp_sn_vs_sc_plot$avg_log2FC_sn< (-1*min_fc) & agrp_sn_vs_sc_plot$avg_log2FC_sc  < (-1*min_fc) & 
                               agrp_sn_vs_sc_plot$p_val_adj_sn < sc_p_cut & agrp_sn_vs_sc_plot$p_val_adj_sc < sc_p_cut] = "down in fasting - both"

agrp_sn_vs_sc_plot$regulated = factor(agrp_sn_vs_sc_plot$regulated,levels = c( "not - both","down in fasting - one", "up in fasting - one","down in fasting - both","up in fasting - both"))

## calculate cor for text:
cor(agrp_sn_vs_sc_plot$avg_log2FC_sc[agrp_sn_vs_sc_plot$regulated !=  "not - both"],agrp_sn_vs_sc_plot$avg_log2FC_sn[agrp_sn_vs_sc_plot$regulated !=  "not - both"],use = "pairwise.complete.obs")
cor.test(agrp_sn_vs_sc_plot$avg_log2FC_sc[agrp_sn_vs_sc_plot$regulated !=  "not - both"],agrp_sn_vs_sc_plot$avg_log2FC_sn[agrp_sn_vs_sc_plot$regulated !=  "not - both"],
         use = "pairwise.complete.obs",method="pearson",alternative = "two.sided")

# make plot for figure:
agrp_Campvell_vs_sn_plot = ggplot(agrp_sn_vs_sc_plot,aes(x=avg_log2FC_sc,y=avg_log2FC_sn))+
  # geom_point(alpha=0.6,aes(color=regulated))+
  scattermore::geom_scattermore(alpha=0.6,aes(color=regulated),pixels = c(2048,2048),pointsize = 5)+
  geom_smooth(method="lm",color="grey60")+
  scale_color_manual(values = c("up in fasting - both" = fasting_color,"not - both" = bg_col,"down in fasting - both" =adlib_color,
                                "up in fasting - one" = "#a3bbc2", "down in fasting - one" = "#ae9bc2"))+
  xlab("log2FC sc-seq")+ylab("log2FC sn-seq")+
  # cowplot::theme_cowplot()+
  theme_bw()+
  theme(text=element_text(size=text_size))
agrp_Campvell_vs_sn_plot

# agrp_Campvell_vs_sn_plot_r = scUtils::rasterize_ggplot(agrp_Campvell_vs_sn_plot,pixel_raster = rasterize_pixels,pointsize = rasterize_point_size)
# agrp_Campvell_vs_sn_plot_r
# 
# agrp_Campvell_vs_sn_plot_r + geom_smooth(method="lm",color="grey60")

# save
ggsave(filename = paste0(results_path_figure5,"agrp_Campvell_vs_sn_plot.png"),
       plot = agrp_Campvell_vs_sn_plot, "png",dpi=450,width=300,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"agrp_Campvell_vs_sn_plot.pdf"),
       plot = agrp_Campvell_vs_sn_plot, "pdf",dpi=450,width=300,height = 200,units="mm")

#source
data.table::fwrite(agrp_sn_vs_sc_plot,paste0(results_path_figure5,"source_figure5_e_agrp_sn_vs_campbell_DEG.txt"),sep = "\t")

##########
### Example genes for last section of Figure
##########

# Zbtb16, Fam107b, Vgf, Sv2c, 1700016P03Rik
genes_to_use =c("Zbtb16", "Fam107b", "Vgf", "Sv2c")

## Option 1: just a nice feature plot of a gene of choice
# Idents(dowsett_subset) <- "predicted_K31_named"
# example_gene_plot = Seurat::FeaturePlot(dowsett_subset,"1700016P03Rik",split.by = "Diet",label = TRUE,label.size = 3.5,repel = TRUE,pt.size = 0.4,
#                                         keep.scale="feature",order = TRUE,combine = FALSE)
# example_gene_plot[[1]] = example_gene_plot[[1]]+NoAxes()+theme(panel.border = element_blank())+ggplot2::ggtitle("adlib")#+scale_color_gradient(low = "lightgrey",high = "#8c390a")
# example_gene_plot[[2]] = example_gene_plot[[2]]+NoAxes()+theme(panel.border = element_blank())+ggplot2::ggtitle("fasted")#+scale_color_gradient(low = "lightgrey",high = "#8c390a")
# example_gene_plot = cowplot::plot_grid(plotlist = example_gene_plot)
# example_gene_plot

## Option 2: some violin plot
# use all other cells and top 5 from enrichement plus pomc plus cck
# "C286-99: Lepr.Agrp.GABA-3","C286-114: Sox14.Unc13c.Lef1.GABA-2","C286-122: Tafa4.Ppp1r17.GABA-2","C286-45: Rai14.Hmcn2.vGLUT-2","C286-75: Ttr.Pomc.vGLUT-4","C286-70: Anxa2.Pomc.vGLUT-4","other celltypes"
classes_to_use = c("C286-178","C286-149","C286-139","C286-151","C286-75","C286-77")
dowsett_subset@meta.data$custom_annotation = "other celltypes"
dowsett_subset@meta.data$custom_annotation[dowsett_subset@meta.data$C286_modified %in% classes_to_use] = dowsett_subset@meta.data$C286_named[dowsett_subset@meta.data$C286_modified %in% classes_to_use]
dowsett_subset@meta.data$custom_annotation[dowsett_subset@meta.data$custom_annotation %in% c("C286-176: Serpina3n.Npy.Agrp.GABA-4","C286-177: Acvr1c.Npy.Agrp.GABA-4")] = "C286-178: Lepr.Agrp.GABA-4"
dowsett_subset@meta.data$custom_annotation = factor(dowsett_subset@meta.data$custom_annotation,
                                                    levels=rev(c("C286-178: Lepr.Agrp.GABA-4","C286-149: Grp.Ppp1r17.GABA-1","C286-139: Myo5b.Sox14.Lef1.GABA-1",
                                                                 "C286-151: St18.Lhx6.GABA-1","C286-75: Anxa2.Pomc.GLU-5","C286-77: Ttr.Pomc.GLU-5","other celltypes")))

# make plot
Idents(dowsett_subset) <- "custom_annotation"
interesting_genes_violin_plot = Seurat::VlnPlot(dowsett_subset,features = genes_to_use,split.by = "Diet",cols=c(fasting_color,adlib_color),
                                                same.y.lims = TRUE,combine = TRUE,stack = TRUE,adjust = 0.75,flip=TRUE,pt.size = 1) 
interesting_genes_violin_plot = interesting_genes_violin_plot+xlab("")
interesting_genes_violin_plot

# save
ggsave(filename = paste0(results_path_figure5,"interesting_genes_violin_plot.png"),
       plot = interesting_genes_violin_plot, "png",dpi=450,width=280,height = 200,units="mm")
ggsave(filename = paste0(results_path_figure5,"interesting_genes_violin_plot.pdf"),
       plot = interesting_genes_violin_plot, "pdf",dpi=450,width=280,height = 200,units="mm")

# source
source_figure5_f = bind_cols(dowsett_subset@meta.data[,c("Cell_ID","custom_annotation")],Seurat::FetchData(dowsett_subset,genes_to_use))
data.table::fwrite(source_figure5_f,paste0(results_path_figure5,"source_figure5_f_vln_examples.txt"),sep="\t")



