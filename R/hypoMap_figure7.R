##########
### Load & Prepare
##########

results_path_figure7 = "figure_outputs/figure_7/"
system(paste0("mkdir -p ",results_path_figure7))

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
### Load ish data
##########

data_path = "data_inputs/"

#### load ish results
ish_quantification_glp1r = data.table::fread(paste0(data_path,"ish_quantification_glp1r_updated.csv"),data.table = F)

#### load ish results
ish_quantification_pnoc_crabp1 = data.table::fread(paste0(data_path,"ish_quantification_pnoc_crabp1.csv"),data.table = F)

#### load ish results
ish_quantification_pnoc_sst = data.table::fread(paste0(data_path,"ish_quantification_pnoc_sst.csv"),data.table = F)

##########
### Manually add rows for Pomc and Sst without their second marker
##########

# add pomc alone
ish_quantification_glp1r_pomc = ish_quantification_glp1r[ish_quantification_glp1r$Experiment=="Pomc Anxa2",]
ish_quantification_glp1r_pomc$Experiment="Pomc"
ish_quantification_glp1r_pomc$total_1_2=NA
ish_quantification_glp1r_pomc$total_1_2_Glp1r=NA
ish_quantification_glp1r = rbind(ish_quantification_glp1r,ish_quantification_glp1r_pomc)
# add sst alone
ish_quantification_glp1r_sst = ish_quantification_glp1r[ish_quantification_glp1r$Experiment=="Sst Unc13c",]
ish_quantification_glp1r_sst$Experiment="Sst"
ish_quantification_glp1r_sst$total_1_2=NA
ish_quantification_glp1r_sst$total_1_2_Glp1r=NA
ish_quantification_glp1r = rbind(ish_quantification_glp1r,ish_quantification_glp1r_sst)

# add trh alone
ish_quantification_glp1r_trh = ish_quantification_glp1r[ish_quantification_glp1r$Experiment=="Trh Nkx2-4",]
ish_quantification_glp1r_trh$Experiment="Trh"
ish_quantification_glp1r_trh$total_1_2=NA
ish_quantification_glp1r_trh$total_1_2_Glp1r=NA
ish_quantification_glp1r = rbind(ish_quantification_glp1r,ish_quantification_glp1r_trh)

##########
### make percentages Glp1r
##########

# the data contains cell counts for (example in brackets)
# main marker positive cells (Pomc)
# main-marker + second marker positive cells (Pomc+Anxa2)
# main marker + glp1r positive cells (Pomc + Glp1r)
# main marker + second marker + glp1r positive cells (Pomc+ Anxa2 + Glp1r)

# I make pct 1 as: main marker + glp1r  / main marker
# and pct 2 as: main marker + second marker + glp1r / main-marker + second marker

# the if no second marker is available I cop the values from pct 1 into this column: making pct2 the main value for plotting !!!

# make pcts
ish_quantification_glp1r$total_1_pct = ish_quantification_glp1r$total_1_glp1r / ish_quantification_glp1r$total_1 * 100
ish_quantification_glp1r$total_2_pct = ish_quantification_glp1r$total_1_2_Glp1r / ish_quantification_glp1r$total_1_2 * 100
ish_quantification_glp1r$total_2_pct[is.na(ish_quantification_glp1r$total_2_pct)] = ish_quantification_glp1r$total_1_pct[is.na(ish_quantification_glp1r$total_2_pct)]
ish_quantification_glp1r$total_2_pct[is.nan(ish_quantification_glp1r$total_2_pct)] = 0

# Example summary
ish_quantification_glp1r %>% group_by(Experiment) %>% dplyr::summarise(mean(total_2_pct))

##########
### add mouse ids dataframe
##########

## I assume that the first part of Section is the mouse id
ish_quantification_glp1r$Section[ish_quantification_glp1r$Section == "49B2section2"] = "49B2 section2"
ish_quantification_glp1r$mouse_id = sapply(ish_quantification_glp1r$Section,function(x){strsplit(x," |_")[[1]][1]})

ish_quantification_glp1r %>% dplyr::group_by(Experiment,mouse_id) %>% dplyr::summarise(n()) %>% dplyr::count()

#ish_quantification_glp1r_mouse = ish_quantification_glp1r %>% dplyr::group_by(mouse_id,Experiment) %>% 
#  summarise(across(c(total_1,total_1_2,total_1_glp1r,total_1_2_Glp1r), sum,na.rm=TRUE))

# save result to folder:
data.table::fwrite(ish_quantification_glp1r,file=paste0(results_path_figure7,"pct_expressed_cells_glp1r_ISH.txt"),sep="\t")


##########
### Plot dotplots (all)
##########

magically_transform_to_the_right_pcts = function(count_df_ish,start_col=2,sep_char = "/", max_ref_genes = 2){
  
  split_colnames = stringr::str_split(colnames(count_df_ish)[start_col:ncol(count_df_ish)],pattern = sep_char)
  names(split_colnames) = colnames(count_df_ish)[start_col:ncol(count_df_ish)]
  pct_res_list = list()
  for(i in 1:length(split_colnames)){
    current_col =  names(split_colnames)[i]
    current_genes = split_colnames[[i]]
    if(length(current_genes) > 1){
      dennominator_col = paste0(current_genes[1:min(max_ref_genes,(length(current_genes)-1))],collapse = sep_char) 
      dennominator_data = count_df_ish[,dennominator_col]
      nominator_data = count_df_ish[,current_col] 
      current_pct = nominator_data/dennominator_data * 100
      newname = paste0(current_col,"_of_",dennominator_col)
      pct_res_list[[newname]] = current_pct
    }
    
  }
  pct_res_df =  as.data.frame(do.call(cbind,pct_res_list))
  pct_res_df = as.data.frame(cbind(count_df_ish[1:max(1,(start_col-1))],pct_res_df))
  
  return(pct_res_df)
}

# I need to reformat some columns so that magically_transform_to_the_right_pcts actually calulcates the right pcts (e.g. Pnoc+sst of sst and not Pnoc+sst of Pnoc
colnames(ish_quantification_pnoc_crabp1)[colnames(ish_quantification_pnoc_crabp1) == "Pnoc/Crabp1"] = "Crabp1/Pnoc"
colnames(ish_quantification_pnoc_crabp1)[colnames(ish_quantification_pnoc_crabp1) == "Pnoc/Crabp1/Tmem215"] = "Crabp1/Pnoc/Tmem215"
colnames(ish_quantification_pnoc_crabp1)[colnames(ish_quantification_pnoc_crabp1) == "Pnoc/Crabp1/Htr3b"] = "Crabp1/Pnoc/Htr3b"
colnames(ish_quantification_pnoc_crabp1)[colnames(ish_quantification_pnoc_crabp1) == "Pnoc/Crabp1/Tmem215/Htr3b"] = "Crabp1/Pnoc/Tmem215/Htr3b"
#
colnames(ish_quantification_pnoc_sst)[colnames(ish_quantification_pnoc_sst) == "Pnoc/Sst"] = "Sst/Pnoc"
colnames(ish_quantification_pnoc_sst)[colnames(ish_quantification_pnoc_sst) == "Pnoc/Sst/Nts"] = "Sst/Pnoc/Nts"
colnames(ish_quantification_pnoc_sst)[colnames(ish_quantification_pnoc_sst) == "Pnoc/Sst/Unc13c"] = "Sst/Pnoc/Unc13c"
colnames(ish_quantification_pnoc_sst)[colnames(ish_quantification_pnoc_sst) == "Pnoc/Sst/Nts/Unc13"] = "Sst/Pnoc/Nts/Unc13"

pcts_ish_quantification_pnoc_crabp1 = magically_transform_to_the_right_pcts(ish_quantification_pnoc_crabp1[,!colnames(ish_quantification_pnoc_crabp1) %in% "Pnoc/Crabp1-"],start_col=2,sep_char = "/", max_ref_genes = 2)

pcts_ish_quantification_pnoc_sst = magically_transform_to_the_right_pcts(ish_quantification_pnoc_sst[,!colnames(ish_quantification_pnoc_sst) %in% "Pnoc/Sst-"],start_col=2,sep_char = "/", max_ref_genes = 2)

# save
data.table::fwrite(pcts_ish_quantification_pnoc_crabp1,file=paste0(results_path_figure7,"pct_expressed_cells_pnoc_crabp1_ISH.txt"),sep="\t")
data.table::fwrite(pcts_ish_quantification_pnoc_sst,file=paste0(results_path_figure7,"pct_expressed_cells_pnoc_sst_ISH.txt"),sep="\t")

# test to format Glp1r in the same way:

# sep_char = "/"
# pcts_ish_quantification_glp1r_long = ish_quantification_glp1r %>% tidyr::gather(key="totals",value = "count",-Experiment,-Location,-Section)#tidyr::spread(key = Experiment,value = 4:7)
# pcts_ish_quantification_glp1r_long$gene1 = sapply(pcts_ish_quantification_glp1r_long$Experiment,function(x){stringr::str_split(x,pattern = " ")[[1]][1]})
# pcts_ish_quantification_glp1r_long$gene2 = sapply(pcts_ish_quantification_glp1r_long$Experiment,function(x){stringr::str_split(x,pattern = " ")[[1]][2]})
# pcts_ish_quantification_glp1r_long$new_total = NA
# pcts_ish_quantification_glp1r_long$new_total[pcts_ish_quantification_glp1r_long$totals == "total_1"] = paste0(pcts_ish_quantification_glp1r_long$gene1[pcts_ish_quantification_glp1r_long$totals == "total_1"])
# pcts_ish_quantification_glp1r_long$new_total[pcts_ish_quantification_glp1r_long$totals == "total_1_glp1r"] = paste0(pcts_ish_quantification_glp1r_long$gene1[pcts_ish_quantification_glp1r_long$totals == "total_1_glp1r"],sep_char,"Glp1r")
# pcts_ish_quantification_glp1r_long$new_total[pcts_ish_quantification_glp1r_long$totals == "total_1_2"] = paste0(pcts_ish_quantification_glp1r_long$gene1[pcts_ish_quantification_glp1r_long$totals == "total_1_2"],sep_char,pcts_ish_quantification_glp1r_long$gene2[pcts_ish_quantification_glp1r_long$totals == "total_1_2"])
# pcts_ish_quantification_glp1r_long$new_total[pcts_ish_quantification_glp1r_long$totals == "total_1_2_Glp1r"] = paste0(pcts_ish_quantification_glp1r_long$gene1[pcts_ish_quantification_glp1r_long$totals == "total_1_2_Glp1r"],sep_char,pcts_ish_quantification_glp1r_long$gene2[pcts_ish_quantification_glp1r_long$totals == "total_1_2"],sep_char,"Glp1r")
# pcts_ish_quantification_glp1r_long =pcts_ish_quantification_glp1r_long[!is.na(pcts_ish_quantification_glp1r_long$count),]
# pcts_ish_quantification_glp1r_wide = pcts_ish_quantification_glp1r_long %>% tidyr::spread(key = new_total,value = count)
# 
# pcts_ish_quantification_glp1r = magically_transform_to_the_right_pcts(ish_quantification_pnoc_sst[,!colnames(ish_quantification_pnoc_sst) %in% "Pnoc/Sst-"],start_col=2,sep_char = "/", max_ref_genes = 2)
# 
# TODO: need to split by experimt to handle wide dataframes properly


##########
### load exports from above
##########

pct_expressed_cells_pnoc_sst_ISH = data.table::fread("figure_outputs/figure_7/pct_expressed_cells_pnoc_sst_ISH.txt")
pct_expressed_cells_pnoc_crabp1_ISH = data.table::fread("figure_outputs/figure_7/pct_expressed_cells_pnoc_crabp1_ISH.txt")
pct_expressed_cells_glp1r_ISH = data.table::fread("figure_outputs/figure_7/pct_expressed_cells_glp1r_ISH.txt")

##########
### dotplot fun
##########

plot_ish =function(plot_data,ylim=120){
  
  p = ggplot2::ggplot(data = plot_data,
                      mapping = ggplot2::aes(
                        x = comp,
                        y = pct
                      )
  ) +
    ggplot2::stat_summary(
      fun.data= ggplot2::mean_se,
      fun.args = list(mult=1),
      geom= "errorbar",
      width = 0.12,
      color = "red"
    ) +
    #width=error_bar_width,size=error_bar_size)
    ggplot2::stat_summary(
      fun = "mean",
      colour = "red",
      size = 0.5
    ) +
    ggplot2::geom_jitter(
      width = 0.05,
      height = 0,
      alpha = 0.6,
      shape = 16,
      size= 2
    ) +
    ggplot2::theme_minimal() +
    ggplot2::ylim(0,ylim) +
    ggplot2::ylab("Pct. positive cells") +
    scale_y_continuous(breaks=c(0,25,50,75,100))+
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 15) ,
      axis.title.y = ggplot2::element_text(size = 20) ,
      axis.text.x = ggplot2::element_text(size = 12),
    )
  
 return(p) 
}

##########
### Plot dotplots glp1r
##########

plot_input = pct_expressed_cells_glp1r_ISH
# remove comparison we don't include
plot_input = plot_input[!plot_input$Experiment %in% c("Npw Nkx2-4","Trh"),]
# make cols for function
plot_input$comp = plot_input$Experiment
plot_input$pct = plot_input$total_2_pct
# manually order accoring to plots in paper : Pomc -- Sst -- Ghrh -- Tbx19/Anxa2 -- Nkx2-4/Trh -- oxt
#plot_input$comp = factor(plot_input$comp,levels = c("Pomc","Pomc Anxa2","Sst","Sst Unc13c","Ghrh","Tbx19 Anxa2","Trh Nkx2-4","Oxt" ))

# make similar to pnoc
plot_input$comp = gsub(" ","/",plot_input$comp)
plot_input$comp = paste0(plot_input$comp,"/Glp1r_of_",plot_input$comp)
addline_format <- function(x){ gsub('_','\n',x)}
plot_input$comp = addline_format(plot_input$comp)

# manually order accoring to plots in paper : Pomc -- Sst -- Ghrh -- Tbx19/Anxa2 -- Nkx2-4/Trh -- oxt
plot_input$comp = factor(plot_input$comp,levels = c("Pomc/Glp1r\nof\nPomc","Pomc/Anxa2/Glp1r\nof\nPomc/Anxa2","Sst/Glp1r\nof\nSst","Sst/Unc13c/Glp1r\nof\nSst/Unc13c","Ghrh/Glp1r\nof\nGhrh" ,"Tbx19/Anxa2/Glp1r\nof\nTbx19/Anxa2","Trh/Nkx2-4/Glp1r\nof\nTrh/Nkx2-4","Oxt/Glp1r\nof\nOxt" ))


# plot ish 
p_glp1r = plot_ish(plot_input,ylim = 110)

#define test levels:
comp_test = list(c("Pomc/Glp1r\nof\nPomc","Pomc/Anxa2/Glp1r\nof\nPomc/Anxa2"),
                 c("Sst/Glp1r\nof\nSst","Sst/Unc13c/Glp1r\nof\nSst/Unc13c"))

# with correction:
pairwise_tests = p_glp1r$data %>% 
  rstatix::wilcox_test(formula = pct ~ comp,comparisons = comp_test) %>% 
  rstatix::adjust_pvalue(method = 'none') # hochberg
p_glp1r = p_glp1r + ggpubr::stat_pvalue_manual(
  data = pairwise_tests, label = "p.adj", 
  y.position = c(90)
)
# without correction:
#p_glp1r = p_glp1r + ggpubr::stat_compare_means(comparisons = comp_test,method = "wilcox.test",label="p.format" ) 

## show
p_glp1r

ggsave(filename = paste0(results_path_figure7,"statplot_glp1r.png"),
       plot = p_glp1r, "png",dpi=450,width=400,height = 150,units="mm")
ggsave(filename = paste0(results_path_figure7,"statplot_glp1r.pdf"),
       plot = p_glp1r, "pdf",dpi=450,width=400,height = 150,units="mm")

# source data
source_figure7_c_dotplot =p_glp1r$data
data.table::fwrite(source_figure7_c_dotplot,paste0(results_path_figure7,"source_figure7_c_glp1r.txt"),sep="\t")

##########
### Plot dotplots Pnoc 
##########

pcts_ish_quantification_pnoc_crabp1_long =  pct_expressed_cells_pnoc_crabp1_ISH %>% tidyr::gather(key="comp",value="pct",-section)
pcts_ish_quantification_pnoc_sst_long =  pct_expressed_cells_pnoc_sst_ISH %>% tidyr::gather(key="comp",value="pct",-section)

plot_input = bind_rows(pcts_ish_quantification_pnoc_sst_long,pcts_ish_quantification_pnoc_crabp1_long)

addline_format <- function(x){gsub('_','\n',x)}
plot_input$comp = addline_format(plot_input$comp)

plot_input$comp = factor(plot_input$comp,unique(plot_input$comp))

p_pnoc = plot_ish(plot_input)

#define test levels:
comp_test_pnoc = list(c("Sst/Pnoc/Nts\nof\nSst/Pnoc","Sst/Pnoc/Unc13c\nof\nSst/Pnoc"),
                      c("Sst/Pnoc/Unc13c\nof\nSst/Pnoc","Sst/Pnoc/Nts/Unc13\nof\nSst/Pnoc"),
                      c("Sst/Pnoc/Nts\nof\nSst/Pnoc","Sst/Pnoc/Nts/Unc13\nof\nSst/Pnoc"),
                      c("Crabp1/Pnoc/Tmem215\nof\nCrabp1/Pnoc","Crabp1/Pnoc/Htr3b\nof\nCrabp1/Pnoc"),
                      c("Crabp1/Pnoc/Tmem215\nof\nCrabp1/Pnoc","Crabp1/Pnoc/Tmem215/Htr3b\nof\nCrabp1/Pnoc"),
                      c("Crabp1/Pnoc/Htr3b\nof\nCrabp1/Pnoc","Crabp1/Pnoc/Tmem215/Htr3b\nof\nCrabp1/Pnoc"))

# with correction:
pairwise_tests = p_pnoc$data %>% 
  rstatix::wilcox_test(formula = pct ~ comp,comparisons = comp_test_pnoc) %>% 
  rstatix::adjust_pvalue(method = 'hochberg')
p_pnoc = p_pnoc + ggpubr::stat_pvalue_manual(
  data = pairwise_tests, label = "p.adj", 
  y.position = c(106, 111, 116)
)
# without correction:
#p_pnoc = p_pnoc + ggpubr::stat_compare_means(comparisons = comp_test_pnoc,method = "wilcox.test",label="p.format" ) 

# show
p_pnoc

ggsave(filename = paste0(results_path_figure7,"statplot_pnoc.png"),
       plot = p_pnoc, "png",dpi=450,width=400,height = 160,units="mm")
ggsave(filename = paste0(results_path_figure7,"statplot_pnoc.pdf"),
       plot = p_pnoc, "pdf",dpi=450,width=400,height = 160,units="mm")

# source data
source_figure7_f_dotplot =p_pnoc$data
data.table::fwrite(source_figure7_f_dotplot,paste0(results_path_figure7,"source_figure7_f_pnoc.txt"),sep="\t")

##########
### stats
##########
std_err <- function(x) sd(x)/sqrt(length(x))

## glp1r
stat_input = pct_expressed_cells_glp1r_ISH
# remove comparison we don't include
stat_input = stat_input[!stat_input$Experiment %in% c("Npw Nkx2-4","Trh"),]
# make cols for function
stat_input$comp = stringr::str_replace(stat_input$Experiment,pattern = " ",replacement = "/") 
stat_input$pct = stat_input$total_2_pct
stat_input$comp = factor(stat_input$comp,levels = c("Pomc","Pomc/Anxa2","Sst","Sst/Unc13c","Ghrh" ,"Tbx19/Anxa2","Trh/Nkx2-4","Oxt" ))
stat_res = stat_input %>% dplyr::arrange(comp) %>% dplyr::group_by(comp) %>% dplyr::summarise(mean = round(mean(pct),2), std_err = round(std_err(pct),2))
stat_res$pasted = paste0(stat_res$comp,": ",stat_res$mean,"±",stat_res$std_err)

paste0("mean ± s.e.m: ",paste0(stat_res$pasted,collapse = "; "))

stat_res

## pnoc + sst
pnoc_sst_text_df = data.frame(comp = names(colMeans(pct_expressed_cells_pnoc_sst_ISH[,2:ncol(pct_expressed_cells_pnoc_sst_ISH)])),
                                 mean = round(colMeans(pct_expressed_cells_pnoc_sst_ISH[,2:ncol(pct_expressed_cells_pnoc_sst_ISH)]),2),
                                 std_err = round(apply(pct_expressed_cells_pnoc_sst_ISH[,2:ncol(pct_expressed_cells_pnoc_sst_ISH)],2,std_err),2))
pnoc_sst_text_df$text_short = stringr::str_split_fixed(pnoc_sst_text_df$comp,"_of_",n=2)[,1]
pnoc_sst_text_df$pasted = paste0(pnoc_sst_text_df$text_short,": ",pnoc_sst_text_df$mean,"±",pnoc_sst_text_df$std_err)

paste0("mean ± s.e.m: ",paste0(pnoc_sst_text_df$pasted,collapse = "; "))

## pnoc + crabp1
pnoc_crabp1_text_df = data.frame(comp = names(colMeans(pct_expressed_cells_pnoc_crabp1_ISH[,2:ncol(pct_expressed_cells_pnoc_crabp1_ISH)])),
                           mean = round(colMeans(pct_expressed_cells_pnoc_crabp1_ISH[,2:ncol(pct_expressed_cells_pnoc_crabp1_ISH)]),2),
                           std_err = round(apply(pct_expressed_cells_pnoc_crabp1_ISH[,2:ncol(pct_expressed_cells_pnoc_crabp1_ISH)],2,std_err),2))
pnoc_crabp1_text_df$text_short = stringr::str_split_fixed(pnoc_crabp1_text_df$comp,"_of_",n=2)[,1]
pnoc_crabp1_text_df$pasted = paste0(pnoc_crabp1_text_df$text_short,": ",pnoc_crabp1_text_df$mean,"±",pnoc_crabp1_text_df$std_err)

paste0("mean ± s.e.m: ",paste0(pnoc_crabp1_text_df$pasted,collapse = "; "))
# comp
# "Pomc/Glp1r\nof\nPomc","Pomc/Anxa2/Glp1r\nof\nPomc/Anxa2"
# comp_test = list(c("Pomc/Glp1r\nof\nPomc","Pomc/Anxa2/Glp1r\nof\nPomc/Anxa2"),c("Sst/Glp1r\nof\nSst","Sst/Unc13c/Glp1r\nof\nSst/Unc13c"))
# p_glp1r + ggpubr::stat_compare_means(comparisons = comp_test,method = "wilcox.test",label="p.format" ) 
