##########
### Load & Prepare
##########

require(tidyverse)
require(data.table)
require(ggplot2)

results_path = "figure_outputs/figure_6/"
system(paste0("mkdir -p ",results_path))

# path with output files
data_path = "data_inputs/"

#### load ish results
ish_quantification = data.table::fread(paste0(data_path,"ish_quantification_glp1r_updated.csv"))

### TODO:
# Paul has made the final version of these plots !! (what exactly was his input ?)

##########
### Manually add rows for Pomc and Sst without their second marker
##########

# add pomc alone
ish_quantification_pomc = ish_quantification[ish_quantification$Experiment=="Pomc Anxa2",]
ish_quantification_pomc$Experiment="Pomc"
ish_quantification_pomc$total_1_2=NA
ish_quantification_pomc$total_1_2_Glp1r=NA
ish_quantification = rbind(ish_quantification,ish_quantification_pomc)
# add sst alone
ish_quantification_sst = ish_quantification[ish_quantification$Experiment=="Sst Unc13c",]
ish_quantification_sst$Experiment="Sst"
ish_quantification_sst$total_1_2=NA
ish_quantification_sst$total_1_2_Glp1r=NA
ish_quantification = rbind(ish_quantification,ish_quantification_sst)

# add trh alone
ish_quantification_trh = ish_quantification[ish_quantification$Experiment=="Trh Nkx2-4",]
ish_quantification_trh$Experiment="Trh"
ish_quantification_trh$total_1_2=NA
ish_quantification_trh$total_1_2_Glp1r=NA
ish_quantification = rbind(ish_quantification,ish_quantification_trh)

##########
### make percentages
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
ish_quantification$total_1_pct = ish_quantification$total_1_glp1r / ish_quantification$total_1 * 100
ish_quantification$total_2_pct = ish_quantification$total_1_2_Glp1r / ish_quantification$total_1_2 * 100
ish_quantification$total_2_pct[is.na(ish_quantification$total_2_pct)] = ish_quantification$total_1_pct[is.na(ish_quantification$total_2_pct)]
ish_quantification$total_2_pct[is.nan(ish_quantification$total_2_pct)] = 0

# Example summary
ish_quantification %>% group_by(Experiment) %>% dplyr::summarise(mean(total_2_pct))

##########
### add mouse ids dataframe
##########

## I assume that the first part of Section is the mouse id
ish_quantification$Section[ish_quantification$Section == "49B2section2"] = "49B2 section2"
ish_quantification$mouse_id = sapply(ish_quantification$Section,function(x){strsplit(x," |_")[[1]][1]})

ish_quantification %>% dplyr::group_by(Experiment,mouse_id) %>% dplyr::summarise(n()) %>% dplyr::count()

#ish_quantification_mouse = ish_quantification %>% dplyr::group_by(mouse_id,Experiment) %>% 
#  summarise(across(c(total_1,total_1_2,total_1_glp1r,total_1_2_Glp1r), sum,na.rm=TRUE))

# save result to folder:
data.table::fwrite(ish_quantification,file=paste0(results_path,"pct_expressed_cells_clusters_ISH.txt"),sep="\t")

##########
### Plot dotplots (all)
##########

ish_quantification$x = paste0(ish_quantification$Experiment,"_",ish_quantification$Location)
## all
ggplot(ish_quantification, aes(x=x, y=total_2_pct)) + 
  geom_dotplot(binaxis='y', stackdir='center',dotsize=0.3,fill="grey40",color="grey40") + 
  stat_summary(fun.data= mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) +  # mean_sdl
  stat_summary(fun=mean, geom="point", color="black",size=3)

##########
### Plot function
##########

# simple function to keep plotting succinct.
# Careful: Partly uses global environment vars!
# and also hardcodes labels etc. within function

# function to plot and export:
plot_quant = function(group1,regions = c("rostral","caudal","ARC","VMH","PVH")){
  ish_quantification_current= ish_quantification[ish_quantification$Experiment %in% group1 & ish_quantification$Location %in% c(regions),]
  p_quant = ggplot(ish_quantification_current, aes(x=Experiment, y=total_2_pct)) + 
    stat_summary(fun.data= mean_se, fun.args = list(mult=1),geom="errorbar", color="grey40", width=error_bar_width,size=error_bar_size) + # mean_sdl
    # stat_summary(fun=mean, geom="point", color="grey40",size=point_size_mean)+
    stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - error_bar_width / 4, yend=..y..),size=error_bar_size, color="grey40")+
    stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + error_bar_width / 4, yend=..y..),size=error_bar_size, color="grey40")+
    geom_jitter(size=point_size,shape=21,fill="black",color="darkred",height = 0,width = jitter_width) + 
    ylim(0,100) +ylab("Pct. Glp1r-positive")+ # xlab("Region") +
    theme(text = element_text(size=text_size),panel.background = element_rect(fill = "white"),axis.title.x = element_blank())#+ggtitle(group)
  p_quant
}

##########
### Plot dotplots (per experiment)
##########

# plot params:
text_size = 25
point_size=5
error_bar_width = 0.2
error_bar_size=1.25
jitter_width = 0.02
unique(ish_quantification$Experiment)

### make plots

#pomc
group = c("Pomc","Pomc Anxa2")
p_pomc = plot_quant(group)
p_pomc
ggsave(filename = paste0(results_path,"Pomc","_dotplot.pdf"),
       plot = p_pomc, "pdf",dpi=600,width=880,height = 880,units="mm")

# ghrh
group = c("Ghrh")
p_Ghrh = plot_quant(group)
p_Ghrh
ggsave(filename = paste0(results_path,"Ghrh","_dotplot.pdf"),
       plot = p_Ghrh, "pdf",dpi=600,width=880,height = 880,units="mm")

# oxt
group = c("Oxt")
p_Oxt = plot_quant(group)
p_Oxt
ggsave(filename = paste0(results_path,"Oxt","_dotplot.pdf"),
       plot = p_Oxt, "pdf",dpi=600,width=880,height = 880,units="mm")

# sst
group = c("Sst","Sst Unc13c")
p_Sst = plot_quant(group)
p_Sst
ggsave(filename = paste0(results_path,"Sst","_dotplot.pdf"),
       plot = p_Sst, "pdf",dpi=600,width=880,height = 880,units="mm")

# tbx19
group = "Tbx19 Anxa2"
p_Tbx19 = plot_quant(group)
p_Tbx19
ggsave(filename = paste0(results_path,"Tbx19","_dotplot.pdf"),
       plot = p_Tbx19, "pdf",dpi=600,width=880,height = 880,units="mm")

# trh
group = "Trh Nkx2-4"
p_Trh = plot_quant(group)
p_Trh
ggsave(filename = paste0(results_path,"Trh","_dotplot.pdf"),
       plot = p_Trh, "pdf",dpi=600,width=880,height = 880,units="mm")

# npw
group = "Npw Nkx2-4"
p_Npw = plot_quant(group)
p_Npw
ggsave(filename = paste0(results_path,"Npw","_dotplot.pdf"),
       plot = p_Npw, "pdf",dpi=600,width=880,height = 880,units="mm")


##########
### Save objects for Paul
##########

# I am saving the data as a list and then as an rds with a single object.
export_list = list(ish_quantification = ish_quantification,p_pomc = p_pomc,p_Ghrh=p_Ghrh, p_Oxt = p_Oxt, p_Npw = p_Npw, p_Trh = p_Trh, p_Tbx19 = p_Tbx19, p_Sst =p_Sst)

# test
export_list$p_Tbx19

# save
saveRDS(export_list,paste0(results_path,"Figure_4_data.rds"))
