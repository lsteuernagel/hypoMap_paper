
##########
### Prepare gene input
##########

results_path_extended_figure6 = "figure_outputs/figure_extended_6/"

library(tidyverse)

#setwd("/Users/lsteuernagel/Documents/projects/hypoMap_paper/")

### Be careful: Check transcript not only gene to get accurate TSS!
background_genes_all = data.table::fread(paste0("figure_outputs/figure_4/","per_gene_correlations.txt"),data.table=F)
ieg_genes = data.table::fread("data_inputs/immediate_early_genes_wuetal.txt")$gene

# make subsample for refernce during enrichemnt
nsample = 1000
set.seed(123)
background_genes_sample = sample(x=background_genes_all$gene[!background_genes_all$gene %in% c("1700016P03Rik",ieg_genes)],nsample)

# define genes we need to retrive promoters for:
genes_of_interest = c("1700016P03Rik",ieg_genes[ieg_genes %in% background_genes_all$gene],background_genes_sample)
genes_of_interest = unique(genes_of_interest)

##########
###  Retrieve refTSS_v3.3 promoter positions
##########

# these files are not included in this repo ---> download and adjust filepath

### load files
# https://www.sciencedirect.com/science/article/pii/S0022283619302530
# http://reftss.clst.riken.jp/datafiles/current/mouse/

# load TSS positions
refTSS_v3_ids = data.table::fread("/Users/lsteuernagel/Downloads/refTSS_v3.3_mouse.ids.list",data.table=FALSE,header=FALSE)
# need to format and extract ids plus positions
refTSS_v3_ids_formatted = refTSS_v3_ids %>% dplyr::mutate(
  tss_id = stringr::str_extract(V2,pattern= "mm_[0-9]+\\.[0-9]+"),
  chromosome= stringr::str_extract(V1,pattern= "chr[0-9]+"),
  start = as.numeric(stringr::str_remove(stringr::str_extract(V1,":[0-9]+"),pattern=":|-")),
  end = as.numeric(stringr::str_remove(stringr::str_extract(V1,"-[0-9]+"),pattern=":|-")),
  transcription_start_site = (start+end)/2,
  strand = stringr::str_extract(V2,pattern= "\\+|\\-")
) %>% dplyr::select(-V1,-V2)
#refTSS_v3_ids_formatted$chromosome = stringr::str_extract(refTSS_v3_ids_formatted$V1,pattern= "chr[0-10]+")
# load gene to TSS id file
refTSS_v3_gene_annotations = data.table::fread("/Users/lsteuernagel/Downloads/refTSS_v3.3_mouse_annotation.txt",data.table=FALSE)
colnames(refTSS_v3_gene_annotations)[1] = "tss_id"
# extract TTS for genes of interest 
refTSS_v3_gene_annotations_selected = refTSS_v3_gene_annotations %>% dplyr::select(tss_id,Transcript_name,gene =Gene_symbol) %>%
  dplyr::filter(gene %in% genes_of_interest)
# add coordinates
selected_gene_TSS = dplyr::left_join(refTSS_v3_gene_annotations_selected,refTSS_v3_ids_formatted,by=c("tss_id"="tss_id"))

# reduce overlapping promoters 
min_diff = 500
selected_gene_TSS = selected_gene_TSS %>% dplyr::group_by(gene) %>% dplyr::mutate(
  tss_diff = transcription_start_site - dplyr::lag(transcription_start_site , n=1)- min_diff,
  tss_diff = case_when(!is.na(tss_diff) ~ tss_diff, TRUE ~ 0)
) %>% dplyr::filter(tss_diff >=0)

# rename into mouse_gene_info
mouse_gene_info = selected_gene_TSS %>% dplyr::select(gene, Transcript_name,tss_id, chrom = chromosome, transcription_start_site, strand) %>%
  na.omit()

##########
###  Use BSgenome to retrieve sequence information from the mouse genome (here: promoters)
##########

# This retrieves all DNA sequences for the promoters defined above

# we retrieve + strand only (genes on - are flipped) --> TFBS tools can then still query both strands!

library("BSgenome.Mmusculus.UCSC.mm10")
library(Rsamtools)
library(BSgenome)
library(stringi)
library(dplyr)

promoter_upstream = 800
promoter_downstream = 200

for(i in 1:nrow(mouse_gene_info)){
  if(i%%50==0) print(i)
  if(mouse_gene_info$strand[i]=="+"){
    promotor_seq = as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, names =mouse_gene_info$chrom[i], 
                                       start=mouse_gene_info$transcription_start_site[i]-(promoter_upstream-1), 
                                       end=mouse_gene_info$transcription_start_site[i]+(promoter_downstream+1),
                                       strand="+"))
  }else{
    promotor_seq = as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, mouse_gene_info$chrom[i],
                                       start=mouse_gene_info$transcription_start_site[i]-(promoter_downstream+1),
                                       end=mouse_gene_info$transcription_start_site[i]+(promoter_upstream-1),
                                       strand="+")) # I flip the start and end when it is the minus strand but i still retrieve the plus strand sequence --> then query both strands with TFBS tools!
  }
  mouse_gene_info$promoter_seq[i] = promotor_seq
}


##########
###  Use TFBS tools to find TF sites
##########

# We query JASPAR motifs for CREB1 (and SRF?)

# load TFBSTools
library(TFBSTools)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)

# load the JASPAR motif database
library(JASPAR2020)

# extract motifs corresponding to vertebrates transcription factors
pwm_library = TFBSTools::getMatrixSet(
  JASPAR2020,
  opts=list(
    all_versions =TRUE,
    collection = 'CORE',
    matrixtype = 'PWM',
    tax_group = "vertebrates"
  ))

# subset library if wanted ---> CREB1: MA0018 , MA0083: SRF , CREB3: MA0638, CREB-related factors: MA0840
selected_pwms = names(pwm_library)[grepl("MA0018|MA0638|MA0840|MA0083",names(pwm_library))]
# pwm_library_to_use = pwm_library
pwm_library_to_use = pwm_library[selected_pwms]

#represents the quantile between the minimal and the maximal possible value from the PWM
min_score = "80%"

## find selected pwms in all selected genes
all_stats_list = list()
for(i in 1:nrow(mouse_gene_info)){
  if(i%%50==0) print(i)
  # subject <- Biostrings::DNAString("GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG")
  current_sequence <- Biostrings::DNAString(mouse_gene_info$promoter_seq[i])
  # find motifs --> search both strands
  sitesetList <- TFBSTools::searchSeq(x = pwm_library_to_use, current_sequence, seqname=mouse_gene_info$Transcript_name[i],min.score=min_score, strand="*")
  for(siteSet_name in names(sitesetList)){
    temp = as.data.frame(writeGFF3(sitesetList[[siteSet_name]]))
    if(nrow(temp)>0){
      temp$relScore = relScore(sitesetList[[siteSet_name]])
      temp$siteSet_name = siteSet_name
      all_stats_list[[paste0(i,siteSet_name)]] = temp
    }
  }
  
}
all_stats = do.call(rbind,all_stats_list)

## add gene info
gene_info_short = mouse_gene_info[,c("tss_id","Transcript_name","gene")] %>%
  dplyr::group_by(gene) %>% dplyr::add_count(name="n_transcripts") 
all_stats = dplyr::full_join(all_stats,gene_info_short,by=c("seqname"="Transcript_name"))
## add pwm info
pwm_mapping=data.frame(pwm_id = names(name(pwm_library_to_use)),pwm_name=name(pwm_library_to_use))
all_stats = dplyr::full_join(all_stats,pwm_mapping,by=c("siteSet_name"="pwm_id"))

#data.table::fwrite(all_stats,"data_inputs/all_pwm_enrichment_result")
#all_stats = data.table::fread("table_outputs/all_pwm_enrichment_result",data.table=F)

##########
###  Summarize and visualize results
##########

# colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

text_size = 20
repel_label_size =6

library(ggplot2)
library(cowplot)


## how many per promoter/transcript/gene
pwm_per_gene = all_stats %>% dplyr::group_by(gene, siteSet_name) %>% dplyr::add_count(name="total_pwms_in_gene") %>% 
  dplyr::mutate(sum_score = sum(relScore)) %>%
  dplyr::distinct(gene, siteSet_name,.keep_all=TRUE) %>% 
  dplyr::select(gene, siteSet_name, pwm_name, total_pwms_in_gene, n_transcripts,sum_score) %>%
  dplyr::mutate(sum_score_per_transcript = sum_score / n_transcripts) %>%
  dplyr::mutate(pwms_per_transcript = total_pwms_in_gene / n_transcripts)

# add classification
pwm_per_gene$class = "Background"
pwm_per_gene$class[pwm_per_gene$gene %in% c("1700016P03Rik",ieg_genes[ieg_genes %in% background_genes_all$gene])] = "IEGs"

## save:
pwms_of_interest = c("MA0018.1","MA0018.2")
data.table::fwrite(pwm_per_gene[pwm_per_gene$siteSet_name %in% pwms_of_interest,c("gene", "class","siteSet_name","pwm_name","total_pwms_in_gene","n_transcripts","pwms_per_transcript")],
                   paste0(results_path_extended_figure6,"tfbs_detected_per_gene_promoter.csv"))

pwm_per_gene = data.table::fread(paste0(results_path_extended_figure6,"tfbs_detected_per_gene_promoter.csv"),data.table=F)

# for MA0018.1
pwms_of_interest = c("MA0018.1","MA0018.2")
gois= c("Fos","Jund","1700016P03Rik","Jun","Egr1")#,"Btg2","Gem","Btg2","Noct","Sgk1")
plot_list = list()
for(pwm in pwms_of_interest){
  
  # add genes that have 0 TFBSs to density estimation:
  current_pwm_stat = pwm_per_gene[pwm_per_gene$siteSet_name %in%pwm ,]
  add_0_genes = gene_info_short$gene[!gene_info_short$gene %in% current_pwm_stat$gene]
  if(length(add_0_genes)>0){
    add_0_genes = data.frame(gene = add_0_genes, pwms_per_transcript = 0)
    add_0_genes$class = "Background"
    add_0_genes$class[add_0_genes$gene %in% c("1700016P03Rik",ieg_genes[ieg_genes %in% background_genes_all$gene])] = "IEGs"
    current_pwm_stat = bind_rows(current_pwm_stat,add_0_genes)
  }
  
  # plot:
  plot_list[[pwm]] = ggplot2::ggplot(current_pwm_stat,aes(x=pwms_per_transcript,fill=class))+geom_density(alpha=0.5)+
    geom_point(data=current_pwm_stat[current_pwm_stat$gene %in% gois,],mapping = aes(x=pwms_per_transcript,y=0.05),show.legend = FALSE)+
    ggrepel::geom_label_repel(data=current_pwm_stat[current_pwm_stat$gene %in% gois,],mapping = aes(x=pwms_per_transcript,y=0.05,label = gene),angle = 45,max.overlaps=10,size=repel_label_size,show.legend = FALSE)+
    xlab("PWM per non-overlapping promoter")+ylab("Density")+  #"Summed adjusted p-value"
    labs(fill='Gene')+scale_color_manual(values=short_palette)+scale_fill_manual(values=short_palette)+
    ggtitle(paste0("CREB1: ",pwm))+
    theme( text = element_text(size=text_size),
           axis.title.x = element_text(size = text_size),
           panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90"),
           panel.background = element_rect(fill="white"))#+ theme_light()
}

pb  = cowplot::plot_grid(plotlist=plot_list,ncol=1)
pb


#save
ggsave(filename = paste0(results_path_extended_figure6,"creb1_TFBS_density.png"),
       plot = pb, "png",dpi=450,width=300,height = 250,units="mm")
ggsave(filename = paste0(results_path_extended_figure6,"creb1_TFBS_density.pdf"),
       plot = pb, "pdf",dpi=450,width=300,height = 250,units="mm")


## source data
source_ext_figure6b_1 = plot_list[[1]]$data
data.table::fwrite(source_ext_figure6b_1,paste0(results_path_extended_figure6,"source_ext_figure6_b_tfbs_1.txt"),sep="\t")

source_ext_figure6b_2 = plot_list[[2]]$data
data.table::fwrite(source_ext_figure6b_2,paste0(results_path_extended_figure6,"source_ext_figure6_b_tfbs_2.txt"),sep="\t")

#minor