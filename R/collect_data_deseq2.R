##########
### Define paths
##########

small_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/hypoMap_paper/data_inputs/"
large_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/"
# define some relevant files for loading of data
key="hypothalamusMapNeurons_v4"
scHarmonize_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/"

##########
### transform results
##########

#corinna glp1r:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2021-04-corinna-traps/release/dereportr/hypo/deseq_diff/deseq2_diff.csv"
corinna_deseq2_diff_bac_hypo = as.data.frame(fread(paste0(datafile),data.table = F))
# save as glp1r
data.table::fwrite(corinna_deseq2_diff_bac_hypo,paste0(small_data_path,"bacTRAP_deseq2_glp1r.csv"))

#alex pnoc:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2019-02-alex-trap-nfcore/release/trapdiff/bactrap/de.rds"
alex_deseq2_diff_bac_pnoc = as.data.frame(readRDS(paste0(datafile)))
alex_deseq2_diff_bac_pnoc = alex_deseq2_diff_bac_pnoc[alex_deseq2_diff_bac_pnoc$comparison == "ip_input_cd",]
# save as pnoc
data.table::fwrite(alex_deseq2_diff_bac_pnoc,paste0(small_data_path,"bacTRAP_deseq2_pnoc.csv"))

#nasim pomc lepr:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-02-nasim-bactrap/release/deseq/de/lepr_deseq2_diff.xlsx"
nasim_deseq2_diff_bac_pomc_lepr = as.data.frame(readxl::read_xlsx(paste0(datafile),sheet = 1))
# save as pomc/lepr
data.table::fwrite(nasim_deseq2_diff_bac_pomc_lepr,paste0(small_data_path,"bacTRAP_deseq2_pomc_lepr.csv"))

#nasim pomc glp1r:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-02-nasim-bactrap/release/deseq/de/glp1r_deseq2_diff.xlsx"
nasim_deseq2_diff_bac_pomc_glp1r =  as.data.frame(readxl::read_xlsx(paste0(datafile),sheet = 1))
# save as pomc/lepr
data.table::fwrite(nasim_deseq2_diff_bac_pomc_glp1r,paste0(small_data_path,"bacTRAP_deseq2_pomc_glp1r.csv"))

#kasia agrp:
# datafile = "/beegfs/scratch/bruening_scratch/pklemm/2020-05-kasia-bactrap/release/trapdiff/fasted_vs_control/de.rds"
# kasia_bac_agrp=  as.data.frame(readRDS(paste0(datafile)))
# kasia_bac_agrp = kasia_bac_agrp[kasia_bac_agrp$comparison == "ip_input_control",]
# # save as agrp
# data.table::fwrite(kasia_bac_agrp,paste0(small_data_path,"bacTRAP_deseq2_agrp.csv"))

# almu agrp
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2021-08-almu-bactrap/release/trapdiff/fed_ko_control/de.rds"
almu_bac_agrp=  as.data.frame(readRDS(paste0(datafile)))
almu_bac_agrp = almu_bac_agrp[almu_bac_agrp$comparison == "ip_input_control_fed",]
# save as agrp
data.table::fwrite(almu_bac_agrp,paste0(small_data_path,"bacTRAP_deseq2_agrp.csv"))

#alai Pomc agrp:
datafile = "/beegfs/scratch/bruening_scratch/pklemm/2021-05-alai-bactrap/release/trapdiff/dre_cre/de.rds"
alai_bac_pomc =  as.data.frame(readRDS(paste0(datafile)))
alai_bac_pomc = alai_bac_pomc[alai_bac_pomc$comparison =="ip_input_cre",]
# save as pomc
data.table::fwrite(alai_bac_pomc,paste0(small_data_path,"bacTRAP_deseq2_pomc.csv"))



