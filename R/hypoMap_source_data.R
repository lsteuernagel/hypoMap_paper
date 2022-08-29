##########
### Load & Prepare
##########

# I don't want to have excel screw up numbers in scientifc notation:
#options(scipen = 999)
#fwrite_scipen = 999 # I'll leave it one for some sc-seq tables because they are difficult to read else

results_path_source = "source_outputs/"
system(paste0("mkdir -p ",results_path_source))

# load required functions
require(dplyr)
require(ggplot2)
require(Seurat)
#require(openxlsx)
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
### Load source files per figure dir and save as excel
##########


all_figure_dirs = list.dirs("figure_outputs/",recursive = F)
# all_figure_dirs = all_figure_dirs[,c(2,3,8:16)]

for(dir in all_figure_dirs){
  print(dir)
  all_source_files = list.files(dir,pattern = "source_",full.names = TRUE) 
  if(dir == "figure_outputs//figure_extended_1" ){all_source_files = "figure_outputs//figure_extended_1/source_ext_figure1_b_metrics.txt"}
  all_source_files_short = list.files(dir,pattern = "source_",full.names = FALSE) 
  if(dir == "figure_outputs//figure_extended_1" ){all_source_files_short = "source_ext_figure1_b_metrics.txt" }
  source_table_list = list()
  for(i in 1:length(all_source_files)){
    current_table = data.table::fread(paste0(all_source_files[i]),data.table = FALSE) 
    current_name = gsub("source_figure[0-9]+_|source_ext_figure[0-9]+_|\\.txt","",all_source_files_short[i])
    source_table_list[[current_name]] = current_table
  }
  # shorten colnames
  names(source_table_list)[sapply(names(source_table_list),nchar) > 31] = sapply(names(source_table_list)[sapply(names(source_table_list),nchar) > 31],substr,start=1,stop=31)
  # filename
  filename = gsub("figure_outputs//","",dir)
  #write
  WriteXLS::WriteXLS(x = source_table_list,ExcelFileName = paste0(results_path_source,filename,"_source.xlsx"),col.names=TRUE)
}



