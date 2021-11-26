#load data etc
require(tidyverse)
require(data.table)
require(Seurat)
require(ggplot2)

if(file.exists("utils.R")){source("utils.R")}

### Set parameters
global_seed = 123467# seed
map_name = "hypothalamus_neurons_reference"
map_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_objects/"

### Load map
map_seurat_path = paste0(map_path,map_name,".rds")
neuron_map_seurat = readRDS(map_seurat_path)
#neuron_map_seurat@meta.data = cbind(neuron_map_seurat@meta.data ,neuron_map_seurat@misc$other_metdata)

# TODO: update this !??!?
## load full map
project_name = "hypothalamus_full_map"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/harmonization_results/"
seurat_file_name = paste0(project_path,project_name,".h5Seurat")
message(Sys.time(),": Load seurat object.." )
full_map_seurat = SeuratDisk::LoadH5Seurat(seurat_file_name)