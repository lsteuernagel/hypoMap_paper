# hypoMap_paper

Figure creation and final analysis scripts for the HypoMap paper. This code does not contain the full pipeline used to build HypoMap, please see the section below for an overview of associated repositories. 

## Associated Repositories

- [scIntegration](https://github.com/lsteuernagel/scIntegration) is the pipeline for evaluation of scVI hyperparameters to optimize the HypoMap model.

- [scHarmonization](https://github.com/lsteuernagel/scHarmonization) is the pipeline to build HypoMap including all downstream clustering and annotation. 

- [scUtils](https://github.com/lsteuernagel/scUtils) is an R package that provides various helper functions to analyze single-cell data that are used in  scIntegration and scHarmonization.

- [scCoco](https://github.com/lsteuernagel/scCoco) is an R package that provides wrapper functions around the [cocoframer R package](https://github.com/AllenInstitute/cocoframer) and is used for the spatial prediction step in scHarmonization.

- [mapscvi](https://github.com/lsteuernagel/mapscvi) is an R package that allows quick mapping of new single cell data onto the existing HypoMap scVI model.

## HypoMap and NucSeq data

Small data input tables are already included in the data_inputs/ folder.

For large data files (like Seurat objects) that are required please see the data at: [Will be ADDED].

Please also see our [interactive cellxgene viewer](https://www.mrl.ims.medschl.cam.ac.uk)

It will be necessary to adjust the path ('large_data_path') to your local directory in the Figure creation scripts when re-running scripts.

## Overview

R/ contains all scripts for Figure creation.

- R/utility_functions.R: Contains various loading and general utility functions used throughout the figure scripts.
- R/plot_functions: A set of plot functions that are used by various scripts.

figure_outputs/ contains tables generated within the Figure scripts that are useful for analysis or required as source data, as well as all pdf and png files of the plots used in the figures of the paper.

table_outputs/ contains all supplementary tables as csv files as well as the final supplementary file as .xlsx

source_outputs/ contains .xlsx files with the source data for each figure.
