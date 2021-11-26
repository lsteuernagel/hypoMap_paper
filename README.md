# hypoMap_paper
Figure creation and final analysis scripts for the hypoMap paper.

### Required data

Small data tables are already included in the data_inputs/ folder.

For large data files (like Seurat objects) that are required please see the data at: XXX.XX

It might be necessary to adjust the paths to these data to your local directory!

### Overview

R/ contains all scripts for Figure creation.

- R/collect_data.R is a help to gather the data in one place (Might be removed before publication)
- R/load_data.R is sa simple helper that loads the most used standard data (to avoid duplication of this code at the beginning of each script)
- R/plot_functions: A set of plot functions that are used by various scripts.

data_outputs/ contains some of the generated tables within the Figure scripts that are useful for analysis

plot_outputs/ contains pdf and png files of the results.