# hypoMap_paper

Figure creation and final analysis scripts for the hypoMap paper.

### Required data

Small data tables are already included in the data_inputs/ folder.

For large data files (like Seurat objects) that are required please see the data at: XXX.XX

It will be necessary to adjust the path ('large_data_path') to your local directory in the Figure creation scripts!


### Overview

R/ contains all scripts for Figure creation.

- R/collect_data.R is a help to gather the data in one place (Might be removed before publication)
- R/utility_functions.R: Contains various loading and general utility functions used throughout the figure scripts.
- R/plot_functions: A set of plot functions that are used by various scripts.

figure_outputs/ contains some of the generated tables within the Figure scripts that are useful for analysis, as well as all pdf and png files of the plots used in the figures of the paper.

![HypoMap UMAP (from Figure 1)](figure_outputs/figure_1/neuron_HypoMap_annotated_clusters.pdf.png)