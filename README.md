# WES_analysis

* The scripts were used to analyze Whole Exome Sequencing data

* The metadata and data belongs to the Netherlands Cancer Institue

* The scripts belongs to Yejin Park (yejin.park1016@gmail.com)

* Contents


  1) General_analysis.Rmd

Using 'maftools' library to create oncoplot(co-oncoplot) and VAF distribution comparison
  
  2) VAF_plot_basic_functions.R

Basic functions needed in building VAF plots

  3) CNAqc_nonsyno_mut_VAF_plot.Rmd

Compare shared mutatoins between samples on VAF to chromosome map
Functions are sourced from VAF_plot_basic_functions.R

  4) VAF_group_analysis.Rmd

Normalized VAF comparison between treatment responding and non-responding group

  5) FACETS_to_Pyclone_input.py

Originated from https://github.com/vanallenlab/PyCloneTSVGeneration_FACETS_or_TITAN/blob/master/generate_wgs_pyclone_input.py
Using only FACET output and maf files to create PyClone input format

  6) Pyclone_to_citup.Rmd

Use FACET Snakemake output to prepare PyClone input, and later use PyClone output to CITUP and compare the clusters and clonality tree

  7) Timescape_tree_visualization.Rmd

Use Pyclone-vi output to build a pylogenetic tree using CITUP and Timescape, but can be used for PyClone output

  8) cooncoplot_analysis.Rmd

Using 'maftools' library, build co-oncoplot on primary tumor and lymph nodes
