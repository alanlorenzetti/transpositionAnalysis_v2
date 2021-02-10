#!/bin/bash

# alorenzetti 20210209

# description ####
# this is the scaffold script
# running a collection of other scripts
# to analyze an Oxford Nanopore DNA-Seq
# experiment
#
# the main task is to find insertions 
# and excisions of transposable elements
# using a reference genome and a database
# of known insertion sequences

bash ./scripts/01_run_guppy.sh
bash ./scripts/02_run_deepbinner.sh
bash ./scripts/03_run_porechop.sh
bash ./scripts/04_run_sniffles_ref.sh
bash ./scripts/05_run_filtlong_stats.sh
#Rscript ./scripts/06_run_filtlong_stats_ggplot2.R # called inside 05_run_filtlong_stats.sh
bash ./scripts/07_run_tombo.sh
Rscript ./scripts/08_20190826_is_clusters_perbarcode_plots.R
Rscript ./scripts/09_20191004_circlize.R
bash ./scripts/10_run_computeGC_ref20181106.sh
