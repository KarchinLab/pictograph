#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=2G
#$ -l h_vmem=3G
rlib="${HOME}/Library/R/3.10-bioc-release"
Rscript="${HOME}/bin/Rscript"
R_LIBS_USER=$rlib $Rscript "mh_trees.R" 
