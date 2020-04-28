#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pictograph

time Rscript --vanilla $script $opts
