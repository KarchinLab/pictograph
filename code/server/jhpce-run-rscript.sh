#!/bin/bash
#$ -m e
#$ -M lzheng163@gmail.com
#$ -cwd

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pictograph
echo "script = $script"
echo "opts = $opts"

time Rscript --vanilla $script $opts
