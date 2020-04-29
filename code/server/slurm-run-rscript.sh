#!/bin/bash
#SBATCH --mail-user=lzheng163@gmail.com
#SBATCH --mail-type=END

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pictograph

time Rscript --vanilla $script $opts
