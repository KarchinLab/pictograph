#!/bin/bash

burninscript=mh_trees_sim1.R
script=mh_trees_sim1_postAdmat_MC3.R
logdir=log/$script

mkdir -p $logdir

postadmatdir=~/GitHub/pictograph/output/$burninscript
outdir=~/GitHub/pictograph/output/$script

numIter=100000
thin=100

for seed in {1..5}; do
#for seed in 1; do
	jobname=trees_numIter1e+05_thin${thin}_seed$seed
	opts="-s $seed -o $outdir -i $numIter -t $thin -p $postadmatdir/${jobname}_post.admat.tsv"

	sbatch --job-name=$jobname --output=$logdir/${jobname}.out --export=ALL,script=$script,opts="$opts" slurm-run-rscript.sh
done
