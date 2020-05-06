#!/bin/bash

script=tree_plots.R
logdir=log/$script

mkdir -p $logdir

mh_trees_script=mh_trees_sim1_postAdmat_MC3.R
resultsdir=~/GitHub/pictograph/output/${mh_trees_script}

numIter=100000
thin=100

for seed in {1..5}; do
#for seed in 1; do
	name=trees_numIter1e+05_thin100_seed$seed
    treeResults=$resultsdir/${name}.rds
	jobname=$name
	opts="-r $treeResults -o $resultsdir -n $name"

	sbatch --job-name=$jobname --output=$logdir/${mh_trees_script}_${jobname}.out --export=ALL,script=$script,opts="$opts" slurm-run-rscript.sh
done