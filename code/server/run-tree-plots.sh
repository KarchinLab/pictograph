#!/bin/bash

script=tree_plots.R
logdir=log/$script

mkdir -p $logdir

mh_trees_script=mh_trees_sim1_mutate10.R
mh_trees_script=IP30_mh_trees_weightedMutate10.R
#outdir=~/GitHub/pictograph/output/$script
resultsdir=~/GitHub/pictograph/output/${mh_trees_script}

numIter=100000
thin=100

for seed in {1..5}; do
#for seed in 1; do
	name=trees_numIter1e+05_thin100_seed$seed
    treeResults=$resultsdir/${name}.rds
	jobname=$name
	opts="-r $treeResults -o $resultsdir -n $name"

	sbatch -c 8 --job-name=$jobname --output=$logdir/${mh_trees_script}_${jobname}.out --export=ALL,script=$script,opts="$opts" slurm-run-rscript.sh
done
