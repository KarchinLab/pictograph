#!/bin/bash

script=mh_trees_sim1.R
logdir=log/mh_trees_sim1
mkdir -p $logdir

outdir=/home/lzheng/GitHub/pictograph/output/hn07_mh_trees_sim1

numIter=100000

for seed in {1..5}; do
	jobname=numIter${numIter}_seed$seed
	opts="-s $seed -o $outdir -i $numIter"

	sbatch --cpus-per-task=12 --ntasks=1 --job-name=$jobname --output=$logdir/${jobname}.out --export=ALL,script=$script,opts="$opts" slurm-run-rscript.sh
done
