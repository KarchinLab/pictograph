#!/bin/bash

#script=IP30_mh_trees_weightedMutate10.R
script=IP30_mh_trees.R
logdir=log/$script

mkdir -p $logdir

outdir=/home/lzheng/GitHub/pictograph/output/$script

numIter=100000
thin=100

for seed in {1..5}; do
#for seed in 1; do
	jobname=numIter${numIter}_seed${seed}_thin${thin}
	opts="-s $seed -o $outdir -i $numIter -t $thin"

	sbatch -c 8 --job-name=$jobname --output=$logdir/${jobname}.out --export=ALL,script=$script,opts="$opts" slurm-run-rscript.sh
done
