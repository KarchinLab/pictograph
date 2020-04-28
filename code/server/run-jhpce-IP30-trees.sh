#!/bin/bash

script=IP30_mh_trees_weightedMutate10.R
logdir=log/jhpce_$script

mkdir -p $logdir

outdir=~/GitHub/pictograph/output/jhpce_$script

numIter=100000
thin=100

mem=4G

#for seed in {1..5}; do
for seed in 1; do
    jobname=numIter${numIter}_seed${seed}_thin${thin}
    opts="-s $seed -o $outdir -i $numIter -t $thin"

    qsub -cwd -v script=$script,opts="$opts" -l mem_free=$mem,h_vmem=$mem -N $jobname -o $logdir/${jobname}.o -e $logdir/${jobname}.e jhpce-run-rscript.sh
done
