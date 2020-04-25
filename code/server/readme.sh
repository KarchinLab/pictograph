# make conda environment using environment.yml
conda env create -f environment.yml

# tree MH (with weighted mutate, 10^tier) on simulated data -- 5 different chains
mh_trees_sim1_mutate10.R 
./run-trees-sim1.sh
