# make conda environment using environment.yml
conda env create -f environment.yml

# tree MH (with weighted mutate, 10^tier) on simulated data
mh_trees_sim1_mutate10.R 
# tree MH (plain mutate)
mh_trees_sim1.R

# run tree mh --  5 different chains
./run-trees-sim1.sh

# plot tree mh results -- change mh_trees script
./run-tree-plots.sh
