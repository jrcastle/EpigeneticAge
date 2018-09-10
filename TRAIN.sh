#!/bin/sh
#sbatch --job-name=Train --partition=FatComp --mail-type ALL --mail-user jrca253@uky.edu ./TRAIN.sh
Rscript trainModel.R